import sys
import os
import signal
import subprocess
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QVBoxLayout,
    QGridLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QWidget,
    QFileDialog,
    QTextEdit,
)
from PySide6.QtCore import QThread, Signal, Qt
from PySide6.QtGui import QFont


class PipelineThread(QThread):
    log_signal = Signal(str)  # For logging messages
    completed_signal = Signal()  # For notifying when the pipeline is complete
    running_process = None  # Reference to the running process

    def __init__(self, form_data):
        super().__init__()
        self.form_data = form_data
        self.stopped = False

    def convert_to_wsl_path(self, path):
        """Convert a Windows-style path to a WSL-compatible Linux-style path."""
        try:
            result = subprocess.run(
                ["wsl", "wslpath", "-a", path],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True
            )
            return result.stdout.strip()
        except subprocess.CalledProcessError as e:
            self.log_signal.emit(f"⚠️ Error converting path to WSL: {e.stderr.strip()}")
            return path  # Return the original path if conversion fails

    def run(self):
        try:
            # Prepare reference folder path
            ref_folder = self.form_data["ref_folder"]
            if not ref_folder.endswith("/"):
                ref_folder += "/"

            # Prepare results folder path with /DiffExpResults/
            results_folder = self.form_data["results_folder"]
            if not results_folder.endswith("/"):
                results_folder += "/"
        


            # Prepare results subdirectories
            tkd_dir = os.path.join(results_folder, "Diff_Exp_Results/")
            sam_dir = os.path.join(tkd_dir, "sam/")
            bam_dir = os.path.join(tkd_dir, "bam/")
            counts_dir = os.path.join(tkd_dir, "counts/")
            os.makedirs(sam_dir, exist_ok=True)
            os.makedirs(bam_dir, exist_ok=True)
            os.makedirs(counts_dir, exist_ok=True)

            # Derive file prefixes and paths
            fastq_prefix = os.path.basename(self.form_data["fastq1"]).split("_")[0]
            reference_index = os.path.join(
                ref_folder,
                os.path.splitext(os.path.basename(self.form_data["reference_genome"]))[0]
            )
            sam_file = os.path.join(sam_dir, f"{fastq_prefix}_aligned.sam")
            sam_sorted_file = os.path.join(sam_dir, f"{fastq_prefix}_aligned_sorted.sam")
            sorted_bam = os.path.join(bam_dir, f"{fastq_prefix}_aligned.sorted.bam")
            counts_file = os.path.join(counts_dir, f"{fastq_prefix}_gene_expression.txt")

            def run_command(command, description):
                """Run a shell command inside WSL and redirect output to the app."""
                try:
                    self.log_signal.emit(f"🔄 Running: {description}...")

                    # Prepare the command to run inside WSL
                    wsl_command = f"wsl bash -c \"{command}\""

                    # Use platform-specific process handling
                    self.running_process = subprocess.Popen(
                        wsl_command,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=True,
                        text=True,
                        creationflags=subprocess.CREATE_NEW_PROCESS_GROUP
                    )

                    for line in iter(self.running_process.stdout.readline, ""):
                        if self.stopped:
                            self.running_process.terminate()
                            self.log_signal.emit("❌ Process terminated by user.")
                            return
                        self.log_signal.emit(line.strip())

                    for err_line in iter(self.running_process.stderr.readline, ""):
                        self.log_signal.emit(f"⚠️ {err_line.strip()}")

                    self.running_process.wait()
                    if self.running_process.returncode != 0:
                        raise subprocess.CalledProcessError(self.running_process.returncode, command)

                    self.log_signal.emit(f"✅ Completed: {description}.")
                except subprocess.CalledProcessError as e:
                    self.log_signal.emit(f"❌ Error in {description}: {e}")
                    self.stopped = True
                    return

            # Convert all paths to WSL-compatible paths
            ref_genome_path = self.convert_to_wsl_path(self.form_data["reference_genome"])
            ref_index_path = self.convert_to_wsl_path(reference_index)
            fastq1_path = self.convert_to_wsl_path(self.form_data["fastq1"])
            fastq2_path = self.convert_to_wsl_path(self.form_data["fastq2"])
            gtf_file_path = self.convert_to_wsl_path(self.form_data["gtf_file"])
            sam_file_path = self.convert_to_wsl_path(sam_file)
            sorted_sam_path = self.convert_to_wsl_path(sam_sorted_file)
            sorted_bam_path = self.convert_to_wsl_path(sorted_bam)
            counts_file_path = self.convert_to_wsl_path(counts_file)


            # Run the pipeline steps
            run_command(f"hisat2-build -p {self.form_data['threads']} {ref_genome_path} {ref_index_path}", "Building HISAT2 index")
            if self.stopped:
                return

            run_command(f"hisat2 -p {self.form_data['threads']} -x {ref_index_path} -1 {fastq1_path} -2 {fastq2_path} -S {sam_file_path}", "Aligning reads")
            if self.stopped:
                return

            run_command(f"samtools sort -@ {self.form_data['threads']} {sam_file_path} -o {sorted_sam_path}", "Sorting SAM file")
            if self.stopped:
                return

            run_command(f"samtools view -@ {self.form_data['threads']} -S -b {sorted_sam_path} > {sorted_bam_path}", "Converting SAM to BAM")
            if self.stopped:
                return

            run_command(f"samtools index -@ {self.form_data['threads']} {sorted_bam_path}", "Indexing BAM file")
            if self.stopped:
                return

            run_command(f"featureCounts -T {self.form_data['threads']} -p -a {gtf_file_path} -o {counts_file_path} -t {self.form_data['t_option']} -g {self.form_data['g_option']} {sorted_bam_path}", "Quantifying read counts")
            if self.stopped:
                return


            self.log_signal.emit("🎉 Analysis complete.")
            self.completed_signal.emit()

        except Exception as e:
            self.log_signal.emit(f"❌ Error: {str(e)}")
            self.completed_signal.emit()  # Notify MainWindow that the process ended

    def stop(self):
        self.stopped = True
        if self.running_process:
            self.running_process.terminate()
            self.log_signal.emit("❌ Process forcefully terminated.")
        self.running_process = None  # Clear the reference to avoid reusing a terminated process

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("NGS MADE EASY BY - KTMDRVD")
        self.setGeometry(100, 100, 900, 700)
        self.setStyleSheet("background-color: white;")
        self.pipeline_thread = None
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout()

        # Title
        title = QLabel("Transcriptomics Pipeline (Diff. Exp.)")
        title.setFont(QFont("Sans Serif", 20, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title)

        # Grid layout for input fields
        grid_layout = QGridLayout()

        # Input fields
        self.ref_folder = self.create_input_with_browse(grid_layout, "Reference Folder:", 0, True)
        self.results_folder = self.create_input_with_browse(grid_layout, "Results Folder:", 1, True)
        self.fastq1 = self.create_input_with_browse(grid_layout, "FASTQ File 1:", 2, False)
        self.fastq2 = self.create_input_with_browse(grid_layout, "FASTQ File 2:", 3, False)
        self.reference_genome = self.create_input_with_browse(grid_layout, "Reference Genome:", 4, False)
        self.gtf_file = self.create_input_with_browse(grid_layout, "GTF/GFF File:", 5, False)

        # Subheading for parameters
        gtf_subheading = QLabel("Parameters for extraction based on GTF/GFF file: Based on FeatureCount")
        gtf_subheading.setFont(QFont("Sans Serif", 10, QFont.Bold))
        gtf_subheading.setAlignment(Qt.AlignLeft)
        grid_layout.addWidget(gtf_subheading, 6, 0, 1, 2)

        self.t_option = QLineEdit("exon")
        grid_layout.addWidget(QLabel("-t Option:"), 7, 0)
        grid_layout.addWidget(self.t_option, 7, 1)

        self.g_option = QLineEdit("Parent")
        grid_layout.addWidget(QLabel("-g Option:"), 8, 0)
        grid_layout.addWidget(self.g_option, 8, 1)

        self.threads = QLineEdit("40")
        grid_layout.addWidget(QLabel("Threads:"), 9, 0)
        grid_layout.addWidget(self.threads, 9, 1)

        main_layout.addLayout(grid_layout)

        # Log display
        self.log_display = QTextEdit()
        self.log_display.setReadOnly(True)
        main_layout.addWidget(self.log_display)

        # Start button
        self.start_button = QPushButton("Start Analysis")
        self.start_button.clicked.connect(self.toggle_pipeline)
        main_layout.addWidget(self.start_button)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

    def create_input_with_browse(self, grid_layout, label, row, folder):
        lbl = QLabel(label)
        inp = QLineEdit()
        btn = QPushButton("Browse")
        btn.clicked.connect(lambda: inp.setText(QFileDialog.getExistingDirectory() if folder else QFileDialog.getOpenFileName()[0]))
        grid_layout.addWidget(lbl, row, 0)
        grid_layout.addWidget(inp, row, 1)
        grid_layout.addWidget(btn, row, 2)
        return inp

    def toggle_pipeline(self):
        if self.pipeline_thread and self.pipeline_thread.isRunning():
            # Stop the currently running process
            self.pipeline_thread.stop()
            self.start_button.setText("Start Analysis")
        else:
            # Reinitialize and restart the pipeline
            self.pipeline_thread = None  # Clear the old thread reference
            self.start_pipeline()

    def start_pipeline(self):
        self.log_display.clear()
        form_data = {
            "ref_folder": self.ref_folder.text(),
            "results_folder": self.results_folder.text(),
            "fastq1": self.fastq1.text(),
            "fastq2": self.fastq2.text(),
            "reference_genome": self.reference_genome.text(),
            "gtf_file": self.gtf_file.text(),
            "t_option": self.t_option.text(),
            "g_option": self.g_option.text(),
            "threads": int(self.threads.text()),
        }
        # Initialize a new pipeline thread
        self.pipeline_thread = PipelineThread(form_data)
        self.pipeline_thread.stopped = False  # Reset the stopped flag
        self.pipeline_thread.log_signal.connect(self.log_display.append)
        self.pipeline_thread.completed_signal.connect(self.reset_ui_after_completion)
        self.start_button.setText("Stop Analysis")  # Update button text
        self.pipeline_thread.start()

    def reset_ui_after_completion(self):
        self.start_button.setText("Start Analysis")
        self.pipeline_thread = None  # Clear the pipeline thread reference
        self.log_display.append("🎉 Analysis complete.")
        self.log_display.clear()
        self.ref_folder.clear()
        self.results_folder.clear()
        self.fastq1.clear()
        self.fastq2.clear()
        self.reference_genome.clear()
        self.gtf_file.clear()
        self.t_option.setText("exon")
        self.g_option.setText("Parent")
        self.threads.setText("40")



if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
