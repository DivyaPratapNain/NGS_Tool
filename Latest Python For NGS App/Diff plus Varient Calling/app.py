
import sys
import os
import subprocess
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QWidget, QFileDialog, QTextEdit, QSizePolicy , QHBoxLayout
)
from PySide6.QtCore import QThread, Signal, Qt
from PySide6.QtGui import QFont


class PipelineThread(QThread):
    log_signal = Signal(str)
    completed_signal = Signal()
    running_process = None

    def __init__(self, form_data, mode):
        super().__init__()

        

        self.form_data = form_data
        self.mode = mode
        self.stopped = False

    def convert_to_wsl_path(self, path):
        try:
            result = subprocess.run(
                ["wsl", "wslpath", "-a", path],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
            )
            return result.stdout.strip()
        except subprocess.CalledProcessError as e:
            self.log_signal.emit(f"⚠️ Error converting path to WSL: {e.stderr.strip()}")
            return path

    def run(self):
        try:
            ref_folder = self.form_data["ref_folder"]
            if not ref_folder.endswith("/"):
                ref_folder += "/"

            results_folder = self.form_data["results_folder"]
            if not results_folder.endswith("/"):
                results_folder += "/"

            mode_folder = "Diff_Exp_Results/" if self.mode == "diff" else "Variant_Calling_Results/"
            sam_dir = os.path.join(results_folder, mode_folder, "sam/")
            bam_dir = os.path.join(results_folder, mode_folder, "bam/")
            output_dir = os.path.join(results_folder, mode_folder, "counts/" if self.mode == "diff" else "vcf/")
            os.makedirs(sam_dir, exist_ok=True)
            os.makedirs(bam_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)

            fastq_prefix = os.path.basename(self.form_data["fastq1"]).split("_")[0]
            reference_index = os.path.join(
                ref_folder, os.path.splitext(os.path.basename(self.form_data["reference_genome"]))[0]
            )
            sam_file = os.path.join(sam_dir, f"{fastq_prefix}_aligned.sam")
            sam_sorted_file = os.path.join(sam_dir, f"{fastq_prefix}_aligned_sorted.sam")
            sorted_bam = os.path.join(bam_dir, f"{fastq_prefix}_aligned.sorted.bam")
            output_file = os.path.join(output_dir, f"{fastq_prefix}_gene_expression.txt" if self.mode == "diff"
                                       else f"{fastq_prefix}_final_variants.vcf")

            def run_command(command, description):
                try:
                    self.log_signal.emit(f"🔄 Running: {description}...")
                    wsl_command = f'wsl bash -c "{command}"'
                    self.running_process = subprocess.Popen(
                        wsl_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                        shell=True, text=True, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP
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

            # Convert paths
            ref_genome_path = self.convert_to_wsl_path(self.form_data["reference_genome"])
            ref_index_path = self.convert_to_wsl_path(reference_index)
            fastq1_path = self.convert_to_wsl_path(self.form_data["fastq1"])
            fastq2_path = self.convert_to_wsl_path(self.form_data["fastq2"])
            sam_file_path = self.convert_to_wsl_path(sam_file)
            sorted_sam_path = self.convert_to_wsl_path(sam_sorted_file)
            sorted_bam_path = self.convert_to_wsl_path(sorted_bam)

            run_command(f"hisat2-build -p {self.form_data['threads']} {ref_genome_path} {ref_index_path}", "Building HISAT2 index")
            if self.stopped: return

            run_command(f"hisat2 -p {self.form_data['threads']} -x {ref_index_path} -1 {fastq1_path} -2 {fastq2_path} -S {sam_file_path}", "Aligning reads")
            if self.stopped: return

            run_command(f"samtools sort -@ {self.form_data['threads']} {sam_file_path} -o {sorted_sam_path}", "Sorting SAM file")
            if self.stopped: return

            run_command(f"samtools view -@ {self.form_data['threads']} -S -b {sorted_sam_path} > {sorted_bam_path}", "Converting SAM to BAM")
            if self.stopped: return

            run_command(f"samtools index -@ {self.form_data['threads']} {sorted_bam_path}", "Indexing BAM file")
            if self.stopped: return

            if self.mode == "diff":
                gtf_file_path = self.convert_to_wsl_path(self.form_data["gtf_file"])
                output_path = self.convert_to_wsl_path(output_file)
                run_command(f"featureCounts -T {self.form_data['threads']} -p -a {gtf_file_path} -o {output_path} -t {self.form_data['t_option']} -g {self.form_data['g_option']} {sorted_bam_path}", "Quantifying read counts")
            else:
                ref_genome_path = self.convert_to_wsl_path(self.form_data["reference_genome"])  # ✅ use actual FASTA
                raw_bcf = os.path.join(output_dir, f"{fastq_prefix}_raw.bcf")
                raw_vcf = os.path.join(output_dir, f"{fastq_prefix}_variants.vcf")

                # ✅ Convert all to WSL paths
                raw_bcf_wsl = self.convert_to_wsl_path(raw_bcf)
                raw_vcf_wsl = self.convert_to_wsl_path(raw_vcf)
                output_file_wsl = self.convert_to_wsl_path(output_file)

                run_command(f"bcftools mpileup -f {ref_genome_path} {sorted_bam_path} -O b -o {raw_bcf_wsl}", "Generating BCF with mpileup")
                run_command(f"bcftools call --ploidy {self.form_data['ploidy']} -m -v {raw_bcf_wsl} -o {raw_vcf_wsl}", "Calling variants with BCFtools")
                run_command(f"vcfutils.pl varFilter {raw_vcf_wsl} > {output_file_wsl}", "Filtering variants")

            self.log_signal.emit("🎉 Analysis complete.")
            self.completed_signal.emit()
        except Exception as e:
            self.log_signal.emit(f"❌ Error: {str(e)}")
            self.completed_signal.emit()

    def stop(self):
        self.stopped = True
        if self.running_process:
            self.running_process.terminate()
            self.log_signal.emit("❌ Process forcefully terminated.")
        self.running_process = None


class ModeSelector(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setStyleSheet("""
            QLabel {
                font-family: 'Segoe UI';
                font-size: 14px;
                background: transparent;
                border: none;
            }

            QLineEdit {
                padding: 4px;
                border: 1px solid #ccc;
                border-radius: 4px;
            }

            QTextEdit {
                background-color: #f9f9f9;
                border: 1px solid #ccc;
                border-radius: 4px;
                padding: 6px;
            }

            QWidget {
                background-color: qlineargradient(
                    spread:pad, x1:0, y1:0, x2:1, y2:1,
                    stop:0 #f0f2f5, stop:1 #dbe9f4
                );
            }

            QPushButton {
                background-color: #4fa3d1;
                color: white;
                border-radius: 8px;
                padding: 12px 24px;
                font-size: 16px;
                font-weight: bold;
                transition: transform 0.2s ease-in-out;  /* smoother animation */
            }

            QPushButton:hover {
                background-color: #6fc3f7;
                transform: translateY(-5px);  /* smooth move up */
            }
        """)

        self.setWindowTitle("NGS MADE EASY")
        self.setGeometry(200, 200, 1000, 600)

        layout = QVBoxLayout()
        layout.setContentsMargins(50, 100, 50, 100)
        layout.setSpacing(30)
        title = QLabel("NGS Analysis Toolkit")
        title.setFont(QFont("Arial", 24, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)
        subtitle = QLabel("Choose the type of analysis you want to perform")
        subtitle.setFont(QFont("Arial", 16))
        subtitle.setAlignment(Qt.AlignCenter)
        self.diff_btn = QPushButton("🧬 Differential Expression")
        self.var_btn = QPushButton("🧪 Variant Calling")
        self.diff_btn.setFixedSize(300, 50)
        self.var_btn.setFixedSize(300, 50)
        layout.addWidget(title)
        layout.addWidget(subtitle)
        layout.addSpacing(20)
        layout.addWidget(self.diff_btn, alignment=Qt.AlignCenter)
        layout.addWidget(self.var_btn, alignment=Qt.AlignCenter)
        layout.addStretch()


        self.diff_btn.clicked.connect(lambda: self.load_main_window("diff"))
        self.var_btn.clicked.connect(lambda: self.load_main_window("variant"))

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def load_main_window(self, mode):
        self.main = MainWindow(mode)
        self.main.show()
        self.close()


class MainWindow(QMainWindow):
    def __init__(self, mode):
        super().__init__()

        self.setStyleSheet("""
            QWidget {
                background-color: white;
            }

            QLabel {
                font-family: 'Segoe UI';
                font-size: 12px;
                color: #222;
            }

            QLabel#TitleLabel {
                font-size: 20px;
                font-weight: bold;
                color: #007acc;
            }

            QLabel#SectionHeader {
                font-size: 14px;
                font-weight: bold;
                color: #333;
                margin-top: 15px;
            }

            QLineEdit {
                padding: 4px 8px;
                border: 1px solid #ccc;
                border-radius: 4px;
                background-color: #fff;
                font-size: 12px;
                color: #222;
            }

            QTextEdit {
                background-color: #ffffff;
                border: 1px solid #ccc;
                border-radius: 5px;
                padding: 6px;
                font-size: 12px;
                color: #333;
                font-family: Consolas, monospace;
            }

            QPushButton {
                background-color: #007acc;
                color: white;
                border-radius: 6px;
                padding: 8px 16px;
                font-size: 13px;
                font-weight: bold;
                transition: all 0.2s ease-in-out;
            }

            QPushButton:hover {
                background-color: #3498db;
                transform: translateY(-2px);
            }
        """)

        self.mode = mode
        self.setWindowTitle(f"{'Variant Calling' if mode == 'variant' else 'Differential Expression'} Pipeline")
        self.setGeometry(100, 100, 900, 700)
        self.pipeline_thread = None
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout()

        title = QLabel(f"{'Variant Calling' if self.mode == 'variant' else 'Differential Expression'} Pipeline")
        title.setObjectName("TitleLabel")
        title.setAlignment(Qt.AlignCenter)

        info = QLabel()
        info.setWordWrap(True)
        info.setText(
            "This section allows you to configure and run the pipeline. "
            "Please ensure all fields are correctly filled before starting analysis."
        )
        info.setStyleSheet("color: #444; font-size: 13px;")

        main_layout.addWidget(title)
        main_layout.addSpacing(6)
        main_layout.addWidget(info)
        main_layout.addSpacing(10)


        grid_layout = QGridLayout()

        self.ref_folder = self.create_input_with_browse(grid_layout, "Reference Folder:", 0, True)
        self.results_folder = self.create_input_with_browse(grid_layout, "Results Folder:", 1, True)
        self.fastq1 = self.create_input_with_browse(grid_layout, "FASTQ File 1:", 2, False)
        self.fastq2 = self.create_input_with_browse(grid_layout, "FASTQ File 2:", 3, False)
        self.reference_genome = self.create_input_with_browse(grid_layout, "Reference Genome:", 4, False)

        if self.mode == "diff":
            self.gtf_file = self.create_input_with_browse(grid_layout, "GTF/GFF File:", 5, False)
            grid_layout.addWidget(QLabel("-t Option:"), 6, 0)
            self.t_option = QLineEdit("exon")
            grid_layout.addWidget(self.t_option, 6, 1)

            grid_layout.addWidget(QLabel("-g Option:"), 7, 0)
            self.g_option = QLineEdit("Parent")
            grid_layout.addWidget(self.g_option, 7, 1)
        else:
            grid_layout.addWidget(QLabel("Ploidy:"), 5, 0)
            self.ploidy = QLineEdit("1")
            grid_layout.addWidget(self.ploidy, 5, 1)

        self.threads = QLineEdit("40")
        grid_layout.addWidget(QLabel("Threads:"), 8, 0)
        grid_layout.addWidget(self.threads, 8, 1)

        main_layout.addLayout(grid_layout)

        self.log_display = QTextEdit()
      

        self.log_display.setReadOnly(True)
        main_layout.addWidget(self.log_display)

        self.start_button = QPushButton("Start Analysis")
        self.start_button.clicked.connect(self.toggle_pipeline)
        main_layout.addWidget(self.start_button)

        self.back_button = QPushButton("Back to Home")
        self.back_button.clicked.connect(self.go_back_home)
        main_layout.addWidget(self.back_button)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

    def go_back_home(self):
        self.home = ModeSelector()
        self.home.show()
        self.close()

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
            self.pipeline_thread.stop()
            self.start_button.setText("Start Analysis")
        else:
            self.start_pipeline()

    def start_pipeline(self):
        self.log_display.clear()
        form_data = {
            "ref_folder": self.ref_folder.text(),
            "results_folder": self.results_folder.text(),
            "fastq1": self.fastq1.text(),
            "fastq2": self.fastq2.text(),
            "reference_genome": self.reference_genome.text(),
            "threads": int(self.threads.text()),
        }
        if self.mode == "diff":
            form_data.update({
                "gtf_file": self.gtf_file.text(),
                "t_option": self.t_option.text(),
                "g_option": self.g_option.text(),
            })
        else:
            form_data["ploidy"] = self.ploidy.text()

        self.pipeline_thread = PipelineThread(form_data, self.mode)
        self.pipeline_thread.log_signal.connect(self.log_display.append)
        self.pipeline_thread.completed_signal.connect(self.reset_ui_after_completion)
        self.start_button.setText("Stop Analysis")
        self.pipeline_thread.start()

    def reset_ui_after_completion(self):
        self.start_button.setText("Start Analysis")
        self.pipeline_thread = None


if __name__ == "__main__":
    app = QApplication(sys.argv)
    selector = ModeSelector()
    selector.show()
    sys.exit(app.exec())


