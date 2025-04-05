import sys
import os
import subprocess
import signal
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QWidget, QFileDialog, QTextEdit,
    QProgressBar, QMessageBox, QGroupBox, QHBoxLayout
)
from PySide6.QtCore import QThread, Signal, Qt
from PySide6.QtGui import QFont , QIcon ,QMouseEvent
from PySide6.QtWidgets import QGraphicsDropShadowEffect , QLabel



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
                        wsl_command,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=True,
                        text=True,
                        creationflags=subprocess.CREATE_NEW_PROCESS_GROUP # ✅ Set new process group for clean kill
                    )


                    while self.running_process and self.running_process.poll() is None:
                        if self.stopped:
                            self.running_process.send_signal(signal.CTRL_BREAK_EVENT)
                            self.log_signal.emit("❌ Process terminated by user.")
                            return

                        output = self.running_process.stdout.readline()
                        if output:
                            self.log_signal.emit(output.strip())

                    # After process ends
                    if self.running_process:
                        try:
                            for err_line in self.running_process.stderr:
                                self.log_signal.emit(f"⚠️ {err_line.strip()}")

                            self.running_process.wait()
                            if self.running_process.returncode != 0:
                                raise subprocess.CalledProcessError(self.running_process.returncode, command)

                            self.log_signal.emit(f"✅ Completed: {description}.")
                        except Exception as e:
                            self.log_signal.emit(f"⚠️ Error finishing process: {e}")
                except Exception as e:
                    self.log_signal.emit(f"❌ Error: {str(e)}")
                    if not self.stopped:
                        self.completed_signal.emit()


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
            

            elif self.mode == "lncrna":

                annotation_gtf = self.convert_to_wsl_path(self.form_data["annotation"])

                cpc2_local_path = os.path.join(os.path.dirname(__file__), "cpc2", "bin", "CPC2.py")
                if not os.path.exists(cpc2_local_path):
                    self.log_signal.emit("❌ CPC2.py not found in bundled CPC2 directory!")
                    return
                cpc2_script = self.convert_to_wsl_path(cpc2_local_path)



                aligned_bam = os.path.join(bam_dir, f"{fastq_prefix}_Aligned.sortedByCoord.out.bam")
                aligned_bam_wsl = self.convert_to_wsl_path(aligned_bam)

                run_command(
                    f"STAR --runThreadN {self.form_data['threads']} "
                    f"--genomeDir {self.convert_to_wsl_path(ref_folder)} "
                    f"--readFilesIn {fastq1_path} {fastq2_path} "
                    f"--readFilesCommand zcat "
                    f"--outSAMtype BAM SortedByCoordinate "
                    f"--quantMode TranscriptomeSAM "
                    f"--outFileNamePrefix {os.path.join(bam_dir, fastq_prefix + '_')}",
                    "Aligning with STAR"
                )
                if self.stopped: return

                sample_gtf = os.path.join(output_dir, f"{fastq_prefix}.gtf")
                sample_gtf_wsl = self.convert_to_wsl_path(sample_gtf)
                abundance_file = os.path.join(output_dir, f"{fastq_prefix}_abundance.tsv")

                run_command(
                    f"stringtie {aligned_bam_wsl} -G {annotation_gtf} -o {sample_gtf_wsl} "
                    f"-A {self.convert_to_wsl_path(abundance_file)} -l {fastq_prefix}",
                    "Assembling transcripts with StringTie"
                )
                if self.stopped: return

                run_command(
                    f"gffcompare -r {annotation_gtf} -o {os.path.join(output_dir, 'compare')} {sample_gtf_wsl}",
                    "Comparing transcripts with GFFCompare"
                )
                if self.stopped: return

                novel_gtf = os.path.join(output_dir, "novel_candidates.gtf")
                novel_gtf_wsl = self.convert_to_wsl_path(novel_gtf)

                run_command(
                    f"awk '$3==\"transcript\"' {os.path.join(output_dir, 'compare.annotated.gtf')} "
                    f"| grep 'class_code \"[uix]\"' > {novel_gtf_wsl}",
                    "Extracting novel transcript candidates"
                )
                if self.stopped: return

                transcripts_fa = os.path.join(output_dir, "transcripts.fa")
                transcripts_fa_wsl = self.convert_to_wsl_path(transcripts_fa)
                run_command(
                    f"gffread {novel_gtf_wsl} -g {ref_genome_path} -w {transcripts_fa_wsl}",
                    "Extracting transcript sequences"
                )
                if self.stopped: return

                cpc2_output = os.path.join(output_dir, "cpc2_results.txt")
                run_command(
                    f"python3 {cpc2_script} -i {transcripts_fa_wsl} -o {self.convert_to_wsl_path(cpc2_output)}",
                    "Running CPC2 coding potential prediction"
                )
                if self.stopped: return

                final_ids = os.path.join(output_dir, "noncoding_ids.txt")
                final_gtf = os.path.join(output_dir, "final_ncrna.gtf")
                final_gtf_wsl = self.convert_to_wsl_path(final_gtf)

                run_command(
                    f"awk 'NR>1 && $8 < 0.5 {{ print $1 }}' {self.convert_to_wsl_path(cpc2_output)} "
                    f"> {self.convert_to_wsl_path(final_ids)}",
                    "Filtering non-coding transcripts"
                )
                if self.stopped: return

                run_command(
                    f"grep -Ff {self.convert_to_wsl_path(final_ids)} {novel_gtf_wsl} > {final_gtf_wsl}",
                    "Extracting final novel lncRNAs"
                )
                if self.stopped: return

                run_command(
                    f"featureCounts -T {self.form_data['threads']} -p -a {final_gtf_wsl} -o {self.convert_to_wsl_path(output_file)} {aligned_bam_wsl}",
                    "Quantifying novel lncRNA expression"
                )

            else:
                ref_genome_path = self.convert_to_wsl_path(self.form_data["reference_genome"])
                raw_bcf = os.path.join(output_dir, f"{fastq_prefix}_raw.bcf")
                raw_vcf = os.path.join(output_dir, f"{fastq_prefix}_variants.vcf")

                raw_bcf_wsl = self.convert_to_wsl_path(raw_bcf)
                raw_vcf_wsl = self.convert_to_wsl_path(raw_vcf)
                output_file_wsl = self.convert_to_wsl_path(output_file)

                run_command(f"bcftools mpileup -f {ref_genome_path} {sorted_bam_path} -O b -o {raw_bcf_wsl}", "Generating BCF with mpileup")
                if self.stopped: return

                run_command(f"bcftools call --ploidy {self.form_data['ploidy']} -m -v {raw_bcf_wsl} -o {raw_vcf_wsl}", "Calling variants with BCFtools")
                if self.stopped: return

                run_command(f"vcfutils.pl varFilter {raw_vcf_wsl} > {output_file_wsl}", "Filtering variants")


            self.log_signal.emit("🎉 Analysis complete.")
            self.completed_signal.emit()

        except Exception as e:
            self.log_signal.emit(f"❌ Error: {str(e)}")
            if not self.stopped:
                self.completed_signal.emit()
                
        

    def stop(self):
        self.stopped = True
        if self.running_process:
            try:
                self.running_process.send_signal(signal.CTRL_BREAK_EVENT)
                self.log_signal.emit("❌ Process forcefully terminated.")
            except Exception as e:
                self.log_signal.emit(f"⚠️ Failed to terminate process cleanly: {e}")
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
        layout.setSpacing(20)
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
        self.lncrna_btn = QPushButton("🧠 lncRNA Identification")
        self.lncrna_btn.setFixedSize(300, 50)
        self.diff_btn.setFixedSize(300, 50)
        self.var_btn.setFixedSize(300, 50)
        layout.addWidget(title)
        layout.addWidget(subtitle)
        layout.addSpacing(20)
        layout.addWidget(self.diff_btn, alignment=Qt.AlignCenter)
        layout.addWidget(self.var_btn, alignment=Qt.AlignCenter)
        layout.addWidget(self.lncrna_btn, alignment=Qt.AlignCenter)  
        layout.addStretch()

        self.lncrna_btn.clicked.connect(lambda: self.load_main_window("lncrna"))
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
        self.mode = mode
        self.setWindowTitle(f"{'Variant Calling' if mode == 'variant' else 'Differential Expression'} Pipeline")
        self.setGeometry(100, 100, 900, 700)
        self.pipeline_thread = None

        self.is_dark_mode = False
        self.init_ui()
        self.apply_light_mode()

    def apply_light_mode(self):
        self.setStyleSheet("""
            QWidget { background-color: white; }
            QLabel { font-family: 'Segoe UI'; font-size: 12px; color: #222; }
            QLabel#TitleLabel { font-size: 20px; font-weight: bold; color: #007acc; }
            QLineEdit {
                padding: 4px 8px; border: 1px solid #ccc;
                border-radius: 4px; background-color: #fff;
                font-size: 12px; color: #222;
            }
            QTextEdit {
                background-color: #ffffff; border: 1px solid #ccc;
                border-radius: 5px; padding: 6px; font-size: 12px;
                color: #333; font-family: Consolas, monospace;
            }
            QPushButton {
                background-color: #007acc; color: white; border-radius: 6px;
                padding: 8px 16px; font-size: 13px; font-weight: bold;
                transition: all 0.2s ease-in-out;
            }
            QPushButton:hover {
                background-color: #3498db; transform: translateY(-2px);
            }
        """)

    def apply_dark_mode(self):
        self.setStyleSheet("""
            QWidget { background-color: #121212; }
            QLabel { font-family: 'Segoe UI'; font-size: 12px; color: #ddd; }
            QLabel#TitleLabel { font-size: 20px; font-weight: bold; color: #4fc3f7; }
            QLineEdit {
                padding: 4px 8px; border: 1px solid #444;
                border-radius: 4px; background-color: #1e1e1e;
                font-size: 12px; color: #eee;
            }
            QTextEdit {
                background-color: #1e1e1e; border: 1px solid #444;
                border-radius: 5px; padding: 6px; font-size: 12px;
                color: #ccffcc; font-family: Consolas, monospace;
            }
            QPushButton {
                background-color: #3498db; color: white; border-radius: 6px;
                padding: 8px 16px; font-size: 13px; font-weight: bold;
                transition: all 0.2s ease-in-out;
            }
            QPushButton:hover {
                background-color: #5dade2; transform: translateY(-2px);
            }
        """)

    def init_ui(self):
        main_layout = QVBoxLayout()

        self.theme_toggle = QPushButton("🌓")
        self.theme_toggle.setFixedSize(30, 30)
        self.theme_toggle.setToolTip("Toggle Light/Dark Mode")
        self.theme_toggle.setStyleSheet("border: none; background: none; font-size: 16px;")
        self.theme_toggle.clicked.connect(self.toggle_theme)

        # Add to top-right
        theme_layout = QHBoxLayout()
        theme_layout.addStretch()
        theme_layout.addWidget(self.theme_toggle)
        main_layout.addLayout(theme_layout)

        title_text = {
            "variant": "🧪 Variant Calling",
            "diff": "🧬 Differential Expression",
            "lncrna": "🧠 lncRNA Identification & Expression"
        }.get(self.mode, "NGS Pipeline")
        title = QLabel(f"{title_text} Pipeline")


        title.setObjectName("TitleLabel")
        title.setAlignment(Qt.AlignCenter)

        info = QLabel("Configure and run your NGS pipeline below. Make sure all fields are correctly filled!")
        info.setWordWrap(True)

        main_layout.addWidget(title)
        main_layout.addWidget(info)
        main_layout.addSpacing(10)

        # Grouping form layout
        form_group = QGroupBox("📁 Input Configuration")
        form_layout = QGridLayout()

        self.ref_folder = self.create_input_with_browse(form_layout, "Reference Folder:", 0, True, "Select the folder containing reference indices")
        self.results_folder = self.create_input_with_browse(form_layout, "Results Folder:", 1, True, "Folder where output will be saved")
        self.fastq1 = self.create_input_with_browse(form_layout, "FASTQ File 1:", 2, False, "Select the first FASTQ file")
        self.fastq2 = self.create_input_with_browse(form_layout, "FASTQ File 2:", 3, False, "Select the second FASTQ file")
        self.reference_genome = self.create_input_with_browse(form_layout, "Reference Genome:", 4, False, "FASTA file of the reference genome")

        if self.mode == "diff":
            self.gtf_file = self.create_input_with_browse(form_layout, "GTF/GFF File:", 5, False)
            form_layout.addWidget(QLabel("-t Option:"), 6, 0)
            self.t_option = QLineEdit("exon")
            form_layout.addWidget(self.t_option, 6, 1)

            form_layout.addWidget(QLabel("-g Option:"), 7, 0)
            self.g_option = QLineEdit("Parent")
            form_layout.addWidget(self.g_option, 7, 1)

        elif self.mode == "lncrna":
            self.annotation_file = self.create_input_with_browse(form_layout, "Annotation GTF:", 5, False)
            
            self.gtf_file = None  # for compatibility

        else:
            form_layout.addWidget(QLabel("Ploidy:"), 5, 0)
            self.ploidy = QLineEdit("1")
            form_layout.addWidget(self.ploidy, 5, 1)

        self.threads = QLineEdit("40")
        form_layout.addWidget(QLabel("Threads:"), 8, 0)
        form_layout.addWidget(self.threads, 8, 1)

        form_group.setLayout(form_layout)
        main_layout.addWidget(form_group)

        # Log terminal
        self.log_display = QTextEdit()
        self.log_display.setReadOnly(True)
        self.log_display.setMinimumHeight(220)
        main_layout.addWidget(QLabel("🧾 Terminal Output:"))
        main_layout.addWidget(self.log_display)

        # Progress bar
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.hide()
        main_layout.addWidget(self.progress)

        # Buttons
        self.start_button = QPushButton("▶️ Start Analysis")
        self.back_button = QPushButton("🔙 Back to Home")

        self.start_button.clicked.connect(self.toggle_pipeline)
        self.back_button.clicked.connect(self.go_back_home)

        btn_row = QHBoxLayout()
        btn_row.addWidget(self.start_button)
        btn_row.addWidget(self.back_button)
        main_layout.addLayout(btn_row)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

    def toggle_theme(self):
        self.is_dark_mode = not self.is_dark_mode
        if self.is_dark_mode:
            self.apply_dark_mode()
        else:
            self.apply_light_mode()

    def create_input_with_browse(self, layout, label, row, folder, tooltip=""):
        lbl = QLabel(label)
        inp = QLineEdit()
        inp.setToolTip(tooltip)

        icon = QLabel("📁")
        icon.setFixedSize(28, 28)
        icon.setAlignment(Qt.AlignCenter)
        icon.setCursor(Qt.PointingHandCursor)
        icon.setToolTip("Browse")

        icon.setStyleSheet("""
            QLabel {
                font-size: 18px;
            }
            QLabel:hover {
                transform: scale(1.3);
            }
        """)

        # Make QLabel behave like a button
        def browse():
            selected = QFileDialog.getExistingDirectory() if folder else QFileDialog.getOpenFileName()[0]
            if selected:
                inp.setText(selected)

        icon.mousePressEvent = lambda event: browse() if isinstance(event, QMouseEvent) else None

        layout.addWidget(lbl, row, 0)
        layout.addWidget(inp, row, 1)
        layout.addWidget(icon, row, 2)

        return inp


    def toggle_pipeline(self):
        if self.pipeline_thread and self.pipeline_thread.isRunning():
            self.pipeline_thread.stop()
            self.pipeline_thread.wait()
            self.pipeline_thread = None
            self.start_button.setText("▶️ Start Analysis")
            self.progress.hide()  # ✅ hide progress bar when stopped
        else:
            self.start_pipeline()



    def start_pipeline(self):
        self.log_display.clear()
        self.progress.show()

        form_data = {
            "ref_folder": self.ref_folder.text(),
            "results_folder": self.results_folder.text(),
            "fastq1": self.fastq1.text(),
            "fastq2": self.fastq2.text(),
            "reference_genome": self.reference_genome.text(),
            "threads": int(self.threads.text()),
        }

        if self.mode == "lncrna":
            form_data.update({
                "annotation": self.annotation_file.text(),
                
            })
        elif self.mode == "diff":
            form_data.update({
                "gtf_file": self.gtf_file.text(),
                "t_option": self.t_option.text(),
                "g_option": self.g_option.text(),
            })
        else:
            form_data["ploidy"] = self.ploidy.text()

        self.pipeline_thread = PipelineThread(form_data, self.mode)
        self.pipeline_thread.stopped = False  # ✅ Reset stop flag here!
        self.pipeline_thread.log_signal.connect(self.log_display.append)
        self.pipeline_thread.completed_signal.connect(self.reset_ui_after_completion)
        self.start_button.setText("⏹ Stop Analysis")
        self.pipeline_thread.start()

    def reset_ui_after_completion(self):
        self.start_button.setText("▶️ Start Analysis")
        self.progress.hide()  # ✅ stop the progress bar
        self.pipeline_thread = None
        QMessageBox.information(self, "🎉 Done!", "Analysis completed successfully!")


    def go_back_home(self):
        self.home = ModeSelector()
        self.home.show()
        self.close()



if __name__ == "__main__":
    app = QApplication(sys.argv)
    selector = ModeSelector()
    selector.show()
    sys.exit(app.exec())

