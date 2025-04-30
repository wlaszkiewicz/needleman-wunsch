from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QLabel, QLineEdit, QPushButton, QTabWidget,
                             QSpinBox, QFileDialog, QMessageBox, QRadioButton)
from core.algorithm import NeedlemanWunsch
from gui.widgets.results_tab import AlignmentResultsTab


class MainWindow(QMainWindow):
    """Main application window for the Needleman-Wunsch alignment tool.

    Provides a graphical interface for sequence input, parameter configuration,
    alignment execution, and result visualization.

    Attributes:
        nw (NeedlemanWunsch): The alignment algorithm instance
        seq1_input (QLineEdit): Input field for sequence 1
        seq2_input (QLineEdit): Input field for sequence 2
        match_score (QSpinBox): Widget for match score parameter
        mismatch_penalty (QSpinBox): Widget for mismatch penalty parameter
        gap_penalty (QSpinBox): Widget for gap penalty parameter
        results_tab (AlignmentResultsTab): Tab for displaying alignment results
        matrix_tab (MatrixTab): Tab for visualizing the score matrix
    """

    def __init__(self):
        """Initializes the main window with default settings and UI components."""
        super().__init__()

        # Center the window on screen
        screen_geometry = QApplication.desktop().screenGeometry()
        x = (screen_geometry.width() - 1700) // 2
        y = (screen_geometry.height() - 1500) // 2

        self.setWindowTitle("Needleman-Wunsch Alignment")
        self.setGeometry(x, y, 1800, 1600)
        self.nw = NeedlemanWunsch()
        self.init_ui()

    def init_ui(self):
        """Initializes all UI components and layouts."""
        main_widget = QWidget()
        main_layout = QVBoxLayout()

        # Input section
        input_group = QWidget()
        input_layout = QVBoxLayout()

        # Sequence type selection
        type_layout = QHBoxLayout()
        type_layout.addWidget(QLabel("Sequence Type:"))
        self.auto_detect_btn = QRadioButton("Auto-detect")
        self.auto_detect_btn.setChecked(True)
        self.dna_btn = QRadioButton("DNA")
        self.protein_btn = QRadioButton("Protein")
        type_layout.addWidget(self.auto_detect_btn)
        type_layout.addWidget(self.dna_btn)
        type_layout.addWidget(self.protein_btn)
        input_layout.addLayout(type_layout)

        # Sequence input fields
        seq1_layout = QHBoxLayout()
        seq1_layout.addWidget(QLabel("Sequence 1:"))
        self.seq1_input = QLineEdit()
        self.load_seq1_btn = QPushButton("Load FASTA")
        seq1_layout.addWidget(self.seq1_input)
        seq1_layout.addWidget(self.load_seq1_btn)
        input_layout.addLayout(seq1_layout)

        seq2_layout = QHBoxLayout()
        seq2_layout.addWidget(QLabel("Sequence 2:"))
        self.seq2_input = QLineEdit()
        self.load_seq2_btn = QPushButton("Load FASTA")
        seq2_layout.addWidget(self.seq2_input)
        seq2_layout.addWidget(self.load_seq2_btn)
        input_layout.addLayout(seq2_layout)

        # Scoring parameters
        params_layout = QHBoxLayout()
        params_layout.addWidget(QLabel("Match Score:"))
        self.match_score = QSpinBox()
        self.match_score.setRange(-10, 10)
        self.match_score.setValue(1)
        params_layout.addWidget(self.match_score)

        params_layout.addWidget(QLabel("Mismatch Penalty:"))
        self.mismatch_penalty = QSpinBox()
        self.mismatch_penalty.setRange(-10, 0)
        self.mismatch_penalty.setValue(0)
        params_layout.addWidget(self.mismatch_penalty)

        params_layout.addWidget(QLabel("Gap Penalty:"))
        self.gap_penalty = QSpinBox()
        self.gap_penalty.setRange(-10, 0)
        self.gap_penalty.setValue(-1)
        params_layout.addWidget(self.gap_penalty)

        input_layout.addLayout(params_layout)
        input_group.setLayout(input_layout)
        main_layout.addWidget(input_group)

        # Alignment button
        self.align_btn = QPushButton("Align Sequences")
        main_layout.addWidget(self.align_btn)

        # Results tabs
        self.tabs = QTabWidget()

        # Results display tab
        self.results_tab = AlignmentResultsTab()
        self.results_tab.save_btn.clicked.connect(self.save_results)
        self.tabs.addTab(self.results_tab, "Alignment Results")

        # Matrix visualization tab (imported locally to avoid circular imports)
        from gui.widgets.matrix_tab import MatrixTab
        self.matrix_tab = MatrixTab()
        self.tabs.addTab(self.matrix_tab, "Score Matrix")

        main_layout.addWidget(self.tabs)
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)

        # Connect signals to slots
        self.align_btn.clicked.connect(self.perform_alignment)
        self.load_seq1_btn.clicked.connect(lambda: self.load_fasta(1))
        self.load_seq2_btn.clicked.connect(lambda: self.load_fasta(2))

    def load_fasta(self, seq_num):
        """Loads a sequence from a FASTA file into the specified input field.

        Args:
            seq_num (int): Which sequence to load (1 or 2)

        Note:
            Only reads the first sequence from the file if multiple are present.
            Skips lines starting with '>' (header lines).
        """
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Open FASTA File",
            "",
            "FASTA Files (*.fasta *.fa *.txt)"
        )

        if filepath:
            try:
                with open(filepath, 'r') as f:
                    sequence = ''.join(
                        line.strip()
                        for line in f
                        if not line.startswith('>')
                    )
                    if seq_num == 1:
                        self.seq1_input.setText(sequence)
                    else:
                        self.seq2_input.setText(sequence)
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Error",
                    f"Failed to load FASTA file: {str(e)}"
                )

    def get_selected_sequence_type(self):
        """Gets the currently selected sequence type from radio buttons.

        Returns:
            str: 'auto', 'dna', or 'protein'
        """
        if self.auto_detect_btn.isChecked():
            return 'auto'
        if self.dna_btn.isChecked():
            return 'dna'
        return 'protein'

    def validate_input(self):
        """Validates the user input before alignment.

        Returns:
            bool: True if input is valid, False otherwise

        Raises:
            Shows QMessageBox warnings for invalid input
        """
        seq1 = self.seq1_input.text().strip().upper()
        seq2 = self.seq2_input.text().strip().upper()

        if not seq1 or not seq2:
            QMessageBox.warning(
                self,
                "Input Error",
                "Both sequences must be provided"
            )
            return False

        seq_type = self.get_selected_sequence_type()
        valid_chars = {
            'dna': {'A', 'T', 'C', 'G', 'U', 'N', '-'},
            'protein': {
                'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
                'X', '-', 'B', 'Z', 'J'
            },
            'auto': set()
        }.get(seq_type, set())

        if valid_chars:
            for char in seq1 + seq2:
                if char not in valid_chars:
                    QMessageBox.warning(
                        self,
                        "Input Error",
                        f"Invalid character '{char}'"
                    )
                    return False

        return True

    def perform_alignment(self):
        """Executes the sequence alignment and displays results.

        Handles the complete workflow:
        1. Validates input
        2. Configures algorithm parameters
        3. Runs alignment
        4. Displays results in both tabs
        """
        if not self.validate_input():
            return

        # Configure alignment parameters
        self.nw = NeedlemanWunsch(
            match_score=self.match_score.value(),
            mismatch_penalty=self.mismatch_penalty.value(),
            gap_penalty=self.gap_penalty.value()
        )

        try:
            # Set sequences and perform alignment
            self.nw.set_sequences(
                self.seq1_input.text().strip(),
                self.seq2_input.text().strip(),
                self.get_selected_sequence_type()
            )

            results = self.nw.align()
            self.display_results(results)

            # Visualize paths in matrix tab if available
            if self.nw.alignment_paths:
                self.matrix_tab.plot_all_paths(
                    self.nw,
                    self.nw.alignment_paths
                )

            self.tabs.setCurrentIndex(0)

        except ValueError as ve:
            QMessageBox.critical(self, "Input Error", str(ve))
        except Exception as e:
            QMessageBox.critical(
                self,
                "Alignment Error",
                f"An error occurred: {str(e)}"
            )

    def display_results(self, results):
        """Formats and displays alignment results in the results tab.

        Args:
            results (dict): Alignment results from NeedlemanWunsch.align()
        """
        output = [
            "=== Needleman-Wunsch Global Alignment Results ===",
            f"Sequence Type: {self.nw.seq_type.upper()}",
            f"Match Score: {self.nw.match_score}",
            f"Mismatch Penalty: {self.nw.mismatch_penalty}",
            f"Gap Penalty: {self.nw.gap_penalty}",
            "",
            f"Sequence 1 (len={len(self.nw.seq1)}): {self.nw.seq1[:50]}..."
            if len(self.nw.seq1) > 50
            else f"Sequence 1: {self.nw.seq1}",
            f"Sequence 2 (len={len(self.nw.seq2)}): {self.nw.seq2[:50]}..."
            if len(self.nw.seq2) > 50
            else f"Sequence 2: {self.nw.seq2}",
            "",
            f"Alignment Score: {results['score']}",
            f"Total Optimal Paths: {results['total_paths']}"
            if results['total_paths'] < 10
            else "Total Optimal Paths: >=10",
            f"Identity: {results['identity_percentage']:.2f}%",
            f"Gaps: {results['gap_percentage']:.2f}%",
            ""
        ]

        # Format alignment examples
        if results['examples']:
            output.append("=== Example Alignments ===")
            for i, (align1, align2, symbols) in enumerate(results['examples'][:10], 1):
                # Split long alignments into chunks for better readability
                chunk_size = 100
                chunks = [
                    (
                        align1[i:i + chunk_size],
                        align2[i:i + chunk_size],
                        symbols[i:i + chunk_size]
                    )
                    for i in range(0, len(align1), chunk_size)
                ]

                output.append(f"\nAlignment Example {i}:")
                for chunk_num, (chunk1, chunk2, chunk_symbols) in enumerate(chunks, 1):
                    if len(chunks) > 1:
                        output.append(f"\nChunk {chunk_num}:")
                    output.extend([
                        chunk1,
                        chunk_symbols,
                        chunk2,
                    ])

        self.results_tab.results_text.setPlainText('\n'.join(output))
        self.results_tab.save_btn.setEnabled(True)

    def save_results(self):
        """Saves the current alignment results to a text file."""
        if not hasattr(self.nw, 'alignment_paths') or not self.nw.alignment_paths:
            QMessageBox.warning(
                self,
                "Save Error",
                "No alignment results to save"
            )
            return

        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Save Results",
            "",
            "Text Files (*.txt)"
        )

        if filepath:
            try:
                with open(filepath, 'w') as f:
                    f.write(self.results_tab.results_text.toPlainText())
                QMessageBox.information(
                    self,
                    "Success",
                    "Results saved successfully"
                )
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Save Error",
                    f"Failed to save file: {str(e)}"
                )