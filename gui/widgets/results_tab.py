from PyQt5.QtWidgets import QWidget, QVBoxLayout, QScrollArea, QPushButton, QTextEdit


class AlignmentResultsTab(QWidget):
    """Tab widget for displaying Needleman-Wunsch alignment results in a scrollable view.

    Provides a read-only text display for alignment results with:
    - Fixed-width font for proper alignment formatting
    - Scrollable content area for long alignments
    - Export functionality to save results to file

    Attributes:
        results_text (QTextEdit): The main text display widget
        scroll (QScrollArea): Scroll container for the text display
        save_btn (QPushButton): Button to export results to file
    """

    def __init__(self, parent=None):
        """Initializes the results tab with text display and save functionality.

        Args:
            parent (QWidget, optional): Parent widget. Defaults to None.
        """
        super().__init__(parent)
        self.layout = QVBoxLayout(self)

        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setFontFamily("Consolas")


        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.results_text)

        self.save_btn = QPushButton("Save Results to File")

        self.layout.addWidget(self.scroll)
        self.layout.addWidget(self.save_btn)
