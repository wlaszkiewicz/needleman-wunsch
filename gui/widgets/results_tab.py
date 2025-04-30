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

        # Configure the text display area
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)  # Make non-editable
        self.results_text.setFontFamily("Consolas")  # Monospace font for alignment formatting

        # Set up scrolling container
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)  # Allow text area to expand
        self.scroll.setWidget(self.results_text)  # Nest the text display

        # Add export capability
        self.save_btn = QPushButton("Save Results to File")

        # Assemble the layout
        self.layout.addWidget(self.scroll)  # Scrollable text area takes most space
        self.layout.addWidget(self.save_btn)  # Save button at bottom
