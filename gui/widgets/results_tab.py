from PyQt5.QtWidgets import QWidget, QVBoxLayout, QScrollArea, QPushButton, QTextEdit


class AlignmentResultsTab(QWidget):
    def __init__(self, parent=None):
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

