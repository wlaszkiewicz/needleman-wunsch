from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QMessageBox, QFileDialog
from gui.widgets.canvas import ScoreMatrixCanvas


class MatrixTab(QWidget):
    """Tab widget for visualizing and navigating through Needleman-Wunsch score matrix paths.

    Provides interactive visualization of the alignment matrix with navigation controls
    to explore different optimal alignment paths.

    Attributes:
        matrix_canvas (ScoreMatrixCanvas): The visualization canvas
        current_path_index (int): Index of currently displayed path
        prev_btn (QPushButton): Button to show previous path
        next_btn (QPushButton): Button to show next path
        path_label (QLabel): Displays current path position
        save_plot_btn (QPushButton): Button to save visualization
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        self.current_path_index = 0

        self.matrix_canvas = ScoreMatrixCanvas(self)

        nav_container = QWidget()
        nav_container.setFixedHeight(60)

        nav_layout = QHBoxLayout(nav_container)
        nav_layout.setContentsMargins(10, 10, 10, 10)
        nav_layout.addStretch()

        self.prev_btn = QPushButton("◄ Previous Path")
        nav_layout.addWidget(self.prev_btn)
        nav_layout.addSpacing(20)
        self.path_label = QLabel("Path 0 of 0")
        nav_layout.addWidget(self.path_label)
        nav_layout.addSpacing(20)
        self.next_btn = QPushButton("Next Path ►")
        nav_layout.addWidget(self.next_btn)
        nav_layout.addStretch()

        self.prev_btn.clicked.connect(self.show_previous_path)
        self.next_btn.clicked.connect(self.show_next_path)

        self.save_plot_btn = QPushButton("Save Plot as Image")
        self.save_plot_btn.clicked.connect(self.save_plot)

        self.layout.addWidget(nav_container)
        self.layout.addWidget(self.matrix_canvas, 1)
        self.layout.addWidget(self.save_plot_btn)

        self.prev_btn.setEnabled(False)
        self.next_btn.setEnabled(False)

    def show_previous_path(self):
        """Displays the previous optimal path in the sequence.

        Decrements the current path index and updates the display if not at first path.
        """
        if self.current_path_index > 0:
            self.current_path_index -= 1
            self.update_path_display()

    def show_next_path(self):
        """Displays the next optimal path in the sequence.

        Increments the current path index and updates the display if not at last path.
        """
        if self.current_path_index < len(self.matrix_canvas.all_paths) - 1:
            self.current_path_index += 1
            self.update_path_display()

    def update_path_display(self):
        """Updates the visualization to show the current path and updates UI state.

        Handles:
        - Loading the correct path into the visualization
        - Updating the path counter label
        - Enabling/disabling navigation buttons as needed
        """
        if hasattr(self.matrix_canvas, 'all_paths') and self.matrix_canvas.all_paths:
            path = self.matrix_canvas.all_paths[self.current_path_index]
            self.matrix_canvas.plot_matrix(self.matrix_canvas.nw, path)
            self.path_label.setText(
                f"Path {self.current_path_index + 1} of {len(self.matrix_canvas.all_paths)}"
            )

            self.prev_btn.setEnabled(self.current_path_index > 0)
            self.next_btn.setEnabled(
                self.current_path_index < len(self.matrix_canvas.all_paths) - 1
            )

    def plot_all_paths(self, nw, paths):
        """Initializes the tab with a new set of alignment paths.

        Args:
            nw (NeedlemanWunsch): Configured alignment algorithm instance
            paths (list): List of optimal paths to visualize
        """
        self.matrix_canvas.nw = nw
        self.matrix_canvas.all_paths = paths
        self.current_path_index = 0

        if paths:

            self.prev_btn.setEnabled(False)
            self.next_btn.setEnabled(len(paths) > 1)
            self.path_label.setText(f"Path 1 of {len(paths)}")
            self.update_path_display()
        else:

            self.path_label.setText("No paths available")
            self.prev_btn.setEnabled(False)
            self.next_btn.setEnabled(False)

    def save_plot(self):
        """Saves the current matrix visualization to an image file.

        Supports multiple image formats through file dialog selection.
        Shows appropriate message boxes for success/error cases.
        """
        if not hasattr(self.matrix_canvas, 'nw') or not self.matrix_canvas.nw:
            QMessageBox.warning(
                self,
                "Save Error",
                "No matrix plot to save"
            )
            return

        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Save Plot",
            "",
            "PNG Images (*.png);;JPEG Images (*.jpg *.jpeg);;All Files (*)"
        )

        if filepath:
            try:
                self.matrix_canvas.fig.savefig(
                    filepath,
                    dpi=300,
                    bbox_inches='tight'
                )
                QMessageBox.information(
                    self,
                    "Success",
                    "Plot saved successfully"
                )
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Save Error",
                    f"Failed to save plot: {str(e)}"
                )