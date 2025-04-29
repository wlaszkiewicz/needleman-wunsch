from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QMessageBox, QFileDialog

from gui.widgets.canvas import ScoreMatrixCanvas


class MatrixTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        self.current_path_index = 0

        # Create canvas
        self.matrix_canvas = ScoreMatrixCanvas(self)

        # Navigation controls
        nav_layout = QHBoxLayout()
        self.prev_btn = QPushButton("◄ Previous Path")
        self.next_btn = QPushButton("Next Path ►")
        self.path_label = QLabel("Path 0 of 0")

        self.prev_btn.clicked.connect(self.show_previous_path)
        self.next_btn.clicked.connect(self.show_next_path)

        nav_layout.addWidget(self.prev_btn)
        nav_layout.addWidget(self.path_label)
        nav_layout.addWidget(self.next_btn)

        # Save button
        self.save_plot_btn = QPushButton("Save Plot as Image")
        self.save_plot_btn.clicked.connect(self.save_plot)

        # Add widgets to layout
        self.layout.addWidget(self.matrix_canvas)
        self.layout.addLayout(nav_layout)
        self.layout.addWidget(self.save_plot_btn)

        # Disable buttons initially
        self.prev_btn.setEnabled(False)
        self.next_btn.setEnabled(False)

    def show_previous_path(self):
        if self.current_path_index > 0:
            self.current_path_index -= 1
            self.update_path_display()

    def show_next_path(self):
        if self.current_path_index < len(self.matrix_canvas.all_paths) - 1:
            self.current_path_index += 1
            self.update_path_display()

    def update_path_display(self):
        if hasattr(self.matrix_canvas, 'all_paths') and self.matrix_canvas.all_paths:
            path = self.matrix_canvas.all_paths[self.current_path_index]
            self.matrix_canvas.plot_matrix(self.matrix_canvas.nw, path)
            self.path_label.setText(f"Path {self.current_path_index + 1} of {len(self.matrix_canvas.all_paths)}")

            # Update button states
            self.prev_btn.setEnabled(self.current_path_index > 0)
            self.next_btn.setEnabled(self.current_path_index < len(self.matrix_canvas.all_paths) - 1)

    def plot_all_paths(self, nw, paths):
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
        if not hasattr(self.matrix_canvas, 'nw') or not self.matrix_canvas.nw:
            QMessageBox.warning(self, "Save Error", "No matrix plot to save")
            return

        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Save Plot",
            "",
            "PNG Images (*.png);;JPEG Images (*.jpg *.jpeg);;All Files (*)"
        )

        if filepath:
            try:
                self.matrix_canvas.fig.savefig(filepath, dpi=300, bbox_inches='tight')
                QMessageBox.information(self, "Success", "Plot saved successfully")
            except Exception as e:
                QMessageBox.critical(self, "Save Error", f"Failed to save plot: {str(e)}")

