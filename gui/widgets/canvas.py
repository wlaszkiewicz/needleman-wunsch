import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.lines import Line2D

class ScoreMatrixCanvas(FigureCanvas):
    def __init__(self, parent=None, width=8, height=8, dpi=100):
        self.fig, self.ax = plt.subplots(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)
        self.nw = None
        self.path = None
        self.all_paths = []

    def plot_matrix(self, nw, path=None):
        self.nw = nw
        self.path = path or []
        self.ax.clear()

        if not nw or not nw.score_matrix.any():
            return

        seq1 = nw.seq1
        seq2 = nw.seq2
        matrix = nw.score_matrix
        rows, cols = matrix.shape

        # Draw sequence characters
        for i, char in enumerate(seq1):
            self.ax.text(0, i + 2, char, ha='center', va='center',
                         fontsize=14, fontweight='bold')
        for j, char in enumerate(seq2):
            self.ax.text(j + 2, 0, char, ha='center', va='center',
                         fontsize=14, fontweight='bold')

        # Draw matrix cells
        for i in range(rows):
            for j in range(cols):
                if (i, j) in self.path:
                    rect = plt.Rectangle((j + 0.5, i + 0.5), 1, 1,
                                         facecolor='lightblue', alpha=0.6)
                    self.ax.add_patch(rect)

                self.ax.text(j + 1, i + 1, f"{matrix[i, j]}",
                             ha='center', va='center',
                             color='black', fontsize=14)

        # Draw path arrows
        if len(self.path) > 1:
            for k in range(len(self.path) - 1):
                i1, j1 = self.path[k]
                i2, j2 = self.path[k + 1]

                x1, y1 = j1 + 1, i1 + 1
                x2, y2 = j2 + 1, i2 + 1

                dx, dy = x2 - x1, y2 - y1
                arrow_length = 0.7

                if dx != 0 and dy != 0:  # Diagonal
                    self.ax.arrow(x1 + 0.5, y1 + 0.5, dx * arrow_length, dy * arrow_length,
                                  head_width=0.3, head_length=0.2,
                                  fc='#80577e', ec='#80577e')
                elif dy != 0:  # Vertical
                    self.ax.arrow(x1 + 0.5, y1 + 0.5, 0, dy * arrow_length,
                                  head_width=0.3, head_length=0.2,
                                  fc='#424b66', ec='#424b66')
                else:  # Horizontal
                    self.ax.arrow(x1 + 0.5, y1 + 0.5, dx * arrow_length, 0,
                                  head_width=0.3, head_length=0.2,
                                  fc='#424b66', ec='#424b66')

        # Configure axes and appearance
        self.ax.set_xlim(-0.5, cols + 0.5)
        self.ax.set_ylim(rows + 0.5, -0.5)  # Inverted y-axis

        for i in range(rows + 1):
            self.ax.axhline(i - 0.5, color='gray', linewidth=0.5)
        for j in range(cols + 1):
            self.ax.axvline(j - 0.5, color='gray', linewidth=0.5)

        self.ax.set_xticks([])
        self.ax.set_yticks([])

        # Add title and labels
        self.ax.set_title(f'Score Matrix ({nw.seq_type.upper()})',
                          pad=20, fontsize=16, fontweight='bold')
        self.ax.set_xlabel('Sequence 2', fontsize=14, labelpad=10)
        self.ax.set_ylabel('Sequence 1', fontsize=14, labelpad=10)

        # Add legend
        legend_elements = [
            Line2D([0], [0], marker='$\u2198$', color='#80577e', label='Match/Mismatch',
                   markersize=16, linestyle='None'),
            Line2D([0], [0], marker='$\u2193$', color='#424b66', label='Gap in Seq2',
                   markersize=16, linestyle='None'),
            Line2D([0], [0], marker='$\u2192$', color='#424b66', label='Gap in Seq1',
                   markersize=16, linestyle='None')
        ]
        self.ax.legend(handles=legend_elements, loc='upper right',
                       bbox_to_anchor=(1.3, 1), fontsize=16)



        self.fig.tight_layout(rect=[0, 0, 0.98, 1])
        self.draw()
