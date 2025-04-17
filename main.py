import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText
import numpy as np
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from tkinter import font as tkfont


class NeedlemanWunschApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Needleman-Wunsch Sequence Alignment Tool")
        self.root.geometry("1200x600")
        self.root.minsize(1000, 500)

        # Configure styles
        self.style = ttk.Style()
        self.style.configure('TFrame', background='#f0f0f0')
        self.style.configure('TLabel', background='#f0f0f0', font=('Helvetica', 10))
        self.style.configure('TButton', font=('Helvetica', 10), padding=5)
        self.style.configure('Header.TLabel', font=('Helvetica', 12, 'bold'))
        self.style.configure('Highlight.TFrame', background='#e1f5fe')
        self.style.configure('Sequence.TFrame', background='#e8f5e9')
        self.style.configure('Matrix.TFrame', background='#fff3e0')

        # Fonts
        self.bold_font = tkfont.Font(family='Helvetica', size=10, weight='bold')
        self.mono_font = tkfont.Font(family='Courier', size=10)

        # Default parameters
        self.match_score = 1
        self.mismatch_penalty = 0
        self.gap_penalty = -1

        # Initialize sequences
        self.seq1 = ""
        self.seq2 = ""

        # Create UI
        self.create_widgets()

    def create_widgets(self):
        # Main container
        main_container = ttk.Frame(self.root)
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Header
        header_frame = ttk.Frame(main_container, style='Highlight.TFrame')
        header_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(header_frame, text="Needleman-Wunsch Sequence Alignment",
                  style='Header.TLabel').pack(pady=10)

        # Content area
        content_frame = ttk.Frame(main_container)
        content_frame.pack(fill=tk.BOTH, expand=True)

        # Left panel (Input)
        left_panel = ttk.Frame(content_frame)
        left_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 10))

        # Sequence input section
        seq_frame = ttk.LabelFrame(left_panel, text=" Input Sequences ", padding=10)
        seq_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # Sequence 1
        seq1_frame = ttk.Frame(seq_frame, style='Sequence.TFrame')
        seq1_frame.pack(fill=tk.X, pady=5)
        ttk.Label(seq1_frame, text="Sequence 1:").pack(side=tk.LEFT, padx=(0, 5))
        self.seq1_btn = ttk.Button(seq1_frame, text="Load FASTA", command=self.load_fasta1, width=10)
        self.seq1_btn.pack(side=tk.RIGHT)
        self.seq1_text = ScrolledText(seq_frame, height=8, wrap=tk.WORD, font=self.mono_font)
        self.seq1_text.pack(fill=tk.BOTH, expand=True)

        # Sequence 2
        seq2_frame = ttk.Frame(seq_frame, style='Sequence.TFrame')
        seq2_frame.pack(fill=tk.X, pady=(10, 5))
        ttk.Label(seq2_frame, text="Sequence 2:").pack(side=tk.LEFT, padx=(0, 5))
        self.seq2_btn = ttk.Button(seq2_frame, text="Load FASTA", command=self.load_fasta2, width=10)
        self.seq2_btn.pack(side=tk.RIGHT)
        self.seq2_text = ScrolledText(seq_frame, height=8, wrap=tk.WORD, font=self.mono_font)
        self.seq2_text.pack(fill=tk.BOTH, expand=True)

        # Parameters and Actions container
        param_action_container = ttk.Frame(left_panel)
        param_action_container.pack(fill=tk.X, pady=(0, 10))

        # Parameters section
        param_frame = ttk.LabelFrame(param_action_container, text=" Alignment Parameters ", padding=10)
        param_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 5))

        # Parameter controls
        param_grid = ttk.Frame(param_frame)
        param_grid.pack(fill=tk.X)

        ttk.Label(param_grid, text="Match Score:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.E)
        self.match_entry = ttk.Entry(param_grid, width=8)
        self.match_entry.insert(0, "1")
        self.match_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W)

        ttk.Label(param_grid, text="Mismatch Penalty:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.E)
        self.mismatch_entry = ttk.Entry(param_grid, width=8)
        self.mismatch_entry.insert(0, "-1")
        self.mismatch_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W)

        ttk.Label(param_grid, text="Gap Penalty:").grid(row=2, column=0, padx=5, pady=5, sticky=tk.E)
        self.gap_entry = ttk.Entry(param_grid, width=8)
        self.gap_entry.insert(0, "-2")
        self.gap_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W)

        # Action buttons section
        btn_frame = ttk.LabelFrame(param_action_container, text=" Actions ", padding=10)
        btn_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0))

        self.align_btn = ttk.Button(btn_frame, text="Align Sequences", command=self.run_alignment,
                                    style='Accent.TButton')
        self.align_btn.pack(fill=tk.X, pady=5)

        self.clear_btn = ttk.Button(btn_frame, text="Clear All", command=self.clear_all)
        self.clear_btn.pack(fill=tk.X, pady=5)

        # Right panel (Output)
        right_panel = ttk.Frame(content_frame)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Results notebook (tabbed interface)
        self.results_notebook = ttk.Notebook(right_panel)
        self.results_notebook.pack(fill=tk.BOTH, expand=True)

        # Alignment results tab
        align_tab = ttk.Frame(self.results_notebook)
        self.results_notebook.add(align_tab, text="Alignment Results")

        self.result_text = ScrolledText(align_tab, wrap=tk.WORD, font=self.mono_font)
        self.result_text.pack(fill=tk.BOTH, expand=True)

        # Matrix tab
        matrix_tab = ttk.Frame(self.results_notebook)
        self.results_notebook.add(matrix_tab, text="Scoring Matrix")

        self.canvas_frame = ttk.Frame(matrix_tab)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True)

        # Status bar
        self.status_bar = ttk.Frame(main_container, height=25, relief=tk.SUNKEN)
        self.status_bar.pack(fill=tk.X, pady=(10, 0))
        self.status_label = ttk.Label(self.status_bar, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_label.pack(fill=tk.X)

        # Configure grid weights for resizing
        content_frame.columnconfigure(0, weight=1)
        content_frame.columnconfigure(1, weight=1)

        # Add some color to buttons
        self.style.configure('Accent.TButton', foreground='white', background='#4caf50')
        self.style.map('Accent.TButton',
                       background=[('active', '#388e3c'), ('pressed', '#2e7d32')])

    def update_status(self, message):
        self.status_label.config(text=message)
        self.root.update_idletasks()

    def load_fasta1(self):
        self.load_fasta(1)

    def load_fasta2(self):
        self.load_fasta(2)

    def load_fasta(self, seq_num):
        filepath = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fasta;*.fa"), ("All files", "*.*")],
            title=f"Select Sequence {seq_num} FASTA File"
        )
        if not filepath:
            return

        self.update_status(f"Loading Sequence {seq_num} from {os.path.basename(filepath)}...")

        try:
            with open(filepath, 'r') as f:
                content = f.read()
                sequences = []
                current_seq = []
                for line in content.split('\n'):
                    if line.startswith('>'):
                        if current_seq:
                            sequences.append(''.join(current_seq))
                            current_seq = []
                    else:
                        current_seq.append(line.strip())
                if current_seq:
                    sequences.append(''.join(current_seq))

                if not sequences:
                    messagebox.showerror("Error", "No sequence found in the FASTA file")
                    return

                sequence = sequences[0]
                if seq_num == 1:
                    self.seq1_text.delete(1.0, tk.END)
                    self.seq1_text.insert(tk.END, sequence)
                    self.seq1_btn.config(text=f"Loaded: {os.path.basename(filepath)[:15]}...")
                else:
                    self.seq2_text.delete(1.0, tk.END)
                    self.seq2_text.insert(tk.END, sequence)
                    self.seq2_btn.config(text=f"Loaded: {os.path.basename(filepath)[:15]}...")

                self.update_status(f"Sequence {seq_num} loaded successfully")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to load FASTA file: {str(e)}")
            self.update_status("Error loading sequence")

    def clear_all(self):
        self.seq1_text.delete(1.0, tk.END)
        self.seq2_text.delete(1.0, tk.END)
        self.result_text.delete(1.0, tk.END)
        self.match_entry.delete(0, tk.END)
        self.match_entry.insert(0, "1")
        self.mismatch_entry.delete(0, tk.END)
        self.mismatch_entry.insert(0, "-1")
        self.gap_entry.delete(0, tk.END)
        self.gap_entry.insert(0, "-2")
        self.seq1_btn.config(text="Load FASTA")
        self.seq2_btn.config(text="Load FASTA")

        # Clear the canvas if it exists
        for widget in self.canvas_frame.winfo_children():
            widget.destroy()

        self.update_status("Cleared all inputs and results")

    def validate_input(self):
        # Get sequences
        self.seq1 = self.seq1_text.get(1.0, tk.END).strip().upper()
        self.seq2 = self.seq2_text.get(1.0, tk.END).strip().upper()

        if not self.seq1 or not self.seq2:
            messagebox.showerror("Error", "Both sequences must be provided")
            return False

        # Validate sequences (DNA or protein)
        valid_dna = set('ACGTU')
        valid_protein = set('ACDEFGHIKLMNPQRSTVWY')

        is_dna1 = all(c in valid_dna for c in self.seq1)
        is_dna2 = all(c in valid_dna for c in self.seq2)
        is_protein1 = all(c in valid_protein for c in self.seq1)
        is_protein2 = all(c in valid_protein for c in self.seq2)

        if not (is_dna1 or is_protein1) or not (is_dna2 or is_protein2):
            messagebox.showerror("Error", "Sequences must contain only valid DNA or protein characters")
            return False

        if (is_dna1 and not is_dna2) or (is_protein1 and not is_protein2):
            messagebox.showerror("Error", "Both sequences must be of the same type (DNA or protein)")
            return False

        # Validate parameters
        try:
            self.match_score = int(self.match_entry.get())
            self.mismatch_penalty = int(self.mismatch_entry.get())
            self.gap_penalty = int(self.gap_entry.get())
        except ValueError:
            messagebox.showerror("Error", "All parameters must be integers")
            return False

        return True

    def run_alignment(self):
        if not self.validate_input():
            return

        try:
            self.update_status("Running Needleman-Wunsch alignment...")
            self.align_btn.config(state=tk.DISABLED)
            self.root.update_idletasks()

            # Run Needleman-Wunsch algorithm
            score_matrix, traceback_matrix = self.compute_score_matrix()
            alignments = self.traceback(score_matrix, traceback_matrix)

            # Display results
            self.display_results(score_matrix, alignments)

            # Display matrix
            self.display_matrix(score_matrix, traceback_matrix)

            self.update_status(f"Alignment complete. Found {len(alignments)} optimal alignment(s).")
            self.align_btn.config(state=tk.NORMAL)

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during alignment: {str(e)}")
            self.update_status("Alignment failed")
            self.align_btn.config(state=tk.NORMAL)

    def compute_score_matrix(self):
        m = len(self.seq1)
        n = len(self.seq2)

        # Initialize matrices
        score_matrix = np.zeros((m + 1, n + 1), dtype=int)
        traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)  # 0=diag, 1=up, 2=left

        # Initialize first row and column
        for i in range(1, m + 1):
            score_matrix[i, 0] = i * self.gap_penalty
            traceback_matrix[i, 0] = 1  # up

        for j in range(1, n + 1):
            score_matrix[0, j] = j * self.gap_penalty
            traceback_matrix[0, j] = 2  # left

        # Fill the matrices
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = score_matrix[i - 1, j - 1] + (
                    self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_penalty)
                delete = score_matrix[i - 1, j] + self.gap_penalty
                insert = score_matrix[i, j - 1] + self.gap_penalty

                score_matrix[i, j] = max(match, delete, insert)

                # Determine direction(s) for traceback
                traceback_matrix[i, j] = 0  # default to diagonal
                if score_matrix[i, j] == delete:
                    traceback_matrix[i, j] = 1  # up
                elif score_matrix[i, j] == insert:
                    traceback_matrix[i, j] = 2  # left

        return score_matrix, traceback_matrix

    def traceback(self, score_matrix, traceback_matrix):
        alignments = []
        self.traceback_recursive(score_matrix, traceback_matrix,
                                 len(self.seq1), len(self.seq2),
                                 "", "", alignments)
        return alignments

    def traceback_recursive(self, score_matrix, traceback_matrix, i, j, seq1_part, seq2_part, alignments):
        if i == 0 and j == 0:
            alignments.append((seq1_part[::-1], seq2_part[::-1]))
            return

        current_cell = traceback_matrix[i, j]

        # Check if we're at the beginning of a row or column
        if i == 0:
            self.traceback_recursive(score_matrix, traceback_matrix, i, j - 1,
                                     seq1_part + "-", seq2_part + self.seq2[j - 1], alignments)
            return
        elif j == 0:
            self.traceback_recursive(score_matrix, traceback_matrix, i - 1, j,
                                     seq1_part + self.seq1[i - 1], seq2_part + "-", alignments)
            return

        # Check all possible paths (for multiple optimal alignments)
        paths = []
        if current_cell == 0 or (i > 0 and j > 0 and
                                 score_matrix[i, j] == score_matrix[i - 1, j - 1] +
                                 (self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_penalty)):
            paths.append((i - 1, j - 1, self.seq1[i - 1], self.seq2[j - 1]))

        if current_cell == 1 or (i > 0 and
                                 score_matrix[i, j] == score_matrix[i - 1, j] + self.gap_penalty):
            paths.append((i - 1, j, self.seq1[i - 1], "-"))

        if current_cell == 2 or (j > 0 and
                                 score_matrix[i, j] == score_matrix[i, j - 1] + self.gap_penalty):
            paths.append((i, j - 1, "-", self.seq2[j - 1]))

        # Follow all possible paths
        for path in paths:
            self.traceback_recursive(score_matrix, traceback_matrix,
                                     path[0], path[1],
                                     seq1_part + path[2], seq2_part + path[3],
                                     alignments)

    def display_results(self, score_matrix, alignments):
        self.result_text.delete(1.0, tk.END)
        self.result_text.tag_configure('header', font=self.bold_font)
        self.result_text.tag_configure('highlight', foreground='blue')
        self.result_text.tag_configure('match', foreground='green')
        self.result_text.tag_configure('mismatch', foreground='black')
        self.result_text.tag_configure('gap', foreground='red')

        # Display parameters
        self.result_text.insert(tk.END, "=== Alignment Parameters ===\n", 'header')
        self.result_text.insert(tk.END, f"Match score: {self.match_score}\n")
        self.result_text.insert(tk.END, f"Mismatch penalty: {self.mismatch_penalty}\n")
        self.result_text.insert(tk.END, f"Gap penalty: {self.gap_penalty}\n\n")

        # Display sequences info
        self.result_text.insert(tk.END, "=== Sequences ===\n", 'header')
        self.result_text.insert(tk.END, f"Sequence 1 length: {len(self.seq1)}\n")
        self.result_text.insert(tk.END, f"Sequence 2 length: {len(self.seq2)}\n\n")

        # Display alignment score
        alignment_score = score_matrix[len(self.seq1), len(self.seq2)]
        self.result_text.insert(tk.END, f"Alignment score: {alignment_score}\n\n", 'highlight')

        # Display number of optimal alignments
        self.result_text.insert(tk.END, f"Number of optimal alignments found: {len(alignments)}\n\n", 'highlight')

        # Display all optimal alignments
        for alignment_num, (aligned_seq1, aligned_seq2) in enumerate(alignments, 1):
            self.result_text.insert(tk.END, f"=== Alignment {alignment_num} ===\n", 'header')

            alignment_length = len(aligned_seq1)

            # Calculate identity and gap percentage
            matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
            gaps = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == '-' or b == '-')

            identity_percent = (matches / alignment_length) * 100
            gap_percent = (gaps / alignment_length) * 100

            self.result_text.insert(tk.END, f"Alignment length: {alignment_length}\n")
            self.result_text.insert(tk.END, f"Identical positions: {matches} ({identity_percent:.2f}%)\n", 'match')
            self.result_text.insert(tk.END, f"Gaps: {gaps} ({gap_percent:.2f}%)\n", 'gap')
            self.result_text.insert(tk.END, "\n")

            # Display the alignment
            self.result_text.insert(tk.END, "Alignment:\n", 'header')

            # Display in chunks of 100 characters
            chunk_size = 100
            for i in range(0, alignment_length, chunk_size):
                chunk_seq1 = aligned_seq1[i:i + chunk_size]
                chunk_seq2 = aligned_seq2[i:i + chunk_size]

                # Add line numbers
                self.result_text.insert(tk.END, f"Seq1 {i + 1:5}: ")
                self.insert_colored_sequence(chunk_seq1, chunk_seq2, is_seq1=True)
                self.result_text.insert(tk.END, f"\nSeq2 {i + 1:5}: ")
                self.insert_colored_sequence(chunk_seq2, chunk_seq1, is_seq1=False)
                self.result_text.insert(tk.END, "\n")

                # Add match line
                match_line = []
                for a, b in zip(chunk_seq1, chunk_seq2):
                    if a == b:
                        match_line.append("|")
                    elif a == '-' or b == '-':
                        match_line.append(" ")
                    else:
                        match_line.append("·")
                self.result_text.insert(tk.END, "          " + " ".join(match_line) + "\n\n")

            # Add separator between alignments unless it's the last one
            if alignment_num < len(alignments):
                self.result_text.insert(tk.END, "\n" + "=" * 80 + "\n\n")

    def insert_colored_sequence(self, seq, other_seq, is_seq1):
        for i, char in enumerate(seq):
            if char == '-':
                self.result_text.insert(tk.END, char, 'gap')
            elif i < len(other_seq) and char == other_seq[i]:
                self.result_text.insert(tk.END, char, 'match')
            else:
                self.result_text.insert(tk.END, char, 'mismatch')

    def display_matrix(self, score_matrix, traceback_matrix):
        # Clear previous canvas
        for widget in self.canvas_frame.winfo_children():
            widget.destroy()

        # Create figure with better styling
        fig, ax = plt.subplots(figsize=(10, 6))
        fig.subplots_adjust(bottom=0.2, left=0.2)
        fig.patch.set_facecolor('#f5f5f5')
        ax.set_facecolor('#fafafa')

        # Prepare data for display
        rows = [' '] + list(self.seq1)
        cols = [' '] + list(self.seq2)

        # Display the matrix
        ax.set_xticks(np.arange(len(cols)))
        ax.set_yticks(np.arange(len(rows)))
        ax.set_xticklabels(cols)
        ax.set_yticklabels(rows)

        # Rotate the tick labels and set their alignment
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        # Add grid
        ax.grid(which='both', color='#e0e0e0', linestyle='-', linewidth=0.5)

        # Loop over data dimensions and create text annotations
        for i in range(len(rows)):
            for j in range(len(cols)):
                # Add arrow or symbol indicating direction
                direction = ""
                if i > 0 and j > 0:
                    if traceback_matrix[i, j] == 0:
                        direction = "↖"
                    elif traceback_matrix[i, j] == 1:
                        direction = "↑"
                    elif traceback_matrix[i, j] == 2:
                        direction = "←"

                # Color the cell based on direction
                cell_color = '#ffffff'
                if i > 0 and j > 0:
                    if traceback_matrix[i, j] == 0:
                        cell_color = '#e8f5e9'  # light green for diagonal
                    elif traceback_matrix[i, j] == 1:
                        cell_color = '#e3f2fd'  # light blue for up
                    elif traceback_matrix[i, j] == 2:
                        cell_color = '#ffebee'  # light red for left

                # Draw rectangle with color
                rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                     facecolor=cell_color, edgecolor='#e0e0e0')
                ax.add_patch(rect)

                # Add text
                ax.text(j, i, f"{score_matrix[i, j]}\n{direction}",
                        ha="center", va="center", color="black",
                        fontsize=8 if len(rows) < 20 else 6)

        ax.set_title("Scoring Matrix with Traceback Directions", pad=20)
        ax.set_xlabel("Sequence 2", labelpad=10)
        ax.set_ylabel("Sequence 1", labelpad=10)

        # Create canvas
        canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def save_results(self):
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="Save alignment results"
        )

        if not filepath:
            return

        try:
            with open(filepath, 'w') as f:
                # Write parameters
                f.write("=== Alignment Parameters ===\n")
                f.write(f"Match score: {self.match_score}\n")
                f.write(f"Mismatch penalty: {self.mismatch_penalty}\n")
                f.write(f"Gap penalty: {self.gap_penalty}\n\n")

                # Write sequences info
                f.write("=== Sequences ===\n")
                f.write(f"Sequence 1: {self.seq1}\n")
                f.write(f"Sequence 2: {self.seq2}\n")
                f.write(f"Sequence 1 length: {len(self.seq1)}\n")
                f.write(f"Sequence 2 length: {len(self.seq2)}\n\n")

                # Write alignment results
                results = self.result_text.get(1.0, tk.END)
                f.write(results)

            messagebox.showinfo("Success", f"Results saved to {filepath}")
            self.update_status(f"Results saved to {os.path.basename(filepath)}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {str(e)}")
            self.update_status("Error saving results")


if __name__ == "__main__":
    root = tk.Tk()
    app = NeedlemanWunschApp(root)
    root.mainloop()