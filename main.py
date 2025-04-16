import sys
import os
import matplotlib.pyplot as plt
from collections import deque
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox, scrolledtext


def read_fasta(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    sequences = []
    current_seq = ""
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                sequences.append(current_seq)
                current_seq = ""
        else:
            current_seq += line
    if current_seq:
        sequences.append(current_seq)
    return sequences[:2]  # only need first two sequences


def initialize_matrix(n, m, gap):
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        score[i][0] = i * gap
    for j in range(m + 1):
        score[0][j] = j * gap
    return score


def fill_matrix(seq1, seq2, match, mismatch, gap):
    n, m = len(seq1), len(seq2)
    score = initialize_matrix(n, m, gap)
    path = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up = score[i-1][j] + gap
            left = score[i][j-1] + gap
            max_score = max(diag, up, left)
            score[i][j] = max_score
            if diag == max_score: path[i][j].append((i-1, j-1))
            if up == max_score: path[i][j].append((i-1, j))
            if left == max_score: path[i][j].append((i, j-1))
    return score, path


def traceback_all_paths(path, seq1, seq2):
    n, m = len(seq1), len(seq2)
    paths = []
    stack = deque([((n, m), [], [])])
    while stack:
        (i, j), a1, a2 = stack.pop()
        if i == 0 and j == 0:
            paths.append((''.join(reversed(a1)), ''.join(reversed(a2))))
            continue
        for ni, nj in path[i][j]:
            if ni == i - 1 and nj == j - 1:
                stack.append(((ni, nj), a1 + [seq1[ni]], a2 + [seq2[nj]]))
            elif ni == i - 1:
                stack.append(((ni, nj), a1 + [seq1[ni]], a2 + ['-']))
            else:
                stack.append(((ni, nj), a1 + ['-'], a2 + [seq2[nj]]))
    return paths


def calculate_stats(al1, al2):
    matches = sum(c1 == c2 for c1, c2 in zip(al1, al2))
    length = len(al1)
    gaps = al1.count('-') + al2.count('-')
    identity = (matches / length) * 100
    return length, matches, gaps, identity


def plot_score_matrix(score, path, seq1, seq2):
    plt.figure(figsize=(10, 8))
    plt.imshow(score, cmap='Blues', interpolation='none')
    plt.colorbar(label='Score')
    plt.xticks(range(len(seq2)+1), ['-']+list(seq2))
    plt.yticks(range(len(seq1)+1), ['-']+list(seq1))
    plt.title('Score Matrix with Optimal Path')

    # Plot one optimal path
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        plt.plot(j, i, 'ro')
        if (i > 0 and j > 0) and (i-1, j-1) in path[i][j]:
            i -= 1
            j -= 1
        elif i > 0 and (i-1, j) in path[i][j]:
            i -= 1
        elif j > 0 and (i, j-1) in path[i][j]:
            j -= 1

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()


def run_gui():
    def load_fasta():
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")])
        if file_path:
            try:
                seqs = read_fasta(file_path)
                if len(seqs) >= 2:
                    seq1_entry.delete("1.0", tk.END)
                    seq2_entry.delete("1.0", tk.END)
                    seq1_entry.insert(tk.END, seqs[0])
                    seq2_entry.insert(tk.END, seqs[1])
                else:
                    messagebox.showerror("Error", "FASTA must contain at least two sequences.")
            except Exception as e:
                messagebox.showerror("Error", str(e))

    def align_sequences():
        seq1 = seq1_entry.get("1.0", tk.END).strip().upper()
        seq2 = seq2_entry.get("1.0", tk.END).strip().upper()
        try:
            match = int(match_entry.get())
            mismatch = int(mismatch_entry.get())
            gap = int(gap_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Scores and penalties must be integers.")
            return

        score, path = fill_matrix(seq1, seq2, match, mismatch, gap)
        all_paths = traceback_all_paths(path, seq1, seq2)
        best_align1, best_align2 = all_paths[0]

        length, matches, gaps, identity = calculate_stats(best_align1, best_align2)

        result_output.delete("1.0", tk.END)
        result_output.insert(tk.END, f"Best alignment:\n{best_align1}\n{best_align2}\n")
        result_output.insert(tk.END, f"Score: {score[len(seq1)][len(seq2)]}\n")
        result_output.insert(tk.END, f"Length: {length}, Matches: {matches}, Gaps: {gaps}, Identity: {identity:.2f}%\n")

        with open("needleman_wunsch_result.txt", "w") as f:
            f.write(f"Sequence 1: {seq1}\n")
            f.write(f"Sequence 2: {seq2}\n")
            f.write(f"Match: {match}, Mismatch: {mismatch}, Gap: {gap}\n")
            f.write("\nBest alignment:\n")
            f.write(best_align1 + "\n")
            f.write(best_align2 + "\n")
            f.write(f"\nScore: {score[len(seq1)][len(seq2)]}\n")
            f.write(f"Length: {length}, Matches: {matches}, Gaps: {gaps}, Identity: {identity:.2f}%\n")

        plot_score_matrix(score, path, seq1, seq2)

    root = tk.Tk()
    root.title("Needleman-Wunsch Sequence Aligner")

    tk.Label(root, text="Sequence 1:").grid(row=0, column=0, sticky="nw")
    seq1_entry = scrolledtext.ScrolledText(root, height=4, width=50)
    seq1_entry.grid(row=0, column=1)

    tk.Label(root, text="Sequence 2:").grid(row=1, column=0, sticky="nw")
    seq2_entry = scrolledtext.ScrolledText(root, height=4, width=50)
    seq2_entry.grid(row=1, column=1)

    tk.Button(root, text="Load FASTA", command=load_fasta).grid(row=2, column=1, sticky="e")

    tk.Label(root, text="Match Score:").grid(row=3, column=0)
    match_entry = tk.Entry(root)
    match_entry.insert(0, "1")
    match_entry.grid(row=3, column=1, sticky="w")

    tk.Label(root, text="Mismatch Penalty:").grid(row=4, column=0)
    mismatch_entry = tk.Entry(root)
    mismatch_entry.insert(0, "-1")
    mismatch_entry.grid(row=4, column=1, sticky="w")

    tk.Label(root, text="Gap Penalty:").grid(row=5, column=0)
    gap_entry = tk.Entry(root)
    gap_entry.insert(0, "-2")
    gap_entry.grid(row=5, column=1, sticky="w")

    tk.Button(root, text="Align Sequences", command=align_sequences).grid(row=6, column=1, sticky="e")

    tk.Label(root, text="Output:").grid(row=7, column=0, sticky="nw")
    result_output = scrolledtext.ScrolledText(root, height=10, width=70)
    result_output.grid(row=7, column=1)

    root.mainloop()


if __name__ == '__main__':
    run_gui()
