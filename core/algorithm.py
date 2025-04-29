import numpy as np
from .constants import DNA_BASES, PROTEIN_RESIDUES


class NeedlemanWunsch:
    def __init__(self, match_score=1, mismatch_penalty=0, gap_penalty=-1):
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty
        self.seq1 = ""
        self.seq2 = ""
        self.seq_type = None
        self.score_matrix = []
        self.alignment_paths = []
        self.total_optimal_paths = 0
        self.alignment_examples = []
        self.alignment_score = 0

    def set_sequences(self, seq1, seq2, seq_type='auto'):
        self.seq1 = seq1.upper()
        self.seq2 = seq2.upper()
        self.seq_type = seq_type if seq_type != 'auto' else self.detect_sequence_type()

    def detect_sequence_type(self):
        is_dna = all(c in DNA_BASES for c in self.seq1 + self.seq2)
        is_protein = all(c in PROTEIN_RESIDUES for c in self.seq1 + self.seq2)

        if is_dna: return 'dna'
        if is_protein: return 'protein'
        raise ValueError("Sequences contain invalid characters")

    def align(self):
        if not self.seq1 or not self.seq2:
            raise ValueError("Sequences cannot be empty")

        len1, len2 = len(self.seq1), len(self.seq2)
        self.score_matrix = np.zeros((len1 + 1, len2 + 1), dtype=int)

        #first row and column
        for i in range(len1 + 1):
            self.score_matrix[i][0] = i * self.gap_penalty
        for j in range(len2 + 1):
            self.score_matrix[0][j] = j * self.gap_penalty

        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                match = self.score_matrix[i - 1][j - 1] + (
                    self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_penalty
                )
                delete = self.score_matrix[i - 1][j] + self.gap_penalty
                insert = self.score_matrix[i][j - 1] + self.gap_penalty
                self.score_matrix[i][j] = max(match, delete, insert)

        self.alignment_score = self.score_matrix[len1][len2]
        self.find_all_optimal_paths(len1, len2)
        return self.get_results()

    def find_all_optimal_paths(self, i, j, path=None, paths=None, max_paths=10):
        if paths is None:
            paths = []
        if path is None:
            path = []

        #if we've already found 10 paths
        if len(paths) >= max_paths:
            return

        path.append((i, j))

        if i == 0 and j == 0:
            paths.append(path[::-1])
            return

        possible_moves = []
        current_score = self.score_matrix[i][j]

        # Diagonal move (match/mismatch)
        if i > 0 and j > 0:
            diagonal_score = self.score_matrix[i - 1][j - 1]
            match_value = self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_penalty
            if current_score == diagonal_score + match_value:
                possible_moves.append((i - 1, j - 1))

        # Up move (gap in seq2)
        if i > 0 and current_score == self.score_matrix[i - 1][j] + self.gap_penalty:
            possible_moves.append((i - 1, j))

        # Left move (gap in seq1)
        if j > 0 and current_score == self.score_matrix[i][j - 1] + self.gap_penalty:
            possible_moves.append((i, j - 1))

        # Recursively explore all possible moves
        for move in possible_moves:
            # Stop recursion if we've reached the maximum number of paths
            if len(paths) >= max_paths:
                break
            self.find_all_optimal_paths(move[0], move[1], path.copy(), paths, max_paths)

        self.alignment_paths = paths
        self.total_optimal_paths = len(paths)


    def generate_alignment_examples(self, max_examples=10):
        examples = []
        for path in self.alignment_paths[:max_examples]:
            align1, align2, symbols = [], [], []
            for k in range(1, len(path)):
                i1, j1 = path[k - 1]
                i2, j2 = path[k]

                if i2 == i1 + 1 and j2 == j1 + 1:  # Diagonal
                    align1.append(self.seq1[i1])
                    align2.append(self.seq2[j1])
                    symbols.append('|' if self.seq1[i1] == self.seq2[j1] else '*')
                elif i2 == i1 + 1:  # Vertical (gap in seq2)
                    align1.append(self.seq1[i1])
                    align2.append('-')
                    symbols.append(' ')
                else:  # Horizontal (gap in seq1)
                    align1.append('-')
                    align2.append(self.seq2[j1])
                    symbols.append(' ')

            examples.append((
                ''.join(align1),
                ''.join(align2),
                ''.join(symbols)
            ))
        return examples

    def get_results(self):
        self.alignment_examples = self.generate_alignment_examples()

        # Calculate identity and gap percentages from first example
        if not self.alignment_examples:
            return {
                'score': self.alignment_score,
                'total_paths': 0,
                'examples': [],
                'identity_percentage': 0,
                'gap_percentage': 0,
                'score_matrix': self.score_matrix,
                'seq_type': self.seq_type
            }

        example = self.alignment_examples[0]
        aligned_length = len(example[0])
        matches = example[2].count('|')
        gaps = example[0].count('-') + example[1].count('-')

        identity = (matches / aligned_length * 100) if aligned_length > 0 else 0
        gap_percent = (gaps / aligned_length * 100) if aligned_length > 0 else 0

        return {
            'score': self.alignment_score,
            'total_paths': self.total_optimal_paths,
            'examples': self.alignment_examples,
            'identity_percentage': identity,
            'gap_percentage': gap_percent,
            'score_matrix': self.score_matrix,
            'seq_type': self.seq_type
        }
