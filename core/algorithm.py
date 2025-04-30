import numpy as np
from .constants import DNA_BASES, PROTEIN_RESIDUES


class NeedlemanWunsch:
    """Implementation of the Needleman-Wunsch global sequence alignment algorithm.

    This class provides functionality for aligning biological sequences (DNA or protein)
    using dynamic programming with configurable scoring parameters.

    Attributes:
        match_score (int): Score for matching characters (default: 1)
        mismatch_penalty (int): Penalty for mismatches (default: 0)
        gap_penalty (int): Penalty for gaps (default: -1)
        seq1 (str): First sequence to align
        seq2 (str): Second sequence to align
        seq_type (str): Type of sequence ('dna' or 'protein')
        score_matrix (numpy.ndarray): Dynamic programming score matrix
        alignment_paths (list): List of optimal alignment paths
        total_optimal_paths (int): Total number of optimal paths found
        alignment_examples (list): Example alignments generated from paths
        alignment_score (int): Optimal alignment score
    """

    def __init__(self, match_score=1, mismatch_penalty=0, gap_penalty=-1):
        """Initializes the NeedlemanWunsch aligner with scoring parameters.

        Args:
            match_score (int): Reward for matching characters
            mismatch_penalty (int): Penalty for mismatched characters
            gap_penalty (int): Penalty for introducing gaps
        """
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
        """Sets the sequences to be aligned and detects their type.

        Args:
            seq1 (str): First biological sequence
            seq2 (str): Second biological sequence
            seq_type (str): Sequence type ('dna', 'protein', or 'auto' for detection)

        Raises:
            ValueError: If sequences contain invalid characters for their type
        """
        self.seq1 = seq1.upper()
        self.seq2 = seq2.upper()
        self.seq_type = seq_type if seq_type != 'auto' else self.detect_sequence_type()

    def detect_sequence_type(self):
        """Determines whether sequences are DNA or protein based on their characters.

        Returns:
            str: 'dna' if sequences contain only DNA bases, 'protein' if only amino acids

        Raises:
            ValueError: If sequences contain invalid characters for both DNA and protein
        """
        is_dna = all(c in DNA_BASES for c in self.seq1 + self.seq2)
        is_protein = all(c in PROTEIN_RESIDUES for c in self.seq1 + self.seq2)

        if is_dna:
            return 'dna'
        if is_protein:
            return 'protein'
        raise ValueError("Sequences contain invalid characters")

    def align(self):
        """Performs global sequence alignment using the Needleman-Wunsch algorithm.

        The algorithm proceeds in three main steps:
        1. Initialization of the score matrix
        2. Matrix filling with dynamic programming
        3. Traceback to find optimal alignments

        Returns:
            dict: Alignment results including score, paths, and examples

        Raises:
            ValueError: If sequences are empty or not set
        """
        if not self.seq1 or not self.seq2:
            raise ValueError("Sequences cannot be empty")

        len1, len2 = len(self.seq1), len(self.seq2)

        # Initialize score matrix with gap penalties in first row/column
        self.score_matrix = np.zeros((len1 + 1, len2 + 1), dtype=int)
        for i in range(len1 + 1):
            self.score_matrix[i][0] = i * self.gap_penalty
        for j in range(len2 + 1):
            self.score_matrix[0][j] = j * self.gap_penalty

        # Fill the score matrix
        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                # Calculate scores for three possible moves:
                match = self.score_matrix[i - 1][j - 1] + (
                    self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_penalty
                )
                delete = self.score_matrix[i - 1][j] + self.gap_penalty
                insert = self.score_matrix[i][j - 1] + self.gap_penalty

                # Take the maximum score for the current cell
                self.score_matrix[i][j] = max(match, delete, insert)

        self.alignment_score = self.score_matrix[len1][len2]
        self.find_all_optimal_paths(len1, len2)
        return self.get_results()

    def find_all_optimal_paths(self, i, j, path=None, paths=None, max_paths=10):
        """Recursively finds all optimal alignment paths through the score matrix.

        The number of optimal paths can grow exponentially with sequence length. For two
        sequences of length n, there can be O(2^n) optimal paths in the worst case.
        To maintain reasonable performance, we limit the search to max_paths (default: 10).

        Example:
            For sequences of length 100, there could theoretically be up to ~1.26e+30
            optimal paths (Catalan number growth). Limiting to 10 paths makes the traceback
            feasible while still showing representative alignments.

        Args:
            i (int): Current row position in matrix (seq1 index)
            j (int): Current column position in matrix (seq2 index)
            path (list, optional): Current path being explored. Defaults to None.
            paths (list, optional): Accumulated optimal paths. Defaults to None.
            max_paths (int, optional): Maximum paths to find before stopping. Defaults to 10.

        Note:
            The 10-path limit represents a practical tradeoff between:
            1. Showing multiple alignment variants
            2. Avoiding exponential time complexity
            3. Preventing memory overload for long sequences
        """
        if paths is None:
            paths = []
        if path is None:
            path = []

        # Early termination to prevent combinatorial explosion
        if len(paths) >= max_paths:
            return

        path.append((i, j))

        # Base case: origin reached (complete path found)
        if i == 0 and j == 0:
            paths.append(path[::-1])  # Store reversed path
            self.total_optimal_paths = len(paths)
            return

        # Recursive case: explore valid moves
        current_score = self.score_matrix[i][j]
        possible_moves = []

        # 1. Diagonal move (match/mismatch)
        if i > 0 and j > 0:
            diagonal_score = self.score_matrix[i - 1][j - 1]
            match_value = self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_penalty
            if current_score == diagonal_score + match_value:
                possible_moves.append((i - 1, j - 1))

        # 2. Vertical move (gap in seq2)
        if i > 0 and current_score == self.score_matrix[i - 1][j] + self.gap_penalty:
            possible_moves.append((i - 1, j))

        # 3. Horizontal move (gap in seq1)
        if j > 0 and current_score == self.score_matrix[i][j - 1] + self.gap_penalty:
            possible_moves.append((i, j - 1))

         # Recursively explore all valid moves
        for move in possible_moves:
            if len(paths) >= max_paths:
                break
            self.find_all_optimal_paths(move[0], move[1], path.copy(), paths, max_paths)

        self.alignment_paths = paths
        self.total_optimal_paths = len(paths)

    def generate_alignment_examples(self, max_examples=10):
        """Generates alignment strings from the found optimal paths.

        Args:
            max_examples (int): Maximum number of alignments to generate

        Returns:
            list: Tuples of (aligned_seq1, aligned_seq2, match_symbols)
        """
        examples = []
        for path in self.alignment_paths[:max_examples]:
            align1, align2, symbols = [], [], []

            # Reconstruct alignment by following each path
            for k in range(1, len(path)):
                i1, j1 = path[k - 1]
                i2, j2 = path[k]

                if i2 == i1 + 1 and j2 == j1 + 1:  # Diagonal move (match/mismatch)
                    align1.append(self.seq1[i1])
                    align2.append(self.seq2[j1])
                    symbols.append('|' if self.seq1[i1] == self.seq2[j1] else '*')
                elif i2 == i1 + 1:  # Vertical move (gap in seq2)
                    align1.append(self.seq1[i1])
                    align2.append('-')
                    symbols.append(' ')
                else:  # Horizontal move (gap in seq1)
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
        """Compiles and returns all alignment results in a structured format.

        Returns:
            dict: Contains:
                - score (int): Optimal alignment score
                - total_paths (int): Number of optimal paths found
                - examples (list): Alignment examples
                - identity_percentage (float): Percentage of matching positions
                - gap_percentage (float): Percentage of gap positions
                - score_matrix (numpy.ndarray): Full DP matrix
                - seq_type (str): Type of sequences aligned
        """
        self.alignment_examples = self.generate_alignment_examples()

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