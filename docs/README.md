# Needleman-Wunsch Global Sequence Alignment

Python implementation of the Needleman-Wunsch dynamic programming algorithm for global sequence alignment, with a PyQt5 GUI for visualization.

## Features

- DNA and protein sequence support with automatic type detection
- Configurable scoring parameters: match score, mismatch penalty, gap penalty
- Interactive score matrix visualization with highlighted optimal path(s)
- Multiple optimal alignment paths (up to 10) with navigation between them
- Alignment statistics: identity percentage, gap percentage, alignment score
- FASTA file import and results export

## Usage

```bash
pip install -r requirements.txt
python main.py
```

## Tech Stack

Python, PyQt5, NumPy

## Documentation

Full technical report including algorithm details, complexity analysis, and example results is available in [`docs/report.pdf`](docs/report.pdf).
