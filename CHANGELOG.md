# Project Changelog

## [2025-12-21] - Zero-Shot Track Initialization & Refinement

### Added
- **`PETASE_BLUEPRINT.md`**: Created the initial comprehensive blueprint for the Zero-Shot Predictive Track.
- **`CHANGELOG.md`**: Created this file to track detailed project changes.
- **`src/zero_shot/`**: Created directory structure for zero-shot code.
- **`data/`**: Created directory for dataset storage.

### Changed
- **`requirements.txt`**: Added `biopython` to support sequence alignment for handling insertions/deletions (indels).
- **`src/zero_shot/utils.py`**: 
    - Updated `get_mutations` to use `Bio.Align.PairwiseAligner`.
    - Replaced simple length check with global alignment to robustly identify mutations even when sequence lengths differ (indels).
- **`src/zero_shot/score.py`**:
    - **Major Strategy Shift**: Switched from *Masked Marginal Probability* (scoring only mutated positions) to *Pseudo-Log-Likelihood (PLL)* (scoring the entire sequence).
    - **Reason**: The dataset contains indels (length mismatches), making position-specific masking ambiguous and error-prone. PLL provides a robust, length-independent fitness score by summing the conditional log-probabilities of all residues in the sequence.
- **`PETASE_BLUEPRINT.md`**:
    - Updated **Methodology** section to reflect the switch to Pseudo-Log-Likelihood (PLL).
    - Updated **Formula** to the PLL equation: $\text{Score} = \sum_{i=1}^{L} \log P(x_i \mid x_{\setminus i})$.

## [2025-12-24] - Performance Optimization

### Changed
- **`src/zero_shot/score.py`**:
    - **Vectorization**: Rewrote `compute_pseudo_likelihood` to construct the entire masked batch at once instead of iterating loop-by-loop. This reduces Python overhead and maximizes GPU parallelism.
    - **Mixed Precision**: Added support for `float16` (FP16) inference on CUDA devices, doubling throughput and reducing memory footprint.
    - **Batching**: Increased default batch size to 128 (configurable), allowing entire sequences to be processed in a single forward pass on A100 GPUs.
- **`PETASE_BLUEPRINT.md`**: Updated the "Inference Process" to document the vectorized batching strategy.

