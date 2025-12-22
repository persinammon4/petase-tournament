# PETase 2025 Zero-Shot Predictive Track Blueprint

## 1. Title & Objective
**Title**: Zero-Shot Ranking of PETase Variants using ESM-2 Masked Marginal Probability
**Objective**: To rank a set of PETase variants based on their predicted fitness (activity/expression) using the ESM-2 protein language model in a zero-shot manner, without training on labeled data.

## 2. Data Architecture

### Input Data
*   **File**: `predictive-pet-zero-shot-test-2025.csv`
*   **Description**: Contains the target mutant sequences to be ranked.
*   **Schema**:
    *   `sequence`: The full amino acid sequence of the PETase variant.
    *   `activity_1 (μmol [TPA]/min·mg [E])`: Empty in test set (Target).
    *   `activity_2 (μmol [TPA]/min·mg [E])`: Empty in test set (Target).
    *   `expression (mg/mL)`: Empty in test set (Target).

### Reference Data
*   **File**: `pet-2025-wildtype-cds.csv`
*   **Description**: Contains reference sequences.
*   **Wild Type (WT)**: The first sequence in the `Wt AA Sequence` column will be used as the canonical Wild Type for mutation identification.
    *   *WT Sequence Start*: `AADNPYQRGPDPTNASIEAATGPFAVGT...`

## 3. Methodology (The "Alpha")

### Model
*   **Architecture**: ESM-2 (Evolutionary Scale Modeling).
*   **Checkpoint**: `esm2_t33_650M_UR50D` (650M parameters) is recommended as a balance between performance and inference speed. Larger models (3B, 15B) can be used if compute allows.

### Scoring Strategy: Pseudo-Log-Likelihood (PLL)
We use the **Pseudo-Log-Likelihood (PLL)** of the entire sequence. This method is robust to insertions and deletions (indels) because it scores the sequence as a whole rather than relying on fixed-position masking.

### Formula
For a sequence $x$ of length $L$:

$$
\text{PLL}(x) = \sum_{i=1}^{L} \log P(x_i \mid x_{\setminus i})
$$

Where:
*   $x_{\setminus i}$ is the sequence with the position $i$ masked.
*   $P(x_i \mid x_{\setminus i})$ is the probability assigned by the model to the true amino acid $x_i$ at the masked position.

The final ranking score is the **Log-Likelihood Ratio (LLR)**:
$$
\text{Score} = \text{PLL}(Mutant) - \text{PLL}(WT)
$$

### Inference Process
1.  **Mutation Identification**: Use pairwise alignment to identify if the sequence is a mutant or indel (for logging purposes).
2.  **PLL Calculation**:
    *   For a given sequence, iteratively mask every position $i$ from $1$ to $L$.
    *   Compute the log-probability of the true residue at $i$.
    *   Sum these log-probabilities to get the sequence PLL.
3.  **Scoring**: Compute $\text{Score} = \text{PLL}(Mutant) - \text{PLL}(WT)$.
4.  **Aggregation**: Store the score.

## 4. Implementation Plan

### Phase 1: Setup & Data Prep
*   Initialize Python environment with `torch`, `transformers` (Hugging Face), `pandas`, and `biopython`.
*   Load `pet-2025-wildtype-cds.csv` and extract the WT sequence.
*   Load `predictive-pet-zero-shot-test-2025.csv`.
*   Implement a `get_mutations(wt, mutant)` function using `Bio.Align` to handle indels.

### Phase 2: Zero-Shot Scorer
*   Load the ESM-2 model and tokenizer.
*   Implement the `compute_pseudo_likelihood` function.
*   Implement the inference loop:
    *   Compute PLL for WT once.
    *   Iterate through test sequences and compute PLL.
    *   Calculate LLR.
    *   *Optimization*: Batch masking to maximize GPU utilization.

### Phase 3: Submission Generation
*   Sort the dataframe by `score` (Descending order: higher LLR = better predicted fitness).
*   Format the output CSV to match submission requirements (if any specific format is needed beyond the score).
*   Save as `submission_zero_shot.csv`.

## 5. Repository Structure

```text
.
├── PETASE_BLUEPRINT.md
├── data/
│   ├── pet-2025-wildtype-cds.csv
│   └── predictive-pet-zero-shot-test-2025.csv
├── src/
│   └── zero_shot/
│       ├── __init__.py
│       ├── score.py        # Main inference script
│       └── utils.py        # Mutation parsing and data loading
├── notebooks/
│   └── exploration.ipynb   # For testing model loading and small scale inference
└── requirements.txt
```

## 6. Evaluation

*   **Primary Metric**: **NDCG (Normalized Discounted Cumulative Gain)**.
    *   This measures the quality of the ranking, giving more importance to correctly ranking the top-performing variants.
*   **Proxy Metric**: Since we do not have ground truth labels for the test set yet, we rely on the **LLR Score** itself. Higher LLR indicates the model perceives the mutation as more evolutionarily favorable, which correlates with protein stability and function.
