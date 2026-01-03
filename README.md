# PETase Engineering Tournament 2025 🧬

Welcome to the official repository for our team's entry in the **2025 PETase Engineering Tournament**. This project aims to engineer enhanced variants of the *Ideonella sakaiensis* PETase enzyme with improved thermostability and catalytic efficiency for plastic degradation.

## 🏆 Competition Tracks

We are competing in the **Predictive Track**, which consists of two phases:

1.  **Zero-Shot Prediction**: Ranking variants without training data using evolutionary models.
2.  **Supervised Learning**: Training models on experimental data to predict activity and stability.

## 📂 Repository Structure

```text
.
├── data/                   # Dataset files (Wild Type, Test Sets)
├── src/
│   └── zero_shot/          # Zero-Shot inference pipeline
│       ├── score.py        # Main scoring script (ESM-2 PLL)
│       └── utils.py        # Sequence handling and alignment
├── notebooks/              # Exploratory Data Analysis (EDA)
├── PETASE_BLUEPRINT.md     # Detailed technical plan and methodology
├── CHANGELOG.md            # Record of project updates
├── requirements.txt        # Python dependencies
└── README.md               # Project documentation
```

## 🚀 Getting Started

### Prerequisites

*   Python 3.9+
*   CUDA-capable GPU (recommended) or Apple Silicon (MPS)

### Installation

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/persinammon4/petase-tournament.git
    cd petase-tournament
    ```

2.  **Install dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

## 🧠 Methodology: Zero-Shot Track

For the Zero-Shot phase, we utilize the **ESM-2 (Evolutionary Scale Modeling)** protein language model.

### Scoring Strategy: Pseudo-Log-Likelihood (PLL)
To robustly handle insertions and deletions (indels) in the mutant sequences, we calculate the **Pseudo-Log-Likelihood (PLL)** of the entire sequence rather than just scoring mutated positions.

$$
\text{Score} = \text{PLL}(\text{Mutant}) - \text{PLL}(\text{Wild Type})
$$

Where $\text{PLL}(x) = \sum_{i=1}^{L} \log P(x_i \mid x_{\setminus i})$. A higher score indicates that the language model perceives the mutant sequence as more evolutionarily "natural" or fit compared to the wild type.

## 💻 Usage

To run the Zero-Shot scorer and generate predictions:

```bash
python src/zero_shot/score.py
```

This will:
1.  Load the `esm2_t33_650M_UR50D` model.
2.  Compute the PLL for the Wild Type.
3.  Score all variants in `data/predictive-pet-zero-shot-test-2025.csv`.
4.  Save the ranked results to `submission_zero_shot.csv`.

### Activity Predictions in micromoles

Copy of ChatGPT vibe coding conversation for activity predictions at pH 5.5 and pH 9.0:
https://chatgpt.com/share/695887b5-e000-8004-a4c1-5901ef01116e

There are two models for predicting activity one coded in `seq2ph.py` and the other in `seq2ph_with_solubility.py`. Solubility scores are taken from [an academia based model](https://loschmidt.chemi.muni.cz/soluprot).
Originally, solubility was meant to be used to predict expression.

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
*Built for the 2025 PETase Tournament.*
