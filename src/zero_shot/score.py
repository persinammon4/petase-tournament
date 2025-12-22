import torch
import pandas as pd
from tqdm import tqdm
from transformers import EsmForMaskedLM, AutoTokenizer
import sys
import os
import numpy as np

# Add src to path to import utils
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from src.zero_shot.utils import load_wildtype, load_test_data, get_mutations

def compute_pseudo_likelihood(sequence, model, tokenizer, device, batch_size=8):
    """
    Computes the Pseudo-Log-Likelihood (PLL) of a sequence.
    PLL = Sum( log P(x_i | x_{/i}) ) for all i in sequence.
    
    Args:
        sequence: Amino acid string
        model: ESM model
        tokenizer: ESM tokenizer
        device: torch device
        batch_size: Batch size for masking (higher = faster but more VRAM)
    """
    tensor_seq = tokenizer(sequence, return_tensors="pt")["input_ids"].to(device)
    # tensor_seq shape: [1, seq_len + 2] (cls + seq + eos)
    seq_len = tensor_seq.shape[1] - 2 # Actual protein length
    
    # We only mask the protein residues, not CLS (0) or EOS (last)
    # Indices to mask: 1 to seq_len
    mask_indices = list(range(1, seq_len + 1))
    
    log_likelihoods = []
    
    # Process in batches to save memory
    for i in range(0, len(mask_indices), batch_size):
        batch_indices = mask_indices[i : i + batch_size]
        current_batch_size = len(batch_indices)
        
        # Create a batch of inputs
        # Shape: [current_batch_size, total_seq_len]
        batch_input = tensor_seq.repeat(current_batch_size, 1)
        
        # Apply masks
        for batch_idx, seq_idx in enumerate(batch_indices):
            batch_input[batch_idx, seq_idx] = tokenizer.mask_token_id
            
        with torch.no_grad():
            logits = model(batch_input).logits # [batch, seq_len, vocab]
            
        # Extract log_softmax probabilities
        log_probs = torch.log_softmax(logits, dim=-1)
        
        # Get the log-prob of the TRUE token at the masked position
        for batch_idx, seq_idx in enumerate(batch_indices):
            true_token = tensor_seq[0, seq_idx]
            token_log_prob = log_probs[batch_idx, seq_idx, true_token].item()
            log_likelihoods.append(token_log_prob)
            
    return sum(log_likelihoods)

def main():
    # Configuration
    MODEL_NAME = "facebook/esm2_t33_650M_UR50D" # Change to t36_3B if on Colab A100
    
    # Paths
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    data_dir = os.path.join(base_dir, "data")
    wt_file = os.path.join(data_dir, "pet-2025-wildtype-cds.csv")
    test_file = os.path.join(data_dir, "predictive-pet-zero-shot-test-2025.csv")
    output_file = os.path.join(base_dir, "submission_zero_shot.csv")
    
    # Device setup
    if torch.backends.mps.is_available():
        device = torch.device("mps")
        print("Using MPS (Mac GPU)")
    elif torch.cuda.is_available():
        device = torch.device("cuda")
        print("Using CUDA")
    else:
        device = torch.device("cpu")
        print("Using CPU")

    # Load Model
    print(f"Loading model: {MODEL_NAME}...")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmForMaskedLM.from_pretrained(MODEL_NAME).to(device)
    model.eval()

    # Load Data
    wt_seq = load_wildtype(wt_file)
    test_df = load_test_data(test_file)
    
    print(f"Wildtype Length: {len(wt_seq)}")
    print(f"Test Set Size: {len(test_df)}")

    # 1. Compute WT Score (Baseline)
    print("Computing Wildtype Baseline Score...")
    wt_score = compute_pseudo_likelihood(wt_seq, model, tokenizer, device, batch_size=4)
    print(f"WT Score: {wt_score:.4f}")

    # 2. Score Mutants
    results = []
    print("Scoring Mutants...")
    
    # Adjust batch_size based on your VRAM. 4-8 is safe for 650M on 16GB.
    # If on Colab with A100, you can increase to 32 or 64.
    BATCH_SIZE = 4 
    
    for idx, row in tqdm(test_df.iterrows(), total=len(test_df)):
        mut_seq = row['sequence']
        
        # Skip if invalid (NaN)
        if not isinstance(mut_seq, str):
            continue
            
        try:
            # Calculate Mutant Score
            mut_score = compute_pseudo_likelihood(mut_seq, model, tokenizer, device, batch_size=BATCH_SIZE)
            
            # LLR = Mutant_PLL - WT_PLL
            # Positive LLR means Mutant is more likely than WT (often correlated with stability/fitness)
            llr_score = mut_score - wt_score
            
            results.append({
                'sequence': mut_seq,
                'score': llr_score
            })
            
        except Exception as e:
            print(f"Error scoring sequence {idx}: {e}")
            results.append({
                'sequence': mut_seq,
                'score': -999.0 # Flag for failure
            })

    # Save Results
    results_df = pd.DataFrame(results)
    
    # Rank by score (descending)
    results_df = results_df.sort_values('score', ascending=False)
    
    results_df.to_csv(output_file, index=False)
    print(f"Scoring complete. Results saved to {output_file}")

if __name__ == "__main__":
    main()
