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

def compute_pseudo_likelihood(sequence, model, tokenizer, device, batch_size=64):
    """
    Computes the Pseudo-Log-Likelihood (PLL) of a sequence efficiently.
    PLL = Sum( log P(x_i | x_{/i}) ) for all i in sequence.
    
    Optimized for batch processing and GPU throughput.
    
    Args:
        sequence: Amino acid string
        model: ESM model
        tokenizer: ESM tokenizer
        device: torch device
        batch_size: Batch size for inference (set high for A100, e.g., 128+)
    """
    tensor_seq = tokenizer(sequence, return_tensors="pt")["input_ids"].to(device)
    # tensor_seq shape: [1, seq_len + 2] (cls + seq + eos)
    seq_len = tensor_seq.shape[1] - 2 # Actual protein length
    
    # Create inputs: repeat sequence (seq_len) times
    # Shape: [seq_len, total_len]
    # We use repeat instead of expand to ensure memory is allocated for masking
    inputs = tensor_seq.repeat(seq_len, 1)
    
    # Create diagonal mask indices (1 to seq_len)
    # Row i masks position i+1
    mask_indices = torch.arange(1, seq_len + 1, device=device)
    
    # Apply mask: inputs[i, i+1] = mask_token
    inputs[torch.arange(seq_len), mask_indices] = tokenizer.mask_token_id
    
    total_log_likelihood = 0.0
    
    # Process in batches
    for i in range(0, seq_len, batch_size):
        batch_input = inputs[i : i + batch_size]
        current_batch_indices = mask_indices[i : i + batch_size]
        
        with torch.no_grad():
            # Forward pass
            logits = model(batch_input).logits # [batch, seq_len+2, vocab]
            
        # We only care about the logits at the masked positions
        # Gather logits: [batch, vocab] at the specific masked indices
        # logits[row, mask_idx, :]
        masked_logits = logits[
            torch.arange(len(batch_input), device=device), 
            current_batch_indices, 
            :
        ]
        
        # Compute log_softmax over vocab
        log_probs = torch.log_softmax(masked_logits, dim=-1)
        
        # Get the log-prob of the TRUE token
        true_tokens = tensor_seq[0, current_batch_indices]
        
        # Gather the log prob of the true token
        token_log_probs = log_probs[
            torch.arange(len(batch_input), device=device), 
            true_tokens
        ]
        
        total_log_likelihood += token_log_probs.sum().item()
            
    return total_log_likelihood

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
    
    # Enable FP16/BF16 for speed on CUDA
    if device.type == "cuda":
        model.half() # Use float16
        print("Model converted to float16 for speed.")
        
    model.eval()

    # Load Data
    wt_seq = load_wildtype(wt_file)
    test_df = load_test_data(test_file)
    
    print(f"Wildtype Length: {len(wt_seq)}")
    print(f"Test Set Size: {len(test_df)}")

    # 1. Compute WT Score (Baseline)
    print("Computing Wildtype Baseline Score...")
    # Increase batch size for WT (single sequence)
    wt_score = compute_pseudo_likelihood(wt_seq, model, tokenizer, device, batch_size=128)
    print(f"WT Score: {wt_score:.4f}")

    # 2. Score Mutants
    results = []
    print("Scoring Mutants...")
    
    # Adjust batch_size based on your VRAM. 
    # For A100 (40GB/80GB), you can easily use 128, 256, or even full sequence length (e.g. 512).
    # This allows computing the entire PLL in 1 or 2 forward passes.
    BATCH_SIZE = 128 
    
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
