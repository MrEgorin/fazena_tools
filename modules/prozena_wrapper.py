import subprocess
import os
import tempfile
import pandas as pd
import numpy as np

def run_bind_real(protein_seq, ligands, device="cpu"):
    """
    Executes the real BIND inference via subprocess.
    Expects BIND to be in modules/tools/BIND/
    """
    bind_dir = os.path.join(os.getcwd(), "modules", "tools", "BIND")
    bind_script = os.path.join(bind_dir, "bind.py")
    
    # Check if BIND exists
    if not os.path.exists(bind_script):
        return None, "BIND engine not found in modules/tools/BIND"

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create FASTA
        fasta_path = os.path.join(tmpdir, "target.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">Target\n{protein_seq}\n")
        
        # Create SMILES file
        smiles_path = os.path.join(tmpdir, "library.smi")
        with open(smiles_path, "w") as f:
            for s in ligands:
                f.write(f"{s}\n")
        
        output_path = os.path.join(tmpdir, "results.tsv")
        
        # Execute BIND
        # Usage: python bind.py --proteins <fasta> --ligands <smi> --output <tsv> --bind_device <device>
        cmd = [
            "python", bind_script,
            "--proteins", fasta_path,
            "--ligands", smiles_path,
            "--output", output_path,
            "--bind_device", device,
            "--esm_device", device,
            "--batch_size", "1"
        ]
        
        try:
            # Set PYTHONPATH to include BIND dir so internal imports work
            env = os.environ.copy()
            env["PYTHONPATH"] = bind_dir + os.pathsep + env.get("PYTHONPATH", "")
            
            result = subprocess.run(cmd, env=env, capture_output=True, text=True, check=True)
            
            if os.path.exists(output_path):
                # BIND outputs TSV with columns: Input protein, Input SMILES, pKi, pIC50, pKd, pEC50, Logit, Non-binder probability
                df = pd.read_csv(output_path, sep="\t")
                # Rename columns for internal consistency
                df = df.rename(columns={
                    "Input SMILES": "SMILES",
                    "pKi": "pKi",
                    "pKd": "pKd",
                    "Non-binder probability": "Non-Binder Prob"
                })
                # Ensure numeric
                for col in ["pKi", "pKd", "Logit", "Non-Binder Prob"]:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                
                return df.sort_values(by="pKi", ascending=False).reset_index(drop=True), None
            else:
                return None, f"BIND failed to produce output. Error: {result.stderr}"
                
        except subprocess.CalledProcessError as e:
            return None, f"BIND execution failed: {e.stderr}"
        except Exception as e:
            return None, f"Unexpected error: {str(e)}"
