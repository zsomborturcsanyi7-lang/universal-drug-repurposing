import subprocess
import os

def convert_pdb_to_pdbqt(input_pdb, output_pdbqt):
    """
    Converts a PDB file to PDBQT using Open Babel.
    Flags:
    -xr: preserve the original residue and atom names
    -c: add charges (important for Vina)
    """
    if not os.path.exists(input_pdb):
        print(f"Error: {input_pdb} not found.")
        return False

    cmd = [
        "obabel", input_pdb, 
        "-O", output_pdbqt, 
        "-xr", 
        "--partialcharge", "gasteiger"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"Successfully converted {input_pdb} to {output_pdbqt}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Conversion failed: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Error: Open Babel (obabel) not found in system path. Please ensure it is installed.")
        return False

if __name__ == "__main__":
    convert_pdb_to_pdbqt("receptor/receptor_cleaned.pdb", "receptor/receptor_cleaned.pdbqt")
