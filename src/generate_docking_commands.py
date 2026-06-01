import os

def create_config_file(config_path="config.txt"):
    """
    Creates a configuration file for AutoDock Vina.
    """
    config_content = """
receptor = receptor/receptor_cleaned.pdbqt
center_x = -10.7
center_y = 12.4
center_z = 68.8
size_x = 20.0
size_y = 20.0
size_z = 20.0
exhaustiveness = 8
cpu = 0
"""
    with open(config_path, "w") as f:
        f.write(config_content)
    print(f"Configuration file generated: {config_path}")

def generate_batch_script(ligands_dir="ligands", output_dir="docking_results", batch_script="run_batch_docking.bat"):
    """
    Generates a Windows batch script to run Vina for all ligands.
    """
    if not os.path.exists(ligands_dir):
        print(f"Warning: '{ligands_dir}' directory not found. Skipping batch script generation.")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    with open(batch_script, "w") as f:
        f.write("@echo off\n")
        f.write(f"mkdir {output_dir}\n")
        
        ligand_files = [f for f in os.listdir(ligands_dir) if f.endswith(".pdbqt")]
        total_ligands = len(ligand_files)
        
        for i, ligand_file in enumerate(ligand_files, 1):
            ligand_name = os.path.splitext(ligand_file)[0]
            f.write(f"echo [Progress] Docking {ligand_name} ({i}/{total_ligands})...\n")
            f.write(f"vina_1.2.7_win.exe --config config.txt --ligand {ligands_dir}/{ligand_file} --out {output_dir}/{ligand_name}_out.pdbqt > {output_dir}/{ligand_name}_log.txt 2>&1\n")
                
        f.write("echo [Progress] All docking tasks completed.\n")
    
    print(f"Batch script generated: {batch_script}")

if __name__ == "__main__":
    create_config_file()
    generate_batch_script()
