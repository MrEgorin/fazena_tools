import subprocess
import os
from typing import Dict, List

class AutoDocking:
    def __init__(self, vina_path: str):
        self.vina_path = vina_path

    def dock_protein(self, protein_file: str, ligand_file: str, center_x: float, center_y: float, center_z: float,
                    size_x: float, size_y: float, size_z: float, exhaustiveness: int) -> Dict[str, float]:
        """Docking for protein targets using AutoDock Vina GPU."""
        command = [
            self.vina_path,
            '--receptor', protein_file,
            '--ligand', ligand_file,
            '--center_x', str(center_x),
            '--center_y', str(center_y),
            '--center_z', str(center_z),
            '--size_x', str(size_x),
            '--size_y', str(size_y),
            '--size_z', str(size_z),
            '--exhaustiveness', str(exhaustiveness),
            '--gpu',
            '--out', 'output_protein_docking.pdbqt'
        ]
        try:
            subprocess.run(command, check=True)
            return self.parse_pdbqt('output_protein_docking.pdbqt')
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Docking failed: {e}")

    def dock_rna(self, protein_file: str, ligand_file: str, center_x: float, center_y: float, center_z: float,
                 size_x: float, size_y: float, size_z: float, exhaustiveness: int) -> Dict[str, float]:
        """Docking for RNA targets using standard AutoDock Vina."""
        command = [
            self.vina_path,
            '--receptor', protein_file,
            '--ligand', ligand_file,
            '--center_x', str(center_x),
            '--center_y', str(center_y),
            '--center_z', str(center_z),
            '--size_x', str(size_x),
            '--size_y', str(size_y),
            '--size_z', str(size_z),
            '--exhaustiveness', str(exhaustiveness),
            '--out', 'output_rna_docking.pdbqt'
        ]
        try:
            subprocess.run(command, check=True)
            return self.parse_pdbqt('output_rna_docking.pdbqt')
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Docking failed: {e}")

    def parse_pdbqt(self, pdbqt_file: str) -> Dict[str, float]:
        """Parse the output PDBQT file and extract binding energies for the top 5 hit compounds."""
        with open(pdbqt_file, 'r') as f:
            lines = f.readlines()

        results = []
        in_results_section = False
        for line in lines:
            if line.startswith('REMARK VINA RESULT:'):
                in_results_section = True
            elif line.startswith('ENDMDL'):
                in_results_section = False
            elif in_results_section:
                parts = line.split()
                if len(parts) >= 6 and parts[0] == 'MODEL':
                    binding_energy = float(parts[5])
                    results.append(binding_energy)

        # Sort by binding energy (lower is better)
        results.sort()

        # Extract top 5 binding energies
        top_energies = {f'pose_{i+1}': energy for i, energy in enumerate(results[:5])}
        return top_energies

    def auto_configure_grid(self, protein_file: str) -> Dict[str, float]:
        """Automate grid-box configuration generation."""
        # Placeholder logic for grid configuration
        # For demonstration purposes, we'll use hardcoded values
        center_x = 0.0
        center_y = 0.0
        center_z = 0.0
        size_x = 10.0
        size_y = 10.0
        size_z = 10.0
        exhaustiveness = 8

        return {
            'center_x': center_x,
            'center_y': center_y,
            'center_z': center_z,
            'size_x': size_x,
            'size_y': size_y,
            'size_z': size_z,
            'exhaustiveness': exhaustiveness
        }

    def dock(self, protein_file: str, ligand_file: str, target_type: str) -> Dict[str, float]:
        """Automated docking based on target type."""
        grid_config = self.auto_configure_grid(protein_file)
        if target_type == 'protein':
            return self.dock_protein(protein_file, ligand_file, **grid_config)
        elif target_type == 'rna':
            return self.dock_rna(protein_file, ligand_file, **grid_config)
        else:
            raise ValueError("Invalid target type. Supported types: 'protein', 'rna'")
