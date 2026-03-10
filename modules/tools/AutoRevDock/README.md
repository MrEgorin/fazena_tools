# rDock - An Open-Source Toolkit for Scalable Reverse Docking

An open-source Python toolkit that provides a comprehensive molecular docking workflow, automating the entire process from ligand preparation to binding affinity prediction. The workflow integrates two docking engines—AutoDock Vina and idock—and supports both the original Vina scoring function and an enhanced scoring scheme, Vina_SFCT. It is designed for high-throughput reverse docking with improved accuracy and scalability, making it particularly suitable for large-scale virtual screening tasks.

## 🔬 Overview

rDock is a Python-based molecular docking pipeline that streamlines:
- **Ligand Processing**: SMILES → MOL2 → PDBQT conversion
- **Protein Preparation**: PDB fixing and receptor preparation
- **Binding Site Prediction**: Automated pocket detection using P2rank
- **Molecular Docking**: Support for AutoDock Vina and idock
- **Advanced Scoring**: OnionNet-SFCT rescoring for improved accuracy

## ✨ Features

- 🚀 **Multi-core parallel processing** for high-throughput screening
- 🔄 **Flexible workflow modes**: Custom protein processing or prebuilt database usage
- 📊 **Multiple docking engines**: AutoDock Vina and idock support
- 🎯 **Automated binding site prediction** using P2rank
- 📈 **Advanced rescoring** with OnionNet-SFCT for better binding affinity prediction
- 🛠️ **Robust error handling** and progress tracking
- 📋 **Detailed logging** of all operations

## 🔧 Installation

### Prerequisites

```bash
# Required Python packages
pip install openbabel-wheel
pip install pdbfixer
pip install openmm
pip install biopython
pip install pandas
pip install tqdm

# For SFCT rescoring (separate conda environment recommended)
conda create -n sfct python=3.8
conda activate sfct
# Install OnionNet-SFCT dependencies
```

### External Tools Required

1. **AutoDock Tools** (for prepare_ligand and prepare_receptor)
2. **AutoDock Vina** (optional, for Vina docking)
3. **idock** (optional, for idock docking)
4. **P2rank** (for binding site prediction)
5. **OnionNet-SFCT** (for rescoring)

Place these tools in the `Tools/` directory with the following structure:
```
Tools/
├── autodock_vina/bin/
├── idock-2.2.1/bin/Linux/
├── OnionNet-SFCT/
└── external_tools/p2rank_2.4.2/
```

## 🚀 Quick Start

### Basic Usage

```bash
python rDock.py --ligand ligands.txt --protein proteins/ --mode custom \
                --box_x 20 --box_y 20 --box_z 20 --dock_tool vina
```

### Input File Formats

**Ligand file** (`ligands.txt`):
```
ligand1 CCO
ligand2 c1ccccc1
ligand3 CC(=O)OC1=CC=CC=C1C(=O)O
```

**Protein directory**: Contains PDB files for processing

## 📖 Detailed Usage

### Command Line Arguments

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--ligand` | ✓ | Path to ligand file (SMILES format) |
| `--protein` | ✓ | Protein folder or prebuilt database path |
| `--mode` | ✓ | Processing mode: `custom` or `prebuilt` |
| `--box_x/y/z` | ✓ | Docking box dimensions (Å) |
| `--cpu` | ✗ | Number of CPU cores (default: auto) |
| `--dock_tool` | ✗ | Docking software: `vina` or `idock` (default: idock) |

### Workflow Modes

#### Custom Mode
Process raw PDB files through the complete pipeline:
```bash
python rDock.py --ligand compounds.txt --protein raw_pdbs/ --mode custom \
                --box_x 25 --box_y 25 --box_z 25 --cpu 8
```

#### Prebuilt Mode
Use preprocessed protein database:
```bash
python rDock.py --ligand compounds.txt --protein processed_db/ --mode prebuilt \
                --box_x 20 --box_y 20 --box_z 20 --dock_tool vina
```

## 📁 Output Structure

```
project_directory/
├── ligands_mol2/           # MOL2 format ligands
├── ligands_pdbqt/          # PDBQT format ligands
├── proteins_fixed/         # Fixed PDB structures
├── proteins_pdbqt/         # PDBQT format receptors
├── proteins_pocket/        # Binding site predictions
├── vina_results/           # Docking results
└── vina_SFCT/             # Rescoring results
```

## 🔧 Individual Script Usage

### 1. Ligand Processing
```bash
python 1.Ligand_process.py ligands.txt
```

### 2. Protein Processing
```bash
# Fix PDB structures
python 21.Protein_process.py proteins/ fixed_proteins/ 8

# Prepare receptors
python 22.Proten_prepare.py fixed_proteins/ pdbqt_proteins/ 8
```

### 3. Binding Site Prediction
```bash
python 3.site-prediction.py fixed_proteins/ pocket_results/ 8
```

### 4. Molecular Docking
```bash
# AutoDock Vina
python 42.vina.py pdbqt_proteins/ pdbqt_ligands/ center.txt results/ 20 20 20 8

# idock
python 41.idock.py pdbqt_proteins/ pdbqt_ligands/ center.txt results/ 20 20 20 8
```

### 5. Rescoring
```bash
# Vina results rescoring
python 52.vina-Rescoring.py pdbqt_proteins/ vina_results/ sfct_output/ 8

# idock results rescoring
python 51.idock-Rescoring.py pdbqt_proteins/ idock_results/ sfct_output/ 8
```

## 📊 Results Interpretation

### Docking Results
- **PDBQT files**: 3D binding poses
- **Log files**: Binding energies and statistics
- **CSV files**: Combined scoring results

### Scoring Metrics
- **Vina Score**: AutoDock Vina binding affinity (kcal/mol)
- **idock Score**: idock binding energy (kcal/mol)
- **SFCT Score**: OnionNet-SFCT predicted binding affinity
- **Combined Score**: Average of primary and SFCT scores

## ⚠️ Troubleshooting

### Common Issues

1. **Missing dependencies**: Ensure all required tools are installed and in PATH
2. **Memory issues**: Reduce CPU cores for large datasets
3. **SFCT environment**: Make sure conda environment is properly configured
4. **File permissions**: Check read/write permissions for input/output directories

### Performance Tips

- Use SSD storage for better I/O performance
- Optimize CPU usage based on available memory
- Consider splitting large ligand libraries
- Use prebuilt mode for repeated screening

## 📝 Citation

If you use rDock in your research, please cite:

```bibtex
@software{rdock2024,
  title={rDock: Automated Molecular Docking Workflow},
  author={Luo, Qing},
  year={2024},
  institution={Macao Polytechnic University},
  email={p2214931@mpu.edu.mo}
}
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📧 Contact

- **Author**: Qing Luo
- **Email**: p2214931@mpu.edu.mo, sunluoq@163.com
- **Institution**: Macao Polytechnic University

## 🙏 Acknowledgments

- AutoDock Suite developers
- OpenEye Scientific Software
- P2rank development team
- OnionNet-SFCT authors
- OpenMM and PDBFixer developers

---

**Note**: This tool is designed for research purposes. Please ensure you have appropriate licenses for all third-party software components.
