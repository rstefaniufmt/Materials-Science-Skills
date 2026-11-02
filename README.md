# Materials Science Skills 🔬💻

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Contributions Welcome](https://img.shields.io/badge/Contributions-Welcome-brightgreen.svg)](#contributing)
[![Maintainer](https://img.shields.io/badge/Maintainer-rstefaniufmt-orange.svg)](https://github.com/rstefaniufmt)

A curated collection of modular workflows, computational scripts, and AI agent skills designed for **Materials Science**, **Computational Chemistry**, and **Molecular Modeling** research.

---

## 📖 Overview

**Materials-Science-Skills** provides reusable tools, prompt/skill definitions, and automation scripts to streamline common tasks in materials discovery and atomistic simulations. Whether you are setting up Molecular Dynamics (MD) topologies, parsing Density Functional Theory (DFT) outputs, or analyzing structural properties, this repository centralizes best practices and ready-to-use workflows.

## ✨ Key Capabilities & Modules

* **Atomistic & Molecular Dynamics (MD):** System preparation, force field parameterization, solvation/casting workflows, and trajectory post-processing (e.g., GROMACS, LAMMPS, AmberTools).
* **Quantum Chemistry & Electronic Structure:** Input generation, geometry optimization workflows, and output parsing for DFT and semi-empirical packages (e.g., ORCA, Gaussian, Quantum ESPRESSO).
* **Structural & Property Analysis:** Calculation of Radial Distribution Functions ($g(r)$), Mean Squared Displacement ($\text{MSD}$), Root Mean Square Deviation ($\text{RMSD}$), hydrogen bonding networks, and thermodynamic properties.
* **Data Processing & Cheminformatics:** Automated extraction of molecular descriptors, SMILES/PDB/MOL2 file conversions, and publication-ready visualization scripts.
* **AI Agent Skills (Optional/Modular):** Structured instructions (`SKILL.md` / prompts) enabling LLM-based coding assistants to execute domain-specific materials science pipelines accurately.

---

## 📂 Repository Structure

```text
Materials-Science-Skills/
├── skills/                  # Domain-specific skills, prompts, and workflow definitions
├── scripts/                 # Standalone Python and Bash automation scripts
├── templates/               # Input templates for MD/DFT simulation engines
├── examples/                # Sample datasets and benchmark tutorials
└── README.md                # Project documentation
```
*(Note: Adjust the directory tree above to match your exact local folders).*

---

## 🚀 Getting Started

### Prerequisites

Depending on the specific skill or workflow you intend to run, ensure you have the following installed in your environment:

* **Python** ($\ge 3.9$) with scientific libraries (`numpy`, `scipy`, `matplotlib`, `pandas`, `mdanalysis`, `rdkit`, `ase`)
* Standard command-line utilities (`bash`, `awk`, `sed`)
* Relevant simulation engines (e.g., GROMACS, LAMMPS, ORCA, Packmol) if executing full simulation pipelines.

### Installation

1. Clone this repository to your local machine or HPC cluster:
   ```bash
   git clone https://github.com/rstefaniufmt/Materials-Science-Skills.git
   cd Materials-Science-Skills
   ```

2. *(Optional)* If using this repository as a skill library for AI coding assistants (such as Claude Code, Cursor, or Gemini CLI), copy or symlink the target skill directories into your workspace configuration folder.

---

## 🛠️ Usage

1. Browse the repository directories to locate the specific workflow or skill you need.
2. Check the local documentation or header comments inside each subdirectory/script for required parameters and input formats.
3. Run the scripts directly or invoke the skill within your computational workflow.

---

## 🤝 Contributing

Contributions, bug reports, and new materials science workflows are welcome!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/NewMaterialSkill`)
3. Commit your Changes (`git commit -m 'Add new workflow for property calculation'`)
4. Push to the Branch (`git push origin feature/NewMaterialSkill`)
5. Open a Pull Request

---

## 📄 License

Distributed under the **MIT License**. See `LICENSE` for more information.

---

## 📬 Contact & Affiliation

**Prof. Dr. Ricardo Stefani**  
Federal University of Ouro Preto (UFOP) — Brazil  
GitHub: [@rstefaniufmt](https://github.com/rstefaniufmt)
