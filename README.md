# Diamond-Hydrogen
Utility to create diamond - hydrogen interfaces for VASP simulations

## Contents
```
diamond-hydrogen/
├── pyproject.toml
├── README.md
└── src/
    └── diamond_hydrogen/
        ├── __init__.py
        ├── vasp_tools.py
        └── diamond_hydrogen.py
└── trunk/
    └── cd-100.vasp
    └── cd-110.vasp
    └── cd-111.vasp
```

The `trunk` directory contains the starting unit cells, they are single unitcells, so you will need to enlarge them to an appropriate size, see below in [Examples] (#Example-Usage)

## Installation
Get the package and place in a central location
`git clone https://github.com/joegonzalezMSL/Diamond-Hydrogen.git`
```bash
cd Diamond-Hydrogen
pip install -e .
```
The `-e` flag is included so that if you make changes to the source code, they are pushed to the location without the need to re-install everything.

## Verify where the installation went
```bash
pip show diamond-hydrogen
```

This command should show an output like this:
```bash
Name: diamond-hydrogen
Version: 0.1.0
Location: /path/to/.venv/lib/python3.11/site-packages
Entry-points:
  diamond-hydrogen = diamond_hydrogen.diamond_hydrogen:main
```

Test the installation with
```bash
diamond-hydrogen -h
```

You should see the help menu.


## Example Usage

```
# build a structure from the 100 unit cell, keeping all defaults. Only useful if you have a pre-build lattice
diamond-hydrogen trunk/cd-100.vasp

# Replicate unit cell and add vacuum region
diamond-hydrogen --extend 10 --replicate 3 3 4 trunk/cd-111.vasp -o POSCAR_334_H2

# Specify surface type and H2 gas density
diamond-hydrogen trunk/cd-110.vasp -o poscar-110-h2 -r 3 4 3 --extend 10 -s 111 -d 0.3
```
