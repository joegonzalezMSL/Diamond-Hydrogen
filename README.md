# Diamond-Hydrogen
Utility to create diamond - hydrogen interfaces for VASP simulations

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

## Example Usage
python3 diamond_hydrogen.py trunk/cd-111.vasp -o poscar-111-h2 -r 3 4 3 --extend 10 -s 111
