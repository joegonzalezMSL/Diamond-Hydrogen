from . import vasp_tools
import sys
import numpy as np
import argparse
import textwrap

##constants
HH_BOND = 0.74
H2_MIN_MOL_DISTANCE = 2.0
H2_MOLAR_MASS = 2.016
AVOGADRO = 6.022e23


def get_arguments(argv):

	cli = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
           Utility to add Hydrogen to the surface of a structure, create space 
           for adding H2 gas. Used in creating structures that resemble the 
           interface of diamond with DT Ice used in ICF experiments
         -----------------------------------------------------------------
         '''),
    epilog=textwrap.dedent('''\
         examples:
         -----------------------------------------------------------------
            %(prog)s POSCAR
            %(prog)s --extend=10 --replicate 3 3 4 POSCAR -o poscar-334-h2
            %(prog)s POSCAR -d 0.3 -r 3 3 4 -g 1.2  -s 110 -e 7
         '''))

	cli.add_argument(
		"POSCAR",
		type=str,
		help="starting POSCAR file to which hydrogen will be added")
	cli.add_argument(
		"-r","--replicate",
		dest="reps",
		help="replications to apply to the input POSCAR. (default: 1 1 1)",
		nargs=3,
		type=int,
		default=[1,1,1])
	cli.add_argument(
		"-d","--density",
		dest="H2_DENSITY",
		help="density of H2 gas. (default: 0.19 g/cc)",
		default=0.19,
		type=float)
	cli.add_argument(
		"-g","--gap",
		dest="Z_GAP",
		help="Exclusion region for H2 molecules relative to slab surfaces. (default: 2.5 A)",
		type=float,
		default=2.5)
	cli.add_argument(
		"-b","--bond=",
		dest="CH_BOND",
		type=float,
		help="length of C-H bond for surface terminated atoms. (default: 1.1 A)",
		default=1.1)
	cli.add_argument(
		"-o","--output",
		dest="outname",
		type=str,
		help="Name of file containing the new structure. (default: POSCAR_with_H2) ",
		default="POSCAR_with_H2")
	cli.add_argument(
		"-e","--extend=",
		dest="extendZ",
		default=0.0,
		type=float,
		help="Amount to extend Z-lattice vector to contain H2 gas. (default: 0.0)")
	cli.add_argument(
		"-s","--surface",
		dest="surface",
		help="Define the surface that is oriented along Z, only used for the comment line",
		default="111",
		type=str)
	cli.add_argument(
		'--version',
		action='version',
		version='%(prog)s 4.0.')
	args = cli.parse_args()

	return args


def main():

	args = get_arguments(sys.argv[1:])

	CH_BOND = args.CH_BOND
	H2_DENSITY = args.H2_DENSITY
	Z_GAP = args.Z_GAP
	outname = args.outname
	extendZ = args.extendZ
	poscar = args.POSCAR
	reps = args.reps
	surface = args.surface

	


	lattice, positions, species, num_atoms, numbers = vasp_tools.read_vasp(poscar)

	##replicate unitcell if requested
	if reps != [1, 1, 1]:
	    lattice, positions, num_atoms, numbers = vasp_tools.replicate_unitcell(
	        lattice, positions, num_atoms, numbers, reps
	    )


	##extend z-lattice
	lattice[2][2] += extendZ

	nC = num_atoms[0]

	positions = np.array(positions, dtype=float)
	species = list(species)
	num_atoms = list(num_atoms)
	num_atoms = num_atoms + [0] * (len(species) - len(num_atoms))


	a = float(lattice[0][0])
	b = float(lattice[1][1])
	c = float(lattice[2][2])

	print(f"List of parameters for this run:")
	print(f"Diamond surface:             [{surface}]")
	print(f"C-H Bond length:             {CH_BOND:.2f} A")
	print(f"H-H Bond length:             {HH_BOND:.2f} A")
	print(f"H2 molecule distance:        {H2_MIN_MOL_DISTANCE:.2f} A")
	print(f"H2 gas density:              {H2_DENSITY:.2f} g/cc")
	if extendZ > 0.0:
		print(f"Z-lattice increased by:      {extendZ} A")
	if sum(reps) != 3:
		print(f"Replicating original cell:\n    {reps[0]} x {reps[1]} x {reps[2]}")
	print(f"Unit cell dimensions:        {a:.2f} x {b:.2f} x {c:.2f}")
	print(f"File name for new structure: {outname}")

	print(f"======================================================\n")


	## find top/bottom surface C atoms
	z_coords = positions[:, 2]
	max_z = float(np.max(z_coords))
	min_z = float(np.min(z_coords))
	top_atoms = np.where(np.isclose(z_coords, max_z))[0]
	bottom_atoms = np.where(np.isclose(z_coords, min_z))[0]

	print(f"Atoms added to each surface: {len(bottom_atoms)} ")

	##add H terminations (top and bottom)
	added_H = []
	for idx in top_atoms:
	    x, y, z = positions[idx]
	    added_H.append([float(x), float(y), float(z + CH_BOND)])
	## bottom termination (wrap to top of cell so PBC connects)
	for idx in bottom_atoms:
	    x, y, z = positions[idx]
	    added_H.append([float(x), float(y), float(z - CH_BOND + c)])


	positions = np.vstack((positions, np.array(added_H, dtype=float)))


	species.append("H")
	num_atoms.append(0)

	h_index = species.index("H")
	num_atoms[h_index] += len(added_H)

	##compute how many H2 molecules based on requested density
	rho_A3 = H2_DENSITY * 1e-24  # g/Å^3
	n_A3 = rho_A3 * AVOGADRO / H2_MOLAR_MASS  # molecules per Å^3

	##define gas region limits
	z_start = max_z + Z_GAP
	z_end = c - Z_GAP

	if z_start >= z_end:
	    raise ValueError(
	        f"No room for gas region: z_start={z_start:.2f}, z_end={z_end:.2f}, cell height={c:.2f}"
	    )

	V_top = a * b * (z_end - z_start)
	N_H2 = int(round(n_A3 * V_top))

	##random placement with exclusion rule
	rng = np.random.default_rng()
	molecule_centers = []
	max_attempts = max(1000, N_H2 * 500)

	for attempt in range(max_attempts):
	    if len(molecule_centers) >= N_H2:
	        break

	    cand = np.array([
	        rng.uniform(0.0, a),
	        rng.uniform(0.0, b),
	        rng.uniform(z_start, z_end)
	    ], dtype=float)

	    ##must be above top H layer
	    if cand[2] < (max_z + CH_BOND ):
	        continue

	    ##ensure no overlap with existing molecules
	    if all(np.linalg.norm(cand - c2) >= H2_MIN_MOL_DISTANCE for c2 in molecule_centers):
	        molecule_centers.append(cand)


	##build H2 atoms from centers
	added_H2 = []
	for center in molecule_centers:
	    ##random orientation unit vector
	    v = rng.normal(size=3)
	    v /= np.linalg.norm(v)
	    half = (HH_BOND / 2.0) * v
	    h1 = center - half
	    h2 = center + half
	    added_H2.append(h1)
	    added_H2.append(h2)

	if len(added_H2) > 0:
	    positions = np.vstack((positions, np.array(added_H2, dtype=float)))

	nH = len(added_H) + len(added_H2) 
	print(f"H2 molecules added:          {len(molecule_centers)}")
	print(f"New system size:             {nC+nH} -> {nC} C | {nH} H")
	print(f"")

	h_index = species.index("H")
	num_atoms[h_index] += len(added_H2)

	ntypes = [nC,nH]
	tag = surface + " surface terminated with H and " +str(len(added_H2))+" H2 molecules"

	vasp_tools.write_poscar( outname,lattice,ntypes,species,positions,tag )
	print(f"Job complete! ")


if __name__ == "__main__": 
	main()
