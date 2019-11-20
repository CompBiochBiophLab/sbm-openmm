The models class of sbmOpenMM contains three methods for automatic setting up predefined SBM flavors. It works by initializing a system class with the necessary force field parameters, derived from the input files, to set up default models that are detailed next:

- Coarse grained, alpha-carbon (CA), model

The coarse grained method represents the protein system as beads centered at the alpha carbons of each residue in the protein. It uses harmonic potentials to hold the covalent connectivity, geometry and chirality of the beads. Torsional geometries are modeled with a periodic torsion potential. Native contacts are represented through the use of Lennard-Jones potentials that allow to form and break the connections, permitting complete and local unfolding of the structures.

The force field equation is:

The method to create a CA model is to call:

sbmOpenMM.models.getCAModel(input_pdb_file, input_contacts_file)

- All-heavy-atoms (AA) model

The all-atom model represents the protein system with all its heavy atoms (i.e. excluding hydrogens). It uses harmonic potentials to hold the covalent connectivity, geometry and chirality of the protein residues. Periodic torsional potentials are used to maintain dihedral geometries of backbones and side chains. Native contacts are represented through the use of Lennard-Jones potentials that allow to form and break the connections, permitting complete and local unfolding of the structures.

The force field equation is:

- Multi basin model
