The models class of sbmOpenMM contains three methods for automatic setting up predefined SBM potentials. It works by initializing a system class with the necessary force field parameters, derived from the input files, to set up one of the possible models which are detailed next:

- Coarse grained, alpha-carbon (CA), model

The coarse grained method represents the protein system as beads centered at the alpha carbons of each residue in the protein. It uses harmonic potentials to hold the covalent connectivity and geometry of the beads. Torsional geometries are modeled with a periodic torsion potential. Native contacts are represented through the use of Lennard-Jones potentials that allow to form and break non-bonded interactions, permitting complete and local unfolding of the structures.

The force field equations are:

.. math::
	H_A = \sum_{bonds}V_{bond}+\sum_{angles}V_{angle}+\sum_{torsions}V_{torsion}+\sum_{contacts}V_{LJ_{12-10}}+\sum_{non-contacts}V_{LJ_{12}}

.. math::
        V_{bond} = \frac{k_b}{2}(r-r_0)^2

.. math::
        V_{angle} = \frac{k_a}{2}(\theta-\theta_0)^2

.. math::
        V_{torsion} = k_t(1-cos(\phi-\phi_0))+\frac{1}{2}(1-cos(3(\phi-\phi_0))))

.. math::
        V_{LJ_{12-10}} = \epsilon_{c}(5(\frac{\sigma_{ij}}{r})^{12}-6(\frac{\sigma_{ij}}{r})^{10})

.. math::
        V_{LJ_{12}} = \epsilon_{nc}(\frac{\sigma_{ij}}{r})^{12}


Here the default values are :math:`k_b=20000` kJ/(mol mm) , :math:`k_a=40` kJ/(mol rad), :math:`k_t=1.0` kJ/mol , :math:`\epsilon_{c}=1.0` kJ/mol and :math:`\epsilon_{nc}=1.0` kJ/mol. The geometric parameters are set to the calculated structural value in the input structure, with :math:`r_0` the equilibrium bond distance in nanometers, :math:`\theta_0` the equilibrium angle length in radians, :math:`\phi_0` the equilibrium torsional angle in radians and :math:`\sigma_{ij}` the equilibrium contact distance in nanometers. The variable :math:`r` represents, accordingly, the current bond or contact distance in nanometers, :math:`\theta` the current angle length in radians and :math:`\phi` the current torsional angle in radians.

To create a CA model, call:

sbmOpenMM.models.getCAModel(pdb_file, contacts_file)

Here, pdb_file is the path to the PDB format structure of the protein and  contacts_file is the path to the contact file containing only the CA atoms of the system. This last file should be numbered considering the CA atoms consecutively.

- All-heavy-atoms (AA) model

The all-atom model represents the protein system with all its heavy atoms (i.e. excluding hydrogens). It uses harmonic potentials to hold the covalent connectivity, geometry and chirality of the protein residues. Periodic torsional potentials are used to maintain dihedral geometries of backbones and side chains. Native contacts are represented through the use of Lennard-Jones potentials that allow to form and break non-bonded interactions, permitting complete and local unfolding of the structures.

The force field equation is:

The method to create an AA model is:

sbmOpenMM.models.getAllAtomModel(pdb_file, contacts_file)

Here, pdb_file is the path to the PDB format structure of the protein and  contacts_file is the path to the contact file containing only the non-hydorgen atoms of the protein system.

- Multi basin model

The multi basin model automates the creation of a dual basin native contact potential. It receives as input two sbmOpenMM system classes, either CA or AA models, containing two different definitions of native contacts. One of the configurations is defined as the main model and the other is considered as the alternate model. All forcefield parameters different than the native contacts are copied from the main configuration into the multi basin model. Then, the contacts are compared between the input configurations to define the sets of common and unique contacts. Common contacts with equilibrium length distances that differ more than a threshold are defined as dual basin and are assigned a special non-bonded Gaussian potential. The rest of the contacts are treated as single minima and are modeled with a Lennard-Jones (default) or a single basin Gaussian potential. 
The method to create a multi basin model is:

sbmOpenMM.models.getMultiBasinModel(main_model, alternate_configuration=alternate_model)

Here, main_model and alternate_model are initialized sbmOpenMM system classes containing full force field parameter definitions.
