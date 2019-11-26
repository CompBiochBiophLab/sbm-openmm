The models class of sbmOpenMM contains three methods for automatic setting up predefined SBM potentials. It works by initializing a system class with the necessary force field parameters, derived from the input files, to set up one of the possible models which are detailed next:

Coarse grained, alpha-carbon (CA), model
++++++++++++++++++++++++++++++++++++++++

The coarse grained method represents the protein system as beads centered at the alpha carbons of each residue in the protein. It uses harmonic potentials to hold the covalent connectivity and geometry of the beads. Torsional geometries are modeled with a periodic torsion potential. Native contacts are represented through the use of Lennard-Jones potentials that allow to form and break non-bonded interactions, permitting complete and local unfolding of the structures.

To create a CA model, call:

sbmOpenMM.models.getCAModel(pdb_file, contacts_file)

Here, pdb_file is the path to the PDB format structure of the protein and  contacts_file is the path to the contact file containing only the CA atoms of the system. This last file should be numbered considering the CA atoms consecutively.

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
        V_{LJ_{12}} = \epsilon_{nc}(\frac{\sigma_{nc}}{r})^{12}


Here the default values are :math:`k_b=20000\ kJ/(mol \cdot nm^2)`, :math:`k_a=40\ kJ/(mol \cdot rad^2)`, :math:`k_t=1.0\ kJ/mol`, :math:`\epsilon_{c}=1.0\ kJ/mol`, :math:`\epsilon_{nc}=1.0\ kJ/mol` and :math:`\sigma_{nc}=0.4\ nm`. The geometric parameters are set to the calculated structural values in the input structure, with :math:`r_0` the equilibrium bond distance in nanometers, :math:`\theta_0` the equilibrium angle length in radians, :math:`\phi_0` the equilibrium torsional angle in radians and :math:`\sigma_{ij}` the equilibrium contact distance in nanometers. The variable :math:`r` represents, accordingly, the current bond or (non)contact distance in nanometers, :math:`\theta` the current angle length in radians and :math:`\phi` the current torsional angle in radians.

It is possible to use a :math:`V_{LJ_{12-10-6}}` potential for the native contact interactions, defined as:

.. math::
        V_{LJ_{12-10-6}} = \epsilon_{c}(13(\frac{\sigma_{ij}}{r})^{12}-18(\frac{\sigma_{ij}}{r})^{10}+4(\frac{\sigma_{ij}}{r})^{6})

This potential gives a small energy barrier for contact formation/breaking that emulates a "desolvation effect". To use this potential as the native contact energy function, instead of the :math:`V_{LJ_{12-10}}` potential, give the option contact_force ='12-10-6' to the sbmOpenMM.models.getCAModel() method. 
 
Note that even if the units for the force constants are given in real physical units (e.g. :math:`kJ/mol`), this is just to match the variables used by OpenMM. The models are not parametrized to equate this real physical values and comparison with experiments will require further adjustment to the energy unit system. 

All-heavy-atoms (AA) model
++++++++++++++++++++++++++

The all-atom model represents the protein system with all its heavy atoms (i.e. excluding hydrogens). It uses harmonic potentials to hold the covalent connectivity, geometry and chirality of the protein residues. Periodic torsional potentials are used to maintain dihedral geometries of backbones and side chains. Native contacts are represented through the use of Lennard-Jones potentials that allow to form and break non-bonded interactions, permitting complete and local unfolding of the structures.

The method to create an AA model is:

sbmOpenMM.models.getAllAtomModel(pdb_file, contacts_file)

Here, pdb_file is the path to the PDB format structure of the protein and  contacts_file is the path to the contact file containing only the non-hydorgen atoms of the protein system.

The force field equations are:

.. math::
        H_A = \sum_{bonds}V_{bond}+\sum_{angles}V_{angle}+\sum_{torsions}V_{torsion}+\sum_{impropers}V_{improper}+\sum_{planars}V_{planar}+\sum_{contacts}V_{LJ_{12-10}}+\sum_{non-contacts}V_{LJ_{12}} 

.. math::
        V_{bond} = \frac{k_b}{2}(r-r_0)^2

.. math::
        V_{angle} = \frac{k_a}{2}(\theta-\theta_0)^2

.. math::
        V_{torsion} = k_t(1-cos(\phi-\phi_0))+\frac{1}{2}(1-cos(3(\phi-\phi_0))))

.. math::
        V_{improper} = \frac{k_i}{2}(\chi-\chi_{0})^2

.. math::
        V_{planar} = \frac{k_p}{2}(\chi-\chi_{0})^2

.. math::
        V_{LJ_{12-10}} = \epsilon_{c}(5(\frac{\sigma_{ij}}{r})^{12}-6(\frac{\sigma_{ij}}{r})^{10})

.. math::
        V_{LJ_{12}} = \epsilon_{nc}(\frac{\sigma_{nc}}{r})^{12}

Here the default values are :math:`k_b=10000\ kJ/(mol \cdot nm^2)`, :math:`k_a=80\ kJ/(mol \cdot rad^2)`, :math:`k_i=10.0\ kJ/(mol \cdot rad^2)`, :math:`k_p=20.0\ kJ/(mol \cdot rad^2)`, :math:`\epsilon_{nc}=0.1\ kJ/mol` and :math:`\sigma_{nc}=0.25\ nm`. The values of the torsional :math:`k_t` and native energy constant :math:`\epsilon_{c}` are assigned by the following equations:

.. math::
        k_t=N_{atoms}/3N_{torsions}\ (kJ/mol)
.. math::
        k_c=2N_{atoms}/3N_{contacts}\ (kJ/mol)


Where :math:`N_{atoms}` is the total number of atoms in the system, :math:`N_{torsions}` is the total number of proper torsions assigned by the forcefield and :math:`N_{contacts}` is the number of native contacts in the contact file definition. Additionally, the torsional energy constant :math:`k_t` is further divided by classifying the torsions into backbone and sidechain groups. The assignment is carried out as:

.. math::
        k_{t}^{bb}=2k_t/3
.. math::
        k_{t}^{sc}=k_t/3

Here, :math:`k_{t}^{bb}` and :math:`k_{t}^{sc}` are the torsional energy constant for backbone and sidechain torsion groups, respectively. This grouping of torsions into backbone and side chains is the default behaviour of the sbmOpenMM.models.getAllAtomModel() method. It can be disabled by given the option group_by_bb_and_sc=False.

The geometric parameters are set to the calculated structural values in the input structure, with :math:`r_0` the equilibrium bond distance in nanometers, :math:`\theta_0` the equilibrium angle length in radians, :math:`\phi_0` the equilibrium torsional angle in radians, :math:`\chi_0` the equilibrium improper or planar equilibrium angle in radians and :math:`\sigma_{ij}` the equilibrium contact distance in nanometers. The variable :math:`r` represents, accordingly, the current bond or (non)contact distance in nanometers, :math:`\theta` the current angle length in radians, :math:`\phi` the current proper torsional angle in radians and :math:`\chi` the equilibrium improper or planar torsional angles in radians.

Note that even if the units for the force constants are given in real physical units (e.g. :math:`kJ/mol`), this is just to match the variables used by OpenMM. The models are not parametrized to equate this real physical values and comparison with experiments will require further adjustment to the energy unit system.

Multi basin model
+++++++++++++++++

The multi basin model automates the creation of a dual basin native contact potential. It receives as input two sbmOpenMM system classes, either CA or AA models, containing two different definitions of native contacts. One of the configurations is defined as the main model and the other is considered as the alternate model. All forcefield and topology parameters, different than the native contacts, are passed from the main configuration into the multi basin model. Then, the contacts are compared between the input configurations to define the sets of common and unique contacts. Common contacts with equilibrium length distances that differ more than a threshold are defined as dual basin and are assigned a special non-bonded Gaussian potential. The rest of the contacts are treated as single minima and are modeled with a Lennard-Jones (default) or a single basin Gaussian potential.

The multi basin Gaussian potential is defined as:

.. math::
        V_{Multi-basin} = \epsilon_{C}((1+(\frac{r_{ex}}{r})^{12})\prod_{minima}G(r,r_{0}^{\alpha})-1)
.. math::
        G(r,r_{0}^{\alpha}) = 1-exp(\frac{-(r-r_{0}^{\alpha})^2}{2\sigma^2})

.. math::
        \sigma^{2} = \frac{(r_{0}^{\alpha})^2}{50ln(2)} 

Here, :math:`\epsilon_{C}` is the native contact energy constant inherited from the main configuration, :math:`r_{ex}` is the contact excluded volume radius, :math:`r_{0}^{\alpha}` is the equilibrium distance for the :math:`alpha`-th configuration and :math:`r` is the current contact distance. :math:`\sigma` is a parameter that modulates the well amplitude of the :math:`V_{Multi-basin}` energy function. The single and double basin gaussian potential are distinguished by the number of :math:`r_{0}^{\alpha}` parameters given. 

The Lennard Jones contact potential is inherited accordingly from the CA or AA models used to build the multi basin SBM. 

The method to create a multi basin model is:

sbmOpenMM.models.getMultiBasinModel(main_model, alternate_configuration=alternate_model)

Here, main_model and alternate_model are initialized sbmOpenMM system classes containing full force field parameter definitions.
