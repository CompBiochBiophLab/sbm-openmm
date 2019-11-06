## Reference Manual

## sbmOpenMM.Model

### - getAllAtomModel()

​        Creates an All Atom SBM system class object with default initialized parameters.

```python
getAllAtomModel( pdb_file , contact_file, 

                        default_parameters=True, 

                        default_forces=True, 

                        group_by_bb_and_sc=True,

                        create_system=True,

                        minimise=False)
```

​        Initialises a Full Atom SBM OpenMM system class from a PDB coordinates

​        file and a contact file defining the native contacts for the model

| Parameters         | Type            | Details                                                      |
| ------------------ | --------------- | ------------------------------------------------------------ |
| pdb_file           | string          | Name of the input PDB file                                   |
| contact_file       | float           | Name of the input native contact file. The file can be an output from the SMOG program or a two column file defining the atoms |
| default_parameters | boolean(True)   | Wheter to initilize the system with the default parameters for the All Atom model SBM forcefield. |
| default_forces     | boolean(True)   | Whether to add default SBM All Atom forcefield parameters to the model. Set to False if the parameters which to be majoritarily different from the default. |
| group_by_bb_and_sc | boolean (True)  | Wether to classify the torsions into backbone and side-chain  to partition the torsional energy for each torsion in the system. |
| create_system      | boolean (True)  | If True the function will call the createSystemObject() methodto create the system openmm object. If modifications to the default |
| minimise           | boolean (False) | Whether to minimise the system (with default options) if large |

​        Returns

​        \----------------------------------------

​        sbm : sbmOpenMm.system

​            Initializes a sbmOpenMM.system class with default options for 

​            defining an All Atom SBM force field

### - getCAModel()

​        Creates an Carbon-Alpha only SBM system class object with default 

​        initialized parameters. 

```python
getCAModel(pdb_file, contact_file, 

                   default_parameters=True, 

                   default_forces=True, 

                   create_system=True):
```

​      

| Parameters         | Type           | Details                                                      |
| ------------------ | -------------- | ------------------------------------------------------------ |
| pdb_file           | string         | Name of the input PDB file                                   |
| contact_file       | string         | Name of the input native contact file. The file can be an from the SMOG program or a two column file defining the atoms to be paired. |
| default_parameters | boolean (True) | Wheter to initilize the system with the default parameters for the CA model SBM forcefield. |
| default_forces     | boolean (True) | Whether to add default SBM CA forcefield parameters to the model. Set to False if the parameters which to be majoritarily different |
| create_system      | boolean (True) | If True the function will call the createSystemObject() method to create the system openmm object. If modifications to the default forcefield are necessary this option should be given False. |

​        Returns

​        \-------

​        sbm : sbmOpenMm.system

​            Initializes a sbmOpenMM.system class with default options for 

​            defining an All Atom SBM force field.



## sbmOpenMM.System

##### -  removeHydrogens()

​        Remove hydrogens from the input pdb by using a regexpression pattern.

​        Used specially for creating All Atom models.



##### -  getCAlphaOnly()

​        Filter in only Alpha Carbon atoms from the input pdb and updates 

​        the topology object to add new bonds between them.

​        Used specially for creating C-alpha corse-grained models.



##### -  getAtoms()

​        Reads atoms from topology, add them to the main class and sort them 

​        into a dictionary to store their forcefield properties.



##### -   getBonds()

​        Reads bonds from topology, add them to the main class and sort them 

​        into a dictionary to store their forcefield properties.



##### -   getAngles()

​        Reads bonds from topology and creates a list of all possible angles,

​        add them to the main class and sort them into a dictionary to store 

​        their forcefield properties.



##### -  getProperTorsions()

​        Using the created angles by getAngles() creates a list of all possible 

​        proper torsion dihedral angles, filtering out torsions based on 

​        residue-specific rules (only all-atom model). The torsions are then 

​        added to the main class and sorted into a dictionary to store their

​        forcefield properties.



##### -  getImpropers()

​        Create improper torsions based on backbone and sidechain residue-specific 

​        rules, add them to the main class and sort them into a dictionary 

​        to store their forcefield properties. Used specially for simulating

​        All Atom systems.



##### -  getPlanars()

​        Create planar torsions based on backbone and sidechain residue-specific 

​        rules, add them to the main class and sort them into a dictionary 

​        to store their forcefield properties. Used specially for simulating

​        All Atom systems.



##### -   readContactFile()

​        Reads a file containing native contact information and adds them 

​        into the main class. The file can be smog-style (4 columns) or given 

​        as 2 columns, which will be automatically detected.



##### -  setBondParameters()

​        Allows to change the forcefield parameters for bonded terms.



##### -   setAngleParameters()

​        Allows to change the forcefield parameters for angle terms.



##### -  setProperTorsionParameters()

​        Allows to change the forcefield parameters for proper torsion terms.



##### -  setImproperParameters()

​        Allows to change the forcefield parameters for improper torsion  terms.



##### -  setPlanarParameters()

​        Allows to change the forcefield parameters for planar torsion terms.



##### -  setNativeContactParameters()

​        Allows to change the forcefield parameters for native contact terms.



##### -  setParticlesMasses()

​        Allows to change the mass parameter for each atom in the system.



##### -   addHarmonicBondForces()

​        Creates an harmonic bonded force term for each bond in the main

​        class using their defined forcefield parameters.



##### -   addHarmonicAngleForces()

​        Creates an harmonic angle force term for each angle in the main

​        class using their defined forcefield parameters.



##### -  addPeriodicTorsionForces()

​        Creates an periodic torsion force term for each proper torsion in 

​        the main class using their defined forcefield parameters.



##### -   addGeneralPeriodicTorsionForces()

​        Creates an periodic torsion force term for each proper torsion in 

​        the main class using their defined forcefield parameters.



##### -  addHarmonicImproperForces()

​        Creates an harmonic torsion force term for each improper torsion 

​        in the main class using their defined forcefield parameters. Used 

​        specially for simulating All Atom systems. 



##### -  addHarmonicPlanarForces()

​        Creates an harmonic torsion force term for each planar torsion in

​        the main class using their defined forcefield parameters. Used specially

​        for simulating All Atom systems.



##### -   addLJ12_6ContactForces()

​        Creates a 12-6 Lennard-Jones bond potential for each native contact

​        in the main class using their defined forcefield parameters. Used 

​        specially for simulating All Atom systems.



##### -  addLJ12_10ContactForces()

​        Creates a 12-10 Lennard-Jones bond potential for each native contact

​        in the main class using their defined forcefield parameters. Used 

​        specially for simulating Coarse grained Carbon-Alpha systems.



##### -  addGaussianContactForces()

​        Creates a gaussian single and double basin bond potential for each 

​        native contact in the main class using their defined forcefield 

​        parameters. The contacts are recognized according two the number 

​        of parameters given as values in the attribute system.contacts.



##### -  addLJRepulsionForces()

​        Creates a repulsive-only 12 Lennard-Jones non-bonded potential specifying 

​        a exclusion list for bond, angle, torsion, and native contact terms. 



##### -  groupTorsionsbyBBAndSC()

​        Groups proper torsions by backbone and sidechain torsions. Used 

​        exclusively for simulating All Atom systems.



##### -  getAATorsionParameters()

​        Generate default periodic torsion forcefield parameters, for proper

​        torsions, using pre-defined assignment schemes. Used exclusively for 

​        simulating All Atom systems.



##### -  getAANativeContactParameters()

​        Generate default bonded contact forcefield parameters, for native 

​        contacts, using pre-defined assignment schemes. Used exclusively for 

​        simulating All Atom systems.



##### -  createSystemObject()

​        Create OpenMM system object adding particles, masses and forces. 

​        It also groups the added forces into Force-Groups.



##### -  addParticles()

​        Add particles to the main class OpenMM system instance.



##### -  addSystemForces()

​        Add forces to the main class OpenMM system instance. It also save

​        names of the added forces to include in the reporter class.



##### -  dumpPdb()

​        Writes a pdb file of the system in its current state.



##### -  dumpForceFieldData()

​        Writes to a file the parameters of the SBM forcefield.