#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from simtk.openmm.app import *
from simtk.openmm import *
from simtk import unit

from collections import OrderedDict
import numpy as np
import re


# In[ ]:


class system:
    """
    A class containing methods and parameters for generating Structure Based
    Models (SBM) systems to be simulated using the OpenMM interface. It offers 
    flexibility to create default and custom SBM systems and to easily 
    modify their parameters.

    Attributes
    ----------
    pdb_path : string
        Path to the pdb input file
    pdb : simtk.openmm.app.pdbfile.PDBFile
        Object that holds the information of OpenMM PDB parsing method.
    topology : simtk.openmm.app.topology.Topology
        OpenMM topology of the model.
    positions : simtk.unit.quantity.Quantity
        Atomic positions of the model.
    particles_mass : float or list
        Mass of each particle. If float then uniform masses are given to all 
        particles. If list per-particle masses are assigned.
    model_type : string
        String representing the model type: All-atom (AA), alpha-carbon (CA)
        or multi-basin variants (AA-MB, CA-MB).
    atoms : list
        A list of the current atoms in the model. The items are simtk.openmm.app.topology.Atom
        initialised classes.
    n_atoms : int
        Total numer of atoms in the model.
    bonds : collections.OrderedDict
        A dict that uses bonds (2-tuple of simtk.openmm.app.topology.Atom objects) 
        present in the model as keys and their forcefield properties as values.
    n_bonds : int
        Total numer of bonds in the model.
    angles : collections.OrderedDict
        A dict that uses angles (3-tuple of simtk.openmm.app.topology.Atom objects)
        present in the model as keys and their forcefield properties as values.
    n_angles : int
        Total numer of angles in the model.
    torsions : collections.OrderedDict
        A dict that uses proper torsions (4-tuple of simtk.openmm.app.topology.Atom objects) 
        present in the model as keys and their forcefield properties as values.
    n_torsions : int
        Total numer of proper torsions in the model.
    impropers : collections.OrderedDict
        A dict that uses improper torsions (4-tuple of simtk.openmm.app.topology.Atom objects)
        present in the model as keys and their forcefield properties as values.
    n_impropers : int
        Total numer of improper torsions in the model.
    planars : collections.OrderedDict
        A dict that uses planar torsions (4-tuple of simtk.openmm.app.topology.Atom objects)
        present in the model as keys and their forcefield properties as values.
    n_planars : int
        Total numer of planar torsions in the model.
    contacts : collections.OrderedDict
        A dict that uses native contacts (2-tuple of simtk.openmm.app.topology.Atom objects
        present in the model as keys and their forcefield properties as values.
    n_contacts : int
        Total numer of native contacts in the model.
    torsions_group : dict
        A dict that uses proper torsions (4-tuple of simtk.openmm.app.topology.Atom objects)
        present in the model as keys and the number of torsions (int) that share the same middle 
        bond atoms as values.
    torsions_type : dict
        A dict that uses proper torsions (4-tuple of simtk.openmm.app.topology.Atom objects)
        present in the model as keys and, a string representing wheter the torsion 
        is classified as 'backbone' or 'sidechain' as values.
    energy_constant : dict
        A dict that holds the value for the different energy terms parameters
        used by different SBM models and forces.
    harmonicBondForce : simtk.openmm.openmm.HarmonicBondForce
        Stores the OpenMM HarmonicBondForce initialised-class. Implements
        an harmonic bond potential between pairs of particles, that depends
        quadratically on their distance.
    harmonicAngleForce : simtk.openmm.openmm.HarmonicAngleForce
        Stores the OpenMM HarmonicAngleForce initialised-class. Implements
        an harmonic angle potential between trios of particles, that depends
        quadratically on their angle length.
    periodicTorsionForce : simtk.openmm.openmm.CustomTorsionForce
        Stores the OpenMM CustomTorsionForce initialised-class. Implements 
        a force potential that varies periodically with the value of the 
        proper torsion angle.
    generalPeriodicTorsionForce: simtk.openmm.openmm.CustomTorsionForce
        Stores the OpenMM CustomTorsionForce initialised-class. Implements 
        a general force potential that varies periodically with the value of the 
        proper torsion angle.
    harmonicImproperForce : simtk.openmm.openmm.CustomTorsionForce
        Stores the OpenMM CustomTorsionForce initialised-class. Implements 
        a force potential that varies quadratically with the value of the  
        improper torsion angle.
    harmonicPlanarForce : simtk.openmm.openmm.CustomTorsionForce
        Stores the OpenMM CustomTorsionForce initialised-class. Implements 
        a force potential that varies quadratically with the value of the  
        planar torsion angle.
    lj12_6contactForce : simtk.openmm.openmm.CustomBondForce
        Stores the OpenMM CustomBondForce initialised-class to be applied
        to non-bonded interactions between native contact pairs. Implements 
        a lennard-jones potential with exponents 12 and 6 for the repulsive
        and attractive componenets, respectively.
    lj12_10contactForce : simtk.openmm.openmm.CustomBondForce
        Stores the OpenMM CustomBondForce initialised-class to be applied
        to non-bonded interactions between native contact pairs. Implements 
        a lennard-jones potential with exponents 12 and 10 for the repulsive
        and attractive components, respectively.
    lj12_10_6contactForce
        Stores the OpenMM CustomBondForce initialised-class to be applied
        to non-bonded interactions between native contact pairs. Implements 
        a lennard-jones potential with exponents 12 and 10 for the repulsive
        and attractive components, respectively, and an additional 6-exponent 
        term to model a “desolvation penalty” for forming/breaking contacts.
    singleGaussianContactForce : simtk.openmm.openmm.CustomNonbondedForce
        Stores the OpenMM CustomBondForce initialised-class to be applied
        to non-bonded interactions between native contact pairs with single
        minimum. Implements a mixed lennard-jones (repulsive) and gaussian
        potential (attractive) with separate control of equilibrium distance
        and excluded volume.
    doubleGaussianContactForce : simtk.openmm.openmm.CustomNonbondedForce
        Stores the OpenMM CustomBondForce initialised-class to be applied
        to non-bonded interactions between native contact pairs with two 
        minima. Implements a mixed lennard-jones (repulsive) and gaussian
        potential (attractive) with separate control of the equilibrium 
        distances and excluded volume.
    ljRepulsionForce : simtk.openmm.openmm.CustomNonbondedForce
        Stores the OpenMM CustomBondForce initialised-class to be applied
        to the non-bonded interactions between non-native contact pairs. 
        Implements only a repulsive lennard-jones potential with exponent 
        12.
    forceGroups : collections.OrderedDict
        A dict that uses force names as keys and their corresponding force
        as values.
    system : simtk.openmm.openmm.System
        Stores the OpenMM System initialised class. It stores all the forcefield
        information for the SBM model.
    rf_epsilon : float
        Epsilon parameter used in the repulsion force object.
    rf_sigma : float
        Sigma parameter used in the repulsion force object.
    rf_cutoff : float
        Cutoff value used for repulsion force interactions.
    exclusions : list
        List of added exclusions for the repuslion non-bonded term.
    
    Methods
    -------
    removeHydrogens()
        Remove hydrogens from the input pdb by using a regexpression pattern.
        Used specially for creating all atom (AA) models.
    getCAlphaOnly()
        Filter in only alpha carbon atoms from the input pdb and updates 
        the topology object to add new bonds between them. Used specially for 
        creating alpha-carbon (CA) corse-grained models.
    getAtoms()
        Reads atoms from topology, adds them to the main class and sorts them 
        into a dictionary to store their forcefield properties.
    getBonds()
        Reads bonds from topology, adds them to the main class and sorts them 
        into a dictionary to store their forcefield properties.
    getAngles()
        Reads bonds from topology and creates a list of all possible angles,
        adds them to the main class and sorts them into a dictionary to store 
        their forcefield properties.
    getProperTorsions()
        Using the created angles, usually by the getAngles() method, creates 
        a list of all possible proper torsion dihedral angles, filtering out 
        torsions based on residue-specific rules (only for the all-atom model).
        The torsions are then added to the main class and sorted into a dictionary 
        to store their forcefield properties.
    getImpropers()
        Create improper torsions based on backbone and sidechain residue-specific 
        rules (all-atom model only), adds them to the main class and sorts them 
        into a dictionary to store their forcefield properties.
    getPlanars()
        Create planar torsions based on backbone and sidechain residue-specific 
        rules (all-atom model only), adds them to the main class and sorts them 
        into a dictionary to store their forcefield properties.
    readContactFile()
        Reads a file containing native contact information and adds it into the 
        main class. The file can be smog-style (4 columns) or given as 2 columns.
        The format will be automatically detected. 
    setBondParameters()
        Change the forcefield parameters for bonded terms.
    setAngleParameters()
        Change the forcefield parameters for angle terms.
    setProperTorsionParameters()
        Change the forcefield parameters for proper torsion terms.
    setImproperParameters()
        Change the forcefield parameters for improper torsion 
        terms.
    setPlanarParameters()
        Change the forcefield parameters for planar torsion terms.
    setNativeContactParameters()
        Change the forcefield parameters for native contact terms.
    setParticlesMasses()
        Change the mass parameter for each atom in the system.
    addHarmonicBondForces()
        Creates an harmonic bonded force term for each bond in the main
        class using their defined forcefield parameters.
    addHarmonicAngleForces()
        Creates an harmonic angle force term for each angle in the main
        class using their defined forcefield parameters.
    addPeriodicTorsionForces()
        Creates a periodic torsion force term for each proper torsion in 
        the main class using their defined forcefield parameters.
    addGenericPeriodicTorsionForces()
        Creates a periodic torsion force term for each proper torsion in 
        the main class using their defined forcefield parameters.
    addHarmonicImproperForces()
        Creates an harmonic torsion force term for each improper torsion 
        in the main class using their defined forcefield parameters. Used 
        specially for simulating All Atom systems.
    addHarmonicPlanarForces()
        Creates an harmonic torsion force term for each planar torsion in
        the main class using their defined forcefield parameters. Used specially
        for simulating All Atom systems.
    addLJ12_6ContactForces()
        Creates a 12-6 Lennard-Jones bond potential for each native contact
        in the main class using their defined forcefield parameters. Used 
        specially for simulating All Atom systems.
    addLJ12_10ContactForces()
        Creates a 12-10 Lennard-Jones bond potential for each native contact
        in the main class using their defined forcefield parameters. Used 
        specially for simulating coarse-grained alpha-carbon systems.
    addGaussianContactForces()
        Creates a gaussian single and double basin bond potential for each 
        native contact in the main class using their defined forcefield 
        parameters. The contacts are recognized according two the number 
        of parameters given as values in the attribute system. 3-parameters for
        single-basin potential and 4-paramters for dual-basin potentials.
    addLJRepulsionForces()
        Creates a repulsive-only 12 Lennard-Jones non-bonded potential specifying 
        a exclusion list for bond, angle, torsion, and native contact terms. 
    groupTorsionsbyBBAndSC()
        Groups proper torsions by backbone and sidechain torsions. Used exclusively
        for simulating all-atom SBM.
    getAATorsionParameters()
        Generates default periodic torsion forcefield parameters, for proper
        torsions, using pre-defined assignment schemes. Used exclusively for 
        simulating all-atom SBM.
    getAANativeContactParameters()
        Generates default bonded contact forcefield parameters, for native 
        contacts, using pre-defined assignment schemes. Used exclusively for 
        simulating all-atom SBM.
    createSystemObject()
        Creates OpenMM system object adding particles, masses and forces. 
        It also groups the added forces into Force-Groups for the sbmReporter
        class.
    addParticles()
        Add particles to the system OpenMM class instance.
    addSystemForces()
        Add forces to the system OpenMM class instance. It also save
        names for the added forces to include them in the reporter class.
    dumpPdb()
        Writes a PDB file of the system in its current state.
    dumpForceFieldData()
        Writes to a file the parameters of the SBM forcefield.
    loadForcefieldFromFile()
        Loads forcefield parameters from a sbmOpenMM force field file written with
        the dumpForceFieldData() method.
    
    """
    
    def __init__(self, pdb_path, particles_mass=1.0):
        """
        Initialises the SBM OpenMM system class.
        
        Parameters
        ----------
        pdb_path : string
            Name of the input PDB file
        particles_mass : float or list
            mass of all the particles in the system.
            
        Returns
        -------
        None
        """
        
        #Define structure object attributes
        self.pdb_path = pdb_path
        self.pdb = PDBFile(pdb_path)
        self.topology = self.pdb.topology
        self.positions = self.pdb.positions
        self.particles_mass = particles_mass
        self.model_type = None
        
        #Define geometric attributes
        self.atoms = []
        self.n_atoms = None
        self.bonds = OrderedDict()
        self.n_bonds = None
        self.angles = OrderedDict()
        self.n_angles = None
        self.torsions = OrderedDict()
        self.n_torsions = None
        self.impropers = OrderedDict()
        self.n_impropers = None
        self.planars = OrderedDict()
        self.n_planars = None
        self.contacts = OrderedDict()
        self.n_contacts = None
        
        #Define force field attributes
        self.torsions_group = {}
        self.torsions_type = {}
        self.energy_constant = {}
        
        #Define force attributes
        self.harmonicBondForce = None
        self.harmonicAngleForce = None
        self.periodicTorsionForce = None
        self.generalPeriodicTorsionForce = None
        self.harmonicImproperForce = None
        self.harmonicPlanarForce = None
        self.lj12_6contactForce = None
        self.lj12_10contactForce = None
        self.lj12_10_6contactForce = None
        self.singleGaussianContactForce = None
        self.doubleGaussianContactForce = None
        self.ljRepulsionForce = None
        self.rf_sigma = None
        self.rf_epsilon = None
        self.rf_cutoff = None
        self.exclusions = []
        
        self.forceGroups = OrderedDict()
        
        #Initialise an OpenMM system class instance
        self.system = openmm.System()
        
    def removeHydrogens(self):
        """
        Removes all hydrogen atoms in the topology from the SBM system.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #save all hydrogen atoms
        atomsToRemove = []
        _hydrogen = re.compile("[123 ]*H.*")
        for a in self.topology.atoms():
            if _hydrogen.match(a.name):
                atomsToRemove.append(a)
        
        #Remove all hydrogen atoms
        modeller_topology = modeller.Modeller(self.topology, self.positions)
        modeller_topology.delete(atomsToRemove)
        self.topology = modeller_topology.getTopology()
        self.positions = modeller_topology.getPositions()
        self.model_type = 'AA'
        
    def getCAlphaOnly(self):
        """
        Keeps in the SBM system only the alpha carbon atoms from the OpenMM topology.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #save all non C-alpha atoms
        atomsToRemove = []
        oldIndex = []
        for a in self.topology.atoms():
            if a.name != 'CA':
                atomsToRemove.append(a)
        
        #Remove all non C-alpha atoms
        modeller_topology = modeller.Modeller(self.topology, self.positions)
        modeller_topology.delete(atomsToRemove)
        self.topology = modeller_topology.getTopology()
        self.positions = modeller_topology.getPositions()
        
        #Update system atoms
        atoms = list(self.topology.atoms())
        
        #Add bonds between C-alpha atoms of the same chain
        for i in range(1,len(atoms)):
            current_chain = atoms[i].residue.chain
            if i == 1:
                previous_chain = current_chain
                
            #Add bond only when atoms belong to the same chain
            if current_chain == previous_chain:
                self.topology.addBond(atoms[i-1], atoms[i])
            previous_chain = current_chain
        
        self.model_type = 'CA'
            
    def getAtoms(self):
        """
        Adds atoms in the OpenMM topology instance to the sbmOpenMM system class.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Get Atoms From Topology
        atoms = []
        for atom in self.topology.atoms():
            atoms.append(atom)
        
        #Sort atoms by index
        atoms = sorted(atoms, key=lambda x:x.index)
        
        #Add atoms to sbm object
        self.n_atoms = 0
        for atom in atoms:
            self.atoms.append(atom)
            self.n_atoms += 1
            
    def getBonds(self):
        """
        Adds bonds in the OpenMM topology instance to the sbmOpenMM system class.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Get Bonds From Topology
        bonds = []
        for bond in self.topology.bonds():
            if bond[0].index > bond[1].index:
                bonds.append((bond[1], bond[0]))
            else:
                bonds.append((bond[0], bond[1]))
        
        #Sort bonds by index of first atom
        bonds = sorted(bonds, key=lambda x:x[0].index)
        
        #Add bonds to sbm object
        self.n_bonds = 0
        for bond in bonds:
            p1 = self.positions[bond[0].index]
            p2 = self.positions[bond[1].index]
            bond_length = geometry.bond(p1, p2)
            self.bonds[bond] = (bond_length, None)
            self.n_bonds += 1
            
        #Record which atoms are bonded to each other
        self.bondedTo = {}
        for atom in self.topology.atoms():
            self.bondedTo[atom] = []
        for bond in self.topology.bonds():
            self.bondedTo[bond[0]].append(bond[1])
            self.bondedTo[bond[1]].append(bond[0])
            
    def getAngles(self):
        """
        Adds angles to the sbmOpenMM system class based on the bondings in the topology 
        object.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        # Make a list of all unique angles
        uniqueAngles = set()
        for bond in self.bonds:
            for atom in self.bondedTo[bond[0]]:
                if atom != bond[1]:
                    if atom.index < bond[1].index:
                        uniqueAngles.add((atom, bond[0], bond[1]))
                    else:
                        uniqueAngles.add((bond[1], bond[0], atom))
            for atom in self.bondedTo[bond[1]]:
                if atom != bond[0]:
                    if atom.index > bond[0].index:
                        uniqueAngles.add((bond[0], bond[1], atom))
                    else:
                        uniqueAngles.add((atom, bond[1], bond[0]))
                        
        #Sort angles by index of first atom
        uniqueAngles = sorted(list(uniqueAngles), key=lambda x:x[0].index)
        
        #Add angles to sbm object
        self.n_angles = 0
        for angle in uniqueAngles:
            p1 = self.positions[angle[0].index]
            p2 = self.positions[angle[1].index]
            p3 = self.positions[angle[2].index]
            angle_length = geometry.angle(p1, p2, p3)
            self.angles[angle] = (angle_length, None)
            self.n_angles += 1
        
    def getProperTorsions(self):
        """
        Adds proper torsions to the sbmOpenMM system class based on the bonded angles
        in it. Excludes special torsions for rings and backbones.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        # Make a list of all unique proper torsions
        uniqueTorsions = set()
        for angle in self.angles:
            for atom in self.bondedTo[angle[0]]:
                if atom not in angle:
                    if atom.index < angle[2].index:
                        uniqueTorsions.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniqueTorsions.add((angle[2], angle[1], angle[0], atom))
                        
            for atom in self.bondedTo[angle[2]]:
                if atom not in angle:
                    if atom.index > angle[0].index:
                        uniqueTorsions.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniqueTorsions.add((atom, angle[2], angle[1], angle[0]))
        
        #Create dictionary with ring atoms for generating exclusions
        self.ring_exclusions = {}
        self.ring_exclusions['PRO'] = ['CA', 'N', 'CB', 'CD', 'CG']
        self.ring_exclusions['HIS'] = ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
        self.ring_exclusions['TYR'] = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
        self.ring_exclusions['PHE'] = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
        self.ring_exclusions['TRP'] = ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
        
        #Exclude torsions inside rings
        toExclude = set()
        for torsion in uniqueTorsions:
            if list(np.array([t.residue.index for t in torsion]) - torsion[0].residue.index) == [0, 0, 0, 0] or                list(np.array([t.residue.index for t in torsion]) - torsion[0].residue.index) == [0, 1, 1, 1]:
                count = {}
                for a in torsion:
                    try:
                        count[a.residue.name] += 1
                    except:
                        count[a.residue.name] = 1
                for r in count:
                    if count[r] >= 3:
                        if r in list(self.ring_exclusions.keys()):
                            in_ring = 0
                            for a in torsion:
                                if a.name in self.ring_exclusions[r]:
                                    in_ring += 1
                            if in_ring >= 3:
                                toExclude.add(torsion)
                                
        #Exclude specific backbone torsions
        for torsion in uniqueTorsions:
            if list(np.array([t.residue.index for t in torsion]) - torsion[0].residue.index) == [0, 0, 1, 1]: 
                if [t.name for t in torsion] == ['CA', 'C', 'N', 'CA']:
                    toExclude.add(torsion)
                if [t.name for t in torsion] == ['O', 'C', 'N', 'CA']:
                    toExclude.add(torsion)
            if list(np.array([t.residue.index for t in torsion]) - torsion[0].residue.index) == [0, 1, 1, 1]: 
                if [t.residue.name for t in torsion[1:]] == ['PRO','PRO','PRO']:
                    if [t.name for t in torsion[1:]] == ['N', 'CA', 'C']:
                        toExclude.add(torsion)
            if [t.residue.name for t in torsion[2:]] == ['PRO','PRO']:
                if [t.name for t in torsion[2:]] == ['N', 'CD']:
                    toExclude.add(torsion)
        
        #Remove added torsions
        for torsion in toExclude:
            uniqueTorsions.remove(torsion)
        self.uniqueTorsions = uniqueTorsions
        
        #Sort torsions by the first atom index
        uniqueTorsions = sorted(list(uniqueTorsions), key=lambda x:x[0].index)
        
        #Add torsions to sbm object
        self.n_torsions = 0
        for torsion in uniqueTorsions:
            p1 = self.positions[torsion[0].index]
            p2 = self.positions[torsion[1].index]
            p3 = self.positions[torsion[2].index]
            p4 = self.positions[torsion[3].index]
            torsion_angle = geometry.torsion(p1, p2, p3, p4)
            self.torsions[torsion] = (torsion_angle, None)
            self.n_torsions += 1
            
    def getImpropers(self):
        """
        Adds improper torsions to the sbmOpenMM system class to maintain backbone 
        chiralities.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Iterate by every chain in topology to add backbone impropers
        chains = [c for c in self.topology.chains()]
        for chain in chains:
            #Iterate by every residue in chain to add backbone impropers
            residues = [r for r in chain.residues()]
            for i in range(len(residues)-1):
                ad1 = {a.name:a for a in residues[i].atoms()}
                ad2 = {a.name:a for a in residues[i+1].atoms()}

                if residues[i+1].name != 'PRO':
                    self.impropers[(ad1['CA'], ad1['C'], ad2['N'], ad2['CA'])] = None
                    self.impropers[(ad1['O'], ad1['C'], ad2['N'], ad2['CA'])] = None

            #Add impropers to sbm object
            self.n_impropers = 0
            for improper in self.impropers:
                p1 = self.positions[improper[0].index]
                p2 = self.positions[improper[1].index]
                p3 = self.positions[improper[2].index]
                p4 = self.positions[improper[3].index]
                improper_angle = geometry.torsion(p1, p2, p3, p4)
                self.impropers[improper] = (improper_angle, None)
                self.n_impropers += 1
            
    def getPlanars(self):
        """
        Adds planar torsions to the sbmOpenMM system class to mantain side chain and 
        backbone planar arrangements.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Iterate by every chain in topology to add specific planar restraints
        chains = [c for c in self.topology.chains()]
        for chain in chains:
            #Iterate by every residue in chain to add specific planar restraints
            residues = [r for r in chain.residues()]
            for i in range(len(residues)):
                if i < len(residues) - 1:
                    ad1 = {a.name:a for a in residues[i].atoms()}
                    ad2 = {a.name:a for a in residues[i+1].atoms()}
                    
                    rn1 = residues[i].name[:3]
                    rn2 = residues[i+1].name[:3]

                    self.planars[(ad1['O'], ad1['CA'], ad1['C'], ad2['N'])] = None

                    if rn1 != 'GLY':
                        self.planars[(ad1['N'], ad1['C'], ad1['CA'], ad1['CB'])] = None

                else:
                    ad1 = {a.name:a for a in residues[i].atoms()}
                    ad2 = None

                    rn1 = residues[i].name[:3]
                    rn2 = None
                    
                    if rn1 != 'GLY':
                        self.planars[(ad1['N'], ad1['C'], ad1['CA'], ad1['CB'])] = None
                        
                    if 'OXT' in ad1:
                        self.planars[(ad1['CA'], ad1['O'], ad1['C'], ad1['OXT'])] = None

                if rn1 == 'ARG': 
                    self.planars[(ad1['CZ'], ad1['NH1'], ad1['NH2'], ad1['NE'])] = None
                    self.planars[(ad1['CZ'], ad1['NH1'], ad1['NH2'], ad1['CD'])] = None

                elif rn1 == 'ASN':
                    self.planars[(ad1['CB'], ad1['ND2'], ad1['CG'], ad1['OD1'])] = None

                elif rn1 == 'ASP':
                    self.planars[(ad1['CB'], ad1['OD1'], ad1['CG'], ad1['OD2'])] = None

                elif rn1 == 'GLN':
                    self.planars[(ad1['CG'], ad1['NE2'], ad1['CD'], ad1['OE1'])] = None

                elif rn1 == 'GLU':
                    self.planars[(ad1['CG'], ad1['OE1'], ad1['CD'], ad1['OE2'])] = None

                elif rn1 == 'HIS':
                    self.planars[(ad1['ND1'], ad1['CD2'], ad1['CG'], ad1['CB'])] = None
                    self.planars[(ad1['CG'], ad1['NE2'], ad1['ND1'], ad1['CD2'])] = None
                    self.planars[(ad1['ND1'], ad1['CD2'], ad1['CE1'], ad1['CG'])] = None
                    self.planars[(ad1['CE1'], ad1['CG'], ad1['NE2'], ad1['ND1'])] = None
                    self.planars[(ad1['NE2'], ad1['ND1'], ad1['CD2'], ad1['CE1'])] = None
                    self.planars[(ad1['CD2'], ad1['CE1'], ad1['CG'], ad1['NE2'])] = None

                elif rn1 == 'PHE':
                    self.planars[(ad1['CD1'], ad1['CD2'], ad1['CG'], ad1['CB'])] = None
                    self.planars[(ad1['CG'], ad1['CZ'], ad1['CD1'], ad1['CE2'])] = None
                    self.planars[(ad1['CD1'], ad1['CE2'], ad1['CE1'], ad1['CD2'])] = None
                    self.planars[(ad1['CE1'], ad1['CD2'], ad1['CZ'], ad1['CG'])] = None
                    self.planars[(ad1['CZ'], ad1['CG'], ad1['CE2'], ad1['CD1'])] = None
                    self.planars[(ad1['CE2'], ad1['CD1'], ad1['CD2'], ad1['CE1'])] = None
                    self.planars[(ad1['CD2'], ad1['CE1'], ad1['CG'], ad1['CZ'])] = None

                elif rn1 == 'PRO':
                    self.planars[(ad1['CB'], ad1['CA'], ad1['N'], ad1['CD'])] = None
                    self.planars[(ad1['C'], ad1['CA'], ad1['N'], ad1['CD'])] = None

                elif rn1 == 'TRP':
                    self.planars[(ad1['CB'], ad1['CG'], ad1['CD2'], ad1['CE2'])] = None
                    self.planars[(ad1['CG'], ad1['NE1'], ad1['CH2'], ad1['CE3'])] = None
                    self.planars[(ad1['CD1'], ad1['CE2'], ad1['CZ3'], ad1['CD2'])] = None
                    self.planars[(ad1['NE1'], ad1['CZ2'], ad1['CE3'], ad1['CG'])] = None
                    self.planars[(ad1['CE2'], ad1['CH2'], ad1['CD2'], ad1['CD1'])] = None
                    self.planars[(ad1['CZ2'], ad1['CZ3'], ad1['CG'], ad1['NE1'])] = None
                    self.planars[(ad1['CH2'], ad1['CE3'], ad1['CD1'], ad1['CE2'])] = None
                    self.planars[(ad1['CZ3'], ad1['CD2'], ad1['NE1'], ad1['CZ2'])] = None
                    self.planars[(ad1['CE3'], ad1['CG'], ad1['CE2'], ad1['CH2'])] = None
                    self.planars[(ad1['CD2'], ad1['CD1'], ad1['CZ2'], ad1['CZ3'])] = None

                elif rn1 == 'TYR':
                    self.planars[(ad1['CD1'], ad1['CD2'], ad1['CG'], ad1['CB'])] = None
                    self.planars[(ad1['CG'], ad1['CZ'], ad1['CD1'], ad1['CE2'])] = None
                    self.planars[(ad1['CD1'], ad1['CE2'], ad1['CE1'], ad1['CD2'])] = None
                    self.planars[(ad1['CE1'], ad1['CD2'], ad1['CZ'], ad1['CG'])] = None
                    self.planars[(ad1['CZ'], ad1['CG'], ad1['CE2'], ad1['CD1'])] = None
                    self.planars[(ad1['CE2'], ad1['CD1'], ad1['CD2'], ad1['CE1'])] = None
                    self.planars[(ad1['CD2'], ad1['CE1'], ad1['CG'], ad1['CZ'])] = None
                    self.planars[(ad1['CE1'], ad1['CE2'], ad1['CZ'], ad1['OH'])] = None

                elif rn1 == 'LEU':
                    self.planars[(ad1['CB'], ad1['CG'], ad1['CD1'], ad1['CD2'])] = None

                elif rn1 == 'THR':
                    self.planars[(ad1['CA'], ad1['CB'], ad1['OG1'], ad1['CG2'])] = None

                elif rn1 == 'VAL' or rn1 == 'ILE':
                    self.planars[(ad1['CA'], ad1['CB'], ad1['CG1'], ad1['CG2'])] = None

                if rn2 == 'PRO':
                    self.planars[(ad1['C'], ad2['N'], ad2['CA'], ad2['C'])] = None
                    self.planars[(ad1['C'], ad2['N'], ad2['CA'], ad2['CB'])] = None
            
        #Add planars to sbm object
        self.n_planars = 0
        for planar in self.planars:
            p1 = self.positions[planar[0].index]
            p2 = self.positions[planar[1].index]
            p3 = self.positions[planar[2].index]
            p4 = self.positions[planar[3].index]
            planar_angle = geometry.torsion(p1, p2, p3, p4)
            self.planars[planar] = (planar_angle, None)
            self.n_planars += 1
    
    def readContactFile(self, contact_file, shift=1):
        """
        Reads a file to add native contact information to the sbmOpenMM system class.
        The file format can be 'smog' like (i.e. 4 columns) or '2column' like
        (i.e. one column for each atom). This is auto-dected by the method.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Detect which format of contact file is being used (new formats can be added if useful)
        file_type = system._checkContactFileType(contact_file)
        if file_type == 'smog':
            i = 1
            j = 3
        elif file_type == '2column':
            i = 0
            j = 1
        else:
            raise ValueError('Incorrect contact file format! Please check carefully your contact file')
        
        #Read contact file and add contacts to sbm object
        self.n_contacts = 0
        with open(contact_file,'r') as cf:
            for line in cf:
                if not line.startswith('#'):
                    c1 = self.atoms[int(line.split()[i])-shift]
                    c2 = self.atoms[int(line.split()[j])-shift]
                    p1 = self.positions[c1.index]
                    p2 = self.positions[c2.index]
                    contact_length = geometry.bond(p1, p2)
                    self.contacts[(c1,c2)] = (contact_length, None)
                    self.n_contacts += 1
                
    ## Functions for setting force specific parameters ##
    
    def setBondParameters(self, bond_parameters):
        """
        Set the harmonic bond constant force parameters. The input can be 
        a float, to set the same parameter for all force interactions, or 
        a list, to define a unique parameter for each force interaction.
        
        Parameters
        ----------
        bond_parameters : float or list
            Parameter(s) to set up for the harmonic bond forces. 
            
        Returns
        -------
        None
        """
        
        system._setParameters(self.bonds, bond_parameters)
    
    def setAngleParameters(self, angle_parameters):
        """
        Set the harmonic angle constant force parameters. The input can 
        be a float, to set the same parameter for all force interactions, 
        or a list, to define a unique parameter for each force interaction. 
        
        Parameters
        ----------
        angle_parameters : float or list
            Parameter(s) to set up for the harmonic angle forces. 
            
        Returns
        -------
        None
        """
        
        system._setParameters(self.angles, angle_parameters)
                
    def setProperTorsionParameters(self, torsion_parameters):
        """
        Set the periodic torsion constant force parameters. The input can 
        be a float, to set the same parameter for all force interactions, 
        or a list, to define a unique parameter for each force interaction. 
        
        Parameters
        ----------
        torsion_parameters : float or list
            Parameter(s) to set up for the periodic torsion forces. 
            
        Returns
        -------
        None
        """
        
        system._setParameters(self.torsions, torsion_parameters)
                
    def setImproperParameters(self, improper_parameters):
        """
        Set the harmonic torsion constant force parameters for impropers. 
        The input can be a float, to set the same parameter for all force 
        interactions, or a list, to define a unique paremater for each force 
        interaction.
        
        Parameters
        ----------
        improper_parameters : float or list
            Parameter(s) to set up for the harmonic torsion force parameters 
            for impropers.
            
        Returns
        -------
        None
        """
        
        system._setParameters(self.impropers, improper_parameters)
                
    def setPlanarParameters(self, planar_parameters):
        """
        Set the harmonic torsion constant force parameters for planars. 
        The input can be a float, to set the same parameter for all force 
        interactions, or a list, to define a unique parameter for each force 
        interaction.
        
        Parameters
        ----------
        planar_parameters : float or list
            Parameter(s) to set up for the harmonic torsion force parameters 
            for planars.
            
        Returns
        -------
        None
        """
        
        system._setParameters(self.planars, planar_parameters)
        
    def setNativeContactParameters(self, contact_parameters):
        """
        Set the native interactions constant force parameters for native 
        contacts. The input can be a float, to set the same parameter for 
        all force interactions, or a list, to define a unique parameter 
        for each force interaction.
        
        Parameters
        ----------
        contact_parameters : float or list
            Parameter(s) to set up for the native interaction force for 
            native contacts. 
            
        Returns
        -------
        None
        """
        
        system._setParameters(self.contacts, contact_parameters)
        
    def setParticlesMasses(self, particles_mass):
        """
        Set the mass of the particles in the system. The input can be a 
        float, to set the same mass for all particles, or a list, to define 
        a unique mass for each particle.
        
        Parameters
        ----------
        particles_mass : float or list
            Mass(es) values to add for the particles in the sbmOpenMM system class.
            
        Returns
        -------
        None
        """
        
        self.particles_mass = particles_mass
    
    ## Functions for creating force objects with defined parameters ##
    
    def addHarmonicBondForces(self):
        """
        Creates an openmm.HarmonicBondForce() object with the bonds and 
        parameters setted up in the "bonds" dictionary attribute. The force object 
        is stored at the "harmonicBondForce" attribute.
        
        The force parameters must be contained in self.bonds as follows:
        self.bonds is a dictionary:
            - The keys are 2-tuples for two atom items in self.topology.atoms attribute.
            - The values are a 2-tuple of parameters in the following order:
                first  -> bond0 (quantity)
                second -> k (float)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        self.harmonicBondForce = openmm.HarmonicBondForce()
        for bond in self.bonds:
            self.harmonicBondForce.addBond(bond[0].index,
                                           bond[1].index,
                                           self.bonds[bond][0],
                                           self.bonds[bond][1])
    
    def addHarmonicAngleForces(self):
        """
        Creates an openmm.HarmonicAngleForce() object with the angles and 
        parameters setted up in the "angles" attribute. The force object 
        is stored at the "harmonicAngleForce" attribute.
        
        The force parameters must be contained in self.angles as follows:
        self.angles is dictionary:
            - The keys are 3-tuples for three atoms in self.topology.atoms attribute.
            - The values are a 2-tuple of parameters in the following order:
                first  -> angle0 (quantity)
                second -> k (float)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        self.harmonicAngleForce = openmm.HarmonicAngleForce()
        for angle in self.angles:
            self.harmonicAngleForce.addAngle(angle[0].index,
                                            angle[1].index,
                                            angle[2].index,
                                            self.angles[angle][0],
                                            self.angles[angle][1])
            
    def addPeriodicTorsionForces(self):
        """
        Creates an openmm.CustomTorsionForce() object with the torsions 
        and parameters setted up in the "torsions" attribute. The custom 
        torsion force is initilized with the formula: 
        
        energy = k*((1-cos(theta-theta0))+0.5*(1-cos(3*(theta-theta0))))
        
        The force object is stored at the "periodicTorsionForce" attribute.
        
        The force parameters must be contained in self.torsions as follows:
        self.torsions is a dictionary:
            - The keys are tuples for four atoms in self.topology.atoms attribute.
            - The values are a tuple of parameters in the following order:
                first  -> theta0 (quantity)
                second -> k (float)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        energy_function = "k*((1-cos(theta-theta0)+0.5*(1-cos(3*(theta-theta0)))))"
        self.periodicTorsionForce = openmm.CustomTorsionForce(energy_function)
        self.periodicTorsionForce.addPerTorsionParameter('theta0')
        self.periodicTorsionForce.addPerTorsionParameter('k')
        for torsion in self.torsions:
            self.periodicTorsionForce.addTorsion(torsion[0].index,
                                                 torsion[1].index,
                                                 torsion[2].index,
                                                 torsion[3].index,
                                                 (self.torsions[torsion][0],
                                                 self.torsions[torsion][1]))
            
    def addGeneralPeriodicTorsionForces(self):
        """
        Creates an openmm.CustomTorsionForce() object with the torsions 
        and parameters setted up in the "torsions" attribute. The custom 
        torsion force is initilized with the formula: 
        
        energy = k*(1+cos(n*(theta-theta0)))
        
        The force object is stored at the "generalPeriodicTorsionForce" attribute.
        
        The force parameters must be contained in self.torsions as follows:
        self.torsions is a dictionary:
            - The keys are tuples for four atoms in self.topology.atoms attribute.
            - The values are a list with the parameters as items:
              Each item is a tuple containing:
                first  -> theta0 (quantity)
                second -> k (float)
                third -> n (int)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        energy_function = 'k*(1+cos(n*(theta-theta0)))'
        
        self.generalPeriodicTorsionForce = openmm.CustomTorsionForce(energy_function)
        self.generalPeriodicTorsionForce.addPerTorsionParameter('theta0')
        self.generalPeriodicTorsionForce.addPerTorsionParameter('k')
        self.generalPeriodicTorsionForce.addPerTorsionParameter('n')
        
        for torsion in self.torsions:
            for t in self.torsions[torsion]:
                self.generalPeriodicTorsionForce.addTorsion(torsion[0].index,
                                                     torsion[1].index,
                                                     torsion[2].index,
                                                     torsion[3].index,
                                                    (t[0],
                                                     t[1],
                                                     t[2]))
            
    def addHarmonicImproperForces(self):
        """
        Creates an openmm.CustomTorsionForce() object with the torsions 
        and parameters setted up in the "impropers" attribute. The custom 
        torsion force is initilized with the formula: 
        
        energy = 0.5*k*dtheta_torus^2
        dtheta_torus = dtheta - floor(dtheta/(2*pi)+0.5)*(2*pi)
        dtheta = theta - theta0
        pi = np.pi
        
        The force object is stored at the "harmonicImproperForce" attribute.
        
        The force parameters must be contained in self.torsions as follows:
        self.torsions is a dictionary:
            - The keys are tuples for four atoms in self.topology.atoms attribute.
            - The values are a tuple of parameters in the following order:
                first  -> theta0 (quantity)
                second -> k (float)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        energy_function = '0.5*k*dtheta_torus^2;'
        energy_function += 'dtheta_torus = dtheta - floor(dtheta/(2*pi)+0.5)*(2*pi);'
        energy_function += 'dtheta = theta - theta0;'
        energy_function += 'pi = %f;' % np.pi
        
        self.harmonicImproperForce = openmm.CustomTorsionForce(energy_function)
        self.harmonicImproperForce.addPerTorsionParameter('theta0')
        self.harmonicImproperForce.addPerTorsionParameter('k')
        for improper in self.impropers:
            self.harmonicImproperForce.addTorsion(improper[0].index,
                                                 improper[1].index,
                                                 improper[2].index,
                                                 improper[3].index,
                                                 (self.impropers[improper][0],
                                                 self.impropers[improper][1]))
            
    def addHarmonicPlanarForces(self):
        """
        Creates an openmm.CustomTorsionForce() object with the torsions 
        and parameters setted up in the "planars" attribute. The custom 
        torsion force is initilized with the formula: 
        
        energy = 0.5*k*dtheta_torus^2
        dtheta_torus = dtheta - floor(dtheta/(2*pi)+0.5)*(2*pi)
        dtheta = theta - theta0
        pi = np.pi
        
        The force object is stored at the "harmonicPlanarForce" attribute.
        
        The force parameters must be contained in self.torsions as follows:
        self.torsions is a dictionary:
            - The keys are tuples for four atoms in self.topology.atoms attribute.
            - The values are a tuple of parameters in the following order:
                first  -> theta0 (quantity)
                second -> k (float)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        energy_function = '0.5*k*dtheta_torus^2;'
        energy_function += 'dtheta_torus = dtheta - floor(dtheta/(2*pi)+0.5)*(2*pi);'
        energy_function += 'dtheta = theta - theta0;'
        energy_function += 'pi = %f;' % np.pi
        
        self.harmonicPlanarForce = openmm.CustomTorsionForce(energy_function)
        self.harmonicPlanarForce.addPerTorsionParameter('theta0')
        self.harmonicPlanarForce.addPerTorsionParameter('k')
        for planar in self.planars:
            self.harmonicPlanarForce.addTorsion(planar[0].index,
                                              planar[1].index,
                                              planar[2].index,
                                              planar[3].index,
                                              (self.planars[planar][0],
                                              self.planars[planar][1]))
            
    def addLJ12_6ContactForces(self):
        """
        Creates an openmm.CustomBondForce() object with the bonds and 
        parameters setted up in the "contacts" attribute. The custom 
        bond force is initilized with the formula:
        
        energy = epsilon*((sigma/r)^12-2*(sigma/r)^6)
        
        The force object is stored at the "lj12_6contactForce" attribute.
        
        The force parameters must be contained in self.contacts as follows:
        self.contacts is a dictionary:
            - The keys are tuples for two atoms in self.topology.atoms attribute.
            - The values are a tuple of parameters in the following order:
                first  -> sigma (quantity)
                second -> epsilon (float)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Collect the number of parametes given for each contact in the contact set
        parameters_length = set([len(p) for p in self.contacts.values()])
        
        #If any contact has 2 parameters create a LJ potential
        if 2 in parameters_length:
            energy_function = 'epsilon*((sigma/r)^12-2*(sigma/r)^6)'
            self.lj12_6contactForce = CustomBondForce(energy_function)
            self.lj12_6contactForce.addPerBondParameter('sigma')
            self.lj12_6contactForce.addPerBondParameter('epsilon')
            for contact in self.contacts:
                if len(self.contacts[contact]) == 2:
                    self.lj12_6contactForce.addBond(contact[0].index,
                                                    contact[1].index,
                                                    (self.contacts[contact][0],
                                                    self.contacts[contact][1]))

    def addLJ12_10ContactForces(self):
        """
        Creates an openmm.CustomBondForce() object with the bonds and 
        parameters setted up in the "contacts" attribute. The custom 
        bond force is initilized with the formula:
        
        energy = epsilon*(5*(sigma/r)^12-6*(sigma/r)^10)
        
        The force object is stored at the "lj12_10contactForce" attribute.
        
        The force parameters must be contained in self.contacts as follows:
        self.contacts is a dictionary:
            - The keys are tuples for two atoms in self.topology.atoms attribute.
            - The values are a tuple of parameters in the following order:
                first  -> sigma (quantity)
                second -> epsilon (float)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Collect the number of parametes given for each contact in the contact set
        parameters_length = set([len(p) for p in self.contacts.values()])
        
        #If any contact has 2 parameters create a LJ potential
        if 2 in parameters_length:
            energy_function = 'epsilon*(5*(sigma/r)^12-6*(sigma/r)^10)'
            self.lj12_10contactForce = CustomBondForce(energy_function)
            self.lj12_10contactForce.addPerBondParameter('sigma')
            self.lj12_10contactForce.addPerBondParameter('epsilon')
            for contact in self.contacts:
                if len(self.contacts[contact]) == 2:
                    self.lj12_10contactForce.addBond(contact[0].index,
                                                     contact[1].index,
                                                     (self.contacts[contact][0],
                                                     self.contacts[contact][1]))
                    
    def addLJ12_10_6ContactForces(self):
        """
        Creates an openmm.CustomBondForce() object with the bonds and 
        parameters setted up in the "contacts" attribute. The custom 
        bond force is initilized with the formula:
        
        energy = epsilon*(13*(sigma/r)^12-18*(sigma/r)^10+4*(sigma/r)^6)
        
        The force object is stored at the "lj12_10_6contactForce" attribute.
        
        The force parameters must be contained in self.contacts as follows:
        self.contacts is an ordered dictionary:
            - The keys are tuples for two atoms in self.topology.atoms attribute.
            - The values are a tuple of parameters in the following order:
                first  -> sigma (quantity)
                second -> epsilon (float)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Collect the number of parametes given for each contact in the contact set
        parameters_length = set([len(p) for p in self.contacts.values()])
        
        #If any contact has 2 parameters create a LJ potential
        if 2 in parameters_length:
            energy_function = 'epsilon*(13*(sigma/r)^12-18*(sigma/r)^10+4*(sigma/r)^6)'
            self.lj12_10_6contactForce = CustomBondForce(energy_function)
            self.lj12_10_6contactForce.addPerBondParameter('sigma')
            self.lj12_10_6contactForce.addPerBondParameter('epsilon')
            for contact in self.contacts:
                if len(self.contacts[contact]) == 2:
                    self.lj12_10_6contactForce.addBond(contact[0].index,
                                                     contact[1].index,
                                                     (self.contacts[contact][0],
                                                     self.contacts[contact][1]))
            
    def addGaussianContactForces(self):
        """
        Creates an openmm.CustomBondForce() object with the equilibrium 
        distances and parameters setted up in the "contacts" attribute. 
        Two forces are created according to the number of parameters given:
        
        If three parameters are given the custom bond force is initilized 
        with the formula:
        
        energy = epsilon*((1+(r0/r)^12)*(1-exp(-25*log(2)*(r-r0)^2/(r0*r0)))-1)
        
        The force object is stored at the "singleGaussianContactForce" attribute.
        
        The force parameters must be contained in self.contacts as follows:
        self.contacts is a dictionary:
            - The keys are tuples for tww atoms in self.topology.atoms attribute.
            - The values are a tuple of parameters in the following order:
                first  -> rex (float)
                second -> r0 (quantity)
                third -> epsilon (float)
        
        If four parameters are given the custom bond force is initilized 
        with the formula:
        
        energy = epsilon*((1+(r0/r)^12)*(1-exp(-25*log(2)*(r-r0)^2/(r0*r0)))*
                (1-exp(-25*log(2)*(r-r1)^2/(r1*r1)))-1)
        
        The force object is stored at the "doubleGaussianContactForce" attribute.
        
        The force parameters must be contained in self.contacts as follows:
        self.contacts is a dictionary:
            - The keys are tuples for two atoms in self.topology.atoms attribute.
            - The values are a tuple of parameters in the following order:
                first  -> rex (float)
                second -> r0 (quantity)
                second -> r1 (quantity)
                third -> epsilon (float)
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Collect the number of parametes given for each contact in the contact set
        parameters_length = set([len(p) for p in self.contacts.values()])
        
        #If any contact has 3 parameters create a single minimum gaussian potential
        if 3 in parameters_length:
            energy_function = 'epsilon*((1+(rex/r)^12)*'
            energy_function += '(1-exp(-25*log(2)*(r-r0)^2/(r0*r0)))-1)'
            self.singleGaussianContactForce = CustomBondForce(energy_function)
            self.singleGaussianContactForce.addPerBondParameter('rex')
            self.singleGaussianContactForce.addPerBondParameter('r0')
            self.singleGaussianContactForce.addPerBondParameter('epsilon')
            for contact in self.contacts:
                if len(self.contacts[contact]) == 3:
                    self.singleGaussianContactForce.addBond(contact[0].index,
                                                            contact[1].index,
                                                            (self.contacts[contact][0],
                                                            self.contacts[contact][1],
                                                            self.contacts[contact][2]))
                    
        #If any contact has 4 parameters create a double minimum gaussian potential
        if 4 in parameters_length:
            energy_function = 'epsilon*((1+(rex/r)^12)*'
            energy_function += '(1-exp(-25*log(2)*(r-r0)^2/(r0*r0)))*'
            energy_function += '(1-exp(-25*log(2)*(r-r1)^2/(r1*r1)))-1)'
            self.doubleGaussianContactForce = CustomBondForce(energy_function)
            self.doubleGaussianContactForce.addPerBondParameter('rex')
            self.doubleGaussianContactForce.addPerBondParameter('r0')
            self.doubleGaussianContactForce.addPerBondParameter('r1')
            self.doubleGaussianContactForce.addPerBondParameter('epsilon')
            for contact in self.contacts:
                if len(self.contacts[contact]) == 4:
                    self.doubleGaussianContactForce.addBond(contact[0].index,
                                                            contact[1].index,
                                                            (self.contacts[contact][0],
                                                            self.contacts[contact][1],
                                                            self.contacts[contact][2],
                                                            self.contacts[contact][3]))

            
    def addLJRepulsionForces(self, epsilon=None, sigma=None, cutoff=None, bonded_exclusions_index=3):
        """
        Creates an openmm.CustomNonbondedForce() object with the parameters 
        sigma and epsilon given to this method. The custom non-bonded force
        is initilized with the formula:
        
        energy = 'epsilon*(sigma/r)^12; sigma=0.5*(sigma1+sigma2)'
        
        The method adds exclusions for bonded atoms up until 'bonded_exclusions_index' 
        bond orders (default 3) and also for all the pairs defined in the 'contacts' attribute. 
        The force object is stored at the "ljRepulsionForce" attribute.
        
        Parameters
        ----------
        epsilon : float
            Value of the epsilon constant in the energy function. 
        sigma : float or list 
            Value of the sigma constant (in nm) in the energy function. If float the
            same sigma value is used for every particle. If list a unique 
            parameter is given for each particle.
        cutoff : float
            The cutoff distance (in nm) being used for the nonbonded interactions.
        bonded_exclusions_index : int
            Specified number of bonds to consider when adding non-bonded exclusions.
            
        Returns
        -------
        None
        """
        
        self.rf_epsilon = epsilon
        self.rf_sigma = sigma
        self.rf_cutoff = cutoff
        
        if not isinstance(epsilon, float):
            self.rf_epsilon = float(epsilon)
            
        if not isinstance(sigma, float) and not isinstance(sigma, list):
            try:
                self.rf_sigma = float(sigma)
            except:
                self.rf_sigma = list(sigma)
                
        if not isinstance(cutoff, float):
            self.rf_cutoff = float(cutoff)
        
        energy_function = 'epsilon*(sigma/r)^12; sigma=0.5*(sigma1+sigma2)'
        self.ljRepulsionForce = CustomNonbondedForce(energy_function)
        self.ljRepulsionForce.addGlobalParameter('epsilon', self.rf_epsilon)
        self.ljRepulsionForce.addPerParticleParameter('sigma')
        
        if isinstance(self.rf_sigma, float):
            for atom in self.atoms:
                self.ljRepulsionForce.addParticle((self.rf_sigma,))
                
        elif isinstance(self.rf_sigma, list):
            assert self.n_atoms == len(self.rf_sigma)
            for i,atom in enumerate(self.atoms):
                self.ljRepulsionForce.addParticle((self.rf_sigma[i],))
            
        #Add bonded exclusions
        bonded_exclusions = [(i[0].index, i[1].index) for i in self.bonds]
        self.ljRepulsionForce.createExclusionsFromBonds(bonded_exclusions, bonded_exclusions_index)
        self.exclusions = []
        for i in range(self.ljRepulsionForce.getNumExclusions()):
            pair = self.ljRepulsionForce.getExclusionParticles(i)
            self.exclusions.append(pair)
            
        #Add contacts exclusions
        contact_exclusions = [(i[0].index, i[1].index) for i in self.contacts]
        #Remove already excluded pairs (from bonded exclusions)
        #to avoid 'Multiple exclusions are specified' error
        contact_exclusions = [pair for pair in contact_exclusions if pair in self.exclusions]
        self.ljRepulsionForce.createExclusionsFromBonds(contact_exclusions, 1)
        
        #Add all exclusion pairs
        self.exclusions = []
        for i in range(self.ljRepulsionForce.getNumExclusions()):
            pair = self.ljRepulsionForce.getExclusionParticles(i)
            self.exclusions.append(pair)
        
        self.ljRepulsionForce.setCutoffDistance(self.rf_cutoff)
        
    def groupTorsionsbyBBAndSC(self):
        """
        Classifies the torsions into backbone or sidechain, based on the atoms 
        present in each torsion definition.
        
        The classification is stored in the 'torsions_type' attribute.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        
        #Classify atoms in backbone or sidechain atoms
        uniqueAtoms = set()
        for atom in self.topology.atoms():
            uniqueAtoms.add(atom.name)
        self.backboneAtoms = set(['N', 'CA', 'C', 'O', 'OXT, H, 1H, 2H, 3H'])
        self.sidechainAtoms = uniqueAtoms - self.backboneAtoms
        
        #Group torsions by their middle bond atoms
        for torsion in self.torsions:
            if (torsion[1], torsion[2]) not in self.torsions_group:
                    self.torsions_group[(torsion[1], torsion[2])] = 1
            else:
                    self.torsions_group[(torsion[1], torsion[2])] += 1
            if (torsion[2], torsion[1]) not in self.torsions_group:
                    self.torsions_group[(torsion[2], torsion[1])] = 1
            else:
                    self.torsions_group[(torsion[2], torsion[1])] += 1
        
        #Classify torsions in backbone or sidechain torsion according to middle atoms
        for torsion in self.torsions:
            self.torsions_type[(torsion[1], torsion[2])] = 'BB'
            self.torsions_type[(torsion[2], torsion[1])] = 'BB'
            for atom in torsion[1:3]:
                if atom.name in self.sidechainAtoms:
                    self.torsions_type[(torsion[1], torsion[2])] = 'SC'
                    self.torsions_type[(torsion[2], torsion[1])] = 'SC'
                    break
        
        #Count the number of backbone and sidechin torsions with unique middle bond atoms.
        self.n_torsions = int(sum([1 for t in self.torsions_group])/4.0)
    
    ## Functions for generating default parameters for AA default forcefield ##
    
    def getAATorsionParameters(self, group_by_bb_and_sc=True):
        """
        Calculates the torsion force field parameters by distributing the 
        total torsional energy equally among all the torsions. Optionally,
        uses torsional classifications, in backbone or side-chain torsions,
        to group the energy partition. In the last case considers that BB
        torsion energies are twice the value of sidechain torsion energies.
        
        The energy partioned values are stored in the 'energy_constant' 
        attribute. It contains the entries 'torsion' for the equally partitioned
        torsion energies, or 'BB' and 'SC' for group-partioned torsion energies.
        
        Parameters
        ----------
        group_by_bb_and_sc : boolean
            Wheter to partion the torsional energy in backbone and side-chain 
            groups.
            
        Returns
        -------
        torsion_parameters : list
            List of the torsion force constant parameters for each torsion
            in the SBM system.
        """
        
        if group_by_bb_and_sc:
            #Group torsions by BB and SC
            self.groupTorsionsbyBBAndSC()
            #Calculate energy partitions
            self.energy_constant['torsion'] =  self.n_atoms/3.0/(self.n_torsions)
            self.energy_constant['SC'] = self.energy_constant['torsion']/3.0
            self.energy_constant['BB'] = 2.0*self.energy_constant['torsion']/3.0
            
            torsion_parameters = []
            for torsion in self.torsions:
                middle_bond = (torsion[1],torsion[2])
                k = self.energy_constant[self.torsions_type[middle_bond]] / self.torsions_group[middle_bond]
                torsion_parameters.append(k)
                
        else:
            self.energy_constant['torsion'] =  self.n_atoms/3.0/(self.n_torsions)
            torsion_parameters = self.energy_constant['torsion']
            
        return torsion_parameters
    
    def getAANativeContactParameters(self):
        """
        Calculates the contact force field parameters by distributing the 
        total native contact energy equally among all the contacts.
        
        The energy values are stored in the 'energy_constant' attribute, 
        which contains the entry 'C' for the partitioned energy values.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        contact_parameters : list
            List of the native-contact force constant parameters for each 
            native contact in the SBM system.
        """
        
        self.energy_constant['C'] = 2.0*self.n_atoms/3.0/self.n_contacts
        contact_parameters = self.energy_constant['C']
        return contact_parameters
    
    ## Functions for creating OpenMM system object ##
    
    def createSystemObject(self, check_bond_distances=True, minimise=False, force_threshold=10, bond_threshold=0.24):
        """
        Creates an openmm.System() object using the force field parameters
        given to the SBM 'system' class. It adds particles, forces and 
        creates a force group for each force object. Optionally the method 
        can check for large bond distances (default) and minimise the atomic 
        positions if large forces are found in any atom (default False).
        
        Parameters
        ----------
        minimise : boolean
            Wheter to minimise the system if large forces are found.
        check_bond_distances : boolean
            Wheter to check for large bond distances.
        force_threshold : float
            Treshold to check for large forces.
        bond_threshold : float
            Treshold to check for large bond distances.
            
        Returns
        -------
        None
        """
        
        if check_bond_distances:
            #Check for large bond_distances
            self.checkBondDistances(threshold=bond_threshold)
        
        #Add particles to system
        self.addParticles()
        
        #Add created forces into the system
        self.addSystemForces()
        
        #Create force group for each added force
        for i,name in enumerate(self.forceGroups):
            self.forceGroups[name].setForceGroup(i)
        
        #Check for high forces in atoms and minimise the system if necessary
        self.checkLargeForces(minimise=minimise, threshold=force_threshold)
        
    def checkBondDistances(self, threshold=0.24):
        """
        Searches for large bond distances for the atom pairs defined in 
        the 'bonds' attribute. It raises an error when large bonds are found.
        
        Parameters
        ----------
        threshold : float
            Treshold to check for large bond distances.
            
        Returns
        -------
        None
        """
        
        if isinstance(threshold, float):
            threshold = threshold * unit.nanometer
        
        for b in self.bonds:
            if self.bonds[b][0] >= threshold:
                raise ValueError('The bond distance between atoms '+                str(b[0].index+1)+' and '+str(b[1].index+1)+' is larger '+
                'than '+str(threshold)+' nanometers. Please check your input structure.')
                
    def checkLargeForces(self, threshold=1, minimise=False):
        """
        Prints the SBM system energies of the input configuration of the 
        system. It optionally checks for large forces acting upon all 
        particles in the SBM system and iteratively minimises the system
        configuration until no forces larger than a threshold are found.
        
        Parameters
        ----------
        threshold : float
            Treshold to check for large forces.
        minimise : float
            Whether to iteratively minimise the system until all forces are lower or equal to 
            the thresshold value.
            
        Returns
        -------
        None
        """
        
        minimised = False
        
        #Define test simulation to extract forces
        integrator = LangevinIntegrator(1*unit.kelvin, 1/unit.picosecond, 0.0005*unit.picoseconds)
        simulation = Simulation(self.topology, self.system, integrator)
        simulation.context.setPositions(self.positions)
        state = simulation.context.getState(getForces=True, getEnergy=True)
        
        #Print initial state of the system
        print('The Potential Energy of the system is : %s' % state.getPotentialEnergy())
        for i,n in enumerate(self.forceGroups):
            energy = simulation.context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            print('The '+n.replace('Force','Energy')+' is: '+str(energy)+' kj/mol')
        print('')
            
        if minimise:
            #Find if there is an acting force larger than thresshold
            #Minimise the system until forces have converged
            forces = [np.linalg.norm([f[0]._value,f[1]._value,f[2]._value]) for f in  state.getForces()]
            prev_force = None
            tolerance = 10
            
            while np.max(forces) > threshold:

                #Write atom with largest force if not reported before 
                if np.max(forces) != prev_force:
                    atom = self.atoms[np.argmax(forces)]
                    residue = atom.residue
                    print('Large force %.3f kj/(mol nm) found in:' % np.max(forces))
                    print('Atom: %s' % atom.index)
                    print('Name: %s' % atom.name)
                    print('Resiue: %s %s' % (residue.name,residue.index))
                    print('Minimising system with energy tolerance of %.1f kj/mol' % tolerance)
                    print('')

                simulation.minimizeEnergy(tolerance=tolerance*unit.kilojoule/unit.mole)
                minimised = True
                state = simulation.context.getState(getForces=True)
                prev_force = np.max(forces)
                forces = [np.linalg.norm([f.x,f.y,f.z]) for f in  state.getForces()]
                if tolerance > 1:
                    tolerance -= 1
                elif tolerance > 0.1:
                    tolerance -= 0.1
                elif tolerance == 0.1:
                    raise ValueError('The system could no be minimised at the requested convergence\n'+
                                     'Try to increase the force threshold value to achieve convergence.')
                    

            state = simulation.context.getState(getPositions=True, getEnergy=True) 
            print('After minimisation:')
            print('The Potential Energy of the system is : %s' % state.getPotentialEnergy())
            for i,n in enumerate(self.forceGroups):
                energy = simulation.context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                print('The '+n.replace('Force','Energy')+' is: '+str(energy)+' kj/mol')
            print('All forces are less than %.2f kj/mol/nm' % threshold)
            print('Saving minimised positions')
            print('')
            self.positions = state.getPositions()
                
    def addParticles(self):
        """
        Add a particle to the sbmOpenMM system for each atom in it. The mass 
        of each particle is set up with the values in the 'particles_mass' 
        attribute.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        #Set same mass for each atom
        if isinstance(self.particles_mass, float):
            for atom in self.atoms:
                self.system.addParticle(self.particles_mass)
                
        #Set unique masses for each atom
        if isinstance(self.particles_mass, list):
            assert len(self.particles_mass) == len(self.atoms)
            for i in range(len(self.particles_mass)):
                self.system.addParticle(self.particles_mass[i])
    
    def addSystemForces(self):
        """
        Adds generated forces to the sbmOpenMM system, also adding 
        a force group to the 'forceGroups' attribute dictionary.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        if self.harmonicBondForce != None:
            self.system.addForce(self.harmonicBondForce)
            self.forceGroups['Harmonic Bond Energy'] = self.harmonicBondForce
            
        if self.harmonicAngleForce != None:
            self.system.addForce(self.harmonicAngleForce)
            self.forceGroups['Harmonic Angle Energy'] = self.harmonicAngleForce
            
        if self.periodicTorsionForce != None:
            self.system.addForce(self.periodicTorsionForce)
            self.forceGroups['Periodic Torsion Energy'] = self.periodicTorsionForce
            
        if self.generalPeriodicTorsionForce != None:
            self.system.addForce(self.generalPeriodicTorsionForce)
            self.forceGroups['General Periodic Torsion Energy'] = self.generalPeriodicTorsionForce
            
        if self.harmonicImproperForce != None:
            self.system.addForce(self.harmonicImproperForce)
            self.forceGroups['Harmonic Improper Energy'] = self.harmonicImproperForce
            
        if self.harmonicPlanarForce != None:
            self.system.addForce(self.harmonicPlanarForce)
            self.forceGroups['Harmonic Planar Energy'] = self.harmonicPlanarForce
            
        if self.lj12_6contactForce != None:
            self.system.addForce(self.lj12_6contactForce)
            self.forceGroups['LJ 12-6 Contact Energy'] = self.lj12_6contactForce
            
        if self.lj12_10contactForce != None:
            self.system.addForce(self.lj12_10contactForce)
            self.forceGroups['LJ 12-10 Contact Energy'] = self.lj12_10contactForce
            
        if self.lj12_10_6contactForce != None:
            self.system.addForce(self.lj12_10_6contactForce)
            self.forceGroups['LJ 12-10-6 Contact Energy'] = self.lj12_10_6contactForce
            
        if self.singleGaussianContactForce != None:
            self.system.addForce(self.singleGaussianContactForce)
            self.forceGroups['Single Gaussian Contact Energy'] = self.singleGaussianContactForce
            
        if self.doubleGaussianContactForce != None:
            self.system.addForce(self.doubleGaussianContactForce)
            self.forceGroups['Double Gaussian Contact Energy'] = self.doubleGaussianContactForce
            
        if self.ljRepulsionForce != None:
            self.system.addForce(self.ljRepulsionForce)
            self.forceGroups['LJ 12 Repulsion Energy'] = self.ljRepulsionForce
    
    def dumpPdb(self, output_file):
        """
        Writes a PDB file containing the currently defined abmOpenMM system atoms.
        
        Parameters
        ----------
        output_file : string
            name of the PDB output file. 
            
        Returns
        -------
        None
        """
        
        self.pdb.writeFile(self.topology, self.positions, file=open(output_file, 'w'))
        
    def dumpForceFieldData(self, output_file):
        """
        Writes a file containing the current forcefield parameters in the 
        sbmOpenMM system.
        
        Parameters
        ----------
        output_file : string
            name of the output file. 
            
        Returns
        -------
        None
        """
        
        with open(output_file, 'w') as ff:
            
            ff.write('#### SBM Force Field Parameters ####\n')
            ff.write('\n')
            
            if self.bonds != OrderedDict():
                ff.write('[bonds]\n')
                ff.write('# %2s %3s %-6s %-s \t\t self.bonds[bond][1] %12s %12s\n' % ('at1', 'at2', 'b0', 'k', 'at1_name', 'at2_name'))
                for bond in self.bonds:
                    atom1 = bond[0]
                    atom2 = bond[1]
                    res1 = atom1.residue
                    res2 = atom2.residue
                    ff.write('%4s %4s %.4f %s \t# %12s %12s\n' % (atom1.index+1,
                                                                  atom2.index+1,
                                                                  self.bonds[bond][0]._value,
                                                                  self.bonds[bond][1],
                                                                  atom1.name+'_'+res1.name+'_'+str(res1.index+1),
                                                                  atom2.name+'_'+res2.name+'_'+str(res2.index+1)))
                        
            if self.angles != OrderedDict():
                ff.write('\n')
                ff.write('[angles]\n')
                ff.write('# %2s %3s %4s %-6s %-s \t  %12s %12s %12s\n' % ('at1', 'at2', 'at3', 'a0', 'k', 'at1_name', 'at2_name', 'at3_name'))
                for angle in self.angles:
                    atom1 = angle[0]
                    atom2 = angle[1]
                    atom3 = angle[2]
                    res1 = atom1.residue
                    res2 = atom2.residue
                    res3 = atom3.residue
                    ff.write('%4s %4s %4s %.4f %s \t# %12s %12s %12s\n' % (atom1.index+1,
                                                                           atom2.index+1,
                                                                           atom3.index+1,
                                                                           self.angles[angle][0]._value,
                                                                           self.angles[angle][1],
                                                                           atom1.name+'_'+res1.name+'_'+str(res1.index+1),
                                                                           atom2.name+'_'+res2.name+'_'+str(res2.index+1),
                                                                           atom3.name+'_'+res3.name+'_'+str(res3.index+1)))

            if self.torsions != OrderedDict():
                ff.write('\n')
                ff.write('[torsions]\n')
                ff.write('# %2s %3s %4s %4s %-6s %-s \t\t  %12s %12s %12s %12s\n' % ('at1', 'at2', 'at3', 'at4', 't0', 'k', 'at1_name', 'at2_name', 'at3_name', 'at4_name'))
                for torsion in self.torsions:
                    atom1 = torsion[0]
                    atom2 = torsion[1]
                    atom3 = torsion[2]
                    atom4 = torsion[3]
                    res1 = atom1.residue
                    res2 = atom2.residue
                    res3 = atom3.residue
                    res4 = atom3.residue
                    ff.write('%4s %4s %4s %4s %7.4f %.4f \t# %12s %12s %12s %12s\n' % (atom1.index+1,
                                                                                       atom2.index+1,
                                                                                       atom3.index+1,
                                                                                       atom4.index+1,
                                                                                       self.torsions[torsion][0]._value,
                                                                                       self.torsions[torsion][1],
                                                                                       atom1.name+'_'+res1.name+'_'+str(res1.index+1),
                                                                                       atom2.name+'_'+res2.name+'_'+str(res2.index+1),
                                                                                       atom3.name+'_'+res3.name+'_'+str(res3.index+1),
                                                                                       atom4.name+'_'+res4.name+'_'+str(res4.index+1)))
                        
            if self.impropers != OrderedDict():
                ff.write('\n')
                ff.write('[impropers]\n')
                ff.write('# %2s %3s %4s %4s %-6s %-s \t\t  %12s %12s %12s %12s\n' % ('at1', 'at2', 'at3', 'at4', 'i0', 'k', 'at1_name', 'at2_name', 'at3_name', 'at4_name'))
                for improper in self.impropers:
                    atom1 = improper[0]
                    atom2 = improper[1]
                    atom3 = improper[2]
                    atom4 = improper[3]
                    res1 = atom1.residue
                    res2 = atom2.residue
                    res3 = atom3.residue
                    res4 = atom3.residue
                    ff.write('%4s %4s %4s %4s %7.4f %.4f \t# %12s %12s %12s %12s\n' % (atom1.index+1,
                                                                                       atom2.index+1,
                                                                                       atom3.index+1,
                                                                                       atom4.index+1,
                                                                                       self.impropers[improper][0]._value,
                                                                                       self.impropers[improper][1],
                                                                                       atom1.name+'_'+res1.name+'_'+str(res1.index+1),
                                                                                       atom2.name+'_'+res2.name+'_'+str(res2.index+1),
                                                                                       atom3.name+'_'+res3.name+'_'+str(res3.index+1),
                                                                                       atom4.name+'_'+res4.name+'_'+str(res4.index+1)))
                        
            if self.planars != OrderedDict():
                ff.write('\n')
                ff.write('[planars]\n')
                ff.write('# %2s %3s %4s %4s %-6s %-s \t\t  %12s %12s %12s %12s\n' % ('at1', 'at2', 'at3', 'at4', 'p0', 'k', 'at1_name', 'at2_name', 'at3_name', 'at4_name'))
                for planar in self.planars:
                    atom1 = planar[0]
                    atom2 = planar[1]
                    atom3 = planar[2]
                    atom4 = planar[3]
                    res1 = atom1.residue
                    res2 = atom2.residue
                    res3 = atom3.residue
                    res4 = atom3.residue
                    ff.write('%4s %4s %4s %4s %7.4f %.4f \t# %12s %12s %12s %12s\n' % (atom1.index+1,
                                                                                       atom2.index+1,
                                                                                       atom3.index+1,
                                                                                       atom4.index+1,
                                                                                       self.planars[planar][0]._value,
                                                                                       self.planars[planar][1],
                                                                                       atom1.name+'_'+res1.name+'_'+str(res1.index+1),
                                                                                       atom2.name+'_'+res2.name+'_'+str(res2.index+1),
                                                                                       atom3.name+'_'+res3.name+'_'+str(res3.index+1),
                                                                                       atom4.name+'_'+res4.name+'_'+str(res4.index+1)))
                        
            if self.contacts != OrderedDict():
                #Collect the number of parametes given for each contact in the contact set
                parameters_length = set([len(p) for p in self.contacts.values()])
                
                #Write two parameters contacts
                if 2 in parameters_length:
                    ff.write('\n')
                    ff.write('[contacts]\n')
                    ff.write('#Two parameters contacts\n')
                    ff.write('# %2s %3s %-6s %-s \t  %12s %12s\n' % ('at1', 'at2', 'c0', 'epsilon', 'at1_name', 'at2_name'))
                    for contact in self.contacts:
                        if len(self.contacts[contact]) == 2:
                            atom1 = contact[0]
                            atom2 = contact[1]
                            res1 = atom1.residue
                            res2 = atom2.residue
                            ff.write('%4s %4s %.4f %.4f \t# %12s %12s\n' % (atom1.index+1,
                                                                          atom2.index+1,
                                                                          self.contacts[contact][0]._value,
                                                                          self.contacts[contact][1],
                                                                          atom1.name+'_'+res1.name+'_'+str(res1.index+1),
                                                                          atom2.name+'_'+res2.name+'_'+str(res2.index+1)))
                            
                #Write three parameters contacts
                if 3 in parameters_length:
                    ff.write('#Three parameters contacts\n')
                    ff.write('# %2s %3s %-6s %-6s %-s \t  %12s %12s\n' % ('at1', 'at2', 'excvol', 'c0', 'epsilon', 'at1_name', 'at2_name'))
                    for contact in self.contacts:
                        if len(self.contacts[contact]) == 3:
                            atom1 = contact[0]
                            atom2 = contact[1]
                            res1 = atom1.residue
                            res2 = atom2.residue
                            ff.write('%4s %4s %.4f %.4f %.4f \t# %12s %12s\n' % (atom1.index+1,
                                                                          atom2.index+1,
                                                                          self.contacts[contact][0]._value,
                                                                          self.contacts[contact][1]._value,
                                                                          self.contacts[contact][2],
                                                                          atom1.name+'_'+res1.name+'_'+str(res1.index+1),
                                                                          atom2.name+'_'+res2.name+'_'+str(res2.index+1)))
                #Write four parameters contacts
                if 4 in parameters_length:
                    ff.write('#Four parameters contacts\n')
                    ff.write('# %2s %3s %-6s %-6s %-6s %-s \t  %12s %12s\n' % ('at1', 'at2', 'excvol', 'c0', 'c1', 'epsilon', 'at1_name', 'at2_name'))
                    for contact in self.contacts:
                        if len(self.contacts[contact]) == 4:
                            atom1 = contact[0]
                            atom2 = contact[1]
                            res1 = atom1.residue
                            res2 = atom2.residue
                            ff.write('%4s %4s %.4f %.4f %.4f %.4f \t# %12s %12s\n' % (atom1.index+1,
                                                                          atom2.index+1,
                                                                          self.contacts[contact][0]._value,
                                                                          self.contacts[contact][1]._value,
                                                                          self.contacts[contact][2]._value,
                                                                          self.contacts[contact][3],
                                                                          atom1.name+'_'+res1.name+'_'+str(res1.index+1),
                                                                          atom2.name+'_'+res2.name+'_'+str(res2.index+1)))
            
    def loadForcefieldFromFile(self, forcefield_file):
        """
        Loads force field parameters from a force field file written by the 
        dumpForceFieldData() method into the sbmOpenMM system.
        
        Parameters
        ----------
        forcefield_file : string
            path to the force field file.
            
        Returns
        -------
        None
        """
        
        #Open forcefield file
        with open(forcefield_file, 'r') as ff:
            print('Reading Forcefield parameters from file '+forcefield_file+':')
            print('________________________________________'+'_'*len(forcefield_file))
            
            #Initilise ff-file section booleans to 0
            bonds = False
            angles = False
            torsions = False
            impropers = False
            planars = False
            contacts = False
                        
            #Iterate for all lines in forcefield file
            for i,line in enumerate(ff):
                
                #Ignore comment and empty lines.
                if not line.startswith('#') and line.split() != []:
                    
                    #Turn off all sections when a new is being reading.
                    if line.startswith('['):
                        bonds = False
                        angles = False
                        torsions = False
                        impropers = False
                        planars = False
                        contacts = False
                    else:
                        
                        #Remove comments at the end of line
                        line = line[:line.find('#')]
                        ls = line.split()
                        
                        #Reading [bonds] section
                        if bonds == True:
                            if len(ls) > 4:
                                raise ValueError('More than four parameters given in [bonds] section at line '+str(i)+':\n'+line)
                            if len(ls) < 4:
                                raise ValueError('Less than four parameters given in [bonds] section at line '+str(i)+':\n'+line)
                            at1 = self.atoms[int(ls[0])-1]
                            at2 = self.atoms[int(ls[1])-1]
                            bond_length = float(ls[2])*unit.nanometer
                            k = float(ls[3])
                            self.bonds[(at1,at2)] = (bond_length, k)
                            
                        #Reading [angles] section
                        elif angles == True:
                            if len(ls) > 5:
                                raise ValueError('More than five parameters given in [angles] section at line '+str(i)+':\n'+line)
                            if len(ls) < 5:
                                raise ValueError('Less than five parameters given in [angles] section at line '+str(i)+':\n'+line)
                            at1 = self.atoms[int(ls[0])-1]
                            at2 = self.atoms[int(ls[1])-1]
                            at3 = self.atoms[int(ls[2])-1]
                            angle_length = float(ls[3])*unit.radian
                            k = float(ls[4])
                            self.angles[(at1,at2,at3)] = (angle_length, k)
                            
                        #Reading [torsions] section
                        elif torsions == True:
                            if len(ls) > 6:
                                raise ValueError('More than six parameters given in [torsions] section at line '+str(i)+':\n'+line)
                            if len(ls) < 6:
                                raise ValueError('Less than six parameters given in [torsions] section at line '+str(i)+':\n'+line)
                            at1 = self.atoms[int(ls[0])-1]
                            at2 = self.atoms[int(ls[1])-1]
                            at3 = self.atoms[int(ls[2])-1]
                            at4 = self.atoms[int(ls[3])-1]
                            torsion_length = float(ls[4])*unit.radian
                            k = float(ls[5])
                            self.torsions[(at1,at2,at3,at4)] = (torsion_length, k)
                            
                        #Reading [impropers] section
                        elif impropers == True:
                            if len(ls) > 6:
                                raise ValueError('More than six parameters given in [impropers] section at line '+str(i)+':\n'+line)
                            if len(ls) < 6:
                                raise ValueError('Less than six parameters given in [impropers] section at line '+str(i)+':\n'+line)
                            at1 = self.atoms[int(ls[0])-1]
                            at2 = self.atoms[int(ls[1])-1]
                            at3 = self.atoms[int(ls[2])-1]
                            at4 = self.atoms[int(ls[3])-1]
                            improper_length = float(ls[4])*unit.radian
                            k = float(ls[5])
                            self.impropers[(at1,at2,at3,at4)] = (improper_length, k)
                            
                        #Reading [planars] section
                        elif planars == True:
                            if len(ls) > 6:
                                raise ValueError('More than six parameters given in [planars] section at line '+str(i)+':\n'+line)
                            if len(ls) < 6:
                                raise ValueError('Less than six parameters given in [planars] section at line '+str(i)+':\n'+line)
                            at1 = self.atoms[int(ls[0])-1]
                            at2 = self.atoms[int(ls[1])-1]
                            at3 = self.atoms[int(ls[2])-1]
                            at4 = self.atoms[int(ls[3])-1]
                            improper_length = float(ls[4])*unit.radian
                            k = float(ls[5])
                            self.planars[(at1,at2,at3,at4)] = (improper_length, k)
                        
                        #Reading [contacts] section
                        elif contacts == True:
                            if len(ls) > 4:
                                raise ValueError('More than four parameters given in [contacts] section at line '+str(i)+':\n'+line)
                            if len(ls) < 4:
                                raise ValueError('Less than four parameters given in [contacts] section at line '+str(i)+':\n'+line)
                            at1 = self.atoms[int(ls[0])-1]
                            at2 = self.atoms[int(ls[1])-1]
                            contact_length = float(ls[2])*unit.nanometer
                            k = float(ls[3])
                            self.contacts[(at1,at2)] = (contact_length, k)
                            
                    #Select which section is being reading by changing its boolean to 1
                    if line.startswith('[bonds]'):
                        print('Reading bond parameters')
                        self.bonds = OrderedDict()
                        bonds = True
                    elif line.startswith('[angles]'):
                        print('Reading angle parameters')
                        self.angles = OrderedDict()
                        angles = True
                    elif line.startswith('[torsions]'):
                        print('Reading torsion parameters')
                        self.torsions = OrderedDict()
                        torsions = True
                    elif line.startswith('[impropers]'):
                        print('Reading improper parameters')
                        self.impropers = OrderedDict()
                        impropers = True
                    elif line.startswith('[planars]'):
                        print('Reading planar parameters')
                        self.planars = OrderedDict()
                        planars = True
                    elif line.startswith('[contacts]'):
                        print('Reading native contact parameters')
                        self.contacts = OrderedDict()
                        contacts = True
                        
        print('Done reading Forcefield paramters')
        print('')
    
    def setCAMassPerResidueType(self):
        """
        Sets the masses of the alpha carbon atoms to the average mass
        of its amino acid residue. 
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        aa_masses = {'ALA':  71.0, 'ARG': 157.0, 'ASN': 114.0,
                     'ASP': 114.0, 'CYS': 103.0, 'GLU': 128.0,
                     'GLN': 128.0, 'GLY':  57.0, 'HIS': 138.0,
                     'ILE': 113.0, 'LEU': 113.0, 'LYS': 128.0,
                     'MET': 131.0, 'PHE': 147.0, 'PRO':  97.0, 
                     'SER':  87.0, 'THR': 101.0, 'TRP': 186.0, 
                     'TYR': 163.0, 'VAL':  99.0}
        
        masses = []
        
        for r in self.topology.residues():
            if r.name in aa_masses:
                masses.append(aa_masses[r.name])
                
        self.setParticlesMasses(masses)
        
#TODO:     def setMassPerAtomType(self):
        
    ## User-hidden functions ##
    
    def _setParameters(term, parameters):
        """
        General function to set up or change force field parameters.
        
        Parameters
        ----------
        term : dict
            Dictionary object containing the set of degrees of freedom
            (DOF) to set up attributes to (e.g. bonds attribute)
            
        parameters : integer or float or list
            Value(s) for the specific forcefield parameters. If integer
            or float, sets up the same value for all the DOF in terms. If 
            list, sets a unique parameter for each DOF. 
            
        Returns
        -------
        None
        """
        
        if isinstance(parameters, int):
            parameters = float(parameters)
        
        #Set constant parameter for each item in FF term
        if isinstance(parameters, float):
            for item in term:
                term[item] = term[item][:-1]+(parameters,)
                
        #Set unique parameter for each item in FF term
        if isinstance(parameters, list):
            assert len(parameters) == len(list(term.keys()))
            for i, item in enumerate(term):
                term[item] = term[item][:-1]+(parameters[i],)
    
    def _checkContactFileType(contact_file):
        """
        Function to check the format of the input contact file. 
        
        Parameters
        ----------
        contact_file : string
            Name of the input contact file. 
            
        Returns
        -------
        None
        """
        
        with open(contact_file,'r') as cf:
            lengths = [len(l.split()) for l in cf if not l.startswith('#')]
            lengths = list(set(lengths))
            if len(lengths) != 1:
                return False
            elif lengths[0] == 4:
                return 'smog'
            elif lengths[0] == 2:
                return '2column'
            else:
                return False
