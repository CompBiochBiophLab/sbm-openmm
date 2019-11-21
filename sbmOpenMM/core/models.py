#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from simtk.openmm.app import *
from simtk.openmm import *
from simtk import unit

from collections import OrderedDict
import numpy as np

from .system import system


# In[ ]:


class models:
    """
    A class to hold functions for the automated generation of default SBM models.

    Methods
    -------
    getAllAtomModel(pdb_file, contact_file, kwarg**)
        Creates an all atom SBM system class object with default initialized
        parameters.
    getCAModel(pdb_file, contact_file, kwarg**)
        Creates an alpha-carbon only sbmOpenMM system class object with default 
        initialized parameters.
    getMultiBasinModel(main_model, alternate_configuration, kwarg**)
        Creates a multi basin model from two sbmOpenMM.system class instances.
    """
    
    def getAllAtomModel(pdb_file, contact_file, 
                        default_parameters=True, 
                        default_forces=True, 
                        group_by_bb_and_sc=True,
                        create_system=True,
                        minimise=False):
        """
        Initialises a default full-heavy-atom sbmOpenMM system class from a PDB file and a contact
        file defining the native contacts for the model. The system creation steps are:
        
        1) Add the geometrical parameters for the model.
        2) Add the default force field parameters for the model.
        3) Create the default force objects.
        4) Create the OpenMM system class.
        
        The method can be used to generate an initialized sbmOpenMM system class, that only 
        contains the geometrical parameters, by passing the option default_parameters as False.
        This is useful to store the geometrical values of bonds, angles, dihedrals, etc. in 
        order to add custom parameters and forces.
        
        The method can also be created without the initialisation of the forces classes by 
        setting default_forces to False. This allows to load the default forcefield parameters 
        and to modified them before creating the OpenMM force objects.
        
        Finally, the method can be stopped before creating the OpenMM system class using create_system=False.
        
        Parameters
        ----------
        pdb_file : string
            Path to the input PDB file.
        contact_file : string
            Path to the input native contact file. The file can be an output 
            from the SMOG program (4 column contact file) or a two column file 
            defining the atoms to be paired.
        default_parameters : boolean (True)
            Whether to add default SBM All Atom forcefield parameters to the model.
        default_forces : boolean (True)
            Whether to initilize default SBM All Atom force objects. Set to False if
            the parameters will be different from the default ones.
        group_by_bb_and_sc : boolean (True)
            Wether to classify the torsions into backbone and side-chain to partition 
            the torsional energy for each torsion in the system.
        create_system : boolean (True)
            If True the function will call the createSystemObject() method
            to create an OpenMM system object. If modifications to the default 
            forcefield are necessary this option should be given False.
        minimise : boolean (False)
            Whether to minimise the system (with default options) if large
            forces are found.
        
        Returns
        -------
        sbm : sbmOpenMM.system
            Initializes a sbmOpenMM.system class with default options for 
            defining an All Atom SBM force field.
        """
        
        print('Generating AA SBM for PDB file '+pdb_file)
        print('')
        
        #Set up geometric parameters of the model
        print('Setting up geometrical parameters:')
        print('_________________________________')
        
        sbm = system(pdb_file)
        print('Removing hydrogens from topology')
        sbm.removeHydrogens()
        sbm.getAtoms()
        print('Added '+str(sbm.n_atoms)+' atoms')
        sbm.getBonds()
        print('Added '+str(sbm.n_bonds)+' bonds')
        sbm.getAngles()
        print('Added '+str(sbm.n_angles)+' angles')
        sbm.getProperTorsions()
        print('Added '+str(sbm.n_torsions)+' torsions')
        sbm.getImpropers()
        print('Added '+str(sbm.n_impropers)+' impropers')
        sbm.getPlanars()
        print('Added '+str(sbm.n_planars)+' planars')
        if contact_file != None:
            print('Reading contacts from contact file: '+contact_file)
            sbm.readContactFile(contact_file)
        print('Added '+str(sbm.n_contacts)+' native contacts')
        print('')
        
        #Add default parameters to each interaction term
        if default_parameters:
            print('Setting up default forcefield parameters:')
            print('________________________________________')
            print('Adding default bond parameters:')
            sbm.setBondParameters(10000.0)
            print('Adding default angle parameters:')
            sbm.setAngleParameters(80.0)
            #Get default torsion parameters
            if group_by_bb_and_sc:
                print('grouping torsions in backbone and side-chain groups:')
            torsion_parameters = sbm.getAATorsionParameters(group_by_bb_and_sc=group_by_bb_and_sc)
            print('Adding default torsion parameters:')
            sbm.setProperTorsionParameters(torsion_parameters)
            print('Adding default improper parameters:')
            sbm.setImproperParameters(10.0)
            print('Adding default planar parameters:')
            sbm.setPlanarParameters(20.0)
            #Get default contact parameters
            print('Adding default contact parameters:')
            contact_parameters = sbm.getAANativeContactParameters()
            sbm.setNativeContactParameters(contact_parameters)
            print('')

        #Create default system force objects
        if default_parameters and default_forces:
            print('Adding Forces:')
            print('_____________')
            print('Adding Harmonic Bond Forces')
            sbm.addHarmonicBondForces()
            print('Adding Harmonic Angle Forces')
            sbm.addHarmonicAngleForces()
            print('Adding Periodic Torsion Forces')
            sbm.addPeriodicTorsionForces()
            print('Adding Harmonic Improper Forces')
            sbm.addHarmonicImproperForces()
            print('Adding Periodic Planar Forces')
            sbm.addHarmonicPlanarForces()
            print('Adding Lennard Jones 12-6 Forces to native contacts')
            sbm.addLJ12_6ContactForces()
            print('Adding Lennard Jones 12 non-bonded Forces')
            sbm.addLJRepulsionForces(sigma=0.25, epsilon=0.1, cutoff=1.5)
            print('')

        #Generate the system object and add previously generated forces
        if default_parameters and default_forces and create_system:
            print('Creating System Object:')
            print('______________________')
            print('')
            sbm.createSystemObject(minimise=minimise)

        return sbm

    def getCAModel(pdb_file, contact_file, 
                   default_parameters=True, 
                   default_forces=True, 
                   create_system=True,
                   contact_force='12-10'):
        """
        Initialises a coarse-grained, carbon alpha (CA), sbmOpenMM system class 
        from a PDB file and a contact file defining the native contacts for the 
        coarse grained model.
        
        The system creation steps are:
        
        1) Add the geometrical parameters for the model.
        2) Add the default force field parameters for the model.
        3) Create the default force objects.
        4) Create the OpenMM system class.
        
        The method can be used to generate an initialized sbmOpenMM system class, that only 
        contains the geometrical parameters, by passing the option default_parameters as False.
        This is useful to store the geometrical values of bonds, angles, dihedrals, etc. in 
        order to add custom parameters and forces.
        
        The method can also be created without the initialisation of the forces classes by 
        setting default_forces to False. This allows to load the default forcefield parameters 
        and to modified them before creating the OpenMM force objects.
        
        Finally, the method can be stopped before creating the OpenMM system class using create_system=False.
        
        Parameters
        ----------
        pdb_file : string
            Path to the input PDB file.
        contact_file : string
            Path to the input native contact file. The file can be an output 
            from the SMOG program (4 column contact file) or a two column file 
            defining the atoms to be paired.
        default_parameters : boolean (True)
            Whether to add default SBM CA forcefield parameters to the model.
        default_forces : boolean (True)
            Whether to initilize default SBM CA force objects. Set to False if
            the parameters will be different from the default ones.
        create_system : boolean (True)
            If True the function will call the createSystemObject() method
            to create an OpenMM system object. If modifications to the default 
            forcefield are necessary this option should be given False.
        
        Returns
        -------
        sbm : sbmOpenMM.system
            Initialized sbmOpenMM.system class with default options for defining 
            a coarse-grained CA SBM force field.
        """
        
        print('Generating CA SBM for PDB file '+pdb_file)
        print('')
        #Set up geometric parameters of the model
        print('Setting up geometrical parameters:')
        print('_________________________________')
        sbm = system(pdb_file)
        print('Keeping only carbon alpha atoms in topology')
        sbm.getCAlphaOnly()
        sbm.getAtoms()
        print('Added '+str(sbm.n_atoms)+' CA atoms')
        sbm.getBonds()
        print('Added '+str(sbm.n_bonds)+' bonds')
        sbm.getAngles()
        print('Added '+str(sbm.n_angles)+' angles')
        sbm.getProperTorsions()
        print('Added '+str(sbm.n_torsions)+' torsions')
        if contact_file != None:
            print('Reading contacts from contact file: '+contact_file)
            sbm.readContactFile(contact_file)
        print('Added '+str(sbm.n_contacts)+' native contacts')
        print('')
        
        
        #Add default parameters to each interaction term
        if default_parameters:
            print('Setting up default forcefield parameters:')
            print('________________________________________')
            print('Adding default bond parameters:')
            sbm.setBondParameters(20000.0)
            print('Adding default angle parameters:')
            sbm.setAngleParameters(40.0)
            print('Adding default torsion parameters:')
            sbm.setProperTorsionParameters(1.0)
            print('Adding default contact parameters:')
            sbm.setNativeContactParameters(1.0)
            print('')
            
        #Create default system force objects
        if default_forces:
            print('Adding Forces:')
            print('_____________')
            sbm.addHarmonicBondForces()
            print('Added Harmonic Bond Forces')
            sbm.addHarmonicAngleForces()
            print('Added Harmonic Angle Forces')
            sbm.addPeriodicTorsionForces()
            print('Added Periodic Torsion Forces')
            if contact_force == '12-10':
                sbm.addLJ12_10ContactForces()
                print('Added Lennard Jones 12-10 Forces to native contacts')
            elif contact_force == '12-10-6':
                sbm.addLJ12_10_6ContactForces()
                print('Added Lennard Jones 12-10-6 Forces to native contacts')
            else:
                raise ValueError('Wrong contact_force option, valid options are "12-10" and "12-10-6"')
                
            sbm.addLJRepulsionForces(sigma=0.4, epsilon=1.0, cutoff=1.5)
            print('Added Lennard Jones 12 non-bonded Forces')
            print('')
            
        #Generate the system object and add previously generated forces
        if create_system:
            print('Creating System Object:')
            print('______________________')
            sbm.createSystemObject(minimise=False, check_bond_distances=False)
            print('OpenMM system Object created')
            print('')
            
        return sbm
    
    def getMultiBasinModel(main_model, alternate_configuration, 
                                       double_minima_threshold=0.05, 
                                       excluded_volume_radius=None,
                                       use_lennard_jones=True,
                                       create_system=True):
        """
        Generate a multibasin model from two sbmOpenMM initialized system classes.
        It currently supports only two configurations.
        
        The two configurations (main and alternate) are compared at the level of 
        native contacts to define common and unique contacts. If any common contact
        is defined to have equilibrium distances significantly different between the 
        configurations, a multi basin Gaussian potential is added to explicitely consider 
        both equilibrium distances in the non-bonded energy function. Unique contacts 
        from the alternate configuration are added to the main configuration. All other
        bonded parameters are maintained from the main configuration.
        
        Optionally, an excluded volume term can be given, with the excluded_volume_radius
        option, to control separately the sphere radius from the equilibrium contact distances.
        
        The set of single-basin (or unique) contacts are added, by default, as Lennard-Jones
        potentials (use_lennard_jones=True) or can be added as Gaussian terms (use_lennard_jones=False)
        with separate control of the excluded volume term.
        
        Finally, the method can be stopped before creating the OpenMM system class using
        create_system=False.
        
        Attributes
        ----------
        
        common_contacts : set
            Set containing the common contacts between configurations.
        dual_basin_contacts : set
            Set containing the subset of common contacts that have a dual basin potential
        unique_contacts : dict
            Dictionary containing the sets of unique contacts for each configuration. The keys are 
            integers, 0 for the main configuration, 1 for the alternate configuration. 
        
        Parameters
        ----------
        main_model : sbmOpenMM.system
            Configuration upon which the multi basin forces should be added.
        alternate_configuration : sbmOpenMM.system
            Second configuration to define as a new minima in the native contact force
            object.
        double_minima_threshold : float
            The minimum equilibrium distance difference to consider a native contact having 
            two different equilibrium distances. If shared contacts between the configurations 
            have equilibrium distances that do not differ more than this value, only the main 
            configuration equilibrium distance will be  kept.
        excluded_volume_radius : float
            Radius of the excluded volume term in the Gaussian native-contact energy function.
        use_lennard_jones : boolean (True)
            Whether to use or not Lennard-Jones potential for the single-basin native contact
            set. If False a Gaussian function with separate control of the excluded volume term.
        create_system : boolean (True)
            If True the function will call the createSystemObject() method
            to create an OpenMM system object.
            
        Returns
        -------
        sbm : sbmOpenMM.system
            Initialized sbmOpenMM.system class with added multi basin contact potential.
        """
           
        if not isinstance(main_model, system):
            raise ValueError('Multi basin contact potential only supports 2 configurations.'+
                            'The first, main_model, must be a sbmOpenMM system class \
                             containing the protein system.')
            
        if not isinstance(alternate_configuration, system):
            raise ValueError('Multi basin contact potential only supports 2 configurations.'+
                            'The second, alternate_configuration, must be a sbmOpenMM system class \
                             containing the same protein system in a different conformation.')
        
        assert len(main_model.atoms) == len(alternate_configuration.atoms)
        
        print('Creating Multi Basin Potential')
        print('')
        print('Contact analysis:')
        print('________________')
        
        #Create new sbm system object from main model object
        multi_basin_model = copy.copy(main_model)
        
        #Set contact forcefield parameters to OrderedDict and None
        multi_basin_model.contacts = OrderedDict()
        multi_basin_model.n_contacts = None
        
        #Create attributes to store unique and common contacts
        multi_basin_model.common_contacts = None
        multi_basin_model.dual_basin_contacts = None
        multi_basin_model.unique_contacts = OrderedDict()
        
        #Reset forces and force groups
        multi_basin_model.harmonicBondForce = None
        multi_basin_model.harmonicAngleForce = None
        multi_basin_model.periodicTorsionForce = None
        multi_basin_model.harmonicImproperForce = None
        multi_basin_model.harmonicPlanarForce = None
        multi_basin_model.lj12_6contactForce = None
        multi_basin_model.lj12_10contactForce = None
        multi_basin_model.lj12_10_6contactForce = None
        multi_basin_model.singleGaussianContactForce = None
        multi_basin_model.doubleGaussianContactForce = None
        multi_basin_model.ljRepulsionForce = None
        multi_basin_model.exclusions = []
        multi_basin_model.forceGroups = OrderedDict()
        
        #Initialise new system object
        multi_basin_model.system = openmm.System()
        
        #Save contact sets with atom indexes only.
        contacts = OrderedDict()
        contact_atoms = OrderedDict()
        for i,m in enumerate([main_model]+[alternate_configuration]):
            contact_atoms[i] = {(c[0].index, c[1].index) : (c[0], c[1]) for c in m.contacts}
            contacts[i] = set([(c[0].index, c[1].index) for c in m.contacts ])
            
        #Get sets of common and unique contacts.
        common_contacts = contacts[0].intersection(contacts[1])
        unique_contacts = OrderedDict()
        unique_contacts[0] = contacts[0] - common_contacts
        unique_contacts[1] = contacts[1] - common_contacts
        
        print('Contacts in main configuration: '+str(len(contacts[0])))
        print('Contacts in alternate configuration: '+str(len(contacts[1])))
        print('Common contacts between configurations: '+str(len(common_contacts)))  
        
        #Set attribute to store common and unique contacts set
        multi_basin_model.common_contacts = common_contacts
        multi_basin_model.unique_contacts[0] = unique_contacts[0]
        multi_basin_model.unique_contacts[1] = unique_contacts[1]
        
        #Check common contacts with large differences (larger than double_minima_threshold)
        #in equilibrium position.
        dm_contacts = set()
        for c in common_contacts:
            
            #Get main an alterante contact distances
            atoms = contact_atoms[0][c]
            d0 = main_model.contacts[contact_atoms[0][c]][0]
            d1 = alternate_configuration.contacts[contact_atoms[1][c]][0]
            
            #Set minimum an maximum distance
            r0 = min(d0,d1)
            r1 = max(d0,d1)
            
            #If difference between equilibrium distances is larger than threshold add double minimum parameters.
            if np.absolute(r0 - r1) > double_minima_threshold:
                
                #If excluded_volume_radius not given treat lower wall as the minimum potential distance.
                if excluded_volume_radius == None:
                    multi_basin_model.contacts[atoms] = (r0, r0, r1, main_model.contacts[atoms][-1])
                    
                #If excluded_volume_radius given add separate control of excluded volume.
                else:
                    multi_basin_model.contacts[atoms] = (excluded_volume_radius*unit.nanometer, r0, r1, main_model.contacts[atoms][-1])
                    
                # Add double minima contacts
                dm_contacts.add(c)
            
            #Otherwise just add single basin parameter 
            else:
                #If excluded_volume_radius not given treat lower wall as the minimum potential distance.
                if excluded_volume_radius == None:
                    #Add regular Lennard-Jones parameters : use_lennard_jones=True
                    if use_lennard_jones:
                        multi_basin_model.contacts[atoms] = (d0, main_model.contacts[atoms][-1])
                    #Add single minimum Gaussian : use_lennard_jones=False
                    else:
                        multi_basin_model.contacts[atoms] = (d0, d0, main_model.contacts[atoms][-1])
                #If excluded_volume_radius given add separate control of excluded volume.
                else:
                    multi_basin_model.contacts[atoms] = (excluded_volume_radius*unit.nanometer, d0, main_model.contacts[atoms][1])
        
        print('Common contacts with an equilibrium distance difference larger than %.2f nm: %s' % 
             (double_minima_threshold, len(dm_contacts)))
        
        #Set attribute to store common double minima contacts
        multi_basin_model.dual_basin_contacts = dm_contacts
        
        #Add unique contacts in main configuration
        print('Unique contacts in main configuration: '+str(len(unique_contacts[0])))
        for c in unique_contacts[0]:
            atoms = contact_atoms[0][c]
            d0 = main_model.contacts[contact_atoms[0][c]][0]
            #Add regular Lennard-Jones parameters : use_lennard_jones=True
            if use_lennard_jones:
                multi_basin_model.contacts[atoms] = (d0, main_model.contacts[atoms][1])
            #Add single minimum Gaussian : use_lennard_jones=False
            else:
                multi_basin_model.contacts[atoms] = (d0, d0, main_model.contacts[atoms][1])
                
        #Add unique contacts in alternate configuration
        print('Unique contacts in alternate configuration: '+str(len(unique_contacts[1])))
        for c in unique_contacts[1]:
            atoms = contact_atoms[1][c]
            d1 = alternate_configuration.contacts[contact_atoms[1][c]][0]
            atoms = contact_atoms[1][c]
            #Add regular Lennard-Jones parameters : use_lennard_jones=True
            if use_lennard_jones:
                multi_basin_model.contacts[atoms] = (d1, alternate_configuration.contacts[atoms][1])
            #Add single minimum Gaussian : use_lennard_jones=False
            else:
                multi_basin_model.contacts[atoms] = (d1, d1, alternate_configuration.contacts[atoms][1])
        
        multi_basin_model.n_contacts = len(multi_basin_model.contacts)
        print('Total contacts in multi basin potential: '+str(multi_basin_model.n_contacts))
        print('')
        
        #Add forces to multi_basin_model
        print('Adding Forces:')
        print('_____________')
        
        #Add prexisting forces
        if main_model.harmonicBondForce != None:
            print('Adding Harmonic Bond Forces')
            multi_basin_model.addHarmonicBondForces()
        if main_model.harmonicAngleForce != None:
            print('Adding Harmonic Angle Forces')
            multi_basin_model.addHarmonicAngleForces()
        if main_model.periodicTorsionForce != None:
            print('Adding Periodic Torsion Forces')
            multi_basin_model.addPeriodicTorsionForces()
        if main_model.harmonicImproperForce != None:
            print('Adding Harmonic Improper Forces')
            multi_basin_model.addHarmonicImproperForces()
        if main_model.harmonicPlanarForce != None:
            print('Adding Periodic Planar Forces')
            multi_basin_model.addHarmonicPlanarForces()
        
        if main_model.model_type == 'AA':
            contact_parameters = multi_basin_model.getAANativeContactParameters()
            multi_basin_model.setNativeContactParameters(contact_parameters)
        
        #Add contact forces
        if use_lennard_jones:
            #Add Lennard Jones contacts for AA model if main_model is AA.
            if main_model.lj12_6contactForce != None:
                print('Adding Lennard Jones 12-6 potential to single minimum contacts')
                multi_basin_model.addLJ12_6ContactForces()

            #Add Lennard Jones contacts for CA model if main_model is CA.
            if main_model.lj12_10contactForce != None:
                print('Adding Lennard Jones 12-10 potential to single minimum contacts')
                multi_basin_model.addLJ12_10ContactForces()
        else:
            print('Adding Gaussian potential to single minimum contacts')
            
        #Add gaussian potential
        print('Adding Gaussian potential to double minimum contacts')
        multi_basin_model.addGaussianContactForces()
        
        #Add repulsion forces
        print('Adding repulsion Lennard Jones 12 potential for non native contacts')
        multi_basin_model.addLJRepulsionForces(sigma=multi_basin_model.rf_sigma,
                                               epsilon=multi_basin_model.rf_epsilon,
                                               cutoff=multi_basin_model.rf_cutoff)
        print('')
        print('Creating System Object:')
        print('______________________')
        
        for i,a in enumerate(list(multi_basin_model.topology.atoms())):
            if a != multi_basin_model.atoms[i]:
                print(i,a, multi_basin_model.atoms[i])        
        
        if create_system:
            multi_basin_model.createSystemObject(minimise=False, check_bond_distances=False)
        print('')
        
        if main_model.model_type == 'AA':
            multi_basin_model.model_type = 'AA-MB'
        elif main_model.model_type == 'CA':
            multi_basin_model.model_type = 'CA-MB'
        
        return multi_basin_model

