Structure Based Models (SBMs) are representations of protein systems based on simplifications made over classical Molecular Dynamics (MD) force fields. Their are based on the energy landscape theory of protein folding and the principle of minimal frustration. The models maintain protein structures by focusing on chemical contacts formed at the native protein configuration, ignoring other non-native contacts. This allows for simpler force field definitions which capture essential protein dynamics at a much lower computational expense than traditional MD simulations.

sbmOpenMM is a Python library that offers flexibility to set up SBMs using the MD framework of OpenMM toolkit. It automates the creation of openmm.system classes that contain the necessary force field parameters to run molecular dynamics simulations using a protein structure and a contact map as the only necessary inputs. 

sbmOpenMM is divided in three main classes:

1. geometry
2. models
3. system
   
The first class, geometry, contains methods to calculate the geometrical parameters from the input structures. These parameters are used to define the input conformation as the global minimum configuration in the potential energy function. The second class, models, allows to easily set up predefined SBM models, that encompass coarse grained, all atom and multi basin potentials. The third class, system, is the main class that holds all the methods to define, modify and create SBMs to be simulated with OpenMM.

The library is open-source and offers flexibility to create custom SBMs or to modify the predefined topology based models included in it.
