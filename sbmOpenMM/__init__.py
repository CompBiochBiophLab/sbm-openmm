"""
Python library to run structure based model (SBM) simulations using the OpenMM toolkit


sbmOpenMM is a Python library to run protein structure-based model (SBM) simulations using OpenMM toolkit. The library offers flexibility for creating SBM force fields that can be customised to capture different aspects of protein SBM potential energy exploration.

Considering an input structure, the library automatizes the creation of forces to specify it as the only minimum configuration in the potential energy function. Bonds, angles and torsions are maintained close to their equilibrium configuration, while native contact interactions are allowed to form and break using regular or modified Lennard-Jones potentials. This allows complete and local protein unfolding, restricting the interactions only to the evolutionarily relevant chemical contacts, to explore more thoroughly the relevant configurational space of protein folding and function.

Different granularities for the models can be selected as All-heavy-Atom and alpha-carbon representations. These basic models can also be extended to multi-basin potentials employing more than one input configuration. Here, shared native contacts are modeled with special Gaussian functions to allow for more than one equilibrium distance.

The library offers methods to tailor forcefield parameter for each force term. Combining these basic methods and force implementations, sbmOpenMM offers easy set up of more complex force field definition that can aid in a better exploration of different biophysical phenomena.

sbmOpenMM is divided in three main classes:

    1. geometry

    2. models

    3. system

The first class, geometry, contains methods to calculate the geometrical parameters from the input structures. These parameters are used to define the input conformation as the global minimum configuration in the potential energy function. The second class, models, allows to easily set up predefined SBM models, that encompass coarse grained, all atom and multi basin potentials. The third class, system, is the main class that holds all the methods to define, modify and create SBMs to be simulated with OpenMM.

The library is open-source and offers flexibility to create custom SBMs or to modify the predefined topology based models included in it.
"""

from .core import geometry
from .core import system
from .core import models

from .reporter import sbmReporter
