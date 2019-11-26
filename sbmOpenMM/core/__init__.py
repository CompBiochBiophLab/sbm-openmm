"""
core package of the sbmOpenMM package that contains the main sbmOpenMM classes. 

The sbmOpenMM.core package contains the three sbmOpenMM main classes:

    1. geometry

    2. models

    3. system

The first class, geometry, contains methods to calculate the geometrical parameters from the input structures. These parameters are used to define the input conformation as the global minimum configuration in the potential energy function. The second class, models, allows to easily set up predefined SBM models, that encompass coarse grained, all atom and multi basin potentials. The third class, system, is the main class that holds all the methods to define, modify and create SBMs to be simulated with OpenMM.

"""

from .geometry import geometry
from .models import models
from .system import system
