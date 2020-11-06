"""
reporter package of the sbmOpenMM package that contains the sbmReporter class.

The sbmOpenMM.reproter package contains the sbmReporter class.

sbmReporter is a special class of the OpenMM StateDataReporter class, that additionally
accepts a sbmobject to print the SBM forcefield energies.
"""

from .sbm_reporter import sbmReporter
from .sbm_reporter import readOpenMMReporterFile
