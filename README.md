# sbmOpenMM

## Description

sbmOpenMM is a Python library to run protein structure-based model (SBM) simulations using OpenMM toolkit. The library offers flexibility for creating SBM force fields that can be customised to capture different aspects of protein SBM potential energy exploration.

Considering an input structure, the library automatizes the creation of forces to specify it as the only minimum configuration in the potential energy function. Bonds, angles and torsions are maintained close to their equilibrium configuration, while native contact interactions are allowed to form and break using regular or modified Lennard-Jones potentials. This allows complete and local protein unfolding, restricting the interactions only to the evolutionarily relevant chemical contacts, to explore more thoroughly the relevant configurational space of protein folding and function.

Different granularities for the models can be selected as All-heavy-Atom and alpha-carbon representations. These basic models can also be extended to multi-basin potentials employing more than one input configuration. Here, shared native contacts are modeled with special Gaussian functions to allow for more than one equilibrium distance.

The library offers many methods to tailor forcefield parameters and definitions for each force term. Combining these basic methods and force implementations, sbmOpenMM offers easy set up of more complex force field definition that can aid in a better exploration of different biophysical phenomena.

