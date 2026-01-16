#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from openmm import unit
import numpy as np

class geometry:
    """
    A class to hold methods for calculating geometrical values 
    given sets of atom coordinates.
    """

    def position2Array(position, output_unit):
        """Converts an OpenMM position object quantity into a numpy array.

        Parameters
        ----------
        position : openmm.unit.quantity.Quantity
            Array containing quantity objects [e.g. (x,y,z) array returned
            from positions].
        output_unit : openmm.unit.unit.Unit
            Unit in which to return the items of the array.

        Returns
        -------
        numpy.ndarray
            A numpy array containing the quantity values converted to floats.
        """

        return np.array([c.value_in_unit(output_unit) for c in position])

    def bond(coord1, coord2):
        """Calculate the distance length between two (x,y,z) quantity coordinates.

        Parameters
        ----------
        coord1 : openmm.unit.quantity.Quantity array
            Vector for the first coordinate.
        coord2 : openmm.unit.quantity.Quantity array
            Vector for the second coordinate.

        Returns
        -------
        openmm.unit.quantity.Quantity
            Quantity (value and unit) of the distance length in nanometers.
        """

        coord1 = geometry.position2Array(coord1, unit.nanometer)

        coord2 = geometry.position2Array(coord2, unit.nanometer)

        bond_length = np.linalg.norm(coord2 - coord1)

        return bond_length * unit.nanometer

    def angle(coord1, coord2, coord3):
        """Calculate the angle length between three (x,y,z) quantity coordinates.

        Parameters
        ----------
        coord1 : openmm.unit.quantity.Quantity array
            Vector for the first coordinate.
        coord2 : openmm.unit.quantity.Quantity array
            Vector for the second coordinate.
        coord3 : openmm.unit.quantity.Quantity array
            Vector for the third coordinate.

        Returns
        -------
        openmm.unit.quantity.Quantity
            Quantity (value and unit) of the angle length in radians.
        """

        coord1 = geometry.position2Array(coord1, unit.nanometer)
        coord2 = geometry.position2Array(coord2, unit.nanometer)
        coord3 = geometry.position2Array(coord3, unit.nanometer)

        v1 = coord1 - coord2
        v2 = coord3 - coord2
        cos_theta = np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
        angle = np.arccos(np.clip(cos_theta, -1, 1))

        return angle * unit.radian

    def torsion(coord1, coord2, coord3, coord4):
        """Calculate the torsion angle length between four (x,y,z) quantity 
        coordinates.

        Parameters
        ----------
        coord1 : openmm.unit.quantity.Quantity array
            Vector for the first coordinate.
        coord2 : openmm.unit.quantity.Quantity array
            Vector for the second coordinate.
        coord3 : openmm.unit.quantity.Quantity array
            Vector for the third coordinate.
        coord4 : openmm.unit.quantity.Quantity array
            Vector for the fourth coordinate.

        Returns
        -------
        openmm.unit.quantity.Quantity
            Quantity (value and unit) of the torsion length in radians.
        """

        coord1 = geometry.position2Array(coord1, unit.nanometer)
        coord2 = geometry.position2Array(coord2, unit.nanometer)
        coord3 = geometry.position2Array(coord3, unit.nanometer)
        coord4 = geometry.position2Array(coord4, unit.nanometer)

        v1 = coord2 - coord1
        v2 = coord3 - coord2
        v3 = coord4 - coord3

        c1 = np.cross(v2, v3)
        c2 = np.cross(v1, v2)

        p1 = (v1 * c1).sum(-1)
        p1 *= (v2 * v2).sum(-1) ** 0.5
        p2 = (c1 * c2).sum(-1)

        return np.arctan2(p1, p2) * unit.radian
