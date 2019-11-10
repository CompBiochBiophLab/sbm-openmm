#!/usr/bin/env python
# coding: utf-8

# In[1]:


from simtk.openmm.app.statedatareporter import StateDataReporter


# In[ ]:


class sbmReporter(StateDataReporter):
    """
    A special case of the StateDataReporter class that outputs information about a simulation, 
    such as energy and temperature, etc. to a file. This method specially outputs the sbmOpenMM 
    force group energies inside the passed sbmOpenMM system object.
    
    It is used in the same way as the OpenMM StateDataReporter class, but it takes as additional
    input an instance of the sbmOpenMM object with the option 'sbmObject'.
    """
    
    def __init__(self, file, reportInterval, sbmObject=None, **kwargs):
        """
        Initialises the SBM OpenMM system class.
        
        Parameters
        ----------
        reportInterval : int
            The interval (in time steps) at which to write frames
        sbmObject : sbmOpenMM.system
            The sbmOpenMM system instance to read force groups from.
        **kwargs : openMM StateDataReporter arguments
            
        Returns
        -------
        initialized StateDataReporter class.
        
        """
        super(sbmReporter, self).__init__(file, reportInterval, **kwargs)
        self._sbmObject = sbmObject
        
    def _constructHeaders(self):
        """
        Build headers for the StateDataReporter class. It builds the headers
        for the force groups contained in the sbmOpenMM system instance. 
        
        Parameters
        ----------
        None
            
        Returns
        -------
        headers : list
            List with strings representing the headers to be written to the report file.
        """
        
        headers = super()._constructHeaders()
        if isinstance(self._sbmObject, system):
            for i,n in enumerate(self._sbmObject.forceGroups):
                 headers.append(n+' (kJ/mol)')
                    
        return headers
    
    def _constructReportValues(self, simulation, state):
        """
        Calculates the energies for the force groups in the sbmOpenMM system instance.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        values : list
            List with floats representing the values to be written to the report file.
        """
        
        values = super()._constructReportValues(simulation, state)
        
        if isinstance(self._sbmObject, system):
            for i,n in enumerate(self._sbmObject.forceGroups):
                values.append(simulation.context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
        
        return values

