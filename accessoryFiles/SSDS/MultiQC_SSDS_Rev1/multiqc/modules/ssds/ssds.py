#!/usr/bin/env python

""" MultiQC module to parse output from SSDS reports """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config, BaseMultiqcModule

# Import the ssds submodules
from . import ssstats

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Samtools has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found. """
    
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='SSDS QC reports',
        anchor='ssds', target='SSDS QC',
        href='http://genome.cshlp.org/content/early/2012/03/20/gr.130583.111.full.pdf',
        info=" generates statistics for Single Stranded DNA Sequencing (SSDS).")
        
        # Set up class objects to hold parsed data
        self.sections = list()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()
        
        # Call submodule functions
        n['ssstats'] = ssstats.parse_reports(self)
        if n['ssstats'] > 0:
            log.info("Found {} stats reports".format(n['ssstats']))
        
        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning
        
        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
