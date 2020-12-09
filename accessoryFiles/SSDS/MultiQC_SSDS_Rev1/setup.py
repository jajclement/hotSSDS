#!/usr/bin/env python
"""
MultiQC is a tool to aggregate bioinformatics results across many samples into a single report. It is written in Python and contains modules for a large number of common bioinformatics tools.

You can install MultiQC from PyPI as follows::

    pip install multiqc

Then it's just a case of going to your analysis directory and running the script::

    multiqc .

MultiQC will scan the specified directory (:code:`'.'` is the current dir) and produce a report detailing whatever it finds.

The report is created in :code:`multiqc_report.html` by default. Tab-delimited data files are created in :code:`multiqc_data/` to give easy access for downstream processing.

For more detailed instructions, run :code:`multiqc -h`

See the MultiQC website for documentation and tutorial videos: http://multiqc.info

MultiQC was written by Phil Ewels (http://phil.ewels.co.uk) at SciLifeLab Sweden (http://www.scilifelab.se)
"""

from setuptools import setup, find_packages

version = '0.7dev'

print("""-----------------------------------
 Installing MultiQC version {}
-----------------------------------

""".format(version))

setup(
    name = 'multiqc',
    version = version,
    author = 'Phil Ewels',
    author_email = 'phil.ewels@scilifelab.se',
    description = "Create aggregate bioinformatics analysis reports across many samples and tools",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'http://multiqc.info',
    download_url = 'https://github.com/ewels/MultiQC/releases',
    license = 'GPLv3',
    packages = find_packages(),
    include_package_data = True,
    zip_safe = False,
    scripts = ['scripts/multiqc'],
    install_requires = [
         'markupsafe==1.1.0',
         'kiwisolver==1.1.0',
         'numpy==1.16.6',
         'pyparsing==2.3.1',
         'jinja2==2.8',
         'simplejson==3.17.2',
         'pyyaml==5.3.1',
         'click==7.1.2',
         'matplotlib==2.2.4'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'bamtools = multiqc.modules.bamtools:MultiqcModule',
            'bismark = multiqc.modules.bismark:MultiqcModule',
            'bowtie2 = multiqc.modules.bowtie2:MultiqcModule',
            'bowtie1 = multiqc.modules.bowtie1:MultiqcModule',
            'cutadapt = multiqc.modules.cutadapt:MultiqcModule',
            'fastq_screen = multiqc.modules.fastq_screen:MultiqcModule',
            'fastqc = multiqc.modules.fastqc:MultiqcModule',
            'featureCounts = multiqc.modules.featureCounts:MultiqcModule',
            'hicup = multiqc.modules.hicup:MultiqcModule',
            'kallisto = multiqc.modules.kallisto:MultiqcModule',
            'methylQA = multiqc.modules.methylQA:MultiqcModule',
            'picard = multiqc.modules.picard:MultiqcModule',
            'preseq = multiqc.modules.preseq:MultiqcModule',
            'qualimap = multiqc.modules.qualimap:MultiqcModule',
            'rseqc = multiqc.modules.rseqc:MultiqcModule',
            'salmon = multiqc.modules.salmon:MultiqcModule',
            'samblaster = multiqc.modules.samblaster:MultiqcModule',
            'samtools = multiqc.modules.samtools:MultiqcModule',
            'skewer = multiqc.modules.skewer:MultiqcModule',
            'snpeff = multiqc.modules.snpeff:MultiqcModule',
            'star = multiqc.modules.star:MultiqcModule',
            'ssds = multiqc.modules.ssds:MultiqcModule',
            'tophat = multiqc.modules.tophat:MultiqcModule',
            'trimmomatic = multiqc.modules.trimmomatic:MultiqcModule',
        ],
        'multiqc.templates.v1': [
            'default = multiqc.templates.default',
            'default_dev = multiqc.templates.default_dev',
            'simple = multiqc.templates.simple',
            'geo = multiqc.templates.geo',
        ],
        # 'multiqc.cli_options.v1': [
            # 'my-new-option = myplugin.cli:new_option'
        # ],
        # 'multiqc.hooks.v1': [
            # 'execution_start = myplugin.hooks:execution_start',
            # 'config_loaded = myplugin.hooks:config_loaded',
            # 'before_modules = myplugin.hooks:before_modules',
            # 'after_modules = myplugin.hooks:after_modules',
            # 'execution_finish = myplugin.hooks:execution_finish',
        # ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)

print("""
--------------------------------
 MultiQC installation complete!
--------------------------------
For help in running MultiQC, please see the documentation available
at http://multiqc.info or run: multiqc --help
""")
