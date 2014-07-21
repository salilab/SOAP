.. SOAP documentation master file, created by
   sphinx-quickstart on Wed May 14 11:06:57 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SOAP's documentation!
================================

Installation
------------
1. Install required packages defined in `envConfig.py` (:mod:`detailed instructions <envConfig>`).
2. Copy `envConfig.py` to `env.py` and modify it to work with your setup.

Basic workflow to generate a statistical potential
--------------------------------------------------

1. Preprocess PDB if necessary and generate sequence files. 
2. Prepare decoys. 
3. Run SOAP script to select models. 
4. Run SOAP script to calculate the optimal statistical potential using the best model.
5. Write out the optimal potential in hdf5 or lib format for use in Modeller, IMP or other packages. 

Examples:

.. toctree::
   :maxdepth: 2

   mhc2
   soaplig

Additonal undocumented examples can be found at SOAP/examples. 

Pre-calculated SOAP tables
--------------------------
Pre-calculated SOAP statistical potentials for ranking peptides, loops, proteins, protein interfaces and ligands can be found at : `salilab.org/SOAP <http://salilab.org/SOAP>`_


Modules:
--------

.. toctree::
   :maxdepth: 2

   env
   jobSchedule
   modelSelection
   sequences
   decoys
   feature
   statsTable
   recoveryFucntion
   rankScore
   refineScore   
   benchmark   
   scorer
   sampling
   crossValidate
   loop
   utility
   
Please cite:
------------
Dong GQ, Fan H, Schneidman-Duhovny D, Webb B, Sali A. Optimized atomic statistical potentials: assessment of protein interfaces and loops. Bioinformatics. 2013 Dec 15; 29(24):3158-66.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

