.. SOAP documentation master file, created by
   sphinx-quickstart on Wed May 14 11:06:57 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Calculate SOAP potential for ranking MHC2 peptides.
===================================================

Directory: SOAP/examples/mhc2/

1. Preprocess PDB if necessary and generate sequence files, example file: SOAP/examples/mhc2/pdb_X_2.2A_0.25rfree.pir. Please see the module :mod:`sequences` on how to generate different subsets of PDBs. 

2. Prepare Decoys. The original decoy files generated using modeller, mhc2_original.tar.gz and the prepared decoys, mhc2.tar.gz can be found in the example directory. 

.. literalinclude:: ../examples/mhc2/prepare_mhc2_decoys.py

3. Run SOAP script to select models.

* Select the best models with distance features alone, only considering interface atom pairs.

.. literalinclude:: ../examples/mhc2/select_dist_pot.py 

* Select the best models with distance and accessbility feature, multiple recovery functions, only considering interface atom pairs.

.. literalinclude:: ../examples/mhc2/select_dpp_pot.py

* Select the best models with orientation feature, considering both the interface and the introchain atom pairs.

.. literalinclude:: ../examples/mhc2/select_od_pot.py

4. Run SOAP script to calculate the optimal statistical potential using the best model found in the last step.

.. literalinclude:: ../examples/mhc2/optimize_od_pot.py

5. Write out the optimial potential in hdf5 or lib formart for use in Modeller, IMP or other packages.

.. literalinclude:: ../examples/mhc2/write_soap_mhc2.py
