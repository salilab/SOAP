[![docs](https://readthedocs.org/projects/soapsalilab/badge/)](http://salilab.org/SOAP/doc/)

MOTIVATION:

Statistical potentials have been widely used for modeling whole proteins and their parts (e.g. sidechains and loops) as well as interactions between proteins, nucleic acids and small molecules. Here, we formulate the statistical potentials entirely within a statistical framework, avoiding questionable statistical mechanical assumptions and approximations, including a definition of the reference state.

RESULTS:

We derive a general Bayesian framework for inferring statistically optimized atomic potentials (SOAP) in which the reference state is replaced with data-driven 'recovery' functions. Moreover, we restrain the relative orientation between two covalent bonds instead of a simple distance between two atoms, in an effort to capture orientation-dependent interactions such as hydrogen bonds. To demonstrate this general approach, we computed statistical potentials for protein-protein docking (SOAP-PP) and loop modeling (SOAP-Loop). For docking, a near-native model is within the top 10 scoring models in 40% of the PatchDock benchmark cases, compared with 23 and 27% for the state-of-the-art ZDOCK and FireDock scoring functions, respectively. Similarly, for modeling 12-residue loops in the PLOP benchmark, the average main-chain root mean square deviation of the best scored conformations by SOAP-Loop is 1.5 Å, close to the average root mean square deviation of the best sampled conformations (1.2 Å) and significantly better than that selected by Rosetta (2.1 Å), DFIRE (2.3 Å), DOPE (2.5 Å) and PLOP scoring functions (3.0 Å). Our Bayesian framework may also result in more accurate statistical potentials for additional modeling applications, thus affording better leverage of the experimentally determined protein structures. Availability and implementation: SOAP-PP and SOAP-Loop are available as part of [MODELLER](http://salilab.org/modeller/).

[See more](http://bioinformatics.oxfordjournals.org/content/29/24/3158.long)
