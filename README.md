# FATES
------------------------------

This repository holds the Functionally Assembled Terrestrial Ecosystem Simulator (FATES).  FATES is a numerical terrestrial ecosystem model. Its development and support is primarily supported by the Department of Energy's Office of Science, through the Next Generation Ecosystem Experiment - Tropics ([NGEE-T](https://ngee-tropics.lbl.gov/)) project.

For more information on the FATES model, see our [wiki](https://github.com/NGEET/fates/wiki) and [technical documentation](https://fates-docs.readthedocs.io/en/latest/index.html).


## Important Guides:
------------------------------

[How to Contribute](https://github.com/NGEET/fates/blob/master/CONTRIBUTING.md)

[List of Unsupported or Broken Features](https://github.com/NGEET/fates/wiki/Current-Unsupported-or-Broken-Features)

[Code of Conduct](https://github.com/NGEET/fates/blob/master/CODE_OF_CONDUCT.md)


## Important Note:
------------------------------

**Most users should not need to directly clone this repository.  FATES needs to be run through a host model, and all supported host-models are in charge of cloning and loading the fates software.**

FATES has support to be run via the Energy Exascale Earth System Model (E3SM), the Community Earth System Model (CESM), or its land component, the Community Terrestrial Systems Model (CTSM).

https://github.com/E3SM-Project/E3SM

https://github.com/ESCOMP/cesm
https://github.com/ESCOMP/ctsm


## Important Note About Host-Models and Compatible Branches:
------------------------------------------------------------

The FATES and E3SM teams maintain compatability of the NGEET/FATES master branch with the **E3SM master** branch. When changes to the FATES API force compatability updates with E3SM, there may be some modest lag time.

The FATES team maintains compatability of the NGEET/FATES master branch with the **CTSM fates_next_api** branch.  Since the FATES team uses this branch for its internal testing, this compatability is tightly (immediately) maintained and these two should always be in sync.  However, CTSM master may become out of sync with FATES master for large periods (months) of time.




