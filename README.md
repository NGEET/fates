
![FATES_logo](.github/images/logo_fates_small.png)
------------------------------
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3825473.svg)](https://doi.org/10.5281/zenodo.3825473)

This repository holds the Functionally Assembled Terrestrial Ecosystem Simulator (FATES).  FATES is a numerical terrestrial ecosystem model. Its development and support is primarily supported by the Department of Energy's Office of Science, through the Next Generation Ecosystem Experiment - Tropics ([NGEE-T](https://ngee-tropics.lbl.gov/)) project.

For more information on the FATES model, see our [User's Guide](https://fates-users-guide.readthedocs.io/en/latest/) and [technical documentation](https://fates-docs.readthedocs.io/en/latest/index.html).  

Please submit any questions you may have to the [FATES Github Discussions board](https://github.com/NGEET/fates/discussions).

To receive email updates about forthcoming release tags, regular meeting notifications, and other important announcements, please join the [FATES Google group](https://groups.google.com/g/fates_model).

## Important Guides:
------------------------------

[User's Guide](https://fates-users-guide.readthedocs.io/en/latest/)

[How to Contribute](https://github.com/NGEET/fates/blob/master/CONTRIBUTING.md)

[Table of FATES and Host Land Model API compatability](https://fates-users-guide.readthedocs.io/en/latest/user/Table-of-FATES-API-and-HLM-STATUS.html)

[List of Unsupported or Broken Features](https://fates-users-guide.readthedocs.io/en/latest/user/Current-Unsupported-or-Broken-Features.html)

[Code of Conduct](https://github.com/NGEET/fates/blob/master/CODE_OF_CONDUCT.md)

## Important Note:
------------------------------

**Most users should not need to directly clone this repository.  FATES needs to be run through a host model, and all supported host-models are in charge of cloning and loading the fates software.**

FATES has support to be run via the Energy Exascale Earth System Model (E3SM), the Community Earth System Model (CESM), or its land component, the Community Terrestrial Systems Model (CTSM).

https://github.com/E3SM-Project/E3SM

https://github.com/ESCOMP/cesm

The FATES, E3SM and CTSM teams maintain compatability of the NGEET/FATES master branch with the E3SM master and CTSM master branches respectively. There may be some modest lag time in which the latest commit on the FATES master branch is available to these host land models (HLM) by default.  This is typically correlated with FATES development updates forcing necessary changes to the FATES API.  See the table of [FATES API/HLM compatibility](https://fates-users-guide.readthedocs.io/en/latest/user/Table-of-FATES-API-and-HLM-STATUS.html) for information on which fates tag corresponds to which HLM tag or commit.  
