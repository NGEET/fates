#!/usr/bin/env python3


# Steps
# - import clm landuse-pft data (1/4 degree)
# - import luh2 static data file (1/4 degree)
# - set the mask based on static data file `icwtr` data
# - calculate the percentage of forest, pasture, and other (i.e. range) pft percentages
# - calculate the primary and secondard forest percent using the static data `fstnf`
# primary_secondary_percent = luh2_staticdata.fstnf * forest_pft_percent + (1.- luh2_staticdata.fstnf) * other_pft_percent
# - concatenate all this information together (including mask)
# - add/adjust lat/lon names and time input for regridding
# - regrid using luh2mod.py tooling

