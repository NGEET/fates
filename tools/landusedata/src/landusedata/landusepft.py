import argparse, os, sys
import xarray as xr
import xesmf as xe
from landusedata.landusepftmod import ImportLandusePFTFile
from landusedata.landusepftmod import AddLatLonCoordinates, DefineMask, RenormalizePFTs
from landusedata.luh2mod import ImportStaticLUH2File

def main(args):

    # Open the files
    ds_static  = ImportStaticLUH2File(args.luh2_static_file)
    filelist = [args.clm_surface_file,
                args.clm_luhforest_file,
                args.clm_luhpasture_file,
                args.clm_luhother_file]
    ds_landusepfts = []
    for filename in filelist:
        ds_landusepfts.append(ImportLandusePFTFile(filename))

    # Define the landuse mask based on the static luh2 file
    mask_static = DefineMask(ds_static)

    # Add lat/lon coordinates to the CLM5 landuse data
    for dataset in ds_landusepfts:
        AddLatLonCoordinates(dataset)

    # Calculate the bareground percentage after initializing data array list
    # Normalize the percentages
    percent = []
    percent_bareground = ds_landusepfts[0].PCT_NAT_PFT.isel(natpft=0)
    percent_bareground = (percent_bareground / 100.0) * mask_static
    percent.append(percent_bareground)

    # Renormalize the PCT_NAT_PFT for each dataset using the mask
    for data_array in ds_landusepfts:
        percent.append(RenormalizePFTs(data_array))

    # Calculate the primary and secondary PFT fractions as the forest
    # and nonforest-weighted averages of the forest and other PFT datasets.
    percent[2] = ds_static.fstnf * percent[2] + (1. - ds_static.fstnf) * percent[-1]

    # Note that the list order is:
    # bareground, surface data, primary, pasture, rangeland (other)
    ds_var_names = ['frac_brgnd','frac_csurf','frac_primr','frac_pastr','frac_range']

    # Combine percent data arrays into a new dataset
    ds_percent = xr.Dataset()
    for index,data_array in enumerate(percent):
        ds_percent = ds_percent.merge(data_array.to_dataset(name=ds_var_names[index]))

    # Duplicate the 'primary' data array into a 'secondary' data array.  Eventually
    # this will contain different data from a future CLM landuse x pft update
    ds_percent = ds_percent.merge(percent[2].to_dataset(name='frac_secnd'))

    # Regrid dataset (if necessary)
    # TODO: add option to reuse weights file from luh2 code
    if (not isinstance(args.regrid_target_file,type(None))):
        print('Regridding')
        ds_target = xr.open_dataset(args.regrid_target_file)
        # TODO: reuse (and refactor as necessary) code from luh2 data tool
        ds_target = ds_target.rename_dims(dims_dict={'lsmlat':'lat','lsmlon':'lon'})
        ds_target['lon'] = ds_target.LONGXY.isel(lat=0)
        ds_target['lat'] = ds_target.LATIXY.isel(lon=0)
        regridder = xe.Regridder(ds_percent, ds_target, "conservative_normed")
        ds_regrid = regridder(ds_percent)
        ds_regrid = ds_regrid.rename_dims(dims_dict={'lat':'lsmlat','lon':'lsmlon'})
        output_file = os.path.join(os.getcwd(),args.output)
        # rename the dimensions back to lsmlat/lsmlon?
    else:
        ds_regrid = ds_percent

    # Output dataset to netcdf file
    print('Writing fates landuse x pft dataset to file')
    ds_regrid.to_netcdf(output_file)

if __name__ == "__main__":
    main()
