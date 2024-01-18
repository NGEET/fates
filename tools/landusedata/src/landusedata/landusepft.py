import argparse, os, sys
import xarray as xr
import xesmf as xe

from landusedata.landusepftmod import ImportLandusePFTFile, AddLatLonCoordinates, RenormalizePFTs
from landusedata.utils import ImportLUH2StaticFile, ImportRegridTarget
from landusedata.utils import SetMaskRegridTarget, DefineStaticMask

def main(args):

    # Open the files
    ds_static  = ImportLUH2StaticFile(args.luh2_static_file)
    filelist = [args.clm_surface_file,
                args.clm_luhforest_file,
                args.clm_luhpasture_file,
                args.clm_luhother_file]
    ds_landusepfts = []
    for filename in filelist:
        ds_landusepfts.append(ImportLandusePFTFile(filename))

    # Add lat/lon coordinates to the CLM5 landuse data
    for dataset in ds_landusepfts:
        AddLatLonCoordinates(dataset)

    # Define static luh2 ice/water mask
    mask_icwtr = DefineStaticMask(ds_static)

    # Calculate the bareground percentage after initializing data array list
    # Normalize the percentages
    percent = []
    percent_bareground = ds_landusepfts[0].PCT_NAT_PFT.isel(natpft=0)
    percent_bareground = (percent_bareground / 100.0) * mask_icwtr
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

    # Open and prepare the target dataset
    ds_target = ImportRegridTarget(args.regrid_target_file)

    # Set the mask for the regrid target
    ds_target = SetMaskRegridTarget(ds_target)

    # Create an output dataset to contain individually regridded landuse percent datasets
    ds_output = xr.Dataset()

    # Loop through percentage list and regrid each entry
    for index,data_array in enumerate(percent):

        # Get the name for the new variable
        varname = ds_var_names[index]

        # Convert current percent data array into temporary dataset
        ds_percent = data_array.to_dataset(name=varname)

        # Apply mask for the current dataset
        ds_percent['mask'] = mask_icwtr
        if (varname != 'frac_brgnd'):
            mask_zeropercent = xr.where(ds_percent[varname].sum(dim='natpft') == 0.,0,1)
            ds_percent['mask'] = ds_percent.mask * mask_zeropercent

        # Regrid current dataset
        print('Regridding {}'.format(varname))
        regridder = xe.Regridder(ds_percent, ds_target, "conservative_normed")
        ds_regrid = regridder(ds_percent)

        # Drop mask to avoid conflicts when merging
        if (varname != 'frac_brgnd'):
            # There is no mask currently on the bareground
            ds_regrid = ds_regrid.drop_vars(['mask'])

        # Append the new dataset to the output dataset
        ds_output = ds_output.merge(ds_regrid)

    # Duplicate the 'primary' data array into a 'secondary' data array.  Eventually
    # this will contain different data from a future CLM landuse x pft update
    ds_output['frac_secnd'] = ds_output.frac_primr.copy(deep=True)

    # ds_regrid = ds_regrid.rename_dims(dims_dict={'lat':'lsmlat','lon':'lsmlon'})

    # Output dataset to netcdf file
    print('Writing fates landuse x pft dataset to file')
    output_file = os.path.join(os.getcwd(),args.output)
    ds_output.to_netcdf(output_file)

if __name__ == "__main__":
    main()
