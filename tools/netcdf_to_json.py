import xarray as xr
import json
import numpy as np
import write_json
import argparse
import code  # For development: code.interact(local=dict(globals(), **locals()))
from pathlib import Path
import datetime

def decode_bytes(obj):
    # If it's bytes, decode it to a string
    if isinstance(obj, bytes):
        return obj.decode('utf-8')
    # If it's a list/array of bytes, decode each element
    if isinstance(obj, np.ndarray) and obj.dtype.kind in ['S', 'V']: 
        return [x.decode('utf-8') for x in obj]
    return obj

def fates_nc_to_json_structured(nc_path, json_path):

    time_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    file_name = Path(nc_path).name
    
    # Open the dataset
    ds = xr.open_dataset(nc_path, decode_times=False)
    
    # Initialize our root object with your 3 high-level categories
    root = {
        "attributes": {},
        "dimensions": {},
        "parameters": {}
    }

    # 1. ATTRIBUTES: Extract global history
    # We use .get() to avoid crashing if the attribute is missing
    root["attributes"]["history"] = ds.attrs.get("history")+f"\nConverted from netcdf file:{file_name} to json on {time_str}"

    # 2. DIMENSIONS: Map dimension names to their integer sizes
    for dim_name, size in ds.dims.items():
        root["dimensions"][dim_name] = int(size)

    for var_name, da in ds.data_vars.items():

        val = da.values.flatten()

        # 1. Determine the simple type string
        if np.issubdtype(da.dtype, np.floating):
            dtype_str = "float"
            is_null = np.isnan(val) | (val > 1e30)
            # For floats, replace NaN with None
            data_list = [float(x) if not n else None for x, n in zip(val, is_null)]
            #data_list = [float(x) if np.isfinite(x) else None for x in da.values.flatten()]
        
        elif np.issubdtype(da.dtype, np.integer):
            dtype_str = "int"
            # Check for the common NetCDF _FillValue or missing_value
            fill_value = da.attrs.get('_FillValue') or da.attrs.get('missing_value')
            data_list = [int(x) if x != fill_value else None for x in da.values.flatten()]
        
        elif da.dtype.kind in ['S', 'U', 'O']:
            dtype_str = "string"
            raw_values = da.values
            # For strings, an empty string or a specific fill string becomes None
            def clean_string(s):
                s = s.decode('utf-8').strip() if isinstance(s, bytes) else str(s).strip()
                return s if s != "" else None
        
            if raw_values.ndim == 0:
                data_list = [clean_string(raw_values.item())]
            else:
                data_list = [clean_string(i) for i in raw_values.flatten()]
        else:
            dtype_str = "unknown"
            data_list = da.values.tolist()

        # 2. Handle scalar wrapping (ensure it's a list [val])
        if not da.dims and dtype_str != "string":
            # If the scalar itself was a NaN/FillValue, it's already [None]
            pass 

        # 3. Build the parameter entry
        param_entry = {
            "dtype": dtype_str,
            "dims": list(da.dims) if da.dims else ["scalar"],
            "data": data_list
        }

        # 4. Elevate attributes (and clean them of bytes)
        for attr_name, attr_val in da.attrs.items():
            if attr_name in ['_FillValue', 'missing_value']:
                continue # We don't need these in the JSON if we're using null
            if isinstance(attr_val, bytes):
                param_entry[attr_name] = attr_val.decode('utf-8').strip()
            else:
                param_entry[attr_name] = attr_val

        root["parameters"][var_name] = param_entry
            
    with open(json_path, 'w') as outfile:
    
        write_json.traverse_data(outfile,root)
        print(f'Writing complete')
            

def main():

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument("files", nargs='+', help="Two arguments: the path to the netcdf infile and path to JSON output file.")
    args = parser.parse_args()
    nc_infile = args.files[0]
    json_outfile = args.files[1]
    
    # Execute
    fates_nc_to_json_structured(nc_infile,json_outfile)


if __name__ == "__main__":
    main()
