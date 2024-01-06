import argparse

# Start moving command line arguments to a main package
# The idea here is that the user is simply going to call a couple of top level options
# for either creating the luh2 raw data file or the landusepft data file.  Ideally the user
# will not need to wrap this in a shell script
# This is something we also want to be testable with pytest

def main(argv=None):

    parser = argparse.ArgumentParser(description="FATES landuse data tool")

    # The user should input the landuse_data_output_type.  Do this to simply
    # avoid determining what they want from the inputs.
    parser.add_argument('landuse_type', choices=['luh2','lupft'],
                        help="landuse data file output type to be created, either luh2 timeseries or landuse x pft static mapping")

    # The user needs to input specific file locations
    # Make these grouped based on the landuse_type choice?
    #
    # Let the user override the output?
    #
    # Let the user truncate the timeseries data

    # Parse the arguments
    args = parser.parse_args(argv)

    # Depending on the landuse_type, call the appropriate module
    if args.landuse_type == 'luh2':
        print("calling luh2 code")
    elif args.landuse_type == 'lupft':
        print("calling lupft code")

    return 0

# Gaurd against import time side effects
if __name__ == '__main__':
    raise SystemExit(main())
