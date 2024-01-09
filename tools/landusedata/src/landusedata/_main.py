import argparse

from landusedata.luh2 import luh2print
from landusedata.landusepft import lupftprint

# Start moving command line arguments to a main package
# The idea here is that the user is simply going to call a couple of top level options
# for either creating the luh2 raw data file or the landusepft data file.  Ideally the user
# will not need to wrap this in a shell script
# This is something we also want to be testable with pytest

def main(argv=None):

    # Define top level parser and subparser option for luh2 or landuse x pft data tool commands
    parser = argparse.ArgumentParser(description="FATES landuse data tool")
    subparsers = parser.add_subparsers(required=True, help='subcommand help')


    # The user should input the landuse_data_output_type.  Do this to simply
    # avoid determining what they want from the inputs.
    luh2_parser = subparsers.add_parser('luh2', help='landuse harmonization timeseries data output')
    lupft_parser = subparsers.add_parser('lupft', help='landuse x pft static data map output')
    # parser.add_argument('landuse_type', choices=['luh2','lupft'],
    #                     help="landuse data file output type to be created, either luh2 timeseries or landuse x pft static mapping")

    # Set the default called function for the subparser command
    luh2_parser.set_defaults(func=luh2print)
    lupft_parser.set_defaults(func=lupftprint)

    # Target regrid file
    # Note that this will be the first argument after all the subparser arguments
    parser.add_argument('regrid_target_file', help='target surface data file with desired grid resolution')

    # The user needs to input specific file locations
    # Make these grouped based on the landuse_type choice?
    #
    # Let the user override the output?
    #
    # Let the user truncate the timeseries data
    # TO DO: using the checking function to report back if invalid file input

    # LUH2 subparser arguments
    luh2_parser.add_argument("luh2_state_file",
                         help = "full path of luh2 raw states file")

    # Landuse x pft subparser arguments
    lupft_parser.add_argument("lupft_surf_file",
                         help = "full path of the CLM landuse x pft surface data file")

    # Parse the arguments
    args = parser.parse_args(argv)

    # Call the default function for the given command
    args.func(args)

    # Depending on the landuse_type, call the appropriate module
    # if args.landuse_type == 'luh2':
    #     print("calling luh2 code")
    # elif args.landuse_type == 'lupft':
    #     print("calling lupft code")

    return 0

# Gaurd against import time side effects
if __name__ == '__main__':
    raise SystemExit(main())
