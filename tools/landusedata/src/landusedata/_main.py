import argparse

from landusedata.luh2 import main as luh2main
from landusedata.landusepft import main as lupftmain

def main(argv=None):

    # Define top level parser
    parser = argparse.ArgumentParser(description="FATES landuse data tool")

    # Target regrid file - is there a nice way to share this between subparsers?
    # parser.add_argument('regrid_target_file', help='target surface data file with desired grid resolution')

    # Define subparser option for luh2 or landuse x pft data tool subcommands
    subparsers = parser.add_subparsers(required=True, title='subcommands',
                                       help='landuse data tool subcommand options')
    luh2_parser = subparsers.add_parser('luh2', prog='luh2',
                                        help='generate landuse harmonization timeseries data output')
    lupft_parser = subparsers.add_parser('lupft', prog='lupft',
                                         help='generate landuse x pft static data map output')

    # Set the default called function for the subparser command
    luh2_parser.set_defaults(func=luh2main)
    lupft_parser.set_defaults(func=lupftmain)

    # LUH2 subparser arguments
    luh2_parser.add_argument('regrid_target_file',
                             help='target surface data file with desired grid resolution')
    luh2_parser.add_argument("luh2_static_file",
                             help = "luh2 static data file")
    luh2_parser.add_argument('luh2_states_file',
                             help = "full path of luh2 raw states file")
    luh2_parser.add_argument('luh2_transitions_file',
                             help = "full path of luh2 raw transitions file")
    luh2_parser.add_argument('luh2_management_file',
                             help = "full path of luh2 raw management file")
    luh2_parser.add_argument("-w", "--regridder_weights",
                             default = 'regridder.nc',
                             help = "filename of regridder weights to save")
    luh2_parser.add_argument("-b","--begin",
                             type = int,
                             default = None,
                             help = "beginning of date range of interest")
    luh2_parser.add_argument("-e","--end",
                             type = int,
                             default = None,
                             help = "ending of date range to slice")
    luh2_parser.add_argument("-o","--output",
                             default = 'LUH2_timeseries.nc',
                             help = "output filename")

    # Landuse x pft subparser arguments
    lupft_parser.add_argument('regrid_target_file',
                             help='target surface data file with desired grid resolution')
    lupft_parser.add_argument('luh2_static_file',
                             help = "luh2 static data file")
    lupft_parser.add_argument('clm_luhforest_file',
                              help = "CLM5_current_luhforest_deg025.nc")
    lupft_parser.add_argument('clm_luhpasture_file',
                              help = "CLM5_current_luhpasture_deg025.nc")
    lupft_parser.add_argument('clm_luhother_file',
                              help = "CLM5_current_luhother_deg025.nc")
    lupft_parser.add_argument('clm_surface_file',
                              help = "CLM5_current_surf_deg025.nc")
    lupft_parser.add_argument("-o","--output",
                              default = 'fates_landuse_pft_map.nc',
                              help = "output filename")

    # Parse the arguments
    args = parser.parse_args(argv)

    # Call the default function for the given subcommand
    args.func(args)

    # Return successful completion
    return 0

# Gaurd against import time side effects
if __name__ == '__main__':
    raise SystemExit(main())
