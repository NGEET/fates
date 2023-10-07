
# Landuse x PFT script
# Usage: python landusepft.py -s <luh2_static_file> \
#                             -f <clm_luhforest_file> \
#                             -p <clm_luhpasture_file> \
#                             -o <clm_luhother_file> \
#                             -s <clm_luhsurface_file> \
#                             -O <output_file>

def main():


def CommandLineArgs():

    parser = argparse.ArgumentParser(description="placeholder desc")

    # Required static luh2 data to get the ice/water fraction for masking
    parser.add_argument("-s", "--luh2_static_file",
                        required=True,
                        help = "luh2 static data file")

    parser.add_argument("-f", "--clm_luhforest_file",
                        required=True,
                        help = "CLM5_current_luhforest_deg025.nc")

    parser.add_argument("-p", "--clm_luhpasture_file",
                        required=True,
                        help = "CLM5_current_luhpasture_deg025.nc")

    parser.add_argument("-o", "--clm_luhother_file",
                        required=True,
                        help = "CLM5_current_luhother_deg025.nc")

    parser.add_argument("-s", "--clm_surface_file",
                        required=True,
                        help = "CLM5_current_surf_deg025.nc")

    # Optional output argument
    parser.add_argument("-O","--output",
                        default = 'LUH2_timeseries.nc',
                        help = "output filename")

    args = parser.parse_args()

    return(args)

if __name__ == "__main__":
    main()
