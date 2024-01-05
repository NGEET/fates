import argparse

# Start moving command line arguments to a main package
# The idea here is that the user is simply going to call a couple of top level options
# for either creating the luh2 raw data file or the landusepft data file.  Ideally the user
# will not need to wrap this in a shell script
# This is something we also want to be testable with pytest

def main(argv=None):
    parser = argparse.ArgumentParser(description="FATES landuse data tool")
    parser.add_argument('--output')
    # output_group_argument = parser.add_mutually_exclusive_group(required=True)
    # output_group_argument.add_argument("--luh2_timeseries", help="LUH2 timeseries data file")
    # output_group_argument.add_argument("--landuse_pft", help="landuse by pft static file")

    args = parser.parse_args(argv)

    print(f"data type is {args.output}")

# Gaurd against import time side effects
if __name__ == '__main__':
    raise SystemExit(main())
