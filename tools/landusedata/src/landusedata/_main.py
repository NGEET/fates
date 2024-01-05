import argparse

# Start moving command line arguments to a main package
# The idea here is that the user is simply going to call a couple of top level options
# for either creating the luh2 raw data file or the landusepft data file.  Ideally the user
# will not need to wrap this in a shell script
# This is something we also want to be testable with pytest

def main():
    print("this is main")

# Gaurd against import time side effects
if __name__ == '__main__':
    raise SystemExit(main())
