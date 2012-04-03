# Source this file to make all plots.

# Make reference solution in output directory specified in the Makefile below:
make .output -f Makefile_qref

# Clean up so data and output will be recreated below:
make clean

# Now run code and plot results
make .plots
