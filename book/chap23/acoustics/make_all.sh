
# First run the 1d code in 1drad to provide a reference solution.

cd 1drad
make .output

# Now run the 2d code:

cd ..
make .plots
