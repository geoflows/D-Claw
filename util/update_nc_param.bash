cp Makefile Makefile.bak
cp setplot.py setplot.py.bak
perl -pi -w -e 's/\/valout_geo.f/\/valout_nc_geo.f/g;' Makefile
perl -pi -w -e 's/FFLAGS \?=/FFLAGS ?=-I\/usr\/include  -lnetcdf -lnetcdff/g;' Makefile
perl -pi -w -e 's/return plotdata/plotdata\.format = "netcdf"\t\t\t# use netcdf output\n    return plotdata/g;' setplot.py
