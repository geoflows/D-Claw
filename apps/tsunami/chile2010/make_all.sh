
make new     # needed for mem_storage.mod 
make topo    # create topography files
make .plots
make .htmls
python plot_dart.py  # plot gauge data on DART buoy data plot
