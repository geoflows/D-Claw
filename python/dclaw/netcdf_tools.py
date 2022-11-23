import rioxarray  # activate the rio accessor
import xarray as xr
import numpy as np

_level = "AMR_level"
_qelements = ["h", "hu", "hv", "hm", "pb", "hchi", "delta_a", "eta"]
_qunits = {"h": "meters",
           "hu": "meters squared per second",
           "hv": "meters squared per second",
           "hm": "meters",
           "pb": "newton per meter squared",
           "hchi": "meters",
           "delta_a":"meters",
           "eta": "meters"
           }

time_unit="seconds"
reference_time="model start"
space_unit="meters"

#TODO CHECK REGISTRATION.

def griddata2netcdf(time, X, Y, Q_out, outfile, qlst, write_level, epsg=None):

    # get X and Y coordinates (double check corner referencing)
    #print(X.shape)
    #print(Y.shape)
    #print(Q_out.shape)
    x = X[0, :]
    y = Y[:, 0]
    nj = len(x)
    ni = len(y)

    if 0 in qlst:
        mask = True
        h = Q_out[0, :, :].reshape((ni, nj))
        no_material = h<0.0001
    else:
        mask = False
    # create data_vars dictionary based on what is in qlst and write_level
    data_vars = {}
    for qind in qlst:
        if qind<len(_qelements):
            Qsel = Q_out[qind, :, :].reshape((ni, nj))
            if mask and qind <len(_qelements) -1: # don't do eta.
                Q = np.ma.masked_array(Qsel, no_material).reshape((1,ni, nj))
            else:
                Q = Qsel.reshape((1,ni, nj))
            varname = _qelements[qind]
            data_array_attrs = {"units":_qunits[varname]}
            data_vars[varname] = ([ "time", "y", "x",], Q, data_array_attrs)

    if write_level:
        data_array_attrs = {"units":"None"}
        lev = Q_out[-1, :, :].reshape((1, ni, nj))
        data_vars[_level] = ([ "time", "y", "x",], lev, data_array_attrs)

    ds_attrs = {"description":"D-Claw model output"}
    ds = xr.Dataset(
        data_vars=data_vars,
        coords=dict(
            x=(["x"], x, {"units": space_unit}),
            y=(["y"], y, {"units": space_unit}),
            time=("time", [time], {"units": "seconds"}),
            reference_time=reference_time,
        ),
        attrs=ds_attrs,
        )

    if epsg is not None:
        ds.rio.write_crs(epsg, inplace=True)
        #print(ds.spatial_ref)
    #print(ds)
    ds.to_netcdf(outfile)
