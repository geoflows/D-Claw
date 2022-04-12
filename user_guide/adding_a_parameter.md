# D-Claw How to Add a Parameter

Katy Barnhart

Depending on which module uses the parameter, you'll put it in
clawdata, geodata, or digdata. For the moment lets assume that its
digdata and its used by digclaw_module.

Edit the `__init__` and `write` functions of the class `DigclawInputData`
to define your new parameter, and set its defaults. Recommend setting
defaults so that behavior of default value is the same as behavior before
the parameter was created.

Run `make .data` to test that the variables are read from setrun.py
correctly and placed in setdig.data correctly.

Define the variable and add a line to digclaw_mod.f90 that reads this
parameter in from the setdig.data file and assigns it to the new variable.

Where `use digclaw_module` is stated you now have this variable. To
selectively import write:
`use digclaw_module, only : name_of_var`
