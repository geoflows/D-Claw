# D-Claw Tips for Making Source Files

Katy Barnhart


Some tidbits from Dave

1. Sometimes very very thin layers of h can exist if you init by setting
   b (with settopo) and h (with setqinit, h = q1). This might happen when
   you think it shouldn't because of interpolation to the clawpack grids.

   A better practice is to set b and eta (with setqinit eta = q_last).
   In the places where h = 0 set eta to eta << b. The initialization routines
   will set b and eta to the value of b where eta << b.

2. When setting up a run where you want m to have one value in one place  
   and another value (e.g., clearwater of m=0) in another, use a
   setqinit file for m over the entire domain. Set digdata.m0 to the
   sediment flow value. And add a couple of grid cells of the higher
   sediment concentration around the landslide source. This ensures
   that interpolation in m doesn't create a spurious ring of water
   around the sediment domain.

   Note, setqinit for m takes m not hm.

3. Can use a value of h or eta at a landslide source that doesn't include
   clearwater/ocean, and then also set sea level with that geodata input. 
