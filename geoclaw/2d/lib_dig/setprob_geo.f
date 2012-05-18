c=========================================================================
      subroutine setprob
c=========================================================================

      use geoclaw_module
      use digclaw_module
      use topo_module
      use dtopo_module
      use qinit_module
      use auxinit_module

      implicit double precision (a-h,o-z)


      call set_geo          !# sets basic parameters g and coord system
      call set_tsunami      !# sets parameters specific to tsunamis
      call set_flow_grades  !# sets refinement parameters
      call set_qinit        !# specifies file(s) for initial conditions
      call set_auxinit      !# specifies file(s) for auxiliary variables
      call set_topo         !# specifies topography (bathymetry) files
      call set_dtopo        !# specifies moving topo 'dtopo' from earthquake
      call set_dig
      call set_pinit
      call setregions       !# specifies where refinement is allowed/forced
      call setgauges        !# locations of measuring gauges
      call setfixedgrids    !# specifies output on arbitrary uniform fixed grids

      return
      end
