      subroutine setprob

      use gauges_module

      implicit double precision (a-h,o-z)

      call setdigparams
      call set_gauges('setgauges.data')

      return
      end
