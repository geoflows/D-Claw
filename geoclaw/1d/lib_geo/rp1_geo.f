c =====================================================================
      subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &             fwave,s,amdq,apdq)
c =====================================================================
c
c     # solve Riemann problems for the 1D shallow water equations
c     # with source-term resulting from variable topography b(x,t)
c     #   (h)_t + (u h)_x = 0 
c     #   (uh)_t + (uuh + 0.5*gh**2)_x = -g*h*b_x 
c
c
c     # On input,
c     # ql contains the state vector at the left edge of each cell
c     # qr contains the state vector at the right edge of each cell
c     #
c     # On output, wave contains the fwaves/s,
c     #            s the speeds,
c     #            amdq the  left-going flux difference  A**- \Delta q
c     #            apdq the right-going flux difference  A**+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp1 is called with ql=qr=q.
c
c
c   # This Riemann solver performs and approximate augmented decomposition
c   # The details are described in my PhD thesis:

c   #"Finite Volume Methods and Adaptive Refinement for Tsunami Propagation
c   #   and Inundation." University of Washington, July 2006. 
c   # David L George

c======================================================================
c             last modified 07/27/06
c======================================================================

      implicit none

      integer iter,m,mw,meqn,mwaves,mx,maxmx,mbc,i,maxiter,k

      double precision ql(1-mbc:maxmx+mbc, meqn)
      double precision qr(1-mbc:maxmx+mbc, meqn)
      double precision s(1-mbc:maxmx+mbc, mwaves)
      double precision fwave(1-mbc:maxmx+mbc, meqn, mwaves)
      double precision amdq(1-mbc:maxmx+mbc, meqn)
      double precision apdq(1-mbc:maxmx+mbc, meqn)
      double precision auxl(1-mbc:maxmx+mbc, *)
      double precision auxr(1-mbc:maxmx+mbc, *)

      double precision r(4,4)
      double precision speed(4)
      double precision lambda(4)
      double precision roe(2)
      double precision beta(4)
      double precision del(4)
      double precision A(4,4)

      double precision abs_zero,abs_tol,grav,g
      double precision hL,hR,uL,uR,huR,huL,f1L,f2L,f1R,f2R,bR,bL
      double precision delq1,delq2,delf1,delf2,delphi
      double precision hm,hum,um,u1m,u2m,s1m,s2m,hbar,ubar
      double precision uhat,chat,roe1,roe2,s1R,s2L,sR,sL

      double precision hMax,hMin,se1,se2,heinf
      double precision r41Lbound,r41Ubound,barlam12,tildelam12
      double precision rare1st,rare2st,alpha1,alpha2

      double precision det1,det2,det3,determinant
      double precision Sint,dhdb
      double precision rescheck1,rescheck2

      logical solidwallL,solidwallR,rare1,rare2,singular
      logical sonic,transonic,rarecorrector

      common /geo/ grav,abs_tol

      g=grav
      maxiter=1
      abs_zero= abs_tol

c======================================================================

c======================================================================
c         loop through Riemann problems at each grid cell
c======================================================================

      do 30 i=2-mbc,mx+mbc
c======================================================================
c        Initialize Riemann problem for grid interface
c======================================================================

c         Inform of a bad riemann problem from the start
          if((qr(i-1,1).lt.0.d0).or.(ql(i,1) .lt. 0.d0)) then
             write(*,*) 'Negative input: hl,hr,i=',qr(i-1,1),ql(i,1),i
          endif

c         Set small values to zero to prevent rounding trouble
          if (qr(i-1,1).lt.0.d0) then
              qr(i-1,1)=0.d0
              qr(i-1,2)=0.d0
          endif

          if (ql(i,1).lt.0.d0) then
              ql(i,1)=0.d0
              ql(i,2)=0.d0
          endif

          do mw=1,mwaves
              s(i,mw)=0.d0
              do m=1,meqn
                  fwave(i,m,mw)=0.d0
              enddo
          enddo

          do mw=1,mwaves
              speed(mw)=0.d0
              lambda(mw)=0.d0
              beta(mw)=0.d0
              do m=1,4
                  r(m,mw)=0.d0
              enddo
          enddo

          if (qr(i-1,1).le.abs_tol.and.ql(i,1).le.abs_tol) then
              go to 30
          endif

          hL = qr(i-1,1)
          hR = ql(i,1)
          huL = qr(i-1,2)
          huR = ql(i,2)
          bL = auxr(i-1,1)
          bR = auxl(i,1)

c         Check for dry regions
          if (hR.le.abs_tol) then
              hR=0.d0
              uR=0.d0
          else
              uR=huR/hR
          endif
          if (hL.le.abs_tol) then
              hL=0.d0
              uL=0.d0
          else
              uL=huL/hL
          endif


c=======================================================================
c  Initialize solid-wall problem for a wet/dry interface if necessary
c=======================================================================

          solidwallL=.false.
          solidwallR=.false.

          if (hL.le.abs_tol) then
              if (bL.gt.bR+hR) then
                  call riemann(hR,hR,-uR,uR,hm,s1m,s2m,rare1,rare2,1)
                  if (bL.gt.bR+hm) then
                      solidwallL=.true.
                      hL=hR
                      uL=-uR
                      huL=-huR
                      bL=bR
                  endif
              elseif (bL.gt.bR+hR) then
                  bL=bR+hR
              endif
          endif

          if (hR.le.abs_tol) then
              if (bR.gt.bL+hL) then
                  call riemann(hL,hL,uL,-uL,hm,s1m,s2m,rare1,rare2,1)
                  if (bR.gt.bL+hm) then
                      solidwallR=.true.
                      hR=hL
                      huR=-huL
                      uR=-uL
                      bR=bL
                  endif
              elseif (bR.gt.bL+hL) then
                  bR=bL+hL
              endif
          endif

c=======================================================================
c   set-differences to be decomposed into eigenvectors
c=======================================================================
          f1L=hL*uL
          f2L=0.5d0*g*hL**2 + hL*uL**2
          f1R=hR*uR
          f2R=0.5d0*g*hR**2 + hR*uR**2

          delq1=hR-hL
          delq2=f1R-f1L
          delf1=f1R-f1L
          delf2=f2R-f2L
          delphi=bR-bL

          del(1)=delq1
          del(2)=delq2
          del(3)=delf2
          del(4)=delphi

c=====================================================================
c   Determine Riemann structure and set-up eigenvectors
c=====================================================================
          hMax=max(hR,hL)
          hMin=min(hR,hL)

          maxiter=1
          if ((.not.solidwallL).and.(.not.solidwallR)) then
              call riemann(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,maxiter)
          endif

          sL=uL-sqrt(g*hL)
          sR=uR+sqrt(g*hR)

          uhat=sqrt(g*hL)*uL + sqrt(g*hR)*uR
          uhat=uhat/(sqrt(g*hR)+sqrt(g*hL))
          chat=sqrt(g*0.5d0*(hR+hL))
          roe1=uhat-chat
          roe2=uhat+chat

          lambda(1)=min(sL,roe1) !Einfeldt speed
          lambda(1)=min(lambda(1),s2m) !Modified Eindfeldt speed
          lambda(3)=max(sR,roe2) !Einfeldt speed
          lambda(3)=max(lambda(3),s1m) !Modified Eindfeldt speed
          lambda(2)=.5d0*(lambda(1)+lambda(3)) 
          lambda(4)=0.d0

          se1=lambda(1)
          se2=lambda(3)


c===========determine the corrector wave====
          rarecorrector=.false.
          alpha1=0.5d0
          alpha2=0.9d0
          if (rare1.or.rare2) then
              rare1st=3.d0*(sqrt(g*hL)-sqrt(g*hm))
              rare2st=3.d0*(sqrt(g*hR)-sqrt(g*hm))
              if (rare1.and.se1*s1m.lt.0.d0) alpha1=0.2d0
              if (rare2.and.se2*s2m.lt.0.d0) alpha1=0.2d0
              if (max(rare1st,rare2st).gt.alpha1*(se2-se1).and.
     &              max(rare1st,rare2st).lt.alpha2*(se2-se1)) then
                  rarecorrector=.true.
                  if (rare1st.gt.rare2st) then
                      lambda(2)=s1m
                  else
                      lambda(2)=s2m
                  endif
              endif
          endif
c=============================================

          do mw=1,mwaves
              r(1,mw)=1.d0
              r(2,mw)=lambda(mw)
              r(3,mw)=(lambda(mw))**2
              r(4,mw)=0
          enddo

c     ======Einfeldt depth (positivity)=======
          heinf= (huL-huR+se2*hR-se1*hL)/(se2-se1)
          if (heinf.lt.hMin/5.d0) rarecorrector=.false.
c=============================================
          if (.not.rarecorrector) then
              r(1,2)=0.d0
              r(2,2)=0.d0
              r(3,2)=1.d0
              r(4,2)=0.d0
          endif

c=====Determin the steady state eigenvector r(m,4)==================
c===================================================================
          hbar=0.5d0*(hR+hL)
          ubar=0.5d0*(uR+uL)
          tildelam12= max(0.d0,uL*uR)- g*hbar
          barlam12= ubar**2 - g*hbar

c     =====Positive preserving bounds for r(1,4)==================

          if (se1.lt.-abs_zero.and.se2.gt.abs_zero) then 
              if (delphi.gt.0.d0) then
                  r41Lbound=min((se2-se1)*heinf/(se1*delphi),-1.d0)
                  r41Ubound= -1.d0
              elseif (delphi.lt.0.d0) then
                  r41Lbound=min((se2-se1)*heinf/(se2*delphi),-1.d0)
                  r41Ubound= -1.d0
              else   
                  r41Lbound=0.d0
                  r41Ubound=0.d0
              endif  
          elseif (se1.ge.abs_zero) then
              if (delphi.gt.0.d0) then
                  r41Ubound=(se2-se1)*heinf/(se1*delphi)
                  r41Lbound=0.d0
              elseif (delphi.lt.0.d0) then
                  r41Ubound=-hL/delphi
                  r41Lbound=0.d0
              else   
                  r41Lbound=0.d0
                  r41Ubound=0.d0
              endif
          elseif (se2.le.-abs_zero) then
              if (delphi.gt.0.d0) then
                  r41Ubound=hR/delphi
                  r41Lbound=0.d0
              elseif (delphi.lt.0.d0) then
                  r41Ubound=heinf*(se2-se1)/(se2*delphi)
                  r41Lbound=0.d0
              else   
                  r41Lbound=0.d0
                  r41Ubound=0.d0
              endif
          else
              r41Lbound=0.d0
              r41Ubound=0.d0
          endif

c     ===Source term r(3,4) S^+_- and r(1,4) =============
c     Resonant (near sonic) check for r(1,4)==============

          sonic=.false.
          transonic=.false.
    
          rescheck1=tildelam12*barlam12
          rescheck2= se1*se2 
    
          if ((uL+sqrt(g*hL))*(uR+sqrt(g*hR)).lt.abs_zero) sonic=.true.
          if ((uL-sqrt(g*hL))*(uR-sqrt(g*hR)).lt.abs_zero) sonic=.true.
          if (rescheck1.le.abs_zero) sonic=.true.
          if (rescheck2*barlam12.le.abs_zero) sonic=.true.
          if (rescheck2*tildelam12.le.abs_zero) sonic=.true.
          if (se1.lt.0.d0.and.s1m.gt.0.d0) sonic=.true.
          if (s2m.lt.0.d0.and.se2.gt.0.d0) sonic=.true.
    
          if (sonic) then
              Sint=-g*hbar
              dhdb=0.d0
          else
              Sint=-g*hbar*(tildelam12/barlam12)
              dhdb=-Sint/tildelam12
          endif

          Sint=max(Sint,-g*hMax)
          Sint=min(Sint,-g*hMin)

          dhdb=max(dhdb,r41Lbound)
          dhdb=min(dhdb,r41Ubound)
c     ====================================================
          r(1,4)=dhdb
          r(2,4)=0.d0
          r(3,4)=Sint
          r(4,4)=1.d0
c================================================================
c================================================================   


c==================Subtract off beta(4)*r^4 before solving==========
          do m=1,4
              del(m)=del(m)-r(m,4)*delphi
          enddo
c====================================================================

c==================Solve the 3x3 system for beta's===================
c=======Determine determinant of eigenvector matrix========
          do mw=1,3
              do m=1,3
                  A(m,mw)=r(m,mw)
              enddo
          enddo

          det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
          det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
          det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
          determinant=det1-det2+det3

          singular=.false.
          if (determinant.eq.0.d0) then
              singular=.true.
              write(*,*) 'RP1: **ERROR** Degenerate eigenspace'
              stop
          elseif (dabs(determinant).lt.1.d-99) then
              write(*,*) 'RP1:**WARNING** system is poorly conditioned'
          endif
c========solve for beta(k) using Cramers Rule=================
          do k=1,3
              do mw=1,3
                  do m=1,3
                      A(m,mw)=r(m,mw)
                      A(m,k)=del(m)
                  enddo
              enddo

              det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
              det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
              det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

              beta(k)=(det1-det2+det3)/determinant
          enddo
c==================determine speeds and fwaves ======================

          if (solidwallL) then
              do mw=1,2
                  beta(mw)=0.d0
                  lambda(mw)=0.d0
              enddo
          elseif (solidwallR) then
              do mw=2,3
                  beta(mw)=0.d0
                  lambda(mw)=0.d0
              enddo
          endif

c======================diagnose and debug===========================
          if (.false.) then

          endif
c====================================================================

          beta(4)=0.d0
          do mw=1,mwaves
              s(i,mw)=lambda(mw)
              do m=1,meqn
                  fwave(i,1,mw)=beta(mw)*r(2,mw)
                  fwave(i,2,mw)=beta(mw)*r(3,mw)
              enddo
          enddo
          
 30   continue
c=============End: loop through cell interfaces========================


c============= compute fluctuations====================================
      do m=1,meqn
          do  i=2-mbc, mx+mbc
              amdq(i,m) = 0.d0
              apdq(i,m) = 0.d0
              do  mw=1,mwaves
                  if (s(i,mw).lt.0.d0) then
                      amdq(i,m) = amdq(i,m) + fwave(i,m,mw)
                  elseif (s(i,mw).gt.0.d0) then
                      apdq(i,m) = apdq(i,m) + fwave(i,m,mw)
                  else
                      amdq(i,m) = amdq(i,m) + .5d0*fwave(i,m,mw)
                      apdq(i,m) = apdq(i,m) + .5d0*fwave(i,m,mw)
                  endif
              enddo
          enddo
      enddo

      return
      end


c=============================================================================
      subroutine riemann(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,maxiter)
c=============================================================================

      implicit none

      integer iter,maxiter

      double precision abs_tol,grav,g
      double precision hL,hR,uL,uR,hm,s1m,s2m,u1m,u2m,um,delu
      double precision h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR

      logical rare1,rare2

      common /geo/ grav,abs_tol

      g=grav

c=====================================================================
c     Test for Riemann structure
c=====================================================================
      h_min=min(hR,hL)
      h_max=max(hR,hL)
      delu=uR-uL

      if (h_min.le.0.d0) then
          hm=0.d0
          um=0.d0
          s1m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
          s2m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
          if (hL.le.0.d0) then
              rare2=.true.
              rare1=.false.
          else
              rare1=.true.
              rare2=.false.
          endif
      else

          F_min= delu+2.d0*(sqrt(g*h_min)-sqrt(g*h_max))
          F_max= delu +
     &         (h_max-h_min)*(sqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

          if (F_min.gt.0.d0) then !2-rarefactions
              hm=(1.d0/(16.d0*g))*
     &               max(0.d0,-delu+2.d0*(sqrt(g*hL)+sqrt(g*hR)))**2
              um=sign(1.d0,hm)*(uL+2.d0*(sqrt(g*hL)-sqrt(g*hm)))

              s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
              s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)

              rare1=.true.
              rare2=.true.

          elseif (F_max.le.0.d0) then !2 shocks

c               ===root finding using a Newton iteration on sqrt(h)===
              h0=h_max
              do iter=1,maxiter
                  gL=sqrt(.5d0*g*(1/h0 + 1/hL))
                  gR=sqrt(.5d0*g*(1/h0 + 1/hR))
                  F0=delu+(h0-hL)*gL + (h0-hR)*gR
                  dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+
     &                   gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
                  slope=2.d0*sqrt(h0)*dfdh
                  h0=(sqrt(h0)-F0/slope)**2
              enddo 
c               ======================================================
              hm=h0
              u1m=uL-(hm-hL)*sqrt((.5d0*g)*(1/hm + 1/hL))
              u2m=uR+(hm-hR)*sqrt((.5d0*g)*(1/hm + 1/hR))
              um=.5d0*(u1m+u2m)

              s1m=u1m-sqrt(g*hm)
              s2m=u2m+sqrt(g*hm)

              rare1=.false.
              rare2=.false.

          else !one shock one rarefaction
              h0=h_min
c               ===== root find using secant iteration on sqrt(h)=====
              do iter=1,maxiter
                  F0=delu + 2.d0*(sqrt(g*h0)-sqrt(g*h_max))
     &                  + (h0-h_min)*sqrt(.5d0*g*(1/h0+1/h_min))
                  slope=(F_max-F0)/(h_max-h_min)
                  h0=h0-F0/slope
              enddo
c               =====================================================
              hm=h0
              if (hL.gt.hR) then
                  um=uL+2.d0*sqrt(g*hL)-2.d0*sqrt(g*hm)
                  s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
                  s2m=uL+2.d0*sqrt(g*hL)-sqrt(g*hm)
                  rare1=.true.
                  rare2=.false.
              else
                  s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)
                  s1m=uR-2.d0*sqrt(g*hR)+sqrt(g*hm)
                  um=uR-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hm)
                  rare2=.true.
                  rare1=.false.
              endif
          endif
      endif


c=====================================================================
c   RIEMANN STRUCTURE KNOWN
c=====================================================================
      return

      end
