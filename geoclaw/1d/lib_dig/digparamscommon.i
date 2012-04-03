      common /physicalparams/  rho_s, rho_f, grav, phi_bed,
     &                      phi_int, delta, kappita, mu, alpha,
     &                      m_crit, m_min, c1, c2, m0, dudx_eps,dry_tol,
     &                      smallu,phys_tol, m_eqn0

      common /flumeparams/ theta1,theta2,flumerad,flumelen,hopperlen,
     &                                 hmax,hoppertop,hopperangle

      common /auxindices/ i_S,i_rho,i_tanpsi,i_D,i_tau,
     &                  i_kappa,i_phi,i_sigbed,i_theta,i_kperm,i_alpha,
     &                  i_topo

      common /initparams/ init_htype, init_ptype, p_initialized

      common /startparams/ pimin,init_pmax_ratio,init_ptf
