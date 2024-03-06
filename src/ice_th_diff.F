      SUBROUTINE ice_th_diff(nlay_s,nlay_i,kideb,kiut,numofday)

      use var_thermo_vertical

#if   defined THERMO_LIM
! The original LIM stuff
      use ice_thermo_lim
#elif defined THERMO_FV
! Huwald and Tremblay model
      use ice_thermodynamic_FV
#elif defined THERMO_FE
! FD home-made
      use ice_thermodynamic_FE
#endif
        !!------------------------------------------------------------------
        !!                ***         ROUTINE ice_th_diff       ***
        !! ** Purpose :
        !!   This routine determines the time evolution of snow and sea-ice 
        !!   temperature profiles.
        !! ** Method  :
        !!       This is done by solving the heat equation diffusion with
        !!       a Neumann boundary condition at the surface and a Dirichlet one
        !!       at the bottom. Solar radiation is partially absorbed into the ice.
        !!       The specific heat and thermal conductivities depend on ice salinity
        !!       and temperature to take into account brine pocket melting. The 
        !!       numerical
        !!       scheme is an iterative Crank-Nicolson on a non-uniform multilayer grid 
        !!       in the ice and snow system.
        !!       The successive steps of this routine are
        !!       Vertical grid
        !!           1.  Thermal conductivity at the interfaces of the ice layers
        !!           2.  Internal absorbed radiation
        !!           3.  Scale factors due to non-uniform grid
        !!           4.  Kappa factors
        !!           Then iterative procedure begins
        !!           5.  specific heat in the ice
        !!           6.  eta factors
        !!           7.  surface flux computation
        !!           8.  tridiagonal system terms
        !!           9.  solving the tridiagonal system with Gauss elimination
        !!           Iterative procedure ends according to a criterion on evolution
        !!           of temperature
        !!
        !! ** Arguments :
        !!           kideb , kiut : Starting and ending points on which the 
        !!                         the computation is applied
        !!
        !! ** Inputs / Ouputs : (global commons)
        !!           surface temperature : t_su_b
        !!           ice/snow temperatures   : t_i_b, t_s_b
        !!           ice salinities          : s_i_b
        !!           number of layers in the ice/snow: nlay_i, nlay_s
        !!           total ice/snow thickness : ht_i_b, ht_s_b
        !!
        !! ** External : 
        !!
        !! ** References :
        !!
        !! ** History :
        !!           (02-2003) Martin Vancoppenolle, Louvain-la-Neuve, Belgium
        !!

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'
      LOGICAL :: ln_write=.false.
      LOGICAL,SAVE :: start=.true.

      DO ji = kideb, kiut

        nlsno = nlay_s
        nlice = nlay_i
        if (start) call allocate_thermo_1d
        start=.false.

        do layer = 1, nlay_i
          ti(layer) = t_i_b(ji,layer)
          si(layer) = s_i_b(ji,layer)
          qi(layer) = q_i_b(ji,layer)
          dzi(layer) = deltaz_i_phy(layer)
          swradab_i(layer) = radab_phy_i(layer)
        enddo

        do layer = 1, nlay_s
          ts(layer) = t_s_b(ji,layer)
          qs(layer) = q_s_b(ji,layer)
          dzs(layer) = deltaz_s_phy(layer)
          swradab_s(layer)=radab_s(layer)
        enddo

        hi  = ht_i_b(ji)
        hs  = ht_s_b(ji)
        tsu = t_su_b(ji)
        tbo = t_bo_b(ji)

        fsens = fscbqb(ji)
        flat = fltbqb(ji)
        tair = tp0-10.d0
        qair = 1d-3
        uair = 1.d0
        oceflx = oce_flx
        dwnlw = ratbqb(ji)
        swrad = fsolgb(ji)
        fac_transmi = ab(ji)
        pres   = 1d5
        snowfall = hnpbqb(ji) / 86400.d0 ! convert to m/s because Martin has in in m/day
        seasal = oce_sal
        snosub = parsub
        frac_sni = frtr_si_phy
        ith_cond = thcon_i_swi

! FD if first time, set the inter ice-snow temperature to the mean
#if defined THERMO_FV || defined THERMO_FE
        if (first) then
            ts(nlsno+1)=0.5d0*(ti(1)+ts(nlsno))
            first=.false.
# if defined THERMO_FE
            call init_sigma_ice_FE
# elif defined THERMO_FV
            call init_sigma_ice_FV
# endif
        endif
#endif
!-------------------------------------------------
!
! call to module subroutine
!
#if   defined THERMO_LIM
      call ice_thermo_diff (ddtb,ln_write,numout)
! FD Here should be call to ice salinity diffusion routine
      sinew=si ! default, no change
      call ice_thermo_dh   (ddtb,ln_write,numout)
      call ice_thermo_remap(ddtb,ln_write,numout)

#elif defined THERMO_FV
      call ice_thermo(ddtb)
#elif defined THERMO_FE
      call ice_thermo_FE(ddtb)
#endif
!#
!
!-------------------------------------------------

        do layer = 1, nlay_i
          t_i_b(ji,layer) = ti(layer)
          s_i_b(ji,layer) = si(layer)
          q_i_b(ji,layer) = qi(layer)
          deltaz_i_phy(layer) = dzi(layer)
        enddo

! brine info
        do layer = 1, nlay_i
          rayleigh(layer) = brine_r(layer)
          brvolum (layer) = brine_v(layer)
          brflush (layer) = brine_u(layer)
        enddo

        do layer = 1, nlay_s
          t_s_b(ji,layer) = ts(layer)
          q_s_b(ji,layer) = qs(layer)
          deltaz_s_phy(layer) = dzs(layer)
        enddo

        ht_i_b(ji) = hi
        ht_s_b(ji) = hs
        t_su_b(ji) = tsu
        t_bo_b(ji) = tbo

        dh_s_tot(ji) = dhs
        dh_i_surf(ji) = dhi_surf
        dh_i_bott(ji) = dhi_bot
        dh_snowice(ji) = dh_sni

      fbbqb(ji) = oce_flx       ! Oceanic  heat flux and precipitations
      fcsb(ji)  = fsens     ! sensible and latent heat flux
      fleb(ji)  = flat
      qsfcb(ji) = q0
      tfs(ji)   = tmelt

      fratsb(ji)= netlw
      fc_su(ji)   = fcsu 
      fc_bo_i(ji) = fcbo
      fprec     = fprecip
      e_skel    = eskel
      qsummer   = massmelt
      fsbp      = fsalt
      s_i_new   = si_acc_new
      fsnic     = fcons_sni

      DO layer = 0, nlay_s
         fc_s(ji,layer) = hfc_s(layer)
      END DO
      DO layer = 0, nlay_i
         fc_i(ji,layer) = hfc_i(layer)
      END DO

      DO layer = 1, nlay_i
         z_i_phy(layer) = zi(layer)
      ENDDO
      DO layer = 1, nlay_s
         z_s_phy(layer) = zs(layer)
      ENDDO

      enddo ! ji = kideb, kiut
!
!------------------------------------------------------------------------------
! End of ice_th_diff
      END SUBROUTINE
