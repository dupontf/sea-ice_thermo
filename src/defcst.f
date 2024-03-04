      SUBROUTINE defcst(n99)

        !!------------------------------------------------------------------
        !!                ***  ROUTINE defcst       ***
        !! ** Purpose :
        !!           Defines constants of the model
        !! ** Method  :
        !!           Definitions !!! 
        !!
        !! ** Arguments :
        !!           n99
        !!
        !! ** Inputs / Ouputs : (global commons)
        !!
        !! ** External : 
        !!
        !! ** References : Vancoppenolle et al., JGR 2007
        !!
        !! ** History :
        !!       (1) CLIO, Goosse and Fichefet, JGR, 1999.
        !!       (2) LIM-1D, Vancoppenolle et al., JGR, 2007.
        !!       (3) BIO-LIM, Martin Vancoppenolle, 2008
        !!
        !!------------------------------------------------------------------
        !! * Arguments

      INCLUDE 'type.com'
      INCLUDE 'const.com'
      INCLUDE 'para.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'dynami.com'
      INCLUDE 'reper.com'

      ! name of the experiment
      CHARACTER(len=8) exp_id

      ! Formats for reading the initial salinity profile
      CHARACTER(len=1) zc1
      CHARACTER(len=2) zc2
      CHARACTER(len=7) zformat1
      CHARACTER(len=7) zformat2
      CHARACTER(len=8) zformat3
! FD additions
      DOUBLE PRECISION latitude
!
!-----------------------------------------------------------------------
!  1 ) Lecture des parametres du run
!-----------------------------------------------------------------------
!
      WRITE(numout,*) ' * defcst : '
      WRITE(numout,*) ' ~~~~~~~~~~ '

      OPEN(unit=10,file='run.param',status='old')

      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*) exp_id     ! name of the experiment
      READ(10,*)
      READ(10,*) ddtb       ! time step
      READ(10,*)
      READ(10,*) nstart     ! number of the first iteration
      READ(10,*)
      READ(10,*) nend       ! number of the last iteration
      READ(10,*)
      READ(10,*) nyear1     ! initial year
      READ(10,*)
      READ(10,*) nfr_out    ! output frequency
      READ(10,*)
      READ(10,*) latitude   ! latitude
      READ(10,*) 

      ! to remove in the long run
      zlatz     = latitude
      premjour  = nstart
      dts(kmax) = ddtb

      nitrun    = nend - nstart + 1

      CLOSE(10)
!
!-----------------------------------------------------------------------
!  2 ) Mathematical constants
!-----------------------------------------------------------------------
!
      zero    = 0.d0              ! zero
      one     = 1.d0              ! one
      cstmin  = -1.d20            ! used in defgrid
      cstmax  =  1.d20            ! used in defgrid
 
      pi     = 4.0 * atan(one)
      radian = pi / 180.0
      omega  = 2.0 * pi / 86164.0 ! coriolis factor

!------------------------------------------------------------------------------
!  3 ) Physical constants
!------------------------------------------------------------------------------

      !-----------------------
      ! Fundamental constants
      !-----------------------
      gpes   = 9.80d0             ! gravity
      stefan = 5.6697d-08         ! stefan-boltzmann constant
      vkarmn = 0.40d0             ! von karmann constant
      cevap  = 2.5d+06            ! heat transfer coefficient for latent heat... check!!!
      zemise = 0.97d0             ! R, sea water emissivity (remove also)

      !---------------
      ! Ocean physics
      !---------------
      tpw     = 273.16d0          ! water triple point
      rho0    = 1025.0            ! ocean mean density
      cpw     = 3.99d+03          ! seawater specific heat
      oce_sal = 34.0

      visc_br = 1.79e-3           ! dynamic viscosity of water at 0C 
      beta_ocs= 0.78237           ! regulates density changes due to changes in salinity

      !-----------------
      ! Sea ice physics
      !-----------------
 
      !-- thermal properties
      tfsn   = 273.16d0           ! snow melting point
      tfsg   = 273.16d0           ! sea ice melting point
      xkn    = 0.31d0             ! 0.31d0 ISPOL Olivier Lecomte communication ! ref value 0.31d0 ! snow thermal conductivity
      xkg    = 2.034d0            ! pure ice thermal conductivity
      rhog   = 917.0              ! sea ice density
      rhon   = 330.d0             ! 355.d0 ISPOL Olivier Lecomte communication ! ref value 330.0  ! snow density
      cpg    = 2.062d+03          ! sea ice specific heat

      xlgm   = 3.335d+05          ! massive latent heat, ice
      xlgn   = 3.335d+05          !    "       "     " , snow
      xsn    = 2.834d+06          ! sublimation latent heat
      tmut   = 0.054d0            ! rate between seawater freezing point and
      betak1 = 0.09d0             ! first th.cond. constant
      betak2 = 0.011d0            ! second th.cond. constant
      emig   = 0.99d0             ! surface emissivity

      e_tres = 0.05

      !---------------------------
      ! Model physical parameters
      !---------------------------
      OPEN(unit=25,file='icephys.param', status='old')

      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*)
      READ(25,*) n_i              ! number of layers in the ice
      READ(25,*)
      READ(25,*) n_s              ! number of layers in the snow
      READ(25,*)
      READ(25,*) parsub           ! switch for sublimation or not
      READ(25,*)
      READ(25,*) tabq_ano         ! Prescribed air temperature anomalies
      READ(25,*)
      READ(25,*) flu_beta         ! fraction of meltwater percolating
      READ(25,*)
      READ(25,*) flu_bvtr         ! permeability threshold
      READ(25,*)
      READ(25,*) rad_inot_s       ! inot in snow
      READ(25,*)
      READ(25,*) rad_inot_i       ! inot in ice
      READ(25,*)
      READ(25,*) rad_kappa_s      ! attenuation coefficient in snow
      READ(25,*)
      READ(25,*) rad_kappa_i      ! attenuation coefficient in ice
      READ(25,*)
      READ(25,*) frtr_si_phy      ! fractionation coeff in snow ice
      READ(25,*)
      READ(25,*) d_br_mol         ! molecular diffusivity of brine
      READ(25,*)
      READ(25,*) d_br_tur         ! turbulent diffusivity of brine
      READ(25,*)
      READ(25,*) ra_c             ! critical rayleigh number over which convection starts
      READ(25,*)
      READ(25,*) ra_smooth        ! coefficient to smooth the hyperbolic tangential for ra_c
      READ(25,*)
      READ(25,'(a2)') gravdr      ! type of gravity drainage
      READ(25,*)
      READ(25,*) delta_cw         ! Cox and weeks gravity drainage parameter
      READ(25,*)
      READ(25,*) thcon_i_swi      ! conductivity formula switch
      READ(25,*)

      CLOSE(25)

!------------------------------------------------------------------------------
      RETURN
      END
