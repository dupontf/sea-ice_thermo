module ice_thermodynamic_FV

  use var_thermo_vertical

implicit none

!------------------------------------------------------------------------
!     Some model constants and parameters
!------------------------------------------------------------------------
real(8), dimension(0:maxlay) :: dzzi,dzzs


contains

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
subroutine init_sigma_ice_FV
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 implicit none
! locals
 integer k

! parameters for ice thermodynamics and heat fluxes
  real(8),save :: tzero  = 273.15_8       ! tzero C => K

!------------------------------------------------------------------------
!     Some integer constants
!------------------------------------------------------------------------

      ni=   nlice
      ns=ni+nlsno
      if (ns > maxlay ) then
       write(*,'(a,i4,a,i4)') ' ni',ni,' ns',ns
       write(*,*) 'please redefine the number of ice/snow levels in ice_thermoynamic.f90'
       stop
      endif

!---------------------------------------------------------------------
! sin distribution for ice
! zi varies from 0 (bottom) to 1 (surface)
!---------------------------------------------------------------------

! pi=atan(1.0_8)*4.0_8
! do k=0,maxlay-1
!    zi(k)=0.5_8*(1.0_8+sin(pi*(dble(k)/dble(maxlay-1)-0.5_8)))
! enddo

!---------------------------------------------------------------------
! uniform distribution for ice
! zi varies from 0 (bottom) to 1 (surface)
!---------------------------------------------------------------------

!for ice
 do k=0,ni-1
    zi(k) = dble(k)/dble(nlice)
 enddo
   zi(ni) = 1.0_8
 do k=ni+1,ns-1
    zi(k) = 1.0_8 + dble(k-ni)/dble(nlsno)
 enddo
   zi(ns) = 2.0_8
! for snow
   zs(0:ni) = 0.0_8
 do k=ni+1,ns-1
    zs(k) = dble(k-ni)/dble(nlsno)
 enddo
   zs(ns) = 1.0_8

! cell oriented coordinate for ice
   cs(1,1:ni) = zi(0:ni-1) ! bottom z
   cs(2,1:ni) = zi(1:ni  ) ! top z
! cell oriented coordinate for snow
   cs(1,ni+1:ns) = zs(ni:ns-1) ! bottom z
   cs(2,ni+1:ns) = zs(ni+1:ns) ! top z

! add a new variable 1-zs for velocity calculation (only non zero for the ice layers)
   os(:,:) = 1.0_8 - cs(:,:)

!for ice
 dzzi(0:ns)=0.0_8
 do k=1,ni
    dzi(k) = 1.0_8/dble(nlice)
    dzzi(k) = 1.0_8/dble(nlice)
 enddo
 dzzi(0 ) = 0.5_8/dble(nlice)
 dzzi(ni) = 0.5_8/dble(nlice)

!for snow
 dzzs(0:ns)=0.0_8
 dzzs(ni) = 0.5_8/dble(nlsno)
 do k=ni+1,ns
    dzi(k) = 1.0_8/dble(nlsno)
    dzzs(k) = 1.0_8/dble(nlsno)
 enddo
 dzzs(ns) = 0.5_8/dble(nlsno)

!---------------------------------------------------------------------
! other quantities
!---------------------------------------------------------------------

  ziold =  zi
 dziold = dzi

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
end subroutine init_sigma_ice_FV

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ice_thermodynamic.f90 file
!
! contain subroutines: - ice_thermo

! FD summer 2016: mushy layer physics included

! Rewritten by F. Cyr with few changes:
!         -> in ice_thermo
!           - new param: nuu, a, b, Smax 
!           - new variable: Cbrine, zn, salin
!           - new melting scheme (brine pockets) (not usefull for now!)
!           - new lateral melting scheme and corresponding parameters&variables
!         -> in bottom_ice_creation
!           - new input var: "wind"
!           - new variables: dhi_underice, dhi_openwater, hi_accr
!           - new param Cfa, Cfw & cc
!           - brand new ice_creation scheme!
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!        SUBROUTINE ice_thermo      
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine ice_thermo(dtice)

  implicit none
!arguments
  real(8)       dtice

!locals
!===============================================================================
!     Parameter statement defines layers and integration periode
!===============================================================================
!....&S..1.........2.........3.........4.........5.........6.........7..C......8

      real(8) :: &
        Fbase,&! heat flux at ice base			[W/m2]
        S    ,&! initial salinity profile		[psu]
        T    ,&! initial temp. profile			[C]
        zini(0:maxlay),&! z-values of init interpol tprof 	[m]
        Tsurf           ! surface boundary conditions           [C]

      CHARACTER (len=10) :: iatflx, ocnflx, ocnsal, ocntmp, uiorel, &
                            salinity, bbc, sbc

      real(8) :: uiofix, & ! constant relative ice-ocean velocity		[m/s]
        sigma, &! Stefan Bolzmann constant			[W/m/K4]
        temp0, &! melt temperature of snow and ice		[K]
        tiny,  &! very small number (1d-9)			[-]
        Tdiff    ! temperature difference in Euler step		[C]

      real(8) :: &
     epsilon,&! coefficient of emissivity			[-]
     hicut,&! cut off (minimum) ice  thickness		[m]
     hscut,&! cut off (minimum) snow thickness		[m]
     hslim  ! minimum snow thickness			[m]

      logical   Tsbc           ! Fixed T surface boudnary conditions

      real(8) :: &
     dhidt,&! total rate of change of ice  thickness	[m/s]
     dhsdt,&! total rate of change of snow thickness	[m/s]
     dspdt,&! rate of change of snow thickn. due thin_snow_melting [m/s]
     dssdt,&! rate of change of snow thickn. at the surface [m/s]
     dsbdt,&! rate of change of snow thickn. @ the sno base [m/s]
     disdt,&! rate of change of ice  thickn. at the surface [m/s]
     dibdt,&! rate of change of ice  thickn. @ the ice base [m/s]
     dTlay,&! Tnew-Told of a layer				[C]
     dTmax,&! max temp difference in total ice slab		[C]
     dEin ,&! difference in internal energy			[J/m2]
     Ein0 ,&! initial internal energy			[J/m2]
     Einp ,&! energy input per time step			[J/m2]
     Elay ,&! internal energy in a layer			[J/m2]
     Fcss ,&! cond. heat flux at snow surface		[W/m2]
     Fcsb ,&! cond. heat flux at snow base			[W/m2]
     Fcis ,&! cond. heat flux at ice  surface		[W/m2]
     Fcib ,&! cond. heat flux at ice  base			[W/m2]
     Fcsu ,&! cond. heat flux at surface		[W/m2]
     Fsh  ,&! sensible heat flux at ice surf.		[W/m2]
     Flh  ,&! latent   heat flux at ice surf.		[W/m2]
     Fprec,&! heat flux due to precipitation		[W/m2]
     Frain,&! heat flux due to warm rain			[W/m2]
     Fnet ,&! net atmos. energy flux at surf.		[W/m2]
     Fnet0,&! net atmos. energy flux at surf.		[W/m2]
     Fmlt ,&! heat flux at surf equiv to melt		[W/m2]
     Frad ,&! sw radiation heat flux (surf minus base)	[W/m2]
     gbase,&! total growth  at the ice base			[m]
    mbase,&! total melting at the ice base			[m]
    msurf,&! total melting at the ice surface		[m]
    qice ,&! specific humidity  at ice surf.		[kg/kg]
    Eint ,&! internal energy in all layers		[J/m2]
    Fr(0:ns+1)      ,&! sw rad heat flux at s/i internal level 	[W/m2]
    ki(0:ns+1)      ,&! ice  thermal conductivity			[W/m/C]
    ks(0:ns+1)      ,&! snow thermal conductivity			[W/m/C]
    kki(0:ns+1)     ,&! ice  thermal conductivity			[W/m/C]
    kks(0:ns+1)     ,&! ice  thermal conductivity			[W/m/C]
    em(0:ns+1)      ,&! volumetric enthalpiy	[W/m3]
    em0(0:ns+1)     ,&! old volumetric enthalpiy	[W/m3]
    eb(0:ni+1)      ,&! brine enthalpiy	[W/m3]
    eb0(0:ni+1)     ,&! old brine enthalpiy	[W/m3]
    cp(0:ns+1)      ,&! heat capacity for ice	[W/kg/K^-1]
    dedp(0:ns+1)    ,&! heat capacity for ice	[W/kg/K^-1]
    dsidt(0:ns+1)   ,&! 1st derivate of liquidus salinity relative to temperature
    R(0:ns+1)       ,&! penetrating shortwave radiation		[W/m3]
    sali(0:ns+1,0:1),&! internal sea ice salinity			[psu]
    siold(0:ns+1)   ,&! initial ice salinity			[psu]
    phib(0:ns+1,0:1),&! brine fraction			[1]
    sibr(0:ns+1,0:1),&! brine salinity			[psu]
    Tf(0:ns+1)      ,&! ice freezing temp. of salinity S		[C]
    w(0:ns)         ,&! grid advection velocity			[m/s]
    wb(0:ni)        ,&! convective brine velocity			[m/s]
    rho(0:ns+1)    , & ! density   [kg/m3]
    um(ns)         , &   ! melting speed (m/s)
    ub(ni)         , &   ! brine channel lateral speed (m/s)
    hi_b(0:1)       ,&! ice thickness (m)
    hs_b(0:1)       ,&   ! snow thickness (m)
    wsg(0:ns)       ,&! sign of vertical velocity
    wog(0:ns)       ,&! opposed sign of vertical velocity
    wsb(0:ns)       ,&! sign of vertical brine velocity
    wob(0:ns)         ! opposed sign of vertical brine velocity
      real(8), dimension(ns) :: &
             he, heold, henew

      real(8) temp (0:ns+1,0:1) ! internal sea ice temperature		[C]

! full matrix
      real(8) matj(0:ns+ni+1,0:ns+ni+1) ! jacobian matrix		[C]

! LU Solver START
      real(8) matband(3*(ns+ni+2),ns+ni+2) ! band matrix
      integer ntot,ml,mu,mt,info,bandmax
      integer ipvt(ns+ni+2)
      real(8) work(ns+ni+2)
      real(8) rcond,residual
! LU Solver END

! tri-diagonal solver for temperature only
!      real(8) AT1(0:ns+1), BT1(0:ns+1), CT1(0:ns+1)
      real(8) DT0(0:ns+ni+1), tout(0:ns+ni+1), Rtrans

      real(8) :: sice, &
             latmelt, ftot, sstznew, &
             Tside,deltaT,wlat,rside, &
             sum0, dhi, dho, dq, fwf0
      real(8) :: hminice=1d-3
      real(8) :: energy_bot

  character (len=20) :: Sprofile
  integer i,j,counter,k
  real(8) dzf, es, zssdqw, zrchu1, zrchu2, q0, zref, k0, k1
  real(8), dimension(0:ns+1) :: tiold, tinew
  real(8) elays0
  logical :: debug=.false.
  logical :: extra_debug=.false.
  real(8) :: hsold, tsuold, sfallold
  logical :: thin_snow_active
  real(8) :: snow_precip, rain_precip, fwf, &
                      fthin_snow, energy_snow_melt, sume2, em_thin_snow
! parameters for ice thermodynamics and heat fluxes
  real(8),save :: tzero  = 273.15_8       ! tzero C => K

  real(8),save :: latvap = 2.5e6_8        ! J/kg latent heat of vaporization
  real(8),save :: latsub = 2.834e6_8      ! J/kg latent heat of sublimation
  logical :: adv_upwind = .false.
  logical :: adv_upwind_br = .false.
  logical :: adv_br = .false.
  logical, dimension(ns) :: internal_melt
  integer, dimension(ns) :: im
  integer jm,jk,nm,js
  real(8) :: sumum, sumsal, dsaldt, maxvel
  real(8) liquifrac_bot, minfraliq
  real(8) :: &
      salinis	=   1.0_8, &    ! salinity at ice surf			[psu]
      salinib	=  10.0_8       ! salinity at ice base			[psu]

! permeability arrays
  real(8), dimension(ni) :: perm, & ! permeability at each ice layer (center)
                              ra, & ! local Raileigh number (cell center)
                             zzp    ! mid-cell z height (m)
  logical lstop_minfraliq, lstop_heat, lstop_salt
  real(8) sumhmelt ! sum of all interior melting velocity
  character, external :: real_to_char

  lstop_minfraliq = .false.
  lstop_heat      = .false.
  lstop_salt      = .false.

! FD debug
!  debug=.true.
!  extra_debug=.true.
!  adv_br = .true.
  adv_upwind_br = .true.
  adv_upwind = .true.
      hicut=1.e-2_8
      bandmax=3*(ns+ni+2) ! need to defined as in the variable declaration

! initialize to false the internal melt
      internal_melt(:) = .false.
      um(:)=0.0_8
      ub(:)=0.0_8
      wb(:)=0.0_8
      fsbr = 0.0_8
      fhbr = 0.0_8

!------------------------------------------------------------------------
!     condition on ice thickness
!     no thermo below 5 cm
!------------------------------------------------------------------------

!      if (hi.lt.hicut) return

!------------------------------------------------------------------------
!     Some logical constants
!------------------------------------------------------------------------

      salinity  = 'yes'          ! = no dynamic salinity
      bbc       = 'fixT'        ! bottom boundqry condition

!------------------------------------------------------------------------
!     Some physical constants
!------------------------------------------------------------------------

      epsilon=   0.990_8! snow/ice emissivity			[-]

      sigma= 5.670e-8_8 ! Stefan-Boltzmann constant	       [W/m2/K4]
      temp0= 273.160_8  ! freezing point temp of fresh water	[K]

      uiofix=   0.000_8 ! constant relative ice-ocean velocity	[m/s]
      liquifrac_bot = 0.4_8
      iliquid_cond = 2
      ith_cond=1
      sibr(0,0) = func_liqu_sa(tbo - temp0)
      salinib = liquifrac_bot * sibr(0,0)

! FD debug test constant salinity
!      salinib = 4.0_8
!      si(:) = 4.0_8
!      liquifrac_bot = salinib / sibr(0,0)

! FD debug
!dtice=3600.0_8
!hs=0.0_8
!ts(:)=tsu
!oceflx=0.0_8
!swradab_s=0.0_8
!swradab_i=0.0_8
!snowfall=0.0_8
!flat=0.0_8
!fsens=0.0_8
!swrad=400.0_8
!dwnlw=300.0_8
!hs=0.0_8
!ti=tbo
!ts=tbo
!tsu=temp0
      Fbase	=   oceflx	! conductive heat flux at ice base 	[W/m2]


!because of potential change in steric height due to the presence of the ice, 
!we assume that the ice is not impacting the volume of water, hence
!the flux of volume between the ice and the ocean is set to zero,
!only the salt flux is impacted as a pseudo-volume flux

! check that this is actually the case in test_bio.f90 and forcing.F90

!-----------------------------------------------------------------
! ice creation in case of near freezing ocean temperature
!-----------------------------------------------------------------


!------------------------------------------------------------------------
!     Some constants for the numerics
!------------------------------------------------------------------------

      tiny	=   1.0e-9_8	! very small number			[-]

      hslim     =   5.0e-4_8    ! limit for computing temperature in snow
      Tdiff     =   1.0e-12_8   ! temperature tolerance in Euler step	[C]
      thin_snow_active =.false. ! reroute some energy to sublimate/melt snow if true
      fwf = 0.0_8
! FD debug
!      hslim     =   0.005_8    ! limit for computing temperature in snow
!      hslim     =   0.01_8    ! limit for computing temperature in snow

!---------------------------------------------------------------------
! other quantities
!---------------------------------------------------------------------

  zi =  ziold
 dzi = dziold

!---------------------------------------------------------------------
! conversion from LIM3
!---------------------------------------------------------------------
      tiold(0)=tbo - temp0
      do j=1,nlice
        tiold(j)=ti(nlice-j+1)  - temp0
      enddo
      do j=1,nlsno
        tiold(j+ni)=ts(nlsno-j+1)  - temp0
      enddo
      Tsurf = tsu - temp0
      tiold(ns+1) = Tsurf

!------------------------------------------------------------------------
!     Euler method initialization
!------------------------------------------------------------------------

      dTlay	=   0.0_8	! temp diff in a layer	[C]
      dTmax	=   0.0_8	! max temp difference	[C]
      
!-----------------------------------------------------------------------
!     Initial salinity profile
!-----------------------------------------------------------------------

      siold(0)=salinib ! but will require some prognostic variable too such as si(nlice+1)
      do j=1,nlice
        siold(j)=si(nlice-j+1)
      enddo
      DO j=ni+1,ns+1
         siold(j) = 0.0_8
      ENDDO
      sali(:,0) = siold(:) ! keep a trace of initial salinity
      sali(:,1) = siold(:) ! keep a trace of initial salinity

      hsold = hs
      tsuold = tsu
      sfallold = snowfall
      fthin_snow = 0.0_8
      dspdt = 0.0_8

1000 continue
!------------------------------------------------------------------------
! init variable
!------------------------------------------------------------------------

      hs_b(0:1)=hs
      hi_b(0:1)=hi
      heold(1:ni)=hi*dzi(1:ni)
      heold(ni+1:ns)=hs*dzi(ni+1:ns)
      he=heold
      henew=he


      tocn  = tiold(0)
      temp(:,0)=tiold
      temp(:,1)=tiold

!------------------------------------------------------------------------
! density
!------------------------------------------------------------------------

      rho(0)=rhoice
      rho(1:ni)=rhoice
      rho(ni+1:ns+1)=rhosno
      
!------------------------------------------------------------------------
!     Vertical (grid) velocities
!------------------------------------------------------------------------

      DO j=0,ns
         w(j) = 0.0_8		! initial grid advection
         wsg(j) = 0.0_8         ! sign of w
         wog(j) = 1.0_8         ! sign of w
      ENDDO
      wsb(:) = 0.0_8
      wob(:) = 0.0_8

!------------------------------------------------------------------------
!     Set initial T-S dependent material properties
!------------------------------------------------------------------------

      DO j=0,ns+1
            Tf(j)     = 0.0_8
            ki(j)     = func_ki(sali(j,0),tiold(j))
            sibr(j,:) = func_liqu_sa(tiold(j))
            phib(j,:) = 0.0_8
            if (sibr(j,0) > 0.0_8) phib(j,:) = sali(j,0) / sibr(j,0)
            em(j)     = func_El_mush(tiold(j),phib(j,0))
            cp(j)     = func_cp_mush(tiold(j),phib(j,0))
            dedp(j)   = func_dedp(tiold(j),phib(j,0))
            dsidt(j)  = func_liqu_sadt(tiold(j))
      ENDDO
      em0(0:ns+1) = em(0:ns+1)

! brine enthalpy
      DO j=0,ni
         eb(j) = cp_wat * tiold(j)
      ENDDO
      eb0(0:ni) = eb(0:ni)

!------------------------------------------------------------------------
!     Set conductivity
!------------------------------------------------------------------------

      DO j=0,ns+1					! snow thermal conductivity
         IF (j .LT. ni) THEN
            ks(j) = 0.0_8
         ELSE
            ks(j) = cond_sno
         ENDIF
      ENDDO
      kki(0)=ki(0)/ (hi_b(0) * dzzi(0))
      DO j=1,ni-1
         kki(j) = (ki(j+1)+ki(j)) * 0.5_8 / (hi_b(0) * dzzi(j))
      ENDDO
      kki(ni) = ki(ni) / (hi_b(0) * dzzi(ni))
      kks(0:ns)=0.0_8
      kks(ni) = ks(ni)*ki(ni)/max(tiny, ki(ni)*dzzs(ni)*hs_b(0) + ks(ni)*dzzi(ni)*hi_b(0) )
      DO j=ni+1,ns-1
         kks(j) = (ks(j+1)+ks(j)) * 0.5_8 / (hs_b(0) * dzzs(j))
      ENDDO
      kks(ns) = ks(ns) / (hs_b(0) * dzzs(ns))
      kki(ni:ns)=kks(ni:ns)

!------------------------------------------------------------------------
!     Set permeability for brine advection 
!------------------------------------------------------------------------

if (adv_br) then
 i_perm_eff=1
 i_perm_for=1
 i_Ra=2 ! GN2013
      call ice_permeability(ni,phib(1,0),perm,dzi)
      do k=1,ni
        zzp(k) = 0.5_8 * (zi(k-1)+zi(k)) * hi
        call ice_rayleigh_local(sibr(k,0),seasal,perm(k),ra(k),zzp(k))
      enddo
      call vertical_brine_velocity_GN(ni,Ra,heold,ub,wb)
! FD debug
if (debug) then
 maxvel = abs(wb(0)) / heold(1)
 do j=1,ni
    maxvel = max(maxvel, abs(wb(j)) / heold(j), abs(ub(j)) / heold(j) )
 enddo
 write(*,*) 'Courant number for brine advection', maxvel * dtice
endif
 if (extra_debug) then
 write(*,*) 'vertical brine velocity'
 write(*,'(300(e10.3))') ub(1:ni)
 write(*,'(300(e10.3))') wb(0:ni)
! write(*,'(300(e10.3))') he(1:ni)/dtice
 endif
endif

!------------------------------------------------------------------------
!     add flushing meltwater above a certain threshold
!------------------------------------------------------------------------

      call reject_brine_velocity(ni,phib(1,0),heold,um,dtice)
      sumhmelt = 0.0_8
      do j=1,ni
         sumhmelt = sumhmelt + um(j)
      enddo

!------------------------------------------------------------------------
!     Initial internal energy in ice (1:ni-1) and snow (ni+1:ns-1)
!------------------------------------------------------------------------

      Einp = 0.0_8				! E input [J/m2]
      Ein0 = 0.0_8

      DO j=1,ns
         Elay = em(j)*rho(j)*heold(j)
         Ein0 = Ein0+Elay
      ENDDO

      dEin = 0.0_8				! diff intl E [J/m2]

      sumsal=0.0_8
      DO j=1,ni
         sumsal = sumsal + heold(j) * sali(j,0)
      ENDDO
      dsaldt= 0.0_8

      ! net longwave radiative flux
      netlw = emi*(dwnlw - stefa*tsu*tsu*tsu*tsu)
      ! sensible and latent heat flux
!      fsens   =  -fsens
! FD debug
if (extra_debug) then
 write(*,*) 'fsens,flat',fsens,MIN( flat , 0.0_8 ),dwnlw,netlw
endif
      flat   =  MIN( flat , 0.0_8 ) ! always negative, as precip 
                                           ! energy already added

      ! derivative of the surface atmospheric net flux
          dzf    =  4.0_8*emi*stefa*tsu*tsu*tsu

      ! surface atmospheric net flux
      Fnet0 =  fac_transmi * swrad + netlw + fsens + flat
      Fnet  = Fnet0

!---------------------------------------------------------------------
! accumulation at the surface
!---------------------------------------------------------------------
     if ( tair < tzero ) then
      snow_precip = snowfall
      rain_precip = 0.0_8
     else
      snow_precip = 0.0_8
      rain_precip = snowfall * rhosno / rhowat
     endif

!---------------------------------------------------------------------
! need to get rid of snow if too thin during melting/sublimating
! thin snow case
!---------------------------------------------------------------------
  if ( hs_b(0) <= hslim ) then

    fthin_snow = 0.0_8
    DO j=ni+1,ns
       Elay = em0(j)*rho(j)*heold(j)
       fthin_snow = fthin_snow + Elay
    ENDDO
    fthin_snow = - fthin_snow / dtice
    do j=ni+1,ns
      heold(j)=0.0_8
    enddo
    hs_b(0)=0.0_8
    thin_snow_active=.true.
    tiold(ni+1:ns) = tiold(ns+1)
    em_thin_snow = func_el(Tf(ns+1),tiold(ns+1))
    Tf(ni+1:ns+1) = Tf(ni)
    sali(ni+1:ns+1,0) = sali(ni,0)
    sali(ni+1:ns+1,1) = sali(ni,1)
    rho(ni+1:ns+1)=rhoice
    kki(ni) = ki(ni) / (hi_b(0) * dzzi(ni))
    cp(ni+1)   = func_cp_mush(tiold(ni+1),phib(ni+1,0))
    dedp(ni+1) = func_dedp(tiold(ni+1),phib(ni+1,0))

  endif ! end thin snow
!---------------------------------------------------------------------

      Fnet0 = Fnet0 - fthin_snow
      Fnet = Fnet - fthin_snow

!------------------------------------------------------------------------
!     Update conductive heat fluxes
!------------------------------------------------------------------------

      
      Fcis = -kki(ni) * ( temp(ns+1,1) - temp(ni  ,1))
      Fcib = -kki(0 ) * ( temp(   1,1) - temp(   0,1))

      IF ( .not.thin_snow_active ) THEN
         Fcss = -kki(ns) * ( temp(ns+1,1) - temp(ns  ,1) )
         Fcsb = -kki(ni) * ( temp(ni+1,1) - temp(ni  ,1) )
      ELSE
         Fcss = Fcis
         Fcsb = 0.0_8
      ENDIF
      energy_bot = em(0)

!-----------------------------------------------------------------------
!     ICE/SNOW case
!-----------------------------------------------------------------------


      Tsbc = .FALSE.
      counter = 0

!-----------------------------------------------------------------------
!     Euler loop
!-----------------------------------------------------------------------

4000  CONTINUE

!-----------------------------------------------------------------------
!     Update atm
!-----------------------------------------------------------------------

      Fnet = Fnet0 - dzf * (temp(ns+1,1)-tiold(ns+1)) ! adjust net surface flux
      fwf = fwf + rain_precip ! add liquid precipitation to the freshwater flux to ocean

! FD debug
if (extra_debug) write(*,*) 'total E/dt, Fnet',Ein0/dtice,Fnet+oceflx
! FD         oceflx= rhoo*cpo*Coi*ABS(uio)*(tocn-temp(0,1))
!-----------------------------------------------------------------------
!     Update snow and ice thickness
!-----------------------------------------------------------------------

      dssdt = snow_precip + dspdt
      dsbdt = 0.0_8		! melt of snow base not allowed

      disdt = 0.0_8

         ! snow depth evolution due to precipitation


! FD: if snow present
      if ( .not.thin_snow_active ) then
         ! snow depth evolution due to melt superposed on precip

         IF (temp(ns+1,1)+tiny .GE. Tf(ns+1) .AND. Fnet+Fcss .GT. 0.0_8) THEN
! FD the code cannot deal with concomittent snow accumulation and melt
            Tsbc = .true.
            dssdt = MIN((Fnet+Fcss)/rhosno/em(ns+1),0.0_8)
            snow_precip = 0.0_8
         ENDIF
      else

      ! update ice surface

         IF (temp(ns+1,1)+tiny .GE. Tf(ns+1) .AND. Fnet+Fcis .GT. 0.0_8) THEN
            Tsbc = .true.
            disdt = MIN((Fnet+Fcis)/rhoice/em(ni),0.0_8)
         ENDIF
      endif

! bottom growth or melt
      dibdt	= - (oceflx-Fcib) / rhoice / energy_bot
! FD debug
if (extra_debug) write(*,*) 'bottom budget',oceflx,Fcib,energy_bot

      dhsdt	= dssdt - dsbdt
      dhidt	= disdt - dibdt - sumhmelt ! original

      hs_b(1)	= hs_b(0) + dhsdt*dtice

!---------------------------------------------------------------------
! need to get rid of snow if wmesh too large
!---------------------------------------------------------------------

   if ( hs_b(1) < 0.0_8 .and. .not.thin_snow_active ) then
      sume2 = 0.0_8
      do j=ni+1,ns
        sume2 = sume2 + em0(j) * rho(j) * heold(j)
      enddo
      energy_snow_melt = -sume2
      fthin_snow = energy_snow_melt / dtice
! FD debug
if (extra_debug) write(*,*) 'switch to thin_snow_active',hs_b(1),fthin_snow,Fnet
      fwf = fwf + hs_b(0) / dtice / rhowat * rhosno
      Tsbc = .false.
      thin_snow_active = .true.
      snow_precip = 0.0_8
      em_thin_snow = 0.0_8
      hs_b(1) = 0.0_8
      dspdt = -hs_b(0) / dtice
      dssdt = dspdt
      Fnet0 = Fnet0 - fthin_snow
      Fnet = Fnet - fthin_snow
      kki(ni) = ki(ni) / (hi_b(0) * dzzi(ni))
      sali(ni+1:ns+1,0) = sali(ni,0)
      sali(ni+1:ns+1,1) = 0.0_8
      temp(ni+1:ns+1,1) = temp(ni,1)
      phib(ni+1:ns+1,0) = 0.0_8
      phib(ni+1:ns+1,1) = 0.0_8
      sibr(ni+1:ns+1,0) = 0.0_8
      sibr(ni+1:ns+1,1) = 0.0_8
      em(ni+1:ns+1)     = em(ni)
      ! update ice surface

      Fcis = -kki(ni) * ( temp(ns+1,1) - temp(ni  ,1))
         IF (temp(ns+1,1)+tiny .GE. Tf(ns+1) .AND. Fnet+Fcis .GT. 0.0_8) THEN
            Tsbc = .true.
            disdt = MIN((Fnet+Fcis)/rhoice/em(ni),0.0_8)
         ENDIF
      dhidt	= disdt - dibdt ! original
      Fcss = Fcis

   endif

      hi_b(1)	= hi_b(0) + dhidt*dtice

      henew(   1:ni  )=dzi(   1:ni  )*hi_b(1)
      henew(ni+1:ns  )=dzi(ni+1:ns  )*hs_b(1)

!-----------------------------------------------------------------------
!     Vertical velocities (grid advection)
!-----------------------------------------------------------------------

      w(0 ) = - dibdt
      w(ni) = - disdt
!      DO j=1,ni-1
!         w(j)	= cs(2,j) * w(ni) + os(2,j) * w(0)
!      ENDDO
      w(ns) = - dssdt
!      DO j=ni+1,ns-1
!         w(j)	= cs(2,j) * w(ns)
!      ENDDO
! add contribution due to interior melting
! compute total loss of mass
      sumum = 0.0_8
      DO j=1,ni
         sumum = sumum + um(j)
      ENDDO
! recompute vertical velocity in ice
      DO j=1,ni-1
         w(j)	= w(j-1) - um(j) + dzi(j) * ( sumum + w(ni)-w(0 ) )
      ENDDO
! compute total loss of mass in snow
      sumum = 0.0_8
      DO j=ni+1,ns
         sumum = sumum + um(j)
      ENDDO
! recompute vertical velocity in ice
      DO j=ni+1,ns-1
         w(j)	= w(j-1) - um(j) + dzi(j) * ( sumum + w(ns)-w(ni) )
      ENDDO

! compute coefficient for upstream or centered advection
     if (adv_upwind) then
      DO j=1,ns-1
         if (w(j) .gt. 0.0_8) wsg(j) = 1.0_8
         wog(j) = 1.0_8 - wsg(j)
      ENDDO
     else ! centered
      DO j=1,ns-1
         wsg(j) = 0.5_8
         wog(j) = 0.5_8
      ENDDO
     endif

     if (adv_upwind_br) then
! always upwind for brine advection
      DO j=1,ni-1
         if (wb(j) .gt. 0.0_8) wsb(j) = 1.0_8
         wob(j) = 1.0_8 - wsb(j)
      ENDDO
     else ! centered
      DO j=1,ni-1
         wsb(j) = 0.5_8
         wob(j) = 0.5_8
      ENDDO
     endif

! bottom and top cases: always upwind
     j=0
         if (w(j) .gt. 0.0_8) wsg(j) = 1.0_8
!         wsg(j) = 1.0_8 ! always up
         wog(j) = 1.0_8 - wsg(j)

!brine
         if (wb(j) .gt. 0.0_8) wsb(j) = 1.0_8
         wob(j) = 1.0_8 - wsb(j)
     j=ns
!         wsg(j) = 0.0_8
!         if (w(j) .gt. 0.0_8) wsg(j) = 1.0_8
!         wog(j) = 1.0_8 - wsg(j)
! for compability with previous versions, when snow is melting we still use the surface value, and if ice is melting we use the upwind value.
         wsg(j) = 0.0_8
         wog(j) = 1.0_8
     j=ni
         wsg(j) = 1.0_8
         wog(j) = 0.0_8

!-----------------------------------------------------------------------
!     Absorbed SW radiation energy per unit volume (R)       
!     1/dz Int_z^z+dz R dz = Rtrans(z+dz) - Rtrans(z)
!-----------------------------------------------------------------------

      if ( .not.thin_snow_active ) then
        DO j=ns,ni+1,-1
          R(j)   = swradab_s(ns-j+1)
        ENDDO
      else
         R(ni+1:ns)=0.0_8
      endif
      R(ni)=0.0_8
      DO j=ni,1,-1
         R(j)   = swradab_i(ni-j+1)
      ENDDO
      R(0)=0.0_8

!-----------------------------------------------------------------------
!     fix bottom liquid fraction in some way
!-----------------------------------------------------------------------

      j=0
        phib(j,1) = liquifrac_bot ! fixed liquidus fraction at bottom if ingoing
        sibr(j,1) = func_liqu_sa(temp(j,1))
        sali(j,1) = phib(j,1) * sibr(j,1)


!-----------------------------------------------------------------------
!     matrix coefficients [matj] {x} = {D}
!-----------------------------------------------------------------------

       if ( thin_snow_active ) then
          ntot=2*ni+2
          Fcsu = Fcis
          ns = ni ! for simplification, ns becomes the total number of active ice(+snow) layers
       else
          ntot=ni+ns+2
          Fcsu = Fcss
       endif
! FD debug
if (extra_debug) write(*,*) 'total unknowns',ntot,thin_snow_active
       matj(0:ntot-1,0:ntot-1) = 0.0_8


! the unknown is w(0)
         k0 = kki(0)
         DT0(0)    = - w(0) * rho(0) * energy_bot + oceflx - Fcib ! line 0 is for w(0)
         matj(0,0) =   rho(0) * energy_bot       ! w increment
         matj(1,0) = - k0                        ! dFcib/dT increment

! implicit term in temperature for BC bottom (derivation of energy_bot)
         matj(1,0) =  matj(1,0) &
                     + w(0) * rho(0) * cp(1) * wog(0)

! liquid fraction terms for BC bottom (derivation of energy_bot)
      j=1
         jm=j+ns+1
         matj(jm ,0) = matj(jm ,0) &
                     + w(0) * rho(0) * dedp(j) * wog(0)

! for normal cells
      DO j=1,ns
         k0 = kki(j-1)
         k1 = kki(j  )
         DT0(j) =  R(j) + ( & ! radiation
                          - rho(j)*(henew(j)*em(j)-heold(j)*em0(j)))/dtice & ! time tendency energy
                        + k0 * ( temp(j-1,1) - temp(j,1) ) &
                        + k1 * ( temp(j+1,1) - temp(j,1) )
         if (j>1) &
         matj(j-1,j) = -k0
         matj(j  ,j) =  k0 + k1 + rho(j) * cp(j) * henew(j) / dtice
         if (j<ns .or. .not.Tsbc) &
         matj(j+1,j) = -k1
      ENDDO

! top condition, assuming at this stage that the free variable is temperature (i.e., no melt)
      j=ns+1
         k0 = kki(j-1)
         matj(j-1,j) = - k0
         matj(j  ,j) =   k0 + dzf
         DT0(j)      =   Fnet + Fcsu

! liquid fraction terms
      DO j=1,ni
         jm=j+ns+1
         matj(jm ,j) =          rho(j) * dedp(j) * henew(j) / dtice
      ENDDO

! reset terms depending on velocity (j=0 and/or j=ns+1)
! Newton terms for change in thickness in rho h E: 
!    d(rho h E)/d(w0)= rho E dh/dw0 = rho E ds; d(h)/d(ws)=-ds
      DO j=1,ni
         matj(0   ,j) =   rho(j) * em(j) * dzi(j)
      ENDDO
      if ( ni==ns .and. Tsbc ) then
       DO j=1,ni
         matj(ns+1,j) = - rho(j) * em(j) * dzi(j)
       ENDDO
      endif
      if ( ns>ni .and. Tsbc ) then
       DO j=ni+1,ns
         matj(ns+1,j) = - rho(j) * em(j) * dzi(j)
       ENDDO
      endif

      j=ns+1
      IF (Tsbc) THEN
! solve for delta w(ns)
         DT0(j)      = + w(j-1) * rho(j-1) * em(j-1) * wsg(j-1) &
                       + w(j-1) * rho(j  ) * em(j  ) * wog(j-1) &
                       + Fnet + Fcsu
         matj(j  ,j) = -          rho(j-1) * em(j-1) * wsg(j-1) &           ! w increment
                       -          rho(j  ) * em(j  ) * wog(j-1)             ! w increment
         matj(j-1,j) = - k0 &                                               ! dFcss/dT increment
                       - w(j-1) * rho(j-1) * cp(j-1) * wsg(j-1)             ! + variation due to top cell temp variation
      ENDIF

! liquid fraction terms for BC top (derivation of energy_bot)
      j=ns+1
      IF (Tsbc .and. ni==ns) THEN
         jm=j+ns ! +1-1
         matj(jm ,j) = matj(jm ,j) &
                       - w(j-1) * rho(j-1) * dedp(j-1) * wsg(j-1)             ! + variation due to top cell temp variation
      ENDIF

!-----------------------------------------------------------------------
!     add transport to the equations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! lower transport in ice
!-----------------------------------------------------------------------

          do j=1,ns
            DT0(j)      = DT0(j)      + rho(j-1) * w(j-1) * em(j-1) * wsg(j-1)
            DT0(j)      = DT0(j)      + rho(j  ) * w(j-1) * em(j  ) * wog(j-1)
          enddo

          j=0
            matj(j+1,j) = matj(j+1,j) + rho(j  ) * w(j  ) * cp(j  ) * wog(j  ) ! condition at bottom in case of upwind advection
          j=1
            matj(j  ,j) = matj(j  ,j) - rho(j  ) * w(j-1) * cp(j  ) * wog(j-1)
! starts at j=2 (not j=1) because we know temperature is not the unknown at j-1=0 and so requires some special treatments
          do j=2,ns
            matj(j-1,j) = matj(j-1,j) - rho(j-1) * w(j-1) * cp(j-1) * wsg(j-1)
            matj(j  ,j) = matj(j  ,j) - rho(j  ) * w(j-1) * cp(j  ) * wog(j-1)
          enddo

! metric terms, newton terms for change in thickness d(rho w E)/d(w0)
! bottom velocity
           do j=1,ni
            matj(0   ,j) = matj(0   ,j) - rho(j-1) * em(j-1) * os(1,j) * wsg(j-1)
            matj(0   ,j) = matj(0   ,j) - rho(j  ) * em(j  ) * os(1,j) * wog(j-1)
           enddo

! top velocity
          if ( Tsbc ) then ! top velocity: ice or snow?
             if (ns==ni) then ! definitily ice
                js=2
             else             ! else must be snow
                js=ni+1
             endif
           do j=js,ns
            matj(ns+1,j) = matj(ns+1,j) - rho(j-1) * em(j-1) * cs(1,j) * wsg(j-1)
            matj(ns+1,j) = matj(ns+1,j) - rho(j  ) * em(j  ) * cs(1,j) * wog(j-1)
           enddo
          endif

! liquid fraction terms
          do j=2,ni
            jm=j+ns+1
            matj(jm-1,j) = matj(jm-1,j) - rho(j-1) * w(j-1) * dedp(j-1) * wsg(j-1)
          enddo
          do j=1,ni
            jm=j+ns+1
            matj(jm  ,j) = matj(jm  ,j) - rho(j  ) * w(j-1) * dedp(j  ) * wog(j-1)
          enddo
   
!-----------------------------------------------------------------------
! upper transport in ice
!-----------------------------------------------------------------------

          do j=1,ns
            matj(j  ,j) = matj(j  ,j) + rho(j  ) * w(j  ) * cp(j  ) * wsg(j)
            DT0(j)      = DT0(j)      - rho(j  ) * w(j  ) * em(j  ) * wsg(j)
            DT0(j)      = DT0(j)      - rho(j+1) * w(j  ) * em(j+1) * wog(j)
          enddo
          do j=1,ns-1
            matj(j+1,j) = matj(j+1,j) + rho(j+1) * w(j  ) * cp(j+1) * wog(j)
          enddo
          j=ns ! covers case of snow fall (w<0)
          if ( .not.Tsbc ) &
            matj(j+1,j) = matj(j+1,j) + rho(j+1) * w(j  ) * cp(j+1) * wog(j)

! metric terms, newton terms for change in thickness
          do j=1,ni-1
            matj(0   ,j) = matj(0   ,j) + rho(j  ) * em(j  ) * os(2,j) * wsg(j)
            matj(0   ,j) = matj(0   ,j) + rho(j+1) * em(j+1) * os(2,j) * wog(j)
          enddo
          if ( Tsbc ) then ! top velocity: ice or snow?
             if (ns==ni) then ! definitily ice
                js=1
             else             ! else must be snow
                js=ni+1
             endif
           do j=js,ns
            matj(ns+1,j) = matj(ns+1,j) + rho(j  ) * em(j  ) * cs(2,j) * wsg(j)
            matj(ns+1,j) = matj(ns+1,j) + rho(j+1) * em(j+1) * cs(2,j) * wog(j)
           enddo
          endif

! liquid fraction terms
          do j=1,ni
            jm=j+ns+1
            matj(jm  ,j) = matj(jm  ,j) + rho(j  ) * w(j  ) * dedp(j  ) * wsg(j)
          enddo
          do j=1,ni-1
            jm=j+ns+1
            matj(jm+1,j) = matj(jm+1,j) + rho(j+1) * w(j  ) * dedp(j+1) * wog(j)
          enddo


!-----------------------------------------------------------------------
! add terms for salinity
!-----------------------------------------------------------------------

        do j=1,ni
         DT0(j ) =  DT0(j ) &
                       + wb(j-1) * rho(j-1) * eb(j-1) * wsb(j-1)   & ! upward brine heat lower transport
                       + wb(j-1) * rho(j  ) * eb(j  ) * wob(j-1)   & ! lower transport
                       - wb(j  ) * rho(j  ) * eb(j  ) * wsb(j  )   & ! upper transport
                       - wb(j  ) * rho(j+1) * eb(j+1) * wob(j  )   & ! upper transport
                       - ub(j  ) * rho(j  ) * eb(j  )              & ! cell flushing
                       - um(j  ) * rho(j  ) * eb(j  )                ! cell flushing
           matj(j  ,j) = matj(j  ,j)                               & ! brine heat advection
                       - wb(j-1) * rho(j  ) * cp_wat  * wob(j-1)   & ! lower transport
                       + wb(j  ) * rho(j  ) * cp_wat  * wsb(j  )   & ! upper transport
                       + ub(j  ) * rho(j  ) * cp_wat               & ! right hand where drainage occurs
                       + um(j  ) * rho(j  ) * cp_wat                 ! right hand where meltwater flushing occurs
           if (j>1) &
           matj(j-1,j) = matj(j-1,j)                  &
                       - wb(j-1) * rho(j-1) * cp_wat  * wsb(j-1)     ! lower transport
           if (j<ni) &
           matj(j+1,j) = matj(j+1,j)                  &
                       + wb(j  ) * rho(j+1) * cp_wat  * wob(j  )     ! upper transport
         jm=ns+1+j
         DT0(jm) =  - (henew(j)*sali(j,1)-heold(j)*sali(j,0)) / dtice & ! time tendency energy
                       + w (j-1) * sali(j-1,1) * wsg(j-1)   & ! lower transport
                       + w (j-1) * sali(j  ,1) * wog(j-1)   & ! lower transport
                       - w (j  ) * sali(j  ,1) * wsg(j  )   & ! upper transport
                       - w (j  ) * sali(j+1,1) * wog(j  )   & ! upper transport
                       - um(j  ) * sibr(j  ,1)                ! melting flushing
         DT0(jm) =  DT0(jm)                                 & ! time tendency energy
                       + wb(j-1) * sibr(j-1,1) * wsb(j-1)   & ! upward brine lower transport
                       + wb(j-1) * sibr(j  ,1) * wob(j-1)   & ! lower transport
                       - wb(j  ) * sibr(j  ,1) * wsb(j  )   & ! upper transport
                       - wb(j  ) * sibr(j+1,1) * wob(j  )   & ! upper transport
                       - ub(j  ) * sibr(j  ,1)                ! flushing of brine into brine channels due to upward brine transport
           matj(jm  ,jm) = henew(j) * sibr(j,1) / dtice            ! diagonal term
           matj(jm  ,jm) = matj(jm  ,jm)                  &
                       - w (j-1) * sibr(j  ,1) * wog(j-1) &  ! lower transport
                       + w (j  ) * sibr(j  ,1) * wsg(j  )    ! upper transport
           if (j>1) &
           matj(jm-1,jm) = matj(jm-1,jm)                  &
                       - w (j-1) * sibr(j-1,1) * wsg(j-1)    ! lower transport
           if (j<ni) &
           matj(jm+1,jm) = matj(jm+1,jm)                  &
                       + w (j  ) * sibr(j+1,1) * wog(j  )    ! upper transport
! temperature effect on salinity
           matj(j   ,jm) =                                               &
                        henew(j) * dsidt(j    ) * phib(j  ,1) / dtice      ! diagonal term
           if (j>1) &
           matj(j-1 ,jm) = matj(j-1 ,jm)                                 &
                       - wb(j-1) * dsidt(j-1  )               * wsb(j-1) & ! upward brine lower transport
                       - w (j-1) * dsidt(j-1  ) * phib(j-1,1) * wsg(j-1)   ! lower transport
           matj(j   ,jm) = matj(j   ,jm)                                 &
                       - wb(j-1) * dsidt(j    )               * wob(j-1) & ! lower transport
                       + wb(j  ) * dsidt(j    )               * wsb(j  ) & ! upper transport
                       - w (j-1) * dsidt(j    ) * phib(j  ,1) * wog(j-1) & ! lower transport
                       + w (j  ) * dsidt(j    ) * phib(j  ,1) * wsg(j  ) & ! upper transport
                       + ub(j  ) * dsidt(j    )                          & ! flushing of brine into brine channels due to upward brine
                       + um(j  ) * dsidt(j    )                            ! melting flushing
           if (j<ni) &
           matj(j+1 ,jm) = matj(j+1 ,jm)                                 &
                       + wb(j  ) * dsidt(j+1  )               * wob(j  ) & ! upper transport
                       + w (j  ) * dsidt(j+1  ) * phib(j+1,1) * wog(j  )   ! upper transport
        enddo

! reset terms depending on velocity (j=0 and/or j=ns+1)
! Newton terms for change in thickness in rho h E: 
!    d(h S)/d(w0)= S dh/dw0 = S ds; d(h)/d(ws)=-ds
      do j=1,ni
         jm=ns+1+j
         matj(0   ,jm) =   sali(j,1) * dzi(j)
      enddo
      if ( Tsbc .and. ns==ni) then ! definitily ice
       do j=1,ni
         jm=ns+1+j
         matj(ns+1,jm) = - sali(j,1) * dzi(j)
       enddo
      endif

! metric terms, newton terms for change in thickness for lower transport

! bottom velocity
           do j=1,ni
         jm=ns+1+j
            matj(0   ,jm) = matj(0   ,jm) - sali(j-1,1) * os(1,j) * wsg(j-1)
            matj(0   ,jm) = matj(0   ,jm) - sali(j  ,1) * os(1,j) * wog(j-1)
           enddo

! top velocity
          if ( Tsbc .and. ns==ni) then ! definitily ice
           do j=2,ni
         jm=ns+1+j
            matj(ns+1,jm) = matj(ns+1,jm) - sali(j-1,1) * cs(1,j) * wsg(j-1)
            matj(ns+1,jm) = matj(ns+1,jm) - sali(j  ,1) * cs(1,j) * wog(j-1)
           enddo
          endif

! metric terms, newton terms for change in thickness for upper transport
          do j=1,ni-1
         jm=ns+1+j
            matj(0   ,jm) = matj(0   ,jm) + sali(j  ,1) * os(2,j) * wsg(j)
            matj(0   ,jm) = matj(0   ,jm) + sali(j+1,1) * os(2,j) * wog(j)
          enddo
          if ( Tsbc  .and. ns==ni) then ! definitily ice
           do j=1,ni
         jm=ns+1+j
            matj(ns+1,jm) = matj(ns+1,jm) + sali(j  ,1) * cs(2,j) * wsg(j)
            matj(ns+1,jm) = matj(ns+1,jm) + sali(j+1,1) * cs(2,j) * wog(j)
           enddo
          endif

!-----------------------------------------------------------------------
! prepare linear solver
!-----------------------------------------------------------------------

      residual=0.0_8
      do j=1,ns-1
!      do j=1,ntot-1
         residual=residual+abs(dt0(j))
! FD debug
if (extra_debug) write(*,*) 'res',j,dt0(j)
      enddo
      if (debug) write(*,*) 'residual',counter,residual
! FD debug
if (extra_debug) then
 write(*,*) 'mat'
 do j=0,ntot-1
  write(*,'(40(a1))') (real_to_char(matj(k,j)),k=0,ntot-1)!,real_to_char(dt0(j))
 enddo
 do j=0,ntot-1
  write(*,'(40(e8.1,1x))') matj(0:ntot-1,j),dt0(j)
 enddo
 write(*,*) 'before solving'
 write(*,'(I6,300(1x,e9.3))') counter, w(0:ns)
! write(*,'(I6,300(1x,e9.3))') counter, um(1:ns)
 write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns+1,1)
 write(*,'(I6,300(f8.3,1x))') counter, sibr(0:ni,1)
 write(*,'(I6,300(f8.3,1x))') counter, phib(0:ni,1)
 write(*,'(I6,300(f8.3,1x))') counter, sali(0:ni,1)
! write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns+1,0)
!stop
endif
        ml = ntot - 1
        mu = ml
        MT = ML + MU + 1
!        write(*,*) 'bandwidth',MT+ML, 3*(ns+ni+2)
      do i=1,ntot
       do k=1,bandmax
         matband(k,i)=0.0_8
       enddo
      enddo
      do i=1,ntot
       do j=1,ntot
         K = I - J + MT
         matband(k,j)=matj(j-1,i-1)
       enddo
      enddo
!       call DGBCO(matband,bandmax,nn,ML,MU,IPVT,RCOND,work)
!       write(*,*) 'rcond for the current matrix is',rcond
      call DGBFA(matband,bandmax,ntot,ML,MU,IPVT,INFO)
!      write(*,*) 'matrix factorization done'
      tout(0:ntot-1)=dt0(0:ntot-1)
      call DGBSL(matband,bandmax,ntot,ML,MU,IPVT,tout,0)
! FD debug
if (extra_debug) then
 write(*,*) 'solution'
! do j=0,ntot-1
!  write(*,'(e10.3)') tout(j)
! enddo
!if (counter==3) stop
endif
!-----------------------------------------------------------------------
!     Solve for the internal temperature
!-----------------------------------------------------------------------

!     if ( .not.thin_snow_active ) then
!      CALL tridag (AT1,BT1,CT1,DT0,tout,ns+2)
!     else
!      CALL tridag (AT1,BT1,CT1,DT0,tout,ni+2)
!     endif

! always melting or accreting ice to the bottom
     temp(0,1) = temp(0,0)
     w(0)=w(0)+tout(0)

! ice column increment to temperature
      DO j=1,ns
         if (internal_melt(j)) then
            um(j) = um(j) + tout(j)
         else
            temp(j,1) = temp(j,0) + tout(j)
         endif
      ENDDO

      if (Tsbc) then
! melting case
        w(ns) = w(ns) + tout(ns+1)
        temp(ns+1,1) = Tf(ns+1)
      else
! cold surface case, update temperature
        temp(ns+1,1) = temp(ns+1,0) + tout(ns+1)
      endif

         DO j=1,ns
            IF (temp(j,1) .GT. Tf(j)+tiny) THEN
! FD debug
if (debug) then
 write(*,*) 'detecting internal melt at',j,temp(j,1),Tf(j)
endif
!               temp(j,1) = 0.5_8 * (Tf(j) + tiold(j))
               temp(j,1) = Tf(j)
!! in case top melting should start first
!               if (j<ns) temp(j,1) = tiold(j)
            ENDIF
         ENDDO

! check on top temperature
      j=ns+1
      IF (temp(j,1) .GT. Tf(j)+tiny) temp(j,1) = Tf(j)

! any temperature above top (being snow) is overwritten
      temp(ns+2:ni+nlsno+1,1) = temp(ns+1,1)

! liquid fraction
      DO j=1,ni
         jm = j+ns+1
         phib(j,1) = phib(j,1) + tout(jm)
         sibr(j,1) = func_liqu_sa(temp(j,1))
         sali(j,1) = phib(j,1) * sibr(j,1)
         eb(j)     = cp_wat * temp(j,1)
      ENDDO

! From now on, ns retrieves its original sense
      ns=ni+nlsno

!------------------------------------------------------------------------
!     Update conductive heat fluxes
!------------------------------------------------------------------------

      IF ( .not.thin_snow_active ) THEN
         Fcss = -kks(ns) * ( temp(ns+1,1) - temp(ns  ,1) )
         Fcsb = -kks(ni) * ( temp(ni+1,1) - temp(ni  ,1) )
      ELSE
         Fcss = 0.0_8
         Fcsb = 0.0_8
      ENDIF
      
      Fcis = -kki(ni) * ( temp(ns+1,1) - temp(ni  ,1))  ! only valid for no-snow case
      Fcib = -kki(0 ) * ( temp(   1,1) - temp(   0,1))

      Fnet = Fnet0 - dzf * (temp(ns+1,1)-tiold(ns+1))

! FD debug
if (extra_debug) then
 write(*,'(I6,300(1x,e9.3))') counter, w(0:ns)
! write(*,'(I6,300(1x,e9.3))') counter, um(1:ns)
 write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns+1,1)
 write(*,'(I6,300(f8.3,1x))') counter, sibr(0:ni,1)
 write(*,'(I6,300(f8.3,1x))') counter, phib(0:ni,1)
 write(*,'(I6,300(f8.3,1x))') counter, sali(0:ni,1)
! write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns+1,0)
! stop
endif

     if ( .not.thin_snow_active ) then
      IF (( temp(ns+1,1)+tiny .GE. Tf(ns+1)) .AND. (Tsbc .EQV. .FALSE.)) THEN
         Tsbc = .TRUE. 
      ENDIF
      IF ( temp(ns+1,1)+tiny .GE. Tf(ns+1) .AND. -Fnet-Fcss .GT. 0.0_8 ) THEN
         Tsbc = .FALSE.
         temp(ns+1,1)=Tf(ns+1)-tiny
     ENDIF
     else
      IF ((temp(ni+1,1)+tiny .GE. Tf(ni+1)) .AND. (Tsbc .EQV. .FALSE.)) THEN
         Tsbc = .TRUE. 
      ENDIF
!      IF ( temp(ni+1,1)+tiny .GE. Tf(ni+1) .AND. -Fnet-Fcis .GT. 0.0_8 .and. w(ni)<=0.0_8) THEN
      IF ( temp(ni+1,1)+tiny .GE. Tf(ni+1) .AND. -Fnet-Fcis .GT. 0.0_8) THEN
         Tsbc = .FALSE.
         temp(ni+1:ns+1,1)=Tf(ni+1)-tiny
      ENDIF
     endif
! FD debug
if (extra_debug) then
 if ( .not.thin_snow_active ) then
  write(*,*) 'SBCond snow',Tsbc,Fnet,Fcss, temp(ns+1,1)
 else
  write(*,*) 'SBCond ice',Tsbc,Fnet,Fcis,kki(ni), temp(ns+1,1), temp(ni,1)
endif
endif

!-----------------------------------------------------------------------
!     Update T-S dependent ice parameters
!     Lf, cp are diagnostics for output
!-----------------------------------------------------------------------
      
!      DO j=0,ni+1
!         ki(j)   = func_ki(sali(j,1),temp(j,1))
!      ENDDO
      DO j=0,ns+1
!         Tf(j)   = Tfreeze1(sali(j,1))
         em(j)   = func_El_mush(temp(j,1),phib(j,1))
         cp(j)   = func_cp_mush(temp(j,1),phib(j,1))
         dedp(j) = func_dedp(temp(j,1),phib(j,1))
         dsidt(j)= func_liqu_sadt(temp(j,1))
      ENDDO
      j=0
      energy_bot = + em(j  ) * wsg(j  ) &
                   + em(j+1) * wog(j  )

!-----------------------------------------------------------------------
!     Euler step
!-----------------------------------------------------------------------
      
      DO j=0,ns+1
         dTlay = ABS(temp(j,1)-temp(j,0))	! T diff in a layer
         dTmax = MAX(dTmax,dTlay)		! max diff in s/ice slab
      ENDDO
      IF (dTmax .GT. Tdiff .and. counter < 100 ) THEN
         counter=counter+1
         DO j=0,ns+1
            temp(j,0) = temp(j,1)		! update temperature
         ENDDO
         dTmax = 0.0_8
         GOTO 4000				! recalculate T profile
      ENDIF

!-----------------------------------------------------------------------
!     Compute internal energy
!-----------------------------------------------------------------------

      Frad=0.0_8
      DO j=1,ns
         Frad	= Frad + R(j)
      ENDDO
      
!-----------------------------------------------------------------------
!     Update energy flux due to snow precipitation
!-----------------------------------------------------------------------
! assumes that the temperature at ns+1 is the surface temperature (independent of presence of snow or not)
! and reset at temperature at ns as in the old snow layer
      
      if ( .not.thin_snow_active ) then
          Fprec = snow_precip * rhosno * em(ns+1)
      else
          Fprec = snow_precip * rhosno * em_thin_snow
          do j=ni+1,ns
             temp(j,1) = tiold(j)
             sali(j,1) = 0.0_8
             Tf(j) = 0.0_8 ! Tfreeze1(sali(j,1))
             rho(j) = rhosno
          enddo
      endif

! fluxes associated with brine
      j=1
      Fsbr =         + wb(j-1) * sibr(j-1,1) * wsb(j-1)   & ! upward brine lower transport
                     + wb(j-1) * sibr(j  ,1) * wob(j-1)   & ! lower transport
                     + w (j-1) * sali(j-1,1) * wsg(j-1)   & ! lower transport
                     + w (j-1) * sali(j  ,1) * wog(j-1)     ! lower transport
      Fhbr =         + wb(j-1) * rho(j-1) * eb(j-1) * wsb(j-1)   & ! upward brine heat lower transport
                     + wb(j-1) * rho(j  ) * eb(j  ) * wob(j-1)     ! lower transport
      DO j=1,ni
         Fsbr = Fsbr - sibr(j,1) * ub(j) &
                     - sibr(j,1) * um(j)
         Fhbr = Fhbr - ub(j  ) * rho(j  ) * eb(j  )
         Fhbr = Fhbr - um(j  ) * rho(j  ) * eb(j  )
      ENDDO
      j=ni ! loss from top melting
         Fsbr = Fsbr &
                     - w (j) * sali(j  ,1) * wsg(j)   & ! upper transport
                     - w (j) * sali(j+1,1) * wog(j)     ! upper transport

      dsaldt = sumsal
      sumsal = 0.0_8
      DO j=1,ni
         sumsal = sumsal + henew(j) * sali(j,1)
      ENDDO
      dsaldt = ( sumsal - dsaldt ) / dtice

! FD debug
if (debug) write(*,*) 'salinity',Fsbr,dsaldt

      Einp = (Fnet+oceflx+Frad+Fprec+Fthin_snow+Fhbr)*dtice	! E input [J/m2]

      
      Eint = 0.0_8

      DO j=1,ns
         Elay = func_El_mush(temp(j,1),phib(j,1))*henew(j)*rho(j)
         Eint = Eint+Elay
      ENDDO
     
      dEin = Eint-Ein0		! diff intl E [J/m2]
! FD debug
if (debug) write(*,*) 'energy',Einp/dtice,dEin/dtice

!     WARNING NI WAS INVENTED FOR TEMPERATURE POINTS WHICH ALSO EXISTS
!     AT THE INTERFACE OF THE ICE AND SNOW, AT THE ICE BASE AND ICE SURFACE
!     WHEN IT COMES TO INTERNAL ENERGY, THIS ONLY EXISTS IN AT THE GRID CENTER
!     AND FOR THIS REASON A NEW PARAMETER SHOULD BE DEFINED LT=LI+LS
!     THIS WOULD MAKE THE LOOP BELOW MUCH CLEANER AND EASIER TO FOLLOW
!     DO LATER IN THE FINAL CLEAN UP STAGE
!     ALSO SHOULDN'T EINT(I) BE INITIALIZE TO ZERO BEFORE THIS LOOP?
!     IT SEEMS TO WORK LIKE THIS BUT WOULD BETTER
!     SALI SHOULD BE DIMENSIONNED TO INCLUDE THE SNOW AS WELL... SNOW CAN 
!     SALINITY EVENTUALLY, BUT MORE IMPORTANTLY IT MAKES THE USE OF THE
!     FUNCTION EL EASIER TO HANDLE

!-----------------------------------------------------------------------
!     Cumulative growth/melt
!-----------------------------------------------------------------------

      msurf = msurf-disdt*dtice
      
      IF (dibdt .LT. 0.0_8) gbase = gbase-dibdt*dtice
      IF (dibdt .GT. 0.0_8) mbase = mbase+dibdt*dtice
      
!-----------------------------------------------------------------------
!     Reset prognostic variables
!-----------------------------------------------------------------------
      
! FD debug
if (extra_debug) then
 write(*,*) 'snow',hs_b,-w(ns)
 write(*,*) 'hice',hi_b,-w(ni),w(0)
 write(*,'(I6,300(1x,e9.3))') counter, w(0:ns)
 write(*,'(I6,300(1x,e9.3))') counter, um(1:ni)
 write(*,'(I6,300(1x,e9.3))') counter, wb(0:ni)
 write(*,'(I6,300(1x,e9.3))') counter, ub(1:ni)
 write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns+1,1)
 write(*,'(I6,300(f8.3,1x))') counter, phib(0:ni,1)
 write(*,'(I6,300(f8.3,1x))') counter, sali(0:ni,1)
 write(*,'(I6,300(f8.3,1x))') counter, sibr(0:ni,1)
endif
! FD debug if energy not conserved
      minfraliq=minval(phib(1:ns,1))
      IF (minfraliq < 0.0_8) THEN
         lstop_minfraliq = .true.
         write(*,*) 'liquidus fraction negative, stop!'
      ENDIF
      IF (abs(dsaldt-Fsbr) > 1e-9) THEN
         lstop_salt = .true.
         write(*,*) 'salt not conserved'
      ENDIF
      IF (abs(Einp-dEin)/dtice > 1e-3) THEN
         lstop_heat = .true.
         write(*,*) 'heat not conserved'
      ENDIF
      IF ( lstop_minfraliq .or. lstop_salt .or. lstop_heat .or. counter > 10 ) THEN
!      IF ( counter > 10 ) THEN
       OPEN(1,file='restart.dat')
       WRITE(1,*) nlice,nlsno,ith_cond,dtice
       WRITE(1,*) ti(1:nlice),ts(1:nlsno),tsuold,tbo
       WRITE(1,*) si(1:nlice)
       WRITE(1,*) hi,hsold
       WRITE(1,*) sfallold,dwnlw,tsuold,tair,qair,uair,swrad,oceflx,pres
       WRITE(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
       WRITE(1,*) fsens,flat
       CLOSE(1)
       STOP 'something bad happened, dumping restart and exiting.'
      ENDIF

      IF (hi_b(1) .LT. hicut) THEN
         hi_b(1) = hicut
         DO j=0,ni
            temp(j,1) = tocn
         ENDDO
      ENDIF
      
      hs_b(0) = hs_b(1)
      hi_b(0) = hi_b(1)

! from here on, ns retrieves its meaning of top snow
      ns=ni+nlsno
      
      DO j=0,ns
         sali(j,0) = sali(j,1)
         temp(j,0) = temp(j,1)
      ENDDO
      

! FD prepare for callback

!    write(*,'(5(f10.3,1x),I6)') tocn,temp(0,1),temp(ns,1),tair,hice,counter

! init variable
      hs = hs_b(1)
      hi = hi_b(1)

!
!------------------------------------------------------------------------------|
!  5) Formation of snow-ice                                                    |
!------------------------------------------------------------------------------|
!
      ! When snow load excesses Archimede's limit, snow-ice interface goes down
      ! under sea-level, flooding of seawater transforms snow into ice
      ! dh_snowice is positive for the ice

      dh_sni = MAX( 0.0_8 , ( rhosno * hs + (rhoice - rhowat ) * hi) / ( rhosno + rhowat - rhoice ) )

      if (dh_sni > 0.0_8 ) then
       hi  = hi + dh_sni
       hs  = hs - dh_sni
       if (debug) then
        WRITE(*,*) ' dh_snowice : ', dh_sni
        WRITE(*,*) ' ht_s_b : ', hs
        WRITE(*,*) ' ht_i_b : ', hi
       endif
      endif

! conversion from LIM3
      tbo=temp(0,1) + temp0
      do j=1,nlice
        ti(nlice-j+1)=temp(j,1) + temp0
      enddo
      do j=1,nlsno
        ts(nlsno-j+1)=temp(j+ni,1) + temp0
      enddo
      tsu=temp(ns+1,1) + temp0

! call back for LIM3 of thickness (needed for radiative transfer)
      do j=1,nlice
        dzi(nlice-j+1) = henew(j)
      enddo
      do j=1,nlsno
        dzs(nlsno-j+1) = henew(j+ni)
      enddo

! call back for bulk salinity and brine variables (only needed for output)
      do j=1,nlice
        si(nlice-j+1)=sali(j,1)
        brine_v(nlice-j+1) = phib(j,1) ! brine volume
        brine_u(nlice-j+1) =   ub(j)   ! brine flushing velocity
        brine_r(nlice-j+1) =   ra(j)   ! brine Rayleigh number
      enddo
      dhi_bot = dibdt
      dhi_surf = disdt - sumhmelt
      dhs = dhsdt

! FD debug
!stop

      end subroutine ice_thermo

!************************************************************************
!     Subroutine tridag: Solve tridiagonal system of equations (Ax=r).
!       a : lower  diagonal element of A
!       b : middle diagonal element of A
!       c : upper  diagonal element of A
!       r : right hand side of the equation (inhomogenities)
!       u : solution x to the system of equation
!       n : dimension of the matrix A (n x n)
!************************************************************************

      SUBROUTINE tridag(a,b,c,r,u,n)
      implicit none
      INTEGER n,NMAX
      real(8) a(n),b(n),c(n),r(n),u(n)
      PARAMETER(NMAX=230)
      INTEGER j
      real(8) bet,gam(NMAX)
      if(b(1).eq.0.0_8) stop 'tridag: check BC'
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j)*gam(j)
         if(bet.eq.0.0_8) stop 'tridag: failed'
         u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
      END SUBROUTINE tridag




end module ice_thermodynamic_FV

!-----------------------------------------------------------------------
!     Fraction of sw radiation penetrating the surface
!-----------------------------------------------------------------------

      FUNCTION func_i0(hs,hscut)

      implicit none
      real(8)	func_i0, hs, hscut
      real(8) :: i0snow    =   0.080_8 ! SW fraction penetrating snow surface
      real(8) :: i0ice     =   0.170_8	! SW fraction penetrating ice  surface

      IF (hs .GT. hscut) THEN
         func_i0 = i0snow
      ELSE
         func_i0 = i0ice
      ENDIF

      END FUNCTION func_i0

!-----------------------------------------------------------------------
!     Saturation vapor pressure
!-----------------------------------------------------------------------

      FUNCTION esat(T)

      implicit none
      real(8)	esat, T, coef1, coef2 	! T=[C]
      real(8) :: temp0	= 273.160_8	! freezing point temp of fresh water	[K]

      coef1	= 21.87456_8			! coeff. over ice
      coef2	=  7.66000_8			! coeff. over ice
      esat	=  6.11e2_8*EXP(MIN(coef1*T/(T+temp0-coef2),10.0_8))

      END FUNCTION esat


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

      function real_to_char(r)
      implicit none
      ! arguments
      real(8) r
      character(len=1) real_to_char

      real_to_char = '*'
      if (abs(r) < 1e-11) real_to_char = '.'

      end function real_to_char
