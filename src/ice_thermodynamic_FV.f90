module ice_thermodynamic_FV

  use var_thermo_vertical

implicit none

!------------------------------------------------------------------------
!     Some model constants and parameters
!------------------------------------------------------------------------
double precision :: &
      salinis	=   1.000d+0, &    ! salinity at ice surf			[psu]
      salinib	=   4.000d+0, &    ! salinity at ice base			[psu]
! FD test      salnice	=  10.000d+0, &    ! salinity of newly formed ice		[psu]
      salnice	=   4.000d+0       ! salinity of newly formed ice		[psu]
double precision, dimension(0:maxlay) :: dzzi,dzzs


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
  double precision,save :: tzero  = 273.15d0       ! tzero C => K

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

! pi=atan(1.d0)*4.d0
! do k=0,maxlay-1
!    zi(k)=0.5d0*(1.d0+sin(pi*(dble(k)/dble(maxlay-1)-0.5d0)))
! enddo

!---------------------------------------------------------------------
! uniform distribution for ice
! zi varies from 0 (bottom) to 1 (surface)
!---------------------------------------------------------------------

 do k=0,ni-1
    zi(k) = dble(k)/dble(nlice)
 enddo
   zi(ni) = 1.d0
 do k=ni+1,ns-1
    zi(k) = 1.d0 + dble(k-ni)/dble(nlsno)
 enddo
   zi(ns) = 2.d0
   zs(0:ni)=0.d0
 do k=ni+1,ns-1
    zs(k) = dble(k-ni)/dble(nlsno)
 enddo
   zs(ns) = 1.d0

 dzzi(0:ns)=0.d0
 do k=1,ni
    dzi(k) = 1.d0/dble(nlice)
    dzzi(k) = 1.d0/dble(nlice)
 enddo
 dzzi(0 ) = 0.5d0/dble(nlice)
 dzzi(ni) = 0.5d0/dble(nlice)
 dzzs(0:ns)=0.d0
 dzzs(ni) = 0.5d0/dble(nlsno)
 do k=ni+1,ns
    dzi(k) = 1.d0/dble(nlsno)
    dzzs(k) = 1.d0/dble(nlsno)
 enddo
 dzzs(ns) = 0.5d0/dble(nlsno)

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
  double precision       dtice

!locals
!===============================================================================
!     Parameter statement defines layers and integration periode
!===============================================================================
!....&S..1.........2.........3.........4.........5.........6.........7..C......8

      DOUBLE PRECISION :: &
        Fbase,&! heat flux at ice base			[W/m2]
        S    ,&! initial salinity profile		[psu]
        T    ,&! initial temp. profile			[C]
        zini(0:maxlay),&! z-values of init interpol tprof 	[m]
        Tsurf           ! surface boundary conditions           [C]

      CHARACTER (len=10) :: iatflx, ocnflx, ocnsal, ocntmp, uiorel, &
                            salinity, bbc, sbc

      double precision :: uiofix, & ! constant relative ice-ocean velocity		[m/s]
        sigma, &! Stefan Bolzmann constant			[W/m/K4]
        temp0, &! melt temperature of snow and ice		[K]
        tiny,  &! very small number (1d-9)			[-]
        Tdiff    ! temperature difference in Euler step		[C]

      DOUBLE PRECISION :: &
     epsilon,&! coefficient of emissivity			[-]
     hicut,&! cut off (minimum) ice  thickness		[m]
     hscut,&! cut off (minimum) snow thickness		[m]
     hslim  ! minimum snow thickness			[m]

      logical   Tsbc           ! Fixed T surface boudnary conditions

      DOUBLE PRECISION :: &
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
    Fr(0:maxlay)    ,&! sw rad heat flux at s/i internal level 	[W/m2]
    ki(0:maxlay)    ,&! ice  thermal conductivity			[W/m/C]
    ks(0:maxlay)    ,&! snow thermal conductivity			[W/m/C]
    kki(0:maxlay)    ,&! ice  thermal conductivity			[W/m/C]
    kks(0:maxlay)    ,&! ice  thermal conductivity			[W/m/C]
    em(0:maxlay)    ,&! volumetric enthalpiy	[W/m3]
    em0(0:maxlay)    ,&! old volumetric enthalpiy	[W/m3]
    qm(0:maxlay)    ,&! volumetric energy of melt (rho*Lf(S,T))	[W/m3]
    cp(0:maxlay)    ,&! heat capacity for ice	[W/kg/K^-1]
    R(0:maxlay)    ,&! penetrating shortwave radiation		[W/m3]
    sali(0:maxlay,0:1),&! internal sea ice salinity			[psu]
    Tf(0:maxlay)    ,&! ice freezing temp. of salinity S		[C]
    w(0:maxlay)  ,&! grid advection velocity			[1/s]
    hi_b(0:1)     ,&! ice thickness (m)
    hs_b(0:1)       ! snow thickness (m)
  double precision, dimension(maxlay) :: &
             he, heold, henew

      double precision &
                temp (0:maxlay,0:1) ! internal sea ice temperature		[C]

      double precision &
                matj(0:maxlay,0:maxlay) ! jacobian matrix		[C]

! LU Solver START
      double precision &
                matband(2*maxlay,maxlay) ! band matrix
      integer ntot,ml,mu,mt,info,bandmax
      integer ipvt(maxlay)
      double precision work(maxlay),rcond,residual
! LU Solver END

      DOUBLE PRECISION AT1(0:maxlay), BT1(0:maxlay), CT1(0:maxlay), DT0(0:maxlay), &
      tout(0:maxlay), hsmemo, Rtrans

      double precision :: sice, &
             latmelt, ftot, sstznew, &
             Tside,deltaT,wlat,rside, &
             sum0, dhi, dho, dq, fwf0
      double precision :: hminice=1d-3
      double precision :: energy_bot

  character (len=20) :: Sprofile
  integer i,j,counter,k
  double precision dzf, es, zssdqw, zrchu1, zrchu2, q0, zref, k0, k1
  double precision, dimension(0:maxlay) :: tiold, tinew
  double precision elays0
  logical :: debug=.false.
  logical :: extra_debug=.false.
  double precision :: hsold, tsuold, sfallold
  logical :: thin_snow_active
  double precision :: snow_precip, rain_precip, fwf, &
                      fthin_snow, energy_snow_melt, sume2, em_thin_snow
! parameters for ice thermodynamics and heat fluxes
  double precision,save :: tzero  = 273.15d0       ! tzero C => K

  double precision,save :: latvap = 2.5d6          ! J/kg latent heat of vaporization
  double precision,save :: latsub = 2.834d6        ! J/kg latent heat of sublimation
  logical :: adv_upwind=.false.
! FD debug
 debug=.true.
 extra_debug=.true.
      hicut=1.d-2

!------------------------------------------------------------------------
!     condition on ice thickness
!     no thermo below 5 cm
!------------------------------------------------------------------------

!      if (hi.lt.hicut) return

!------------------------------------------------------------------------
!     Some logical constants
!------------------------------------------------------------------------

      salinity  = 'no'          ! = no dynamic salinity
      bbc       = 'fixT'        ! bottom boundqry condition

!------------------------------------------------------------------------
!     Some physical constants
!------------------------------------------------------------------------

      epsilon=   0.990d+0! snow/ice emissivity			[-]

      sigma=   5.670d-8! Stefan-Boltzmann constant	       [W/m2/K4]
      temp0= 273.160d+0! freezing point temp of fresh water	[K]

      uiofix=   0.000d+0! constant relative ice-ocean velocity	[m/s]


! FD debug
!hs=0.d0
!ts(:)=tsu
!oceflx=0.d0
!swradab_s=0.d0
!swradab_i=0.d0
!snowfall=0.d0
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

      tiny	=   1.000d-9	! very small number			[-]

      hslim     =   0.0005d0    ! limit for computing temperature in snow
      Tdiff     =   1.000d-12   ! temperature tolerance in Euler step	[C]
      thin_snow_active =.false. ! reroute some energy to sublimate/melt snow if true
      fwf = 0.d0
! FD debug
!      hslim     =   0.005d0    ! limit for computing temperature in snow
!      hslim     =   0.01d0    ! limit for computing temperature in snow

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

      dTlay	=   0.000d+0	! temp diff in a layer	[C]
      dTmax	=   0.000d+0	! max temp difference	[C]
      
!-----------------------------------------------------------------------
!     Initial salinity profile
!-----------------------------------------------------------------------

      sali(0,0)=salnice
      do j=1,nlice
        sali(j,0)=si(nlice-j+1)
      enddo
      DO j=ni+1,ns+1
         sali(j,0) = 0d0
      ENDDO

      hsold = hs
      tsuold = tsu
      sfallold = snowfall
      fthin_snow = 0.d0
      dspdt = 0.d0

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
!     Vertical (grid) velocities
!------------------------------------------------------------------------

      DO j=0,ns
         w(j) = 0d0		! initial grid advection
      ENDDO

!------------------------------------------------------------------------
!     Set initial T-S dependent material properties
!------------------------------------------------------------------------

      DO j=0,ns+1
            Tf(j)   = Tfreeze1(sali(j,0))
            ki(j)   = func_ki(sali(j,0),tiold(j))
            qm(j)   = func_qm(Tf(j),tiold(j))
            em(j)   = func_el(Tf(j),tiold(j))
            cp(j)   = func_cph(Tf(j),tiold(j))
      ENDDO
      em0(0:ns+1) = em(0:ns+1)

      DO j=0,ns+1					! snow thermal conductivity
         IF (j .LT. ni) THEN
            ks(j) = 0d0
         ELSE
            ks(j) = cond_sno
         ENDIF
      ENDDO
      kki(0)=ki(0)/ (hi_b(0) * dzzi(0))
      DO j=1,ni-1
         kki(j) = (ki(j+1)+ki(j)) * 0.5d0 / (hi_b(0) * dzzi(j))
      ENDDO
      kki(ni) = ki(ni) / (hi_b(0) * dzzi(ni))
      kks(0:ns)=0.d0
      kks(ni) = ks(ni)*ki(ni)/max(tiny, ki(ni)*dzzs(ni)*hs_b(0) + ks(ni)*dzzi(ni)*hi_b(0) )
      DO j=ni+1,ns-1
         kks(j) = (ks(j+1)+ks(j)) * 0.5d0 / (hs_b(0) * dzzs(j))
      ENDDO
      kks(ns) = ks(ns) / (hs_b(0) * dzzs(ns))
      kki(ni+1:ns)=kks(ni+1:ns)

!------------------------------------------------------------------------
!     Initial internal energy in ice (1:ni-1) and snow (ni+1:ns-1)
!------------------------------------------------------------------------

      Einp = 0d0				! E input [J/m2]
      Ein0 = 0d0

      DO j=1,ni
         Elay = em(j)*rhoice*heold(j)
         Ein0 = Ein0+Elay
      ENDDO

      DO j=ni+1,ns
         Elay = em(j)*rhosno*heold(j)
         Ein0 = Ein0+Elay
      ENDDO

      dEin = 0d0				! diff intl E [J/m2]

! FD find correct surface temperature
      IF (hs_b(0) > hslim) THEN
        kki(ni) = kks(ni)
      ENDIF

      ! net longwave radiative flux
      netlw = emi*(dwnlw - stefa*tsu*tsu*tsu*tsu)
      ! sensible and latent heat flux
!      fsens   =  -fsens
! FD debug
if (extra_debug) then
 write(*,*) 'fsens,flat',fsens,MIN( flat , 0.d0 ),dwnlw,netlw
endif
      flat   =  MIN( flat , 0.d0 ) ! always negative, as precip 
                                           ! energy already added

      ! pressure of water vapor saturation (Pa)
      es         =  611.d0*10.d0**(9.5d0*(tsu-273.16d0)/(tsu-7.66d0))
      ! intermediate variable
      zssdqw     =  q0*q0*pres/ &
     &              (0.622d0*es)*log(10.d0)*9.5d0* &
     &              ((273.16d0-7.66d0)/(tsu-7.66d0)**2.d0)
      ! derivative of the surface atmospheric net flux
          dzf    =  4.d0*emi*stefa*tsu*tsu*tsu

      ! surface atmospheric net flux
      Fnet0 =  fac_transmi * swrad + netlw + fsens + flat
      Fnet  = Fnet0

!---------------------------------------------------------------------
! accumulation at the surface
!---------------------------------------------------------------------
     if ( tair < tzero ) then
      snow_precip = snowfall
      rain_precip = 0.d0
     else
      snow_precip = 0.d0
      rain_precip = snowfall * rhosno / rhowat
     endif

!---------------------------------------------------------------------
! need to get rid of snow if too thin during melting/sublimating
! thin snow case
!---------------------------------------------------------------------
  if ( hs_b(0) <= hslim ) then

    fthin_snow = 0.d0
    DO j=ni+1,ns
       Elay = em0(j)*rhosno*heold(j)
       fthin_snow = fthin_snow + Elay
    ENDDO
    fthin_snow = - fthin_snow / dtice
    do j=ni+1,ns
      heold(j)=0.d0
    enddo
    hs_b(0)=0.d0
    thin_snow_active=.true.
    tiold(ni+1:ns) = tiold(ns+1)
    em_thin_snow = func_el(Tf(ns+1),tiold(ns+1))
    Tf(ni+1:ns+1) = Tf(ni)
    sali(ni+1:ns+1,0) = sali(ni,0)

  endif ! end thin snow
!---------------------------------------------------------------------

      Fnet0 = Fnet0 - fthin_snow
      Fnet = Fnet - fthin_snow

!------------------------------------------------------------------------
!     Update conductive heat fluxes
!------------------------------------------------------------------------

      IF ( .not.thin_snow_active ) THEN
         Fcss = -kks(ns) * ( temp(ns+1,1) - temp(ns  ,1) )
         Fcsb = -kks(ni) * ( temp(ni+1,1) - temp(ni  ,1) )
      ELSE
         Fcss = 0d0
         Fcsb = 0d0
      ENDIF
      
      Fcis = -kki(ni) * ( temp(ns+1,1) - temp(ni  ,1))
      Fcib = -kki(0 ) * ( temp(   1,1) - temp(   0,0))

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

! FD         oceflx= rhoo*cpo*Coi*ABS(uio)*(tocn-temp(0,1))
!-----------------------------------------------------------------------
!     Update snow and ice thickness
!-----------------------------------------------------------------------

      dssdt = snow_precip + dspdt
      dsbdt = 0d0		! melt of snow base not allowed

      disdt = 0d0

         ! snow depth evolution due to precipitation


! FD: if snow present
      if ( .not.thin_snow_active ) then
         ! snow depth evolution due to melt superposed on precip

         IF (temp(ns+1,1)+tiny .GE. Tf(ns+1) .AND. -Fnet-Fcss .LT. 0d0) THEN
! FD the code cannot deal with concomittent snow accumulation and melt
            Tsbc = .true.
            dssdt = MIN((-Fnet-Fcss)/rhosno/qm(ns+1),0d0)
            snow_precip = 0.d0
         ENDIF
      else

      ! update ice surface

         IF (temp(ns+1,1)+tiny .GE. Tf(ns+1) .AND. -Fnet-Fcis .LT. 0d0) THEN
            Tsbc = .true.
            disdt = MIN((-Fnet-Fcis)/rhoice/qm(ni),0d0)
         ENDIF
      endif

! bottom growth or melt
      energy_bot = qm(0)
      dibdt	= (oceflx-Fcib) / rhoice / energy_bot

      dhsdt	= dssdt - dsbdt
      dhidt	= disdt - dibdt ! original

      hs_b(1)	= hs_b(0) + dhsdt*dtice

!---------------------------------------------------------------------
! need to get rid of snow if wmesh too large
!---------------------------------------------------------------------

   if ( hs_b(1) < 0.d0 .and. .not.thin_snow_active ) then
      sume2 = 0.d0
      do j=ni+1,ns
        sume2 = sume2 + em0(j) * rhosno * heold(j)
      enddo
      energy_snow_melt = -sume2
      fthin_snow = energy_snow_melt / dtice
! FD debug
if (extra_debug) write(*,*) 'switch to thin_snow_active',hs_b(1),fthin_snow
      fwf = fwf + hs_b(0) / dtice / rhowat * rhosno
      Tsbc = .false.
      thin_snow_active = .true.
      snow_precip = 0.d0
      em_thin_snow = 0.d0
      hs_b(1) = 0.d0
      dspdt = -hs_b(0) / dtice
      dssdt = dspdt
      Fnet0 = Fnet0 - fthin_snow
      Fnet = Fnet - fthin_snow
      kki(ni) = ki(ni) / (hi_b(0) * dzzi(ni))
      Tf(ni+1:ns+1) = Tf(ni)
      sali(ni+1:ns+1,0) = sali(ni,0)
      temp(ni+1:ns+1,1) = Tf(ni)
      qm(ni+1:ns+1)   = func_qm(Tf(ni+1),temp(ni+1,1))
      em(ni+1:ns+1)   = func_el(Tf(ni+1),temp(ni+1,1))
      ! update ice surface

      Fcis = -kki(ni) * ( temp(ns+1,1) - temp(ni  ,1))
         IF (temp(ns+1,1)+tiny .GE. Tf(ns+1) .AND. -Fnet-Fcis .LT. 0d0) THEN
            Tsbc = .true.
            disdt = MIN((-Fnet-Fcis)/rhoice/qm(ni),0d0)
         ENDIF
      dhidt	= disdt - dibdt ! original

   endif

      hi_b(1)	= hi_b(0) + dhidt*dtice

      henew(   1:ni  )=dzi(   1:ni  )*hi_b(1)
      henew(ni+1:ns  )=dzi(ni+1:ns  )*hs_b(1)

!-----------------------------------------------------------------------
!     Vertical velocities (grid advection)
!-----------------------------------------------------------------------

      w(0 ) = - dibdt
      w(ni) = - disdt
      DO j=1,ni-1
         w(j)	= zi(j) * w(ni) + (1d0-zi(j)) * w(0 )
      ENDDO
      w(ns) = - dssdt
      DO j=ni+1,ns-1
         w(j)	= zs(j) * w(ni) + (1d0-zs(j)) * w(ns)
      ENDDO

!-----------------------------------------------------------------------
!     Absorbed SW radiation energy per unit volume (R)       
!     1/dz Int_z^z+dz R dz = Rtrans(z+dz) - Rtrans(z)
!-----------------------------------------------------------------------

      if ( .not.thin_snow_active ) then
        DO j=ns,ni+1,-1
          R(j)   = swradab_s(ns-j+1)
        ENDDO
      else
         R(ni+1:ns)=0.d0
      endif
      R(ni)=0.d0
      DO j=ni,1,-1
         R(j)   = swradab_i(ni-j+1)
      ENDDO
      R(0)=0.d0

!-----------------------------------------------------------------------
!     matrix coefficients [matj] {x} = {D}
!-----------------------------------------------------------------------

       if ( thin_snow_active ) then
          ntot=ni+2
       else
          ntot=ns+2
       endif

         matj(0:ntot-1,0:ntot-1) = 0.d0


         k0 = kki(0)
! the unknown is w(0)
         DT0(0)    = - w(0) * rhoice * energy_bot - oceflx + Fcib ! line 0 is for w(0)
         matj(0,0) = rhoice * energy_bot       ! w increment
         matj(1,0) = k0                        ! dFcib/dT increment

      DO j=1,ni
         k0 = kki(j-1)
         k1 = kki(j  )
         DT0(j) =  R(j) + ( & ! radiation
                          - rhoice*(henew(j)*em(j)-heold(j)*em0(j)))/dtice & ! time tendency energy
                        + k0 * ( temp(j-1,1) - temp(j,1) ) &
                        + k1 * ( temp(j+1,1) - temp(j,1) )
         matj(j-1,j) = -k0
         matj(j  ,j) =  k0 + k1 + rhoice * cp(j) * henew(j) / dtice
         matj(j+1,j) = -k1
! Newton terms for change in thickness
         matj(0   ,j) =   rhoice * em(j) * dzi(j)
         if ( Tsbc .and. thin_snow_active ) &
         matj(ni+1,j) = - rhoice * em(j) * dzi(j)
      ENDDO
! FD debug
write(*,*) 'mass',ni,dt0(ni)

     if ( .not.thin_snow_active ) then

      DO j=ni+1,ns
        k0 = kks(j-1)
        k1 = kks(j)
         matj(j-1,j) = -k0
         matj(j  ,j) =  k0 + k1 + rhosno * cp(j) * henew(j) / dtice
         matj(j+1,j) = -k1
         DT0(j)      =  R(j) + ( & ! radiation
                          - rhosno*(henew(j)*em(j)-heold(j)*em0(j)))/dtice & ! time tendency energy
                        + k0 * ( temp(j-1,1) - temp(j,1) ) &
                        + k1 * ( temp(j+1,1) - temp(j,1) )
! Newton terms for change in thickness
         if ( Tsbc) &
         matj(ns+1,j) = - rhosno * em(j) * dzi(j)
      ENDDO

      j=ns+1
       k0 = kks(j-1)
       matj(j-1,j) = - k0
       matj(j  ,j) =   k0 + dzf
       DT0(j)      =   Fnet + Fcss

      IF (Tsbc) THEN
! solve for delta w(ns)
         DT0(j)      = - w(ns) * rhosno * qm(ns+1) + Fnet + Fcss ! line 0 is for w(0)
         matj(j  ,j) = rhosno * qm(ns+1)                         ! w increment
         matj(j-1,j) = - k0                                      ! dFcss/dT increment
      ENDIF

! FD no snow
     else

       j=ni+1
       k0 = kki(j-1)
       matj(j-1,j) = - k0
       matj(j  ,j) =   k0 + dzf
       DT0(j)      =   Fnet + Fcis

      IF (Tsbc) THEN
         DT0(j)      = - w(ni) * rhoice * qm(ni) + Fnet + Fcis ! line 0 is for w(0)
         matj(j  ,j) = rhoice * qm(ni)                         ! w increment
         matj(j-1,j) = - k0                                    ! dFcis/dT increment
!         matj(j-1,j) = - k0 + w(ni) * rhoice * cp(j-1)         ! dFcis/dT increment + variation due to top cell temp variation
! FD debug
write(*,*) 'em top',em(ni),qm(ni),cp(ni),cp(j-1)
      ENDIF
     endif ! condition on hs_b

!-----------------------------------------------------------------------
!     add transport to the equations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! lower transport in ice
!-----------------------------------------------------------------------

! at bottom, treat as centered but temperature is at some position as w(0) so looks like upwind
         j=1
            DT0(j)      = DT0(j)      + rhoice * w(j-1) * em(j-1)

         if (adv_upwind) then
          do j=2,ni
           if (w(j-1)>0.d0) then
            matj(j-1,j) = matj(j-1,j) - rhoice * w(j-1) * cp(j-1)
            DT0(j)      = DT0(j)      + rhoice * w(j-1) * em(j-1)
           else
            matj(j  ,j) = matj(j  ,j) - rhoice * w(j-1) * cp(j  )
            DT0(j)      = DT0(j)      + rhoice * w(j-1) * em(j  )
           endif
          enddo
         else
          do j=2,ni
            matj(j-1,j) = matj(j-1,j) - rhoice * w(j-1) * cp(j-1) * 0.5d0
            DT0(j)      = DT0(j)      + rhoice * w(j-1) * em(j-1) * 0.5d0
            matj(j  ,j) = matj(j  ,j) - rhoice * w(j-1) * cp(j  ) * 0.5d0
            DT0(j)      = DT0(j)      + rhoice * w(j-1) * em(j  ) * 0.5d0
          enddo
! FD debug
write(*,*) 'adv down',ni,dt0(ni)
         endif

! metric terms, newton terms for change in thickness
         j=1
            matj(0   ,j) = matj(0   ,j) - rhoice * em(j-1) * (1d0-zi(j-1))
         if (adv_upwind) then
          do j=2,ni
           if (w(j-1)>0.d0) then
            matj(0   ,j) = matj(0   ,j) - rhoice * em(j-1) * (1d0-zi(j-1))
            if ( Tsbc .and. thin_snow_active ) &
            matj(ni+1,j) = matj(ni+1,j) - rhoice * em(j-1) *      zi(j-1)
           else
            matj(0   ,j) = matj(0   ,j) - rhoice * em(j  ) * (1d0-zi(j-1))
            if ( Tsbc .and. thin_snow_active ) &
            matj(ni+1,j) = matj(ni+1,j) - rhoice * em(j  ) *      zi(j-1)
           endif
          enddo
         else
          do j=2,ni
            matj(0   ,j) = matj(0   ,j) - rhoice * em(j-1) * (1d0-zi(j-1)) * 0.5d0
            matj(0   ,j) = matj(0   ,j) - rhoice * em(j  ) * (1d0-zi(j-1)) * 0.5d0
           if ( Tsbc .and. thin_snow_active ) then
            matj(ni+1,j) = matj(ni+1,j) - rhoice * em(j-1) *      zi(j-1)  * 0.5d0
            matj(ni+1,j) = matj(ni+1,j) - rhoice * em(j  ) *      zi(j-1)  * 0.5d0
           endif
          enddo
         endif
   
!-----------------------------------------------------------------------
! upper transport in ice
!-----------------------------------------------------------------------

         if (adv_upwind) then
          do j=1,ni-1
           if (w(j  )>0.d0) then
            matj(j  ,j) = matj(j  ,j) + rhoice * w(j  ) * cp(j  )
            DT0(j)      = DT0(j)      - rhoice * w(j  ) * em(j  )
           else
            matj(j+1,j) = matj(j+1,j) + rhoice * w(j  ) * cp(j+1)
            DT0(j)      = DT0(j)      - rhoice * w(j  ) * em(j+1)
           endif
          enddo
         else
          do j=1,ni-1
            matj(j  ,j) = matj(j  ,j) + rhoice * w(j  ) * cp(j  ) * 0.5d0
            DT0(j)      = DT0(j)      - rhoice * w(j  ) * em(j  ) * 0.5d0
            matj(j+1,j) = matj(j+1,j) + rhoice * w(j  ) * cp(j+1) * 0.5d0
            DT0(j)      = DT0(j)      - rhoice * w(j  ) * em(j+1) * 0.5d0
          enddo
         endif

! surface treatment assuming melting, looks like an upwind formulation (em=cp_ice*Tf at surface)
! if no melting w=0 anyway
         j=ni
            DT0(j)      = DT0(j)      - rhoice * w(j  ) * em(j)
! FD debug
write(*,*) 'adv up',j,dt0(j),w(j)

! metric terms, newton terms for change in thickness
         if (adv_upwind) then
          do j=1,ni-1
           if (w(j  )>0.d0) then
            matj(0   ,j) = matj(0   ,j) + rhoice * em(j  ) * (1d0-zi(j  ))
            if ( Tsbc .and. thin_snow_active ) &
            matj(ni+1,j) = matj(ni+1,j) + rhoice * em(j  ) *      zi(j  )
           else
            matj(0   ,j) = matj(0   ,j) + rhoice * em(j+1) * (1d0-zi(j  ))
            if ( Tsbc .and. thin_snow_active ) &
            matj(ni+1,j) = matj(ni+1,j) + rhoice * em(j+1) *      zi(j  )
           endif
          enddo
         else
          do j=1,ni-1
            matj(0   ,j) = matj(0   ,j) + rhoice * em(j  ) * (1d0-zi(j  )) * 0.5d0
            matj(0   ,j) = matj(0   ,j) + rhoice * em(j+1) * (1d0-zi(j  )) * 0.5d0
           if ( Tsbc .and. thin_snow_active ) then
            matj(ni+1,j) = matj(ni+1,j) + rhoice * em(j  ) *      zi(j  )  * 0.5d0
            matj(ni+1,j) = matj(ni+1,j) + rhoice * em(j+1) *      zi(j  )  * 0.5d0
           endif
          enddo
         endif
         j=ni ! upper transport in ice+snow is zero, therefore it only remains case of thin_snow and top melt
           if ( Tsbc .and. thin_snow_active ) then
            matj(ni+1,j) = matj(ni+1,j) + rhoice * em(j) *      zi(j  )
           endif


!-----------------------------------------------------------------------
! lower transport in snow
!-----------------------------------------------------------------------
      if ( .not.thin_snow_active ) then

         if (adv_upwind) then
          do j=ni+1,ns
           if (w(j-1)>0.d0) then
            matj(j-1,j) = matj(j-1,j) - rhosno * w(j-1) * cp(j-1)
            DT0(j)      = DT0(j)      + rhosno * w(j-1) * em(j-1)
           else
            matj(j  ,j) = matj(j  ,j) - rhoice * w(j-1) * cp(j  )
            DT0(j)      = DT0(j)      + rhoice * w(j-1) * em(j  )
           endif
          enddo
         else
          do j=ni+1,ns
            matj(j-1,j) = matj(j-1,j) - rhosno * w(j-1) * cp(j-1) * 0.5d0
            DT0(j)      = DT0(j)      + rhosno * w(j-1) * em(j-1) * 0.5d0
            matj(j  ,j) = matj(j  ,j) - rhosno * w(j-1) * cp(j  ) * 0.5d0
            DT0(j)      = DT0(j)      + rhosno * w(j-1) * em(j  ) * 0.5d0
          enddo
         endif


! metric terms, newton terms for change in thickness
        if ( Tsbc ) then
         if (adv_upwind ) then
          do j=ni+1,ns
           if (w(j-1)>0.d0) then
            matj(ns+1,j) = matj(ns+1,j) - rhosno * em(j-1) *      zs(j-1)
           else
            matj(ns+1,j) = matj(ns+1,j) - rhosno * em(j  ) *      zs(j-1)
           endif
          enddo
         else
          do j=ni+1,ns
            matj(ns+1,j) = matj(ns+1,j) - rhosno * em(j-1) *      zs(j-1)  * 0.5d0
            matj(ns+1,j) = matj(ns+1,j) - rhosno * em(j  ) *      zs(j-1)  * 0.5d0
          enddo
         endif
        endif

!-----------------------------------------------------------------------
! upper transport in snow
!-----------------------------------------------------------------------

         if (adv_upwind) then
          do j=ni+1,ns-1
           if (w(j  )>0.d0) then
            matj(j  ,j) = matj(j  ,j) + rhosno * w(j  ) * cp(j  )
            DT0(j)      = DT0(j)      - rhosno * w(j  ) * em(j  )
           else
            matj(j+1,j) = matj(j+1,j) + rhosno * w(j  ) * cp(j+1)
            DT0(j)      = DT0(j)      - rhosno * w(j  ) * em(j+1)
           endif
          enddo
         else
          do j=ni+1,ns-1
            matj(j  ,j) = matj(j  ,j) + rhosno * w(j  ) * cp(j  ) * 0.5d0
            DT0(j)      = DT0(j)      - rhosno * w(j  ) * em(j  ) * 0.5d0
            matj(j+1,j) = matj(j+1,j) + rhosno * w(j  ) * cp(j+1) * 0.5d0
            DT0(j)      = DT0(j)      - rhosno * w(j  ) * em(j+1) * 0.5d0
          enddo
         endif

! force upwind treatment for top surface
        j=ns
            DT0(j)      = DT0(j)      - rhosno * w(j  ) * em(j+1)
        if (.not.Tsbc) then
            matj(j+1,j) = matj(j+1,j) + rhosno * w(j  ) * cp(j+1)
        endif

! metric terms, newton terms for change in thickness
        if ( Tsbc ) then
         if (adv_upwind) then
          do j=ni+1,ns-1
           if (w(j  )>0.d0) then
            matj(ns+1,j) = matj(ns+1,j) + rhosno * em(j  ) *      zs(j  )
           else
            matj(ns+1,j) = matj(ns+1,j) + rhosno * em(j+1) *      zs(j  )
           endif
          enddo
         else
          do j=ni+1,ns-1
            matj(ns+1,j) = matj(ns+1,j) + rhosno * em(j  ) *      zs(j  )  * 0.5d0
            matj(ns+1,j) = matj(ns+1,j) + rhosno * em(j+1) *      zs(j  )  * 0.5d0
          enddo
         endif
         j=ns
            matj(ns+1,j) = matj(ns+1,j) + rhosno * em(j+1) *      zs(j  )
        endif

      endif ! active snow

!-----------------------------------------------------------------------
! prepare linear solver
!-----------------------------------------------------------------------

      residual=0d0
      do j=1,ntot-1
         residual=residual+abs(dt0(j))
! FD debug
if (extra_debug) write(*,*) 'res',j,dt0(j)
      enddo
      if (debug) write(*,*) 'residual',counter,residual
! FD debug
!write(*,*) 'mat'
!do j=0,ntot-1
! write(*,'(6(e10.3,1x),3x,e10.3)') matj(0:ntot-1,j),dt0(j)
!enddo
        ml = ntot - 1
        mu = ml
        MT = ML + MU + 1
!        write(*,*) 'bandwidth',MT, ML
!        bandmax=ML+MT
        bandmax=2*maxlay ! hardcoded for now
      do i=1,ntot
       do k=1,bandmax
         matband(k,i)=0.d0
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
!write(*,*) 'solution'
!do j=0,ntot-1
! write(*,'(e10.3)') tout(j)
!enddo

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
      DO j=1,ni
         temp(j,1) = temp(j,0) + tout(j)
      ENDDO

!treatment of snow
     if ( .not.thin_snow_active ) then

! snow column increment to temperature
      DO j=ni+1,ns
         temp(j,1) = temp(j,0) + tout(j)
      ENDDO
      if (Tsbc) then
! melting case
        w(ns) = w(ns) + tout(ns+1)
        temp(ns+1,1) = Tf(ns+1)
      else
! cold surface case, update temperature
        temp(ns+1,1) = temp(ns+1,0) + tout(ns+1)
      endif

     else ! thin ice case

      if (Tsbc) then
! melting case
        w(ni) = w(ni) + tout(ni+1)
        temp(ni+1,1) = Tf(ni+1)
      else
! cold surface case, update temperature
        temp(ni+1,1) = temp(ni+1,0) + tout(ni+1)
      endif
     endif

         DO j=1,ntot-1
            IF (temp(j,1) .GT. Tf(j)+tiny) THEN
! FD debug
if (debug) then
 write(*,*) 'detecting Temp pb at',j,temp(j,1),Tf(j)
endif
               temp(j,1) = Tf(j)
!               Fnet = Fnet + sum(R(1:ns))
!               Fnet0= Fnet0+ sum(R(1:ns))
!               R(1:ns) = 0.d0
!               kki(j  )=kki(j  )*10.d0
!               kki(j+1)=kki(j+1)*0.1d0
            ENDIF
         ENDDO

! in any case, the temperature of snow is overwritten by that of ice
     if ( thin_snow_active ) then
      temp(ni+2:ns+1,1) = temp(ni+1,1)
     endif

!------------------------------------------------------------------------
!     Update conductive heat fluxes
!------------------------------------------------------------------------

      IF ( .not.thin_snow_active ) THEN
         Fcss = -kks(ns) * ( temp(ns+1,1) - temp(ns  ,1) )
         Fcsb = -kks(ni) * ( temp(ni+1,1) - temp(ni  ,1) )
      ELSE
         Fcss = 0d0
         Fcsb = 0d0
      ENDIF
      
      Fcis = -kki(ni) * ( temp(ns+1,1) - temp(ni  ,1))  ! only valid for no-snow case
      Fcib = -kki(0 ) * ( temp(   1,1) - temp(   0,1))

      Fnet = Fnet0 - dzf * (temp(ns+1,1)-tiold(ns+1))

! FD debug
if (extra_debug) then
 write(*,'(I6,300(1x,e9.3))') counter, w(0:ns)
 write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns+1,1)
! write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns+1,0)
! stop
endif
! FD debug
if (extra_debug) then
 if ( .not.thin_snow_active ) then
  write(*,*) 'SBCond snow',Tsbc,Fnet,Fcss, temp(ns+1,1)
 else
  write(*,*) 'SBCond ice',Tsbc,Fnet,Fcis,kki(ni), temp(ns+1,1), temp(ni,1)
endif
endif

     if ( .not.thin_snow_active ) then
      IF (( temp(ns+1,1)+tiny .GE. Tf(ns+1)) .AND. (Tsbc .EQV. .FALSE.)) THEN
         Tsbc = .TRUE. 
      ENDIF
      IF ( temp(ns+1,1)+tiny .GE. Tf(ns+1) .AND. -Fnet-Fcss .GT. 0d0 ) THEN
         Tsbc = .FALSE.
         temp(ns+1,1)=Tf(ns+1)-tiny
     ENDIF
     else
      IF ((temp(ns+1,1)+tiny .GE. Tf(ns+1)) .AND. (Tsbc .EQV. .FALSE.)) THEN
         Tsbc = .TRUE. 
      ENDIF
      IF ( temp(ns+1,1)+tiny .GE. Tf(ns+1) .AND. -Fnet-Fcis .GT. 0d0 .and. w(ni)==0.d0) THEN
         Tsbc = .FALSE.
         temp(ns+1,1)=Tf(ns+1)-tiny
      ENDIF
     endif

!-----------------------------------------------------------------------
!     Update T-S dependent ice parameters
!     Lf, cp are diagnostics for output
!-----------------------------------------------------------------------
      
      sali(0:ns+1,1)=sali(0:ns+1,0) ! FD not sure what is really done in this code, just for debug

      DO j=0,ni+1
         ki(j)   = func_ki(sali(j,1),temp(j,1))
      ENDDO
      DO j=0,ns+1
         Tf(j)   = Tfreeze1(sali(j,1))
         qm(j)   = func_qm(Tf(j),temp(j,1))
         cp(j)   = func_cph(Tf(j),temp(j,1))
         em(j)   = func_el(Tf(j),temp(j,1))
      ENDDO

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
         dTmax = 0d0
         GOTO 4000				! recalculate T profile
      ENDIF

!-----------------------------------------------------------------------
!     Compute internal energy
!-----------------------------------------------------------------------

      Frad=0.d0
      DO j=ni+1,ns
         Frad	= Frad + R(j)
      ENDDO
      DO j=1,ni
         Frad	= Frad + R(j)
      ENDDO
      
!-----------------------------------------------------------------------
!     Update energy flux due to snow precipitation
!-----------------------------------------------------------------------
      
      if ( .not.thin_snow_active ) then
          Fprec = snow_precip * rhosno * em(ns+1)
      else
          Fprec = snow_precip * rhosno * em_thin_snow
      endif

      Einp = (Fnet+oceflx+Frad+Fprec+Fthin_snow)*dtice	! E input [J/m2]

      
      Eint = 0d0

      DO j=1,ni
         Elay = func_El(Tf(j),temp(j,1))*henew(j)*rhoice
         Eint = Eint+Elay
      ENDDO
     if (thin_snow_active) then
      DO j=ni+1,ns
         Elay = em_thin_snow * henew(j) * rhosno
         Eint = Eint+Elay
      ENDDO
     else
      DO j=ni+1,ns
         Elay = func_El(Tf(j),temp(j,1))*henew(j)*rhosno
         Eint = Eint+Elay
      ENDDO
     endif
     
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
      
      IF (dibdt .LT. 0d0) gbase = gbase-dibdt*dtice
      IF (dibdt .GT. 0d0) mbase = mbase+dibdt*dtice
      
!-----------------------------------------------------------------------
!     Reset prognostic variables
!-----------------------------------------------------------------------
      
! FD debug
if (debug) then
 write(*,*) 'snow',hs_b,-w(ns)
 write(*,*) 'hice',hi_b,-w(ni),w(0)
 write(*,'(I6,300(1x,e9.3))') counter, w(0:ns)
 write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns+1,1)
endif
! FD debug if energy not conserved
      IF ( abs(Einp-dEin)/dtice > 1e-3 ) THEN
       OPEN(1,file='restart.dat')
       WRITE(1,*) nlice,nlsno,ith_cond,dtice
       WRITE(1,*) ti(1:nlice),ts(1:nlsno),tsuold,tbo
       WRITE(1,*) si(1:nlice)
       WRITE(1,*) hi,hsold
       WRITE(1,*) sfallold,dwnlw,tsuold,tair,qair,uair,swrad,oceflx,pres
       WRITE(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
       WRITE(1,*) fsens,flat
       CLOSE(1)
       STOP 'energy not conserved'
      ENDIF

      IF (hi_b(1) .LT. hicut) THEN
         hi_b(1) = hicut
         DO j=0,ni
            temp(j,1) = tocn
         ENDDO
      ENDIF
      
      hs_b(0) = hs_b(1)
      hi_b(0) = hi_b(1)
      
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

      dh_sni = MAX( 0.d0 , ( rhosno * hs + (rhoice - rhowat ) * hi) / ( rhosno + rhowat - rhoice ) )

      if (dh_sni > 0.d0 ) then
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
      DOUBLE PRECISION a(n),b(n),c(n),r(n),u(n)
      PARAMETER(NMAX=230)
      INTEGER j
      DOUBLE PRECISION bet,gam(NMAX)
      if(b(1).eq.0d0) stop 'tridag: check BC'
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j)*gam(j)
         if(bet.eq.0d0) stop 'tridag: failed'
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
      DOUBLE PRECISION	func_i0, hs, hscut
      DOUBLE PRECISION :: i0snow    =   0.080d+0 ! SW fraction penetrating snow surface
      DOUBLE PRECISION :: i0ice     =   0.170d+0	! SW fraction penetrating ice  surface

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
      DOUBLE PRECISION	esat, T, coef1, coef2 	! T=[C]
      DOUBLE PRECISION :: temp0	= 273.160d+0	! freezing point temp of fresh water	[K]

      coef1	= 21.87456d0			! coeff. over ice
      coef2	=  7.66000d0			! coeff. over ice
      esat	=  6.11d2*EXP(MIN(coef1*T/(T+temp0-coef2),10d0))

      END FUNCTION esat

!-----------------------------------------------------------------------
!     Thickness dependent ice salinity (Cox and Weeks, 1974)
!-----------------------------------------------------------------------

      FUNCTION sal(hi)

      implicit none
      DOUBLE PRECISION	sal, hi, hic
      DOUBLE PRECISION :: salinib	=   4.000d+0    ! salinity at ice base			[psu]

      hic = 0.57d0

      IF (hi .LT. hic) THEN
         sal = 14.24d0 - 19.39d0 * hi
      ELSE
         sal = salinib
      ENDIF

      END FUNCTION sal



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------


