module ice_thermodynamic_FV

  use var_thermo_vertical

implicit none

contains

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

      integer :: ni, ns

      DOUBLE PRECISION :: &
        Fbase,&! heat flux at ice base			[W/m2]
        Fsen ,&! sensible heat flux			[W/m2]
        Flat ,&! latent heat flux			[W/m2]
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
     dzic,&! delta z for ice				[-]
     dzsc,&! delta z for snow				[-]
     epsilon,&! coefficient of emissivity			[-]
     hinit,&! initial ice thickness				[m]
     hicut,&! cut off (minimum) ice  thickness		[m]
     hscut,&! cut off (minimum) snow thickness		[m]
     hslim,&! minimum snow thickness			[m]
     i0ice,  &! fraction of sw rad penetrating the surface	[-]
     i0snow, &! fraction of sw rad penetrating the surface	[-]
     ki0,&! thermal heat conductivity of ice		[W/m/C]
     ks0,&! thermal heat conductivity of snow		[W/m/C]
     kappai,&! extinction coefficient of ice			[1/m]
     kappas,&! extinction coefficient of snow		[1/m]
     Lf0,&! specific latent heat of fusion of ice/snow 	[J/kg]
     Lsub,&! specific latent heat of sublim of ice/snow	[J/kg]
     Lvap,&! specific latent heat of vapori of ice/snow	[J/kg]
     mu,&! empirical constant relating S and T		[C/psu]
     salinis,&! salinity of surf. sea ice			[psu]
     salinib,&! salinity of basal sea ice			[psu]
     salnice,&! salinity of new basal ice			[psu]
     zic(0:maxlay),&! transformed z-levels in the ice  (C-grid)	[-]
     zsc(0:maxlay),&! transformed z-levels in the snow (C-grid)	[-]
     zib(0:maxlay),&! transformed z-levels in the ice  (B-grid)	[-]
     zsb(0:maxlay) ! transformed z-levels in the snow (B-grid)	[-]

      logical   Tsbc           ! Fixed T surface boudnary conditions

      DOUBLE PRECISION :: &
     dhidt,&! total rate of change of ice  thickness	[m/s]
     dhsdt,&! total rate of change of snow thickness	[m/s]
     dspdt,&! rate of change of snow thickn. due to precip. [m/s]
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
    lh(0:maxlay,0:1),&! bulk latent heat				[J*m/kg]
    em(0:maxlay)    ,&! volumetric enthalpiy	[W/m3]
    qm(0:maxlay)    ,&! volumetric energy of melt (rho*Lf(S,T))	[W/m3]
    cp(0:maxlay)    ,&! heat capacity for ice	[W/kg/K^-1]
    R(0:maxlay)    ,&! penetrating shortwave radiation		[W/m3]
    sh(0:maxlay,0:1),&! bulk specific heat			[J/kg/C]
    sali(0:maxlay,0:1),&! internal sea ice salinity			[psu]
    Tf(0:maxlay)    ,&! ice freezing temp. of salinity S		[C]
    w(1:maxlay-1)  ,&! grid advection velocity			[1/s]
    hi_b(0:1)     ,&! ice thickness (m)
    hs_b(0:1)       ! snow thickness (m)

     double precision :: &
     func_i0,&
     esat,&
     sal

      double precision &
                temp (0:maxlay,0:1) ! internal sea ice temperature		[C]

      DOUBLE PRECISION AT1(0:maxlay), BT1(0:maxlay), CT1(0:maxlay), DT0(0:maxlay), &
      C1i, C2i, C3i, C4i,    C1s, C2s, C3s, C4s, &
      OT1, PT1, QT1,&
      tout(0:maxlay), hsmemo, Rtrans

      double precision :: sice, &
             latmelt, ftot, sstznew, &
             Tside,deltaT,wlat,rside, &
             sum0, dhi, dho, dq, fwf0
      double precision :: hminice=1d-3
      double precision :: energy_bot

  character (len=20) :: Sprofile
  integer i,j,counter
  double precision dzf, es, zssdqw, zrchu1, zrchu2, q0, zref
  double precision, dimension(0:maxlay) :: tiold, tinew
  double precision elays0
  logical :: debug=.false.
  logical :: extra_debug=.false.
! FD debug
 debug=.true.
! extra_debug=.true.

!------------------------------------------------------------------------
!     Some integer constants
!------------------------------------------------------------------------

      ni=   nlice+1
      ns=ni+nlsno+1 
      if (ns > maxlay ) then
       write(*,'(a,i4,a,i4)') ' ni',ni,' ns',ns
       write(*,*) 'please redefine the number of ice/snow levels in ice_thermoynamic.f90'
       stop
      endif

!------------------------------------------------------------------------
!     Some logical constants
!------------------------------------------------------------------------

      salinity  = 'no'          ! = no dynamic salinity
      bbc       = 'fixT'        ! bottom boundqry condition
! FD debug
!      bbc       = 'flux'        ! bottom boundqry condition
      sbc       = 'flux'        ! surface boundqry condition

!------------------------------------------------------------------------
!     Some physical constants
!------------------------------------------------------------------------

      Lf0=   3.335d+5! spec. latent heat of fusion, ice/snow [J/kg]
      Lsub=   2.834d+6! spec. latent heat of sublim, ice/snow [J/kg]
      Lvap=   2.501d+6! spec. latent heat of vapori, ice/snow [J/kg]

      epsilon=   0.990d+0! snow/ice emissivity			[-]

      ki0=   2.034d+0! thermal conductivity of fresh ice	[W/m/C]
      ks0=   0.500d+0! thermal conductivity of snow		[W/m/C]

      sigma=   5.670d-8! Stefan-Boltzmann constant	       [W/m2/K4]
      temp0= 273.160d+0! freezing point temp of fresh water	[K]

      uiofix=   0.000d+0! constant relative ice-ocean velocity	[m/s]


!------------------------------------------------------------------------
!     Some model constants and parameters
!------------------------------------------------------------------------

      salinis	=   1.000d+0    ! salinity at ice surf			[psu]
      salinib	=   4.000d+0    ! salinity at ice base			[psu]
! FD test      salnice	=  10.000d+0    ! salinity of newly formed ice		[psu]
      salnice	=  4.000d+0    ! salinity of newly formed ice		[psu]

      Fbase	=   oceflx	! conductive heat flux at ice base 	[W/m2]

      dzic       = 1d0/DBLE(nlice)	! delta z for ice			[-]
      dzsc       = 1d0/DBLE(nlsno)	! delta z for snow			[-]


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

      Tdiff	=   1.000d-8	! temperature tolerance in Euler step	[C]
      tiny	=   1.000d-9	! very small number			[-]

      hicut	=   0.010d+0	! minimum ice  thickness		[m]
      hscut	=   0.050d+0	! threshold snow thickness		[m]
      hslim	=   0.005d+0	! minimum snow thickness		[m]
! FD debug
hscut=0.01d0
hscut=0.001d0
hslim=0.0005d0
Tdiff=1d-12

!------------------------------------------------------------------------
!     z-coordinates of T and S,  C-grid (non-dimensional)
!------------------------------------------------------------------------

! init variable
      hs_b(0:1)=hs
      hi_b(0:1)=hi

! conversion from LIM3
      temp(0,:)=tbo - temp0
      do j=1,nlice
        temp(j,:)=ti(nlice-j+1)  - temp0
      enddo
      temp(ni,:)=ts(nlsno+1)  - temp0
      do j=1,nlsno
        temp(j+ni,:)=ts(nlsno-j+1)  - temp0
      enddo
      Tsurf = tsu - temp0
      temp(ns,:) = Tsurf
      tocn  = temp(0,0)

!------------------------------------------------------------------------
!     z-coordinates of vertical temperature, B-grid (non-dimensional)
!------------------------------------------------------------------------

      zic=0.d0
      zsc=0.d0
      DO j=0,ns
         IF (j .EQ. 0) THEN
            zic(j) = 0d0
         ELSEIF (j .EQ. 1) THEN
            zic(j) = dzic/2d0
         ELSEIF (j .GT. 1 .AND. j .LT. ni) THEN
            zic(j) = dzic*(j-0.5d0)
         ELSEIF (j .EQ. ni) THEN
            zic(j) = 1d0
            zsc(j) = 0d0
         ELSEIF (j .EQ. ni+1) THEN
            zsc(j) = dzsc/2d0
         ELSEIF (j .GT. ni+1 .AND. j .LT. ns) THEN
            zsc(j) = dzsc*(j-ni-0.5d0)
         ELSEIF (j .EQ. ns) THEN
            zsc(j) = 1d0
         ENDIF
      ENDDO
      tiold=temp(:,0)

!------------------------------------------------------------------------
!     z-coordinates of vertical velocity, B-grid (non-dimensional)
!------------------------------------------------------------------------

      DO j=1,ns-1
         IF (j .LT. ni) THEN
            zib(j) = dzic*(j-1)
         ELSEIF (j .EQ. ni) THEN
            zib(j) = 1d0
            zsb(j) = 0d0
         ELSEIF (j .GT. ni) THEN
            zsb(j) = dzsc*(j-ni)
         ENDIF
      ENDDO


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
      sali(ni,0)= salinis
         DO j=ni+1,ns
            sali(j,0) = 0d0
         ENDDO

      sali(:,1)=sali(:,0) ! FD not sure what is really done in this code, just for debug

1000 continue
!------------------------------------------------------------------------
!     Vertical (grid) velocities
!------------------------------------------------------------------------

      DO j=1,ns-1
         w(j) = 0d0		! initial grid advection
      ENDDO

!------------------------------------------------------------------------
!     Set initial T-S dependent material properties
!------------------------------------------------------------------------

      DO j=0,ns
            Tf(j)   = Tfreeze1(sali(j,0))
            ki(j)   = func_ki(sali(j,0),temp(j,0))
            qm(j)   = func_qm(Tf(j),temp(j,0))
            em(j)   = func_el(Tf(j),temp(j,0))
            cp(j)   = func_cp(Tf(j),temp(j,0),tiold(j))
            sh(j,0) = func_sh(Tf(j),temp(j,0))
         IF (j .LE. ni) THEN
            lh(j,0) = func_lh(Tf(j),temp(j,0))*hi*dzic
         ELSEIF (j .GT. ni) THEN
            lh(j,0) = func_lh(Tf(j),temp(j,0))*hs*dzsc
         ENDIF
      ENDDO

!     WARNING LH, SH SHOULD NOT EVEN BE DEFINED AT THE NODES, THIS IS TROUBLE

      DO j=0,ns
         sh(j,1) = sh(j,0)
         lh(j,1) = lh(j,0)
      ENDDO

      DO j=0,ns					! snow thermal conductivity
         IF (j .LT. ni) THEN
            ks(j) = 0d0
         ELSEIF (j .GE. ni) THEN
            ks(j) = cond_sno
         ENDIF
      ENDDO

!------------------------------------------------------------------------
!     Initial internal energy in ice (1:ni-1) and snow (ni+1:ns-1)
!------------------------------------------------------------------------

      Einp = 0d0				! E input [J/m2]
      Ein0 = 0d0

      DO j=1,ni-1
         Elay = func_El(Tf(j),temp(j,0))*rhoice*hi_b(0)*dzic
         Ein0 = Ein0+Elay
      ENDDO

      DO j=ni+1,ns-1
         if (hs_b(0)>hslim) then
            Elay = func_El(Tf(j),temp(j,0))*rhosno*hs_b(0)*dzsc
         ELSE
            Elay = 0d0
         ENDIF
         Ein0 = Ein0+Elay
      ENDDO

      dEin = 0d0				! diff intl E [J/m2]

!
!------------------------------------------------------------------------------|
!  7) surface flux computation                                                 |
!------------------------------------------------------------------------------|
!
! FD assume only snow
!     if ( hs_b(0) > hscut ) then
     if ( hs_b(0) > hslim ) then
            Fprec = snowfall*rhosno*(-mlfus+cp_ice*temp(ns,1))
     else
            Fprec = 0.d0
     endif
            Frain = 0d0

! FD find correct surface temperature
      IF (hs_b(0) > hslim) THEN
        tsu = min(tsu,tf(ns)+temp0)
      ELSE
        tsu = min(tsu,tf(ni)+temp0)
      ENDIF

      ! net longwave radiative flux
      netlw = emi*(dwnlw - stefa*tsu*tsu*tsu*tsu)
      ! sensible and latent heat flux
      CALL flx(hi,tsu,tair,qair,fsens, &
     &         flat,q0,zrchu1,zrchu2,uair,zref)
      fsens   =  -fsens
      flat   =  MIN( -flat , 0.d0 ) ! always negative, as precip 
                                           ! energy already added

      ! pressure of water vapor saturation (Pa)
      es         =  611.d0*10.d0**(9.5d0*(tsu-273.16d0)/(tsu-7.66d0))
      ! intermediate variable
      zssdqw     =  q0*q0*pres/ &
     &              (0.622d0*es)*log(10.d0)*9.5d0* &
     &              ((273.16d0-7.66d0)/(tsu-7.66d0)**2.d0)
      ! derivative of the surface atmospheric net flux
          dzf    =  4.d0*emi*stefa*tsu*tsu*tsu + zrchu1+zrchu2*zssdqw

      ! surface atmospheric net flux
      Fnet0 =  fac_transmi * swrad + netlw + fsens + flat
      Fnet  = Fnet0

!------------------------------------------------------------------------
!     Initial conductive heat fluxes
!------------------------------------------------------------------------

!     IF (hs_b(0) .GT. hscut) THEN
     IF (hs_b(0) .GT. hslim) THEN
!FD add case for nlsno=1
      if ( nlsno == 1 ) then
       Fcss	= -ks(ns) &
                  *( 2d0*temp(ns,0)-2d0*temp(ns-1,0))/hs_b(0)

       Fcsb	= -ks(ni) &
                  *(-2d0*temp(ni,0)+2d0*temp(ni+1,0))/hs_b(0)
      else
       Fcss	= -ks(ns) &
                  *( 8d0*temp(ns,0)-9d0*temp(ns-1,0)+temp(ns-2,0)) &
                  /3d0/dzsc/hs_b(0)

       Fcsb	= -ks(ni) &
                  *(-8d0*temp(ni,0)+9d0*temp(ni+1,0)-temp(ni+2,0)) &
                  /3d0/dzsc/hs_b(0)
      endif
!     else if ( hs_b(0) <=hscut .and. hs_b(0) > 0.d0) then ! FD thin snow
!      Fcss	= -ks(ns) &
!                  *(temp(ns,0)-temp(ni,0)) &
!                  /hs_b(0)
!
!      Fcsb	= -ks(ns) &
!                  *(temp(ns,0)-temp(ni,0)) &
!                  /hs_b(0)
      ELSE
         Fcss	= 0d0
         Fcsb	= 0d0
      ENDIF
      
      Fcis	= -ki(ni) &
                 *(8d0*temp(ni,0)-9d0*temp(ni-1,0)+temp(ni-2,0))/3d0/dzic/hi_b(0)
      
      Fcib	= -ki(0) &
                 *(-8d0*temp(0,0)+9d0*temp(1,0)-temp(2,0))/3d0/dzic/hi_b(0)

      Fmlt = Fnet-Fcis

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
!     Update snow and ice thickness
!-----------------------------------------------------------------------

      dssdt = 0d0
      disdt = 0d0

      dsbdt = 0d0		! melt of snow base not allowed
         ! snow depth evolution due to precipitation

      dssdt = snowfall

! FD: if snow present
      if (hs_b(0) > hslim) then
         ! snow depth evolution due to melt superposed on precip

         IF (temp(ns,1)+tiny .GE. Tf(ns) .AND. -Fnet-Fcss .LT. 0d0) THEN
! FD            dssdt = dssdt+MIN((-Fnet-Fcss)/rhosno/qm(ns-1),0d0)
! FD the code cannot deal with concomittent snow accumulation and melt
            dssdt = MIN((-Fnet-Fcss)/rhosno/qm(ns-1),0d0)
            snowfall =0.d0
         ENDIF
      else

      ! update ice surface

         IF (temp(ni,1)+tiny .GE. Tf(ni) .AND. -Fnet-Fcis .LT. 0d0) THEN
           disdt = MIN((-Fnet-Fcis)/rhoice/qm(ni-1),0d0)
         ENDIF
      endif
!c     WARNING WE HAVE TO CHECK THE CASE IF SALNICE NE SALI(0,1)
!c     NOTE Hendrik is coming up with a better way for this.
!c     Using salnice ~= salinib does not conserve energy, 
!c     the idea is that oceflx should compensate

      IF (oceflx-Fcib .LT. 0d0) THEN		! for ice growth: new-ice salin.
         qm(0) = func_qm(Tf(0),temp(0,1))
         energy_bot = qm(0)
      ELSEIF (oceflx-Fcib .GT. 0d0) THEN	! qm_mean_first_layer for melt
         energy_bot = qm(1)
      ENDIF
! FD generalized bottom
         heatflx_bot = oceflx-Fcib
         dibdt	= heatflx_bot / rhoice / energy_bot

      dhsdt	= dssdt - dsbdt
      dhidt	= disdt - dibdt ! original

      hs_b(1)	= hs_b(0) + dhsdt*dtice

! FD need to do something when hs<0
      heatflx_bot = 0.d0
      if ( hs_b(1) < hslim .and. hs_b(0) > hslim) then
         hs_b(0:1)=0.d0
         temp(ni+1:ns,0:1)=Tf(ns)
         goto 1000
! how much energy is left for melting ice? FD
!         heatflx_bot = -qm(ns-1)*hs_b(0)*rhos/dtice
!         Fnet = Fnet + heatflx_bot
!         disdt = MIN((-Fnet-Fcis)/rhoi/qm(ni-1),0d0)
      endif

      hi_b(1)	= hi_b(0) + dhidt*dtice

! FD      fwf = fwf0 + ( dhsdt*rhos + dhidt*rhoi ) / rhoo

! FD should remove the following lines in future
      DO j=0,ni
         Tf(j)   = Tfreeze1(sali(j,1))
         lh(j,1) = func_lh(Tf(j),temp(j,1))*hi_b(1)*dzic
      ENDDO
      DO j=ni+1,ns
         Tf(j)   = Tfreeze1(sali(j,1))
         lh(j,1) = func_lh(Tf(j),temp(j,1))*hs_b(1)*dzsc
      ENDDO

!-----------------------------------------------------------------------
!     Vertical velocities (grid advection)
!-----------------------------------------------------------------------

      DO j=1,ni
         w(j)	= -zib(j)*disdt-(1d0-zib(j))*dibdt
      ENDDO
! FD modif
!     if (hs_b(0) > hscut ) then
     if (hs_b(0) > hslim ) then
      DO j=ni+1,ns-1
         w(j)	= -zsb(j)*dssdt-(1d0-zsb(j))*dsbdt
      ENDDO
     endif

!-----------------------------------------------------------------------
!     Define constant coefficients
!-----------------------------------------------------------------------

      C1i	= dtice/2d0/rhoice/hi_b(1)/dzic
      C2i	= dtice
      C3i	= dtice*hi_b(1)*dzic/rhoice
      C4i	= 3d0*hi_b(1)*dzic

!     if (hs_b(0) > hscut ) then
     if (hs_b(0) > hslim ) then
      C1s	= dtice/2d0/rhosno/hs_b(1)/dzsc
      C2s	= dtice
      C3s	= dtice*hs_b(1)*dzsc/rhosno
      C4s	= 3d0*hs_b(1)*dzsc
      if (nlsno ==1 ) c1s = dtice/rhosno/hs_b(1)
      if (nlsno ==1 ) c4s = hs_b(1)
!     else if ( hs_b(0) <= hscut .and. hs_b(0) > hslim ) then
!      C4s	= hs_b(1)
     endif

!-----------------------------------------------------------------------
!     Absorbed SW radiation energy per unit volume (R)       
!     1/dz Int_z^z+dz R dz = Rtrans(z+dz) - Rtrans(z)
!-----------------------------------------------------------------------

!      if (hs_b(1) > hscut ) then
      if (hs_b(1) > hslim ) then
        DO j=ns-1,ni+1,-1
          R(j)   = swradab_s(ns-j) /dzsc/hs_b(1)
        ENDDO
      else
         R(ni+1:ns-1)=0.d0
      endif
      R(ni)=0.d0
      DO j=ni-1,1,-1
         R(j)   = swradab_i(ni-j) /dzic/hi_b(1)
      ENDDO
      R(0)=0.d0

!-----------------------------------------------------------------------
!     Tri-diagonal matrix coefficients [A] {x} = {D}
!-----------------------------------------------------------------------

3000  CONTINUE

         AT1=0.d0; BT1=0.d0; CT1=0.d0; DT0=0.d0;

         AT1(1)	= -C1i/3d0*16d0*ki(0)
         BT1(1)	=  C1i/3d0*(3d0*ki(2)+3d0*ki(1)+18d0*ki(0)) &
                  +cp(1)*hi_b(1)*dzic
         CT1(1)	= -C1i/3d0*(3d0*ki(2)+3d0*ki(1)+2d0*ki(0))
         DT0(1)	=  C3i*R(1)+cp(1)*tiold(1)*hi_b(1)*dzic &
                  -em(1)*(hi_b(1)-hi_b(0))*dzic

      IF (w(1) .GT. 0d0) THEN
         AT1(1) = AT1(1) - C2i * w(1) * cp(0)
         DT0(1) = DT0(1) + C2i * w(1) *(em(0)-cp(0)*tiold(0))
      ELSE
         BT1(1) = BT1(1) - C2i * w(1) * cp(1)
         DT0(1) = DT0(1) + C2i * w(1) *(em(1)-cp(1)*tiold(1))
      ENDIF
      IF ( w(2) .GT. 0d0 ) THEN
         BT1(1) = BT1(1) + C2i * w(2) * cp(1)
         DT0(1) = DT0(1) - C2i * w(2) *(em(1)-cp(1)*tiold(1))
      ELSE
         CT1(1) = CT1(1) + C2i * w(2) * cp(2)
         DT0(1) = DT0(1) - C2i * w(2) *(em(2)-cp(2)*tiold(2))
      ENDIF

      IF (bbc .EQ. 'fixT') THEN
         BT1(0)	= 1d0
         CT1(0)	= 0d0
         DT0(0)	= tocn
      ELSE
         BT1(0) = ( 8d0*ki(0))
         CT1(0) = (-9d0*ki(0))
         OT1	= ( 1d0*ki(0))
         DT0(0) = ( C4i*Fbase )
         
         BT1(0) = OT1*AT1(1)-CT1(1)*BT1(0)
         CT1(0) = OT1*BT1(1)-CT1(1)*CT1(0)
         DT0(0) = OT1*DT0(1)-CT1(1)*DT0(0)
      ENDIF
      
      DO j=2,ni-2
         AT1(j) = -C1i*(ki(j)+ki(j-1))
         BT1(j) =  C1i*(ki(j+1)+2d0*ki(j)+ki(j-1)) &
                  + cp(j)*hi_b(1)*dzic
         CT1(j) = -C1i*(ki(j+1)+ki(j))
         DT0(j) =  C3i*R(j)+cp(j)*tiold(j)*hi_b(1)*dzic &
                  - em(j)*(hi_b(1)-hi_b(0))*dzic

         IF (w(j) .GT. 0d0) THEN
            AT1(j) = AT1(j) - C2i*w(j  ) * cp(j-1)
            DT0(j) = DT0(j) + C2i*w(j  ) *(em(j-1)-cp(j-1)*tiold(j-1))
         ELSE
            BT1(j) = BT1(j) - C2i*w(j  ) * cp(j  )
            DT0(j) = DT0(j) + C2i*w(j  ) *(em(j  )-cp(j  )*tiold(j  ))
         ENDIF
         IF (w(j+1) .GT. 0d0) THEN
            BT1(j) = BT1(j) + C2i*w(j+1) * cp(j  )
            DT0(j) = DT0(j) - C2i*w(j+1) *(em(j  )-cp(j  )*tiold(j  ))
         ELSE
            CT1(j) = CT1(j) + C2i*w(j+1) * cp(j+1)
            DT0(j) = DT0(j) - C2i*w(j+1) *(em(j+1)-cp(j+1)*tiold(j+1))
         ENDIF
      ENDDO

         AT1(ni-1) = -C1i/3d0*(3d0*ki(ni-2)+3d0*ki(ni-1)+2d0*ki(ni))
         BT1(ni-1) =  C1i/3d0*(3d0*ki(ni-2)+3d0*ki(ni-1)+18d0*ki(ni)) &
                  + cp(ni-1)*hi_b(1)*dzic
         CT1(ni-1) = -C1i/3d0*16d0*ki(ni)
         DT0(ni-1) =  C3i*R(ni-1)+cp(ni-1)*tiold(ni-1)*hi_b(1)*dzic &
                  - em(ni-1)*(hi_b(1)-hi_b(0))*dzic

      IF (w(ni-1) .GT. 0d0) THEN
         AT1(ni-1) = AT1(ni-1) - C2i * w(ni-1) * cp(ni-2)
         DT0(ni-1) = DT0(ni-1) + C2i * w(ni-1) *(em(ni-2)-cp(ni-2)*tiold(ni-2))
      ELSE
         BT1(ni-1) = BT1(ni-1) - C2i * w(ni-1) * cp(ni-1)
         DT0(ni-1) = DT0(ni-1) + C2i * w(ni-1) *(em(ni-1)-cp(ni-1)*tiold(ni-1))
      ENDIF
      IF (w(ni) .GT. 0d0) THEN
         BT1(ni-1) = BT1(ni-1) + C2i * w(ni  ) * cp(ni-1)
         DT0(ni-1) = DT0(ni-1) - C2i * w(ni  ) *(em(ni-1)-cp(ni-1)*tiold(ni-1))
      ELSE
         CT1(ni-1) = CT1(ni-1) + C2i * w(ni  ) * cp(ni  )
         DT0(ni-1) = DT0(ni-1) - C2i * w(ni  ) *(em(ni  )-cp(ni  )*tiold(ni  ))
      ENDIF

! FD     if (hs_b(0)>hscut) then
     if (hs_b(0)>hslim) then
      PT1	= (    -ki(ni)/3d0/hi_b(1)/dzic)
      AT1(ni)	= ( 9d0*ki(ni)/3d0/hi_b(1)/dzic)
!FD add case for nlsno=1
      if ( nlsno == 1 ) then
       BT1(ni)	= (-8d0*ki(ni)/3d0/hi_b(1)/dzic)+ &
                  (-2d0*ks(ni)/hs_b(1))
       CT1(ni)	= ( 2d0*ks(ni)/hs_b(1))
       QT1	= ( 0.d0)
      else 
       BT1(ni)	= (-8d0*ki(ni)/3d0/hi_b(1)/dzic)+ &
                  (-8d0*ks(ni)/3d0/hs_b(1)/dzsc)
       CT1(ni)	= ( 9d0*ks(ni)/3d0/hs_b(1)/dzsc)
       QT1	= (-ks(ni)/3d0/hs_b(1)/dzsc)
      endif

      DT0(ni)	= ( 0d0)

!FD add case for nlsno=1
      if ( nlsno == 1 ) then
       AT1(ni+1) = -C1s*(2d0*ks(ni))
       BT1(ni+1) =  C1s*(2d0*ks(ni)+2d0*ks(ns)) &
                  + cp(ni+1)*hs_b(1)*dzsc
       CT1(ni+1) = -C1s*(2d0*ks(ns))
       DT0(ni+1) =  C3s*R(ni+1)+cp(ni+1)*tiold(ni+1)*hs_b(1)*dzsc &
                  - em(ni+1)*(hs_b(1)-hs_b(0))*dzsc
      else
       AT1(ni+1) = -C1s/3d0*16d0*ks(ni)
       BT1(ni+1) =  C1s/3d0*(18d0*ks(ni)+3d0*ks(ni+1)+3d0*ks(ni+2)) &
                  + cp(ni+1)*hs_b(1)*dzsc
       CT1(ni+1) = -C1s/3d0*(3d0*ks(ni+2)+3d0*ks(ni+1)+2d0*ks(ni))
       DT0(ni+1) =  C3s*R(ni+1)+cp(ni+1)*tiold(ni+1)*hs_b(1)*dzsc &
                  - em(ni+1)*(hs_b(1)-hs_b(0))*dzsc
      endif

      IF (w(ni+1) .GT. 0d0) THEN
         BT1(ni+1) = BT1(ni+1) + C2s*w(ni+1) * cp(ni+1)
         DT0(ni+1) = DT0(ni+1) - C2s*w(ni+1) *(em(ni+1)-cp(ni+1)*tiold(ni+1))
      ELSE
         CT1(ni+1) = CT1(ni+1) + C2s*w(ni+1) * cp(ni+2)
         DT0(ni+1) = DT0(ni+1) - C2s*w(ni+1) *(em(ni+2)-cp(ni+2)*tiold(ni+2))
      ENDIF


      AT1(ni) = CT1(ni+1)*(AT1(ni)*AT1(ni-1)-PT1*BT1(ni-1))
      BT1(ni) = CT1(ni+1)*(BT1(ni)*AT1(ni-1)-PT1*CT1(ni-1))- &
                AT1(ni+1)*QT1*AT1(ni-1)
      CT1(ni) = CT1(ni+1)*CT1(ni)*AT1(ni-1)-BT1(ni+1)*QT1*AT1(ni-1)
      DT0(ni) = CT1(ni+1)*(DT0(ni)*AT1(ni-1)-PT1*DT0(ni-1))- &
                DT0(ni+1)*QT1*AT1(ni-1)

      DO j=ni+2,ns-2
            AT1(j) = -C1s*(ks(j)+ks(j-1))
            BT1(j) =  C1s*(ks(j+1)+2d0*ks(j)+ks(j-1))&
                     + cp(j)*hs_b(1)*dzsc
            CT1(j) = -C1s*(ks(j+1)+ks(j))
            DT0(j) =  C3s*R(j)+cp(j)*tiold(j)*hs_b(1)*dzsc &
                     - em(j)*(hs_b(1)-hs_b(0))*dzsc

         IF (w(j-1) .GT. 0d0) THEN
            AT1(j) = AT1(j) - C2s*w(j-1) * cp(j-1)
            DT0(j) = DT0(j) + C2s*w(j-1) *(em(j-1)-cp(j-1)*tiold(j-1))
         ELSE
            BT1(j) = BT1(j) - C2s*w(j-1) * cp(j  )
            DT0(j) = DT0(j) + C2s*w(j-1) *(em(j  )-cp(j  )*tiold(j  ))
         ENDIF
         IF (w(j) .GT. 0d0) THEN
            BT1(j) = BT1(j) + C2s*w(j  ) * cp(j  )
            DT0(j) = DT0(j) - C2s*w(j  ) *(em(j  )-cp(j  )*tiold(j  ))
         ELSE
            CT1(j) = CT1(j) + C2s*w(j  ) * cp(j+1)
            DT0(j) = DT0(j) - C2s*w(j  ) *(em(j+1)-cp(j+1)*tiold(j+1))
         ENDIF
      ENDDO

!FD add case for nlsno=1
      if ( nlsno > 1 ) then
         AT1(ns-1) = -C1s/3d0*(3d0*ks(ns-2)+3d0*ks(ns-1)+2d0*ks(ns))
         BT1(ns-1) =  C1s/3d0*(3d0*ks(ns-2)+3d0*ks(ns-1)+18d0*ks(ns)) &
                     +cp(ns-1)*hs_b(1)*dzsc
         CT1(ns-1) = -C1s/3d0*16d0*ks(ns)
         DT0(ns-1) =  C3s*R(ns-1)+cp(ns-1)*tiold(ns-1)*hs_b(1)*dzsc &
                     -em(ns-1)*(hs_b(1)-hs_b(0))*dzsc

      IF (w(ns-2) .GT. 0d0) THEN
         AT1(ns-1) = AT1(ns-1) - C2s*w(ns-2) * cp(ns-2)
         DT0(ns-1) = DT0(ns-1) + C2s*w(ns-2) *(em(ns-2)-cp(ns-2)*tiold(ns-2))
      ELSE
         BT1(ns-1) = BT1(ns-1) - C2s*w(ns-2) * cp(ns-1)
         DT0(ns-1) = DT0(ns-1) + C2s*w(ns-2) *(em(ns-1)-cp(ns-1)*tiold(ns-1))
      ENDIF
      IF (w(ns-1) .GT. 0d0) THEN
         BT1(ns-1) = BT1(ns-1) + C2s*w(ns-1) * cp(ns-1)
         DT0(ns-1) = DT0(ns-1) - C2s*w(ns-1) *(em(ns-1)-cp(ns-1)*tiold(ns-1))
      ELSE
         CT1(ns-1) = CT1(ns-1) + C2s*w(ns-1) * cp(ns  )
         DT0(ns-1) = DT0(ns-1) - C2s*w(ns-1) *(em(ns  )-cp(ns  )*tiold(ns  ))
      ENDIF
      endif ! end nlsno > 1

!FD add case for nlsno=1
      if ( nlsno == 1 ) then
       PT1	= ( 0d0)
       AT1(ns)	= ( 2d0*ks(ns))
       BT1(ns)	= (-2d0*ks(ns) &
                   -C4s *                     dzf )
       DT0(ns)	=  -C4s * ( Fnet0 + dzf*tiold(ns) )
      else
       PT1	= (-1d0*ks(ns))
       AT1(ns)	= ( 9d0*ks(ns))
       BT1(ns)	= (-8d0*ks(ns) &
                   -C4s *                     dzf )
       DT0(ns)	=  -C4s * ( Fnet0 + dzf*tiold(ns) )
      endif ! end nlsno = 1

      IF (sbc .EQ. 'flux' .AND. Tsbc) THEN
       AT1(ns)	= 0d0
       BT1(ns)	= 1d0
       DT0(ns)	= Tf(ns)
      ELSEIF (sbc .EQ. 'fixT') THEN
       AT1(ns)	= 0d0
       BT1(ns)	= 1d0
       DT0(ns)	= Tsurf
      ELSEIF (temp(ns,1) .LT. Tf(ns) .OR. &
              temp(ns,1)+tiny .GE. Tf(ns) .AND.  &
             -Fnet-Fcss .GT. 0d0        ) THEN
       AT1(ns)	= PT1*BT1(ns-1)-AT1(ns)*AT1(ns-1)
       BT1(ns)	= PT1*CT1(ns-1)-BT1(ns)*AT1(ns-1)
       DT0(ns)	= PT1*DT0(ns-1)-DT0(ns)*AT1(ns-1)
      ENDIF

! FD thin snow
!     else if ( hs_b(0) <=hscut .and. hs_b(0) > hslim ) then ! FD thin snow
!
!      PT1     = -ki(ni)/3d0/hi_b(1)/hi_b(1)/dzic
!      AT1(ni) =  9d0*ki(ni)/3d0/hi_b(1)/hi_b(1)/dzic
!      BT1(ni) = -8d0*ki(ni)/3d0/hi_b(1)/hi_b(1)/dzic &
!                -ks(ni)/hs_b(1)/hs_b(1) * hs_b(1)/hi_b(1)
!      CT1(ni) =  ks(ni)/hs_b(1)/hs_b(1)
!      DT0(ni) =  0d0
!
!      AT1(ni) = AT1(ni)*AT1(ni-1)-PT1*BT1(ni-1)
!      BT1(ni) = BT1(ni)*AT1(ni-1)-PT1*CT1(ni-1)
!      CT1(ni) = CT1(ni)*AT1(ni-1)
!      DT0(ni) = DT0(ni)*AT1(ni-1)-PT1*DT0(ni-1)
!
!       AT1(ni+1) =  ks(ns)
!       BT1(ni+1) = -ks(ns) &
!                   -C4s *                     dzf
!       DT0(ni+1) = -C4s * ( Fnet0 + dzf*tiold(ns) )
!
!      IF (sbc .EQ. 'flux' .AND. Tsbc) THEN
!         AT1(ni+1) = 0d0
!         BT1(ni+1) = 1d0
!         DT0(ni+1) = Tf(ns)
!      ELSEIF (sbc .EQ. 'fixT') THEN
!         AT1(ni+1) = 0d0
!         BT1(ni+1) = 1d0
!         DT0(ni+1) = Tsurf
!      ELSEIF (temp(ns,1) .LT. Tf(ns) .OR. &
!              temp(ns,1)+tiny .GE. Tf(ns) .AND.  &
!              -Fnet-Fcss .GT. 0d0        ) THEN
!         AT1(ni+1) = AT1(ni+1) ! FD seems useless!
!         BT1(ni+1) = BT1(ni+1)
!         DT0(ni+1) = DT0(ni+1)
!      ENDIF

! FD no snow
     else if ( hs_b(0) <= hslim ) then ! FD no snow

         PT1	= (-1d0*ki(ni))
         AT1(ni)= ( 9d0*ki(ni))
         BT1(ni)= (-8d0*ki(ni) &
                   -C4i *                     dzf )
         DT0(ni)=  -C4i * ( Fnet0 + dzf*tiold(ni) )

      IF (sbc .EQ. 'flux' .AND. Tsbc) THEN
        AT1(ni)	= 0d0
        BT1(ni)	= 1d0
        DT0(ni)	= Tf(ni)
      ELSEIF (sbc .EQ. 'fixT') THEN
        AT1(ni)	= 0d0
        BT1(ni)	= 1d0
        DT0(ni)	= Tsurf
      ELSEIF (temp(ni,1) .LT. Tf(ni) .OR. &
              temp(ni,1)+tiny .GE. Tf(ni) .AND.  &
              -Fnet-Fcis .GT. 0d0        ) THEN
        AT1(ni)	= PT1*BT1(ni-1)-AT1(ni)*AT1(ni-1)
        BT1(ni)	= PT1*CT1(ni-1)-BT1(ni)*AT1(ni-1)
        DT0(ni)	= PT1*DT0(ni-1)-DT0(ni)*AT1(ni-1)
      ENDIF
     endif

      

!-----------------------------------------------------------------------
!     Solve for the internal temperature
!-----------------------------------------------------------------------

! FD     if (hs_b(0)>hscut) then
     if (hs_b(0)>hslim) then
      CALL tridag (AT1,BT1,CT1,DT0,tout,ns+1)
!     else if ( hs_b(0) <=hscut .and. hs_b(0) > hslim ) then ! FD thin snow
!      CALL tridag (AT1,BT1,CT1,DT0,tout,ni+2)
     else if ( hs_b(0) <= hslim ) then ! FD no snow
      CALL tridag (AT1,BT1,CT1,DT0,tout,ni+1)
     endif

      DO j=0,ni
         temp(j,1) = tout(j)
      ENDDO
!FD is snow layers
! FD     if (hs_b(0)>hscut) then
     if (hs_b(0)>hslim) then
      DO j=ni+1,ns
         temp(j,1) = tout(j)
      ENDDO
!     else if ( hs_b(0) <=hscut .and. hs_b(0) > hslim ) then ! FD thin snow
!         j=ni+1
!         temp(j,1) = tout(j)
!      DO j=ni+2,ns
!         temp(j,1) = temp(ni+1,1)
!      ENDDO
     else if ( hs_b(0) <= hslim ) then ! FD no snow
      DO j=ni+1,ns
         temp(j,1) = temp(ni,1)
      ENDDO
     endif

! FD debug
if (extra_debug) then
 write(*,'(I6,300(1x,e9.3))') counter, w(1:ns-1)
 write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns,1)
endif

!FD is snow layers
! FD     if (hs_b(0)>hscut) then
     if (hs_b(0)>hslim) then
      IF ((temp(ns,1) .GT. Tf(ns)) .AND. (Tsbc .EQV. .FALSE.)) THEN
         Tsbc = .TRUE. 
      ENDIF
      IF ( temp(ns,1)+tiny .GE. Tf(ns) .AND. -Fnet-Fcss .GT. 0d0 ) THEN
         Tsbc = .FALSE.
         temp(ns,1)=Tf(ns)-1d-2
      ENDIF
!     else if ( hs_b(0) <=hscut .and. hs_b(0) > hslim ) then ! FD thin snow
!      IF ((temp(ni+1,1) .GT. Tf(ns)) .AND. (Tsbc .EQV. .FALSE.)) THEN
!         Tsbc = .TRUE. 
!      ENDIF
!      IF ( temp(ni+1,1)+tiny .GE. Tf(ns) .AND. -Fnet-Fcss .GT. 0d0 ) THEN
!         Tsbc = .FALSE.
!         temp(ni+1,1)=Tf(ns)-1d-2
!      ENDIF
     else if ( hs_b(0) <= hslim ) then ! FD no snow
      IF ((temp(ni,1) .GT. Tf(ni)) .AND. (Tsbc .EQV. .FALSE.)) THEN
         Tsbc = .TRUE. 
      ENDIF
      IF ( temp(ni,1)+tiny .GE. Tf(ni) .AND. -Fnet-Fcss .GT. 0d0 ) THEN
         Tsbc = .FALSE.
         temp(ni,1)=Tf(ni)-1d-2
      ENDIF
     endif

         DO j=0,ns
            IF (temp(j,1) .GT. Tf(j)) THEN
               temp(j,1) = Tf(j)
            ENDIF
         ENDDO

!-----------------------------------------------------------------------
!     Update T-S dependent ice parameters
!     Lf, cp are diagnostics for output
!-----------------------------------------------------------------------
      
      DO j=0,ni
         Tf(j)   = Tfreeze1(sali(j,1))
         ki(j)   = func_ki(sali(j,1),temp(j,1))
         qm(j)   = func_qm(Tf(j),temp(j,1))
         cp(j)   = func_cp(Tf(j),temp(j,1),tiold(j))
         sh(j,1) = func_sh(Tf(j),temp(j,1))
         lh(j,1) = func_lh(Tf(j),temp(j,1))*hi_b(1)*dzic
      ENDDO
      DO j=ni+1,ns
         Tf(j)   = Tfreeze1(sali(j,1))
         qm(j)   = func_qm(Tf(j),temp(j,1))
         cp(j)   = func_cp(Tf(j),temp(j,1),tiold(j))
         sh(j,1) = func_sh(Tf(j),temp(j,1))
         lh(j,1) = func_lh(Tf(j),temp(j,1))*hs_b(1)*dzsc
      ENDDO

!-----------------------------------------------------------------------

!     Update conductive heat fluxes
!-----------------------------------------------------------------------
      
      Fcis	= -ki(ni) &
                *(8d0*temp(ni,1)-9d0*temp(ni-1,1)+temp(ni-2,1)) &
                /3d0/dzic/hi_b(1)
      
      Fcib	= -ki(0) &
                *(-8d0*temp(0,1)+9d0*temp(1,1)-temp(2,1)) &
                /3d0/dzic/hi_b(1)
      
!FD is snow layers
! FD     if (hs_b(0)>hscut) then
     if (hs_b(0)>hslim) then
!FD add case for nlsno=1
      if ( nlsno == 1 ) then
       Fcss	= -ks(ns) &
                  *( 2d0*temp(ns,1)-2d0*temp(ns-1,1))/hs_b(1)

       Fcsb	= -ks(ni) &
                  *(-2d0*temp(ni,1)+2d0*temp(ni+1,1))/hs_b(1)
      else
       Fcss	= -ks(ns) &
                  *(8d0*temp(ns,1)-9d0*temp(ns-1,1)+temp(ns-2,1)) &
                  /3d0/dzsc/hs_b(1)

       Fcsb	= -ks(ni) &
                  *(-8d0*temp(ni,1)+9d0*temp(ni+1,1) -temp(ni+2,1)) &
                  /3d0/dzsc/hs_b(1)
      endif
!     else if ( hs_b(0) <=hscut .and. hs_b(0) > 0.d0) then ! FD thin snow
!      Fcss	= -ks(ns) &
!                  *(temp(ns,1)-temp(ni,1)) &
!                  /hs_b(1)
!
!      Fcsb	= -ks(ni) &
!                  *(temp(ns,1)-temp(ni,1)) &
!                  /hs_b(1)
     endif

!-----------------------------------------------------------------------
!     Update atm and ocn heat fluxes
!-----------------------------------------------------------------------
      
! FD assume only snow
! FD     if ( hs_b(0) > hscut ) then
     if ( hs_b(0) > hslim ) then
            Fprec = snowfall*rhosno*(-mlfus+cp_ice*temp(ns,1))
     else
            Fprec = 0.d0
     endif
            Frain = 0d0
     if ( hs_b(0) > hslim ) then
            Fnet = Fnet0 - dzf * (temp(ns,1)-tiold(ns))
     else ! FD no snow
            Fnet = Fnet0 - dzf * (temp(ni,1)-tiold(ni))
     endif

      IF (hs_b(0) > hslim) THEN
        tsu = min(tsu,tf(ns)+temp0)
      ELSE
        tsu = min(tsu,tf(ni)+temp0)
      ENDIF

      Fmlt = Fnet-Fcss

! FD         oceflx= rhoo*cpo*Coi*ABS(uio)*(tocn-temp(0,1))

!-----------------------------------------------------------------------
!     Euler step
!-----------------------------------------------------------------------
      
      DO j=0,ns
         dTlay = ABS(temp(j,1)-temp(j,0))	! T diff in a layer
         dTmax = MAX(dTmax,dTlay)		! max diff in s/ice slab
      ENDDO
      IF (dTmax .GT. Tdiff .and. counter < 50 ) THEN
         counter=counter+1
         DO j=0,ns
            temp(j,0) = temp(j,1)		! update temperature
         ENDDO
         dTmax = 0d0
         GOTO 4000				! recalculate T profile
      ENDIF

!-----------------------------------------------------------------------
!     Compute internal energy
!-----------------------------------------------------------------------

      Frad=0.d0
!      IF ( hs_b(0) > hscut ) THEN
      IF ( hs_b(0) > hslim ) THEN
       DO j=ni+1,ns-1
         Frad	= Frad + R(j) * dzsc * hs_b(1) 
       ENDDO
      ENDIF
      DO j=1,ni
         Frad	= Frad + R(j) * dzic * hi_b(1) 
      ENDDO
      
      Einp = (Fnet+oceflx+Frad+Fprec)*dtice	! E input [J/m2]

      
      Eint = 0d0

      DO j=1,ni-1
         Elay = func_El(Tf(j),temp(j,1))*hi_b(1)*rhoice*dzic
         Eint = Eint+Elay
      ENDDO
! FD debug     if (hs_b(0)>hscut) then
     if (hs_b(0) > hslim) then
       DO j=ni+1,ns-1
         Elay = func_El(Tf(j),temp(j,1))*hs_b(1)*rhosno*dzsc
         Eint = Eint+Elay
      ENDDO
     endif
     
      dEin = Eint-Ein0		! diff intl E [J/m2]
! FD debug
if (debug) write(*,*) 'energy',Einp,dEin

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
 write(*,*) 'snow',hs_b,-w(ns-1)
 write(*,*) 'hice',hi_b,-w(ni),w(1)
 write(*,'(I6,300(1x,e9.3))') counter, w(1:ns-1)
 write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns,1)
endif
! FD debug if energy not conserved
!      IF ( abs(Einp-dEin) > 1e-3 ) THEN
!       OPEN(1,file='restart.dat')
!       WRITE(1,*) nlice,nlsno,ith_cond
!       WRITE(1,*) ti(1:nlice),ts(1:nlsno+1),tsu,tbo
!       WRITE(1,*) si(1:nlice)
!       WRITE(1,*) hi,hs
!       WRITE(1,*) snowfall,dwnlw,tsu,tair,qair,uair,swrad,oceflx
!       WRITE(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
!       CLOSE(1)
!       STOP 'energy not conserved'
!      ENDIF

!stop
      hs_b(0) = hs_b(1)
      hi_b(0) = hi_b(1)
      
      DO j=0,ns
         sali(j,0) = sali(j,1)
         temp(j,0) = temp(j,1)
      ENDDO
      
      IF (hi_b(0) .LT. hicut) THEN
         hi_b(0) = hicut
         DO j=0,ni
            temp(j,0) = tocn
         ENDDO
      ENDIF
      

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
       WRITE(*,*) ' dh_snowice : ', dh_sni
       WRITE(*,*) ' ht_s_b : ', hs
       WRITE(*,*) ' ht_i_b : ', hi
      endif

! conversion from LIM3
      tbo=temp(0,1) + temp0
      do j=1,nlice
        ti(nlice-j+1)=temp(j,1) + temp0
      enddo
      ts(nlsno+1)=temp(ni,1) + temp0
      do j=1,nlsno
        ts(nlsno-j+1)=temp(j+ni,1) + temp0
      enddo
      tsu=temp(ns,1) + temp0

! call back for LIM3 of thickness (needed for radiative transfer)
      do j=1,nlice
        dzi(nlice-j+1) = dzic * hi
      enddo
      do j=1,nlsno
        dzs(nlsno-j+1) = dzsc * hs
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


