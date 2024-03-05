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
double precision, dimension(maxlay) :: dzzi,dzzs


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

      ni=   nlice+1
      ns=ni+nlsno+1 
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

 do k=1,ni-1
    zi(k) = dble(k-1)/dble(nlice)
 enddo
   zi(ni) = 1.d0
 do k=ni+1,ns-2
    zi(k) = 1.d0 + dble(k-ni)/dble(nlsno)
 enddo
   zi(ns-1) = 2.d0
   zs(ni) = 0.d0
 do k=ni+1,ns-2
    zs(k) = dble(k-ni)/dble(nlsno)
 enddo
   zs(ns-1) = 1.d0

 do k=1,ni-1
    dzi(k) = 1.d0/dble(nlice)
    dzzi(k) = 1.d0/dble(nlice)
 enddo
 dzzi(1 ) = 0.5d0/dble(nlice)
 dzzi(ni) = 0.5d0/dble(nlice)
 dzzs(ni) = 0.5d0/dble(nlsno)
 do k=ni+1,ns-1
    dzi(k) = 1.d0/dble(nlsno)
    dzzs(k) = 1.d0/dble(nlsno)
 enddo
 dzzs(ns-1) = 0.5d0/dble(nlsno)

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
     epsilon,&! coefficient of emissivity			[-]
     hicut,&! cut off (minimum) ice  thickness		[m]
     hscut,&! cut off (minimum) snow thickness		[m]
     hslim  ! minimum snow thickness			[m]

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
    kki(0:maxlay)    ,&! ice  thermal conductivity			[W/m/C]
    kks(0:maxlay)    ,&! ice  thermal conductivity			[W/m/C]
    em(0:maxlay)    ,&! volumetric enthalpiy	[W/m3]
    qm(0:maxlay)    ,&! volumetric energy of melt (rho*Lf(S,T))	[W/m3]
    cp(0:maxlay)    ,&! heat capacity for ice	[W/kg/K^-1]
    R(0:maxlay)    ,&! penetrating shortwave radiation		[W/m3]
    sali(0:maxlay,0:1),&! internal sea ice salinity			[psu]
    Tf(0:maxlay)    ,&! ice freezing temp. of salinity S		[C]
    w(1:maxlay-1)  ,&! grid advection velocity			[1/s]
    hi_b(0:1)     ,&! ice thickness (m)
    hs_b(0:1)       ! snow thickness (m)
  double precision, dimension(maxlay) :: &
             he, heold, henew

      double precision &
                temp (0:maxlay,0:1) ! internal sea ice temperature		[C]

      DOUBLE PRECISION AT1(0:maxlay), BT1(0:maxlay), CT1(0:maxlay), DT0(0:maxlay), &
      C2i, C3i, C2s, C3s, &
      tout(0:maxlay), hsmemo, Rtrans

      double precision :: sice, &
             latmelt, ftot, sstznew, &
             Tside,deltaT,wlat,rside, &
             sum0, dhi, dho, dq, fwf0
      double precision :: hminice=1d-3
      double precision :: energy_bot

  character (len=20) :: Sprofile
  integer i,j,counter
  double precision dzf, es, zssdqw, zrchu1, zrchu2, q0, zref, k0, k1
  double precision, dimension(0:maxlay) :: tiold, tinew
  double precision elays0
  logical :: debug=.false.
  logical :: extra_debug=.false.
  double precision :: hsold, tsuold, sfallold
! FD debug
 debug=.true.
! extra_debug=.true.

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

      epsilon=   0.990d+0! snow/ice emissivity			[-]

      sigma=   5.670d-8! Stefan-Boltzmann constant	       [W/m2/K4]
      temp0= 273.160d+0! freezing point temp of fresh water	[K]

      uiofix=   0.000d+0! constant relative ice-ocean velocity	[m/s]



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

      Tdiff	=   1.000d-8	! temperature tolerance in Euler step	[C]
      tiny	=   1.000d-9	! very small number			[-]

      hicut	=   0.010d+0	! minimum ice  thickness		[m]
      hslim	=   0.005d+0	! minimum snow thickness		[m]
! FD debug
hslim=0.0005d0
Tdiff=1d-12

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
      tiold(ni)=ts(nlsno+1)  - temp0
      do j=1,nlsno
        tiold(j+ni)=ts(nlsno-j+1)  - temp0
      enddo
      Tsurf = tsu - temp0
      tiold(ns) = Tsurf

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
! FD debug
!      sali(ni,0)= salinis
      sali(ni,0)= sali(ni-1,0)
         DO j=ni+1,ns
            sali(j,0) = 0d0
         ENDDO

      sali(:,1)=sali(:,0) ! FD not sure what is really done in this code, just for debug

      hsold = hs
      tsuold = tsu
      sfallold = snowfall

1000 continue
!------------------------------------------------------------------------
! init variable
!------------------------------------------------------------------------

      hs_b(0:1)=hs
      hi_b(0:1)=hi
      heold(1:ni-1)=hi*dzi(1:ni-1)
      heold(ni)=0.d0
      heold(ni+1:ns-1)=hs*dzi(ni+1:ns-1)
      he=heold
      henew=he


      tocn  = tiold(0)
      temp(:,0)=tiold
      temp(:,1)=tiold

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
            ki(j)   = func_ki(sali(j,0),tiold(j))
            qm(j)   = func_qm(Tf(j),tiold(j))
            em(j)   = func_el(Tf(j),tiold(j))
            cp(j)   = func_cp(Tf(j),tiold(j),tiold(j))
      ENDDO

      DO j=0,ns					! snow thermal conductivity
         IF (j .LT. ni) THEN
            ks(j) = 0d0
         ELSEIF (j .GE. ni) THEN
            ks(j) = cond_sno
         ENDIF
      ENDDO
      kki(1)=ki(0)/ (hi_b(0) * dzzi(1))
      DO j=2,ni-1
         kki(j) = (ki(j-1)+ki(j)) * 0.5d0 / (hi_b(0) * dzzi(j))
      ENDDO
      kki(ni) = ki(ni) / (hi_b(0) * dzzi(ni))
      kks(ni) = ks(ni) / (hs_b(0) * dzzs(ni))
      DO j=ni+1,ns-2
         kks(j) = (ks(j+1)+ks(j)) * 0.5d0 / (hs_b(0) * dzzs(j))
      ENDDO
      kks(ns-1) = ks(ns) / (hs_b(0) * dzzs(ns-1))

!------------------------------------------------------------------------
!     Initial internal energy in ice (1:ni-1) and snow (ni+1:ns-1)
!------------------------------------------------------------------------

      Einp = 0d0				! E input [J/m2]
      Ein0 = 0d0

      DO j=1,ni-1
         Elay = em(j)*rhoice*heold(j)
         Ein0 = Ein0+Elay
      ENDDO

     IF (hs_b(0)>hslim) then
      DO j=ni+1,ns-1
         Elay = em(j)*rhosno*heold(j)
         Ein0 = Ein0+Elay
      ENDDO
     ENDIF

      dEin = 0d0				! diff intl E [J/m2]

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
!     Update conductive heat fluxes
!------------------------------------------------------------------------

     IF (hs_b(0) .GT. hslim) THEN
         Fcss = -kks(ns-1) * ( temp(ns  ,1) - temp(ns-1,1) )
         Fcsb = -kks(ni  ) * ( temp(ni+1,1) - temp(ni  ,1) )
      ELSE
         Fcss = 0d0
         Fcsb = 0d0
      ENDIF
      
      Fcis = -kki(ni) * ( temp(ni,1) - temp(ni-1,1))
      Fcib = -kki(1 ) * ( temp( 1,1) - temp(   0,0))

!-----------------------------------------------------------------------
!     ICE/SNOW case
!-----------------------------------------------------------------------


      Tsbc = .FALSE.
      counter = 0

!-----------------------------------------------------------------------
!     Euler loop
!-----------------------------------------------------------------------

4000  CONTINUE

      Fmlt = Fnet-Fcis

!-----------------------------------------------------------------------
!     Update atm and ocn heat fluxes
!-----------------------------------------------------------------------
      
     if ( hs_b(0) > hslim ) then
            Fprec = -snowfall*rhosno*qm(ns)
     else
            Fprec = 0.d0
     endif
            Frain = 0d0

     if ( hs_b(0) > hslim ) then
            Fnet = Fnet0 - dzf * (temp(ns,1)-tiold(ns))
     else ! FD no snow
            Fnet = Fnet0 - dzf * (temp(ni,1)-tiold(ni))
     endif
! FD debug
if (extra_debug) then
write(*,*) 'Fnet',Fnet, Fnet0, - dzf * (temp(ni,1)-tiold(ni))
endif

      IF (hs_b(0) > hslim) THEN
        tsu = min(tsu,tf(ns)+temp0)
      ELSE
        tsu = min(tsu,tf(ni)+temp0)
      ENDIF

! FD         oceflx= rhoo*cpo*Coi*ABS(uio)*(tocn-temp(0,1))
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
! FD debug
if (extra_debug) then
write(*,*) 'snow surfacae melt if <0',dssdt,temp(ns,1)+tiny .GE. Tf(ns),Fnet,Fcss
endif
      else

      ! update ice surface

         IF (temp(ni,1)+tiny .GE. Tf(ni) .AND. -Fnet-Fcis .LT. 0d0) THEN
           disdt = MIN((-Fnet-Fcis)/rhoice/qm(ni-1),0d0)
         ENDIF
         IF (counter>10.and.w(ni)>1e-7) THEN ! case of a severe limit cycle reached
           disdt = 0.5d0 * ( disdt - w(ni) ) ! apply a relaxation scheme
         ENDIF
! FD debug
if (extra_debug) then
write(*,*) 'ice surfacae melt if <0',disdt,temp(ni,1)+tiny .GE. Tf(ni),Fnet,Fcis
endif
      endif

      IF (oceflx-Fcib .LT. 0d0) THEN		! for ice growth: new-ice salin.
         qm(0) = func_qm(Tf(0),temp(0,1))
         energy_bot = qm(0)
      ELSEIF (oceflx-Fcib .GT. 0d0) THEN	! qm_mean_first_layer for melt
         energy_bot = qm(1)
      ENDIF
! FD generalized bottom
         dibdt	= (oceflx-Fcib) / rhoice / energy_bot

      dhsdt	= dssdt - dsbdt
      dhidt	= disdt - dibdt ! original

      hs_b(1)	= hs_b(0) + dhsdt*dtice

! FD need to do something when hs reaches the limit
      if ( hs_b(1) < hslim .and. hs_b(0) > hslim) then
         hs=0.d0
         goto 1000
      endif

      hi_b(1)	= hi_b(0) + dhidt*dtice

      henew(1:ni-1)=dzi(1:ni-1)*hi_b(1)
      henew(ni+1:ns-1)=dzi(ni+1:ns-1)*hs_b(1)

! FD      fwf = fwf0 + ( dhsdt*rhos + dhidt*rhoi ) / rhoo

! FD only relevant when salinity varies during the integration
! FD case of solving for varying salinity
      DO j=0,ns
         Tf(j)   = Tfreeze1(sali(j,1))
      ENDDO

!-----------------------------------------------------------------------
!     Vertical velocities (grid advection)
!-----------------------------------------------------------------------

      DO j=1,ni
         w(j)	= -zi(j)*disdt-(1d0-zi(j))*dibdt
      ENDDO
     if (hs_b(0) > hslim ) then
      DO j=ni+1,ns-1
         w(j)	= -(zi(j)-1.d0)*dssdt-(2d0-zi(j))*dsbdt
      ENDDO
     endif

!-----------------------------------------------------------------------
!     Define constant coefficients
!-----------------------------------------------------------------------

      C2i	= dtice*rhoice
      C3i	= dtice

!     if (hs_b(0) > hscut ) then
     if (hs_b(0) > hslim ) then
      C2s	= dtice*rhosno
      C3s	= dtice
     endif

!-----------------------------------------------------------------------
!     Absorbed SW radiation energy per unit volume (R)       
!     1/dz Int_z^z+dz R dz = Rtrans(z+dz) - Rtrans(z)
!-----------------------------------------------------------------------

      if (hs_b(1) > hslim ) then
        DO j=ns-1,ni+1,-1
          R(j)   = swradab_s(ns-j)
        ENDDO
      else
         R(ni+1:ns-1)=0.d0
      endif
      R(ni)=0.d0
      DO j=ni-1,1,-1
         R(j)   = swradab_i(ni-j)
      ENDDO
      R(0)=0.d0

!-----------------------------------------------------------------------
!     Tri-diagonal matrix coefficients [A] {x} = {D}
!-----------------------------------------------------------------------

3000  CONTINUE

         AT1=0.d0; BT1=0.d0; CT1=0.d0; DT0=0.d0;


      IF (bbc .EQ. 'fixT') THEN
         BT1(0)	= 1d0
         CT1(0)	= 0d0
         DT0(0)	= tocn
      ELSE
         k0 = kki(1)
         BT1(0) =   k0
         CT1(0) = - k0
         DT0(0) = Fbase
      ENDIF

      DO j=1,ni-1
         k0 = dtice * kki(j  )
         k1 = dtice * kki(j+1)
         AT1(j) = -k0
         BT1(j) =  k0 + k1 + rhoice*cp(j)*henew(j)
         CT1(j) = -k1
         DT0(j) =  C3i*R(j)+rhoice*cp(j)*tiold(j)*henew(j) &
                          - rhoice*em(j)*(henew(j)-heold(j))  ! metric term to close the energy conservation

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

     if (hs_b(0)>hslim) then
       k0 = kki(ni)
       k1 = kks(ni)
       AT1(ni)	=   k0
       BT1(ni)	= - k0 - k1
       CT1(ni)	=   k1

       DT0(ni)	= ( 0d0)

      DO j=ni+1,ns-1
       k0 = dtice * kks(j-1)
       k1 = dtice * kks(j)
            AT1(j) = -k0
            BT1(j) =  k0 + k1 + rhosno*cp(j)*henew(j)
            CT1(j) = -k1
            DT0(j) =  C3s*R(j)+rhosno*cp(j)*tiold(j)*henew(j) &
                             - rhosno*em(j)*(henew(j)-heold(j))

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

       k1 = kks(ns-1)
       AT1(ns)	=   k1
       BT1(ns)	= - k1 - dzf
       DT0(ns)	=  -( Fnet0 + dzf*tiold(ns) )

      IF (sbc .EQ. 'flux' .AND. Tsbc) THEN
       AT1(ns)	= 0d0
       BT1(ns)	= 1d0
       DT0(ns)	= Tf(ns)
      ELSEIF (sbc .EQ. 'fixT') THEN
       AT1(ns)	= 0d0
       BT1(ns)	= 1d0
       DT0(ns)	= Tsurf
      ENDIF

! FD no snow
     else if ( hs_b(0) <= hslim ) then ! FD no snow

       k1 = kki(ni)
         AT1(ni) =   k1
         BT1(ni) = - k1 - dzf
         DT0(ni) =  -( Fnet0 + dzf*tiold(ni) )

      IF (sbc .EQ. 'flux' .AND. Tsbc) THEN
        AT1(ni)	= 0d0
        BT1(ni)	= 1d0
        DT0(ni)	= Tf(ni)
      ELSEIF (sbc .EQ. 'fixT') THEN
        AT1(ni)	= 0d0
        BT1(ni)	= 1d0
        DT0(ni)	= Tsurf
      ENDIF
     endif ! condition on hs_b

      

!-----------------------------------------------------------------------
!     Solve for the internal temperature
!-----------------------------------------------------------------------

     if (hs_b(0)>hslim) then
      CALL tridag (AT1,BT1,CT1,DT0,tout,ns+1)
     else if ( hs_b(0) <= hslim ) then ! FD no snow
      CALL tridag (AT1,BT1,CT1,DT0,tout,ni+1)
     endif

      DO j=0,ni
         temp(j,1) = tout(j)
      ENDDO
     if (hs_b(0)>hslim) then
      DO j=ni+1,ns
         temp(j,1) = tout(j)
      ENDDO
     else if ( hs_b(0) <= hslim ) then ! FD no snow
      DO j=ni+1,ns
         temp(j,1) = temp(ni,1)
      ENDDO
     endif

         DO j=0,ns
            IF (temp(j,1) .GT. Tf(j)) THEN
               temp(j,1) = Tf(j)
            ENDIF
         ENDDO

!------------------------------------------------------------------------
!     Update conductive heat fluxes
!------------------------------------------------------------------------

     IF (hs_b(0) .GT. hslim) THEN
         Fcss = -kks(ns-1) * ( temp(ns  ,1) - temp(ns-1,1) )
         Fcsb = -kks(ni  ) * ( temp(ni+1,1) - temp(ni  ,1) )
      ELSE
         Fcss = 0d0
         Fcsb = 0d0
      ENDIF
      
      Fcis = -kki(ni) * ( temp(ni,1) - temp(ni-1,1))
      Fcib = -kki(1 ) * ( temp( 1,1) - temp(   0,1))

     if ( hs_b(0) > hslim ) then
            Fnet = Fnet0 - dzf * (temp(ns,1)-tiold(ns))
     else ! FD no snow
            Fnet = Fnet0 - dzf * (temp(ni,1)-tiold(ni))
     endif

! FD debug
if (extra_debug) then
 write(*,'(I6,300(1x,e9.3))') counter, w(1:ns-1)
 write(*,'(I6,300(f8.3,1x))') counter, temp(0:ns,1)
endif

     if (hs_b(0)>hslim) then
! FD debug
if (extra_debug) then
write(*,*) 'condition on snow melt', temp(ns,1)+tiny .GE. Tf(ns), Tsbc,-Fnet-Fcss.gt.0.d0
endif
      IF (( temp(ns,1)+tiny .GE. Tf(ns)) .AND. (Tsbc .EQV. .FALSE.)) THEN
         Tsbc = .TRUE. 
      ENDIF
      IF ( temp(ns,1)+tiny .GE. Tf(ns) .AND. -Fnet-Fcss .GT. 0d0 ) THEN
         Tsbc = .FALSE.
         temp(ns,1)=Tf(ns)-1d-2
     ENDIF
     else if ( hs_b(0) <= hslim ) then ! FD no snow
! FD debug
if (extra_debug) then
write(*,*) 'condition on ice melt', (temp(ni,1)+tiny .GE. Tf(ni)), (Tsbc .EQV. .FALSE.), &
 ((temp(ni,1)+tiny .GE. Tf(ni)) .AND. (Tsbc .EQV. .FALSE.))
endif
      IF ((temp(ni,1)+tiny .GE. Tf(ni)) .AND. (Tsbc .EQV. .FALSE.)) THEN
         Tsbc = .TRUE. 
      ENDIF
! FD debug
if (extra_debug) then
write(*,*) 'condition on ice refreeze',temp(ni,1)+tiny .GE. Tf(ni) , -Fnet-Fcis .GT. 0d0, &
( temp(ni,1)+tiny .GE. Tf(ni) .AND. -Fnet-Fcis .GT. 0d0 .and. w(ni)==0.d0 )
endif
      IF ( temp(ni,1)+tiny .GE. Tf(ni) .AND. -Fnet-Fcis .GT. 0d0 .and. w(ni)==0.d0) THEN
         Tsbc = .FALSE.
         temp(ni,1)=Tf(ni)-1d-2
      ENDIF
     endif

! FD debug
if (extra_debug) then
 if (hs_b(0)>hslim) then
  write(*,*) 'SBCond snow',Tsbc,Fnet,Fcss
 else
  write(*,*) 'SBCond ice',Tsbc,Fnet,Fcis
endif
endif

!-----------------------------------------------------------------------
!     Update T-S dependent ice parameters
!     Lf, cp are diagnostics for output
!-----------------------------------------------------------------------
      
      DO j=0,ni
         ki(j)   = func_ki(sali(j,1),temp(j,1))
      ENDDO
      DO j=0,ns
         Tf(j)   = Tfreeze1(sali(j,1))
         qm(j)   = func_qm(Tf(j),temp(j,1))
         cp(j)   = func_cp(Tf(j),temp(j,1),tiold(j))
      ENDDO

!-----------------------------------------------------------------------
!     Euler step
!-----------------------------------------------------------------------
      
      DO j=0,ns
         dTlay = ABS(temp(j,1)-temp(j,0))	! T diff in a layer
         dTmax = MAX(dTmax,dTlay)		! max diff in s/ice slab
      ENDDO
      IF (dTmax .GT. Tdiff .and. counter < 100 ) THEN
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
      IF ( hs_b(0) > hslim ) THEN
       DO j=ni+1,ns-1
         Frad	= Frad + R(j)
       ENDDO
      ENDIF
      DO j=1,ni
         Frad	= Frad + R(j)
      ENDDO
      
      Einp = (Fnet+oceflx+Frad+Fprec)*dtice	! E input [J/m2]

      
      Eint = 0d0

      DO j=1,ni-1
         Elay = func_El(Tf(j),temp(j,1))*henew(j)*rhoice
         Eint = Eint+Elay
      ENDDO
     if (hs_b(0) > hslim) then
       DO j=ni+1,ns-1
         Elay = func_El(Tf(j),temp(j,1))*henew(j)*rhosno
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
      IF ( abs(Einp-dEin) > 1e-3 ) THEN
       OPEN(1,file='restart.dat')
       WRITE(1,*) nlice,nlsno,ith_cond
       WRITE(1,*) ti(1:nlice),ts(1:nlsno+1),tsuold,tbo
       WRITE(1,*) si(1:nlice)
       WRITE(1,*) hi,hsold
       WRITE(1,*) sfallold,dwnlw,tsuold,tair,qair,uair,swrad,oceflx,pres
       WRITE(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
       CLOSE(1)
       STOP 'energy not conserved'
      ENDIF

! FD debug
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


