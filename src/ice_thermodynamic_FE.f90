module ice_thermodynamic_FE

  use var_thermo_vertical

implicit none

contains

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
subroutine init_sigma_ice_FE
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 implicit none
! locals
 integer k
 double precision, dimension(maxlay) :: qm, em, tme, sali
 double precision :: tt1, tt2, tseafrz

! parameters for ice thermodynamics and heat fluxes
  double precision,save :: tzero  = 273.15d0       ! tzero C => K

!------------------------------------------------------------------------
!     Some integer constants
!------------------------------------------------------------------------

      ni=nlice
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

!---------------------------------------------------------------------
! other quantities
!---------------------------------------------------------------------

 do k=1,ns
   dzi(k)=zi(k)-zi(k-1)
 enddo
  ziold =  zi
 dziold = dzi

!---------------------------------------------------------------------
! set the temperature fields
!---------------------------------------------------------------------

  tseafrz=Tfreeze1(seasal)      ! temperature at which sea water freezes (degC)
  tiold(0) = tseafrz
  tocn     = tseafrz

!---------------------------------------------------------------------
! transfer salinity and compute melting temperature
!---------------------------------------------------------------------

  do k=1,ni
     sali(k)=si(ni-k+1)
  enddo
  DO k=ni+1,ns
     sali(k) = 0d0
  ENDDO

  do k=1,ns
     tme(k)=Tfreeze1(sali(k))
  enddo

!---------------------------------------------------------------------
! recompute temperature profile from centered values
!---------------------------------------------------------------------

  do k=1,ni
       tt1 = ti(ni-k+1)-tzero
       em(k)= tt1
  enddo
  do k=ni+1,ns
       tt1 = ts(ns-k+1)-tzero
       em(k)= tt1
  enddo
  do k=1,ns-1
       tiold(k)=0.5d0*(em(k)+em(k+1))
  enddo
  tiold(ns)=tsu-tzero

  return

  do k=1,ni
       tt1 = ti(ni-k+1)-tzero
       qm(k) = func_qmelt(tme(k),tt1)
       em(k) = func_el(tme(k),tt1)
  enddo
  do k=ni+1,ns
       tt1 = ts(ns-k+1)-tzero
       qm(k) = func_qmelt(tme(k),tt1)
       em(k) = func_el(tme(k),tt1)
  enddo

  call invert_energy_all(ns,em,tiold,sali)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
end subroutine init_sigma_ice_FE
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ice_thermodynamic.f90 file
!
! contain subroutines: - ice_thermo
!                      - bottom_ice_creation
!
! Rewritten by F. Cyr with few changes:
!         -> in bottom_ice_creation
!           - new input var: "wind"
!           - new variables: dhi_underice, dhi_openwater, hi_accr
!           - new param Cfa, Cfw & cc
!           - brand new ice_creation scheme!
!
! Fred Dupont: 2008 11 01: rewrite the code so the ice is formed or melt before SST changes
! Fred Dupont: 2009 08 01: rewrite the code for delta temperature
! Fred Dupont: 2009 08 03: rewrite the code for variable salinity and sigma transfer of energy.
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!        SUBROUTINE ice_thermo      
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine ice_thermo_FE(dtice)

  implicit none


!-----------------------------------------------------------------------------
!                          THERMODYNAMICS
!-----------------------------------------------------------------------------

! specific heat capicity multiplied by mean sea density 
  double precision,save :: cheat=4.2174d3*1.025d3  ! J/K/m3

! parameters for ice thermodynamics and heat fluxes
  double precision,save :: tzero  = 273.15d0       ! tzero C => K

  double precision,save :: latvap = 2.5d6          ! J/kg latent heat of vaporization
  double precision,save :: latsub = 2.834d6        ! J/kg latent heat of sublimation
  double precision,save :: mlfus = 0.334d6        ! J/kg latent heat of fusion


  double precision, parameter :: StefanBoltzman = 5.67d-8

!-----------------------------------------------------------------------------
! arguments
  double precision       dtice
!locals
  double precision  &
             sstz, bshf, Fnet, dzf, fwf, &
             es, zrchu1,zrchu2,zref,zssdqw
  integer    k,l
  integer    lice_top
  double precision, dimension(0:maxlay) :: &
             rhs, a, b, diag, wmesh, wold, a1,b1,c1,d1
  double precision, dimension(maxlay) :: &
             kice, &
             he, heold, henew, &
             re, tme, dhe, &
             em, cm, qm, km, rad, sali
  double precision blhf, &
             submer_is, hdraft, sice, dhi, dhs, dq, &
             m1, m2, k1, theta_ther, eps, topflux, sum0, alb, estef, hmelt, &
             ftot, &
             Cbrine,zn,salin,&
             Tside,deltaT,wlat,rside, &
             dhi_underice, dhi_openwater, dho, hi_accr, capa, efus, tt, sim, &
             qsr, snow_precip, rain_precip, w1, w2, w3, w4, xi, slope, &
             hiold, hsold, fac, maxvel, slop2, &
             hsnew, hinew, tm, sn, sume, heat_melt_to_ocean, &
             energy_ice_growth, tt1, tt2, xis, xib
  integer, dimension(maxlay) :: p_res

  logical fixed_tb,melt_ts,melt_tb,panic
  integer :: iteration, maxiter = 100
  double precision ftot2, sume2, flux2, pi, energy_snow_melt, &
            sublimation_speed, bottom_speed, rom
  logical :: energy_check=.false., thin_snow_active, debug=.false.
  double precision :: tseafrz, sumrad
      double precision :: hslim=0.0005d0, hicut=1d-2
  double precision Fprec, ssnow ! heat flux due to precipitation [W/m2]


!FD debug
energy_check=.true.
debug=.true.
hslim=1d-3
!---------------------------------------------------------------------
! other quantities
!---------------------------------------------------------------------

  zi =  ziold
 dzi = dziold

!---------------------------------------------------------------------
! init some parameters
!---------------------------------------------------------------------

  
  theta_ther=1.d0    ! time-implicit factor
  fixed_tb=.false.   ! fixed bottom ice temperature at sea water freezing temp
  melt_ts=.false.   ! fixed surface ice temperature at melting temp
  melt_tb=.false.   ! fixed bottom ice temperature at melting temp
  eps = 1.d-5        ! reset hice to eps if below
  thin_snow_active =.false. ! reroute some energy to sublimate/melt snow

!---------------------------------------------------------------------
! set the conductivity in material
!---------------------------------------------------------------------

  kice(1:ni)   =cond_ice
  kice(ni+1:ns)=cond_sno

!---------------------------------------------------------------------
! set to zero some integrable quantities
!---------------------------------------------------------------------

  fwf=0.d0
  bshf=0.d0
  heat_melt_to_ocean = 0.d0
  energy_ice_growth  = 0.d0
  energy_snow_melt   = 0.d0
  bottom_speed       = 0.d0
  sublimation_speed  = 0.d0

!---------------------------------------------------------------------
! set the temperature fields
!---------------------------------------------------------------------

  tseafrz = Tfreeze1(seasal)      ! temperature at which sea water freezes (degC)
  tocn    = tseafrz

! FD we are forced to use tiold, because ti/ts are considered 
! FD to be averaged-layer value by the parent program
  tinew = tiold
  timid = tiold

! other local notations

  re(1:ni   )=rhoice
  re(ni+1:ns)=rhosno

!---------------------------------------------------------------------
! transfer salinity and compute melting temperature
!---------------------------------------------------------------------

  do k=1,ni
     sali(k) = si(ni-k+1)
  enddo
  DO k=ni+1,ns
     sali(k) = 0.d0
  ENDDO

  do k=1,ns
     tme(k)=Tfreeze1(sali(k))
  enddo

  do k=1,ni
!   tt = 0.25d0*(tiold(k-1)+tinew(k-1)+tiold(k)+tinew(k))
   tt = 0.5d0*(tiold(k-1)+tiold(k)) ! both LIM3 and Huwald are not doing any fancier, so why not doing the same?
   kice(k) = func_ki(sali(k),tt)      ! this can be move up and out of the iteration loop.
  enddo

!---------------------------------------------------------------------
! do computation if only enough ice
!---------------------------------------------------------------------

 if (hi.gt.hminice) then

! transformation from average thickness (including lead, i.e. volume) to
! thickness in fraction of cell covered by ice

  hsold=hs
  hiold=hi
  heold(   1:ni)=hiold*dzi(   1:ni)
  heold(ni+1:ns)=hsold*dzi(ni+1:ns)
  he   =heold
  henew=he

  wmesh(0:ns)=0.d0
  dhe = 0.d0
  dhs = 0.d0


!---------------------------------------------------------------------
! layer detection
!---------------------------------------------------------------------

  lice_top=ni
  if (hsold > hslim ) then
    lice_top=ns
  endif

!---------------------------------------------------------------------
! set the radiative transfer (already in W/m2)
!---------------------------------------------------------------------

      sumrad=0.d0
      DO k=lice_top,ni+1,-1
         Rad(k) = swradab_s(ns-k+1)
      ENDDO
      DO k=ni,1,-1
         Rad(k) = swradab_i(ni-k+1)
      ENDDO
      do k=1,lice_top
         sumrad = sumrad + rad(k)
      enddo

!---------------------------------------------------------------------
! degrade to one layer if too close to the limit
! needed during lateral consolidation and melting
! the surface temperature might jump
!---------------------------------------------------------------------

     do k=1,ns
       tt1 = 0.5d0*(tiold(k-1)+tiold(k))
       qm(k) = func_qmelt(tme(k),tt1)
       em(k) = func_el(tme(k),tt1)
     enddo

!---------------------------------------------------------------------
! latent heat calculation, sublimation processes and precipitation
!---------------------------------------------------------------------

      tsu    = tiold(lice_top) + tzero
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
      Fnet =  fac_transmi * swrad + netlw + fsens + flat
if (debug) write(*,*) 'Fnet',Fnet,fac_transmi*swrad,netlw,fsens,flat,qair,q0,tair,tsu,uair,dzf

!---------------------------------------------------------------------
! accumulation at the surface
!---------------------------------------------------------------------
     if (tair < tzero ) then
      snow_precip = snowfall
      rain_precip = 0.d0
     else
      snow_precip = 0.d0
      rain_precip = snowfall * rhosno / rhowat
     endif

   fwf = fwf + rain_precip ! add liquid precipitation to the freshwater flux to ocean

!---------------------------------------------------------------------
! need to get rid of snow if too thin during melting/sublimating
! thin snow case
!---------------------------------------------------------------------
  if ( hsold <= hslim ) then

! then 1- compute total energy in snow
!      2- if Fnet>0, we assume that we can use some of the heat flux to melt the thin snow
!      3- if Fnet<=0,we assume that some of the sublimation can be used to sublimate the snow

    sume2 = 0.d0
    do k=ni+1,ns
      sume2 = sume2 + em(k) * re(k) * heold(k)
    enddo

    if ( Fnet > 0.d0 ) then ! use Fnet if positive

      energy_snow_melt = min( - sume2, Fnet * dtice )! energy (negative) required to melt the thin snow (in entirity if Fnet large enough)
      dhs = 0.d0
      if ( hsold > 0.d0 ) then
       dhs = energy_snow_melt / sume2 * hsold ! melted snow thickness (negative)
      endif
      snow_precip = 0.d0
! FD debug
!write(*,*) 'hsold < hslim',hsold,dhs,Fnet,energy_snow_melt/dtice

    else !sublimation can be used too

      sume2 = sume2 - latvap * re(k) * hsold ! the snow needs to be sublimated, so the total required energy of melting is higher!
      dq = max( sume2, flat * dtice * snosub ) ! negative energy required to sublimate the thin snow (in entirity if flat large enough, partial otherwise)
      if ( hsold > 0.d0 ) then
       dhs = snow_precip * dtice - dq / sume2 * hsold
      else
       dhs = snow_precip * dtice
! FD debug
!write(*,*) 'precip on no snow',dhs
      endif
      ! recalculate the top flux
      Flat =  Flat - dq / dtice ! (-dq>0 so this is in effect reducing the latent heat)
      Fnet =  Fnet - dq / dtice ! (-dq>0 so this is in effect reducing the latent heat)

    endif

    thin_snow_active = .true.

! new velocity at the snow-air interface
    wmesh(ns) = - dhs / dtice

!---------------------------------------------------------------------
  endif ! end thin snow
!---------------------------------------------------------------------

! write(*,*) 'sublimation 2',energy_snow_melt,sublimation_speed

! sublimation at the surface
   sublimation_speed = snosub * flat/(re(lice_top)*(-em(lice_top)+latvap))

if (debug) write(*,*) 'sublimation',sublimation_speed*dtice,flat,snow_precip*dtice


!---------------------------------------------------------------------
! sensible heat flux with ocean

!  bshf = cheat*Ch*u_sea*(tiold(0)-sstz)
   bshf = -oceflx


!-----------------------------------------------------------------------
!     Update energy flux due to snow precipitation
!-----------------------------------------------------------------------
      
      if ( .not.thin_snow_active ) then
          tm    = tme(lice_top)
          tt    = tair - tzero
          efus  = cp_ice * ( tt - tm ) - mlfus + cp_wat * tm
          Fprec = snow_precip * rhosno * efus
      else ! FD no snow
          Fprec = snow_precip * rhosno * em(ns)
      endif

! surface flux
  topflux = Fnet

if (debug) write(*,*) 'topflux',topflux,' botflux',-bshf

!---------------------------------------------------------------------
! vertical regridding due to the sublimation/condensation at the top of the ice/snow
! vertical regridding due to the ice formation at the bottom of the ice


  if (lice_top==ns) then
   wmesh(ns) = - sublimation_speed - snow_precip
   wmesh(ni) =0.d0
   wmesh(0 ) = bottom_speed
  else
   wmesh(ni) = - sublimation_speed ! sublimation or melt
   wmesh(0 ) =        bottom_speed ! ice growth or melt
  endif

! apply a linear scaling between bottom and top ice growth/decay
   fac=0.d0
   do k=1,ni-1
     fac = fac + dzi(k)
     wmesh(k)=wmesh(0 ) * (1.d0-fac) + wmesh(ni) * fac
   enddo
   fac=0.d0
   do k=ni+1,ns-1
     fac = fac + dzi(k)
     wmesh(k)=wmesh(ni) * (1.d0-fac) + wmesh(ns) * fac
   enddo
! FD debug
!wmesh=0.d0


  wold=wmesh

!---------------------------------------------------------------------
! check for significant disappearance of ice
! switch back to one layer approximation



     ftot = dtice * ( Fnet - bshf )
     sume=0.d0
     sume2=0.d0
     do k=1,ns
       sume = sume + heold(k) * em(k) * re(k)
       sume2 = sume2 + heold(k) * qm(k) * re(k)
     enddo

  if ( ftot>(-sume) .or. sstz>tme(ns) ) then
     call degrade_to_two_layer_uniform(  &
                lice_top,&
                re,em,tme,heold,tiold,tinew,topflux,bshf,sali,&
                fwf,hinew,hsnew,sstz,dtice,wmesh,theta_ther)
     iteration = 0
     rhs(0:ns) = 0.d0
     flux2     = 0.d0
     he(   1:ni) = hinew*dzi(   1:ni)
     he(ni+1:ns) = hsnew*dzi(ni+1:ns)
     goto 2000 ! finalize and leave
  endif

!---------------------------------------------------------------------
 do iteration = 1, maxiter ! convergence on temperature and wmesh

  rhs=0.d0
  a=0.d0
  b=0.0d0
  diag=0.d0
  cm=0.d0

  do k=1,lice_top
    tt = 0.25d0*(tiold(k-1)+tinew(k-1)+tiold(k)+tinew(k))
    tt1 = 0.5d0*(tiold(k-1)+tiold(k))
    tt2 = 0.5d0*(tinew(k-1)+tinew(k))

    capa = func_cp ( tme(k), tt1, tt2)
    cm(k)= capa ! needed?
    m1=heold(k)*re(k)/3.d0*capa
    m2=heold(k)*re(k)/6.d0*capa
    k1=dtice/heold(k)*kice(k)
    km(k)=k1 ! needed?
    rhs(k-1) =  rhs(k-1) + k1 *( tiold(k)-tiold(k-1)) + 0.5d0 * rad(k) * dtice
    rhs(k)   =  rhs(k)   - k1 *( tiold(k)-tiold(k-1)) + 0.5d0 * rad(k) * dtice

    w1=wmesh(k-1)*re(k)*dtice
    w2=wmesh(k  )*re(k)*dtice
    rhs(k-1) =  rhs(k-1) - w1 * em(k)
    rhs(k  ) =  rhs(k  ) + w2 * em(k)
    w1=w1 * capa * 0.5d0
    w2=w2 * capa * 0.5d0
!    rhs(k-1) =  rhs(k-1) - w1 * ( tinew(k) + tinew(k-1) )
!    rhs(k  ) =  rhs(k  ) + w2 * ( tinew(k) + tinew(k-1) )

    a  (k-1) =  a  (k-1) - k1 + m2 + w1
    diag(k-1)= diag(k-1) + k1 + m1 + w1
    diag(k  )= diag(k  ) + k1 + m1 - w2
    b  (k)   =  b  (k  ) - k1 + m2 - w2
enddo

!---------------------------------------------------------------------
! implicit terms

! implicit terms  at surface
  diag(lice_top) = diag(lice_top) + theta_ther * dtice * dzf

!---------------------------------------------------------------------


!---------------------------------------------------------------------
! add BCs:

  rhs(0       )= rhs(0       ) - dtice * bshf
  rhs(lice_top)= rhs(lice_top) + dtice * Fnet - energy_snow_melt

!---------------------------------------------------------------------
! in case of melting/growth at fixed temperature, the velocity is the unknown
!---------------------------------------------------------------------

  flux2= 0.d0

k=1
      tinew(k-1) = tseafrz ! keep basal temperature at sea temperature (note that it would make sense that the melting occurs at melt temperature)
      sum0 = em(k) * re(k) * dtice
      w1=wmesh(k-1)*re(k)*dtice
      rhs(k-1) =  rhs(k-1) + w1 * em(k) ! add back
      rhs(k)   = rhs(k) - b(k) * ( tinew(k-1) - tiold(k-1) )
      b(k)= 0.d0
k=0
      rhs(k) =  rhs(k) - diag(k) * ( tinew(k) - tiold(k) )
      diag(k)= sum0

    if ( melt_ts ) then
k=lice_top
      w2=wmesh(k  )*re(k)*dtice
      rhs(k  ) =  rhs(k  ) - w2 * em(k) ! add back

      tinew(k) = tme(k) ! keep top surface temperature at melt point
      sum0 = em(k) * re(k) * dtice
      rhs(k) =  rhs(k) - diag(k) * ( tinew(k) - tiold(k) )
      diag(k)= - sum0

      rhs(k-1)   = rhs(k-1) - a(k-1) * ( tinew(k) - tiold(k) )
      a(k-1)= 0.d0
    endif

!---------------------------------------------------------------------
! add BCs for wmesh:

!  efus = cp_ice * ( tme(lice_top) - tiold(lice_top) ) + mlfus * ( 1.d0 - tme(lice_top)/tiold(lice_top) ) - cp_wat * tme(lice_top)
! precipitation energy transfer
  if (lice_top == ns ) then
    flux2= flux2 + Fprec
  endif
! sublimation energy transfer
  if (.not.melt_ts) then
    efus  = em(lice_top)
    flux2 = flux2 + sublimation_speed * re(lice_top) * efus
  endif

! final combination
  rhs (lice_top) = rhs (lice_top) + dtice * flux2

!---------------------------------------------------------------------
! tri-diagonal solver call

!if (debug) then
!   write(*,'(20(1x,e10.3))') a(0:lice_top)
!   write(*,'(20(1x,e10.3))') diag(0:lice_top)
!   write(*,'(20(1x,e10.3))') b(0:lice_top)
!stop
!endif
!   write(*,*) rhs(0:lice_top)
  call trisolverD(lice_top+1,a(0),b(0),diag(0),rhs(0))
!write(*,*) rhs(0:lice_top)

! new temperature profile
  tinew(0:lice_top)=rhs(0:lice_top)+tiold(0:lice_top)
  if (thin_snow_active) tinew(lice_top+1:ns)=tiold(lice_top+1:ns)

  if (melt_ts ) then
k=lice_top
  if (iteration.gt.10) then
     wmesh(k) = 0.5d0 * wmesh(k) + 0.5d0 * rhs(k)
  else
     wmesh(k) = rhs(k)
  endif
     tinew(k) = tme(k)
     rhs  (k) = tinew(k) - tiold(k) ! needed for posterio diagnostic
  endif

! fixed bottom temperature at sea water temperature
k=0
     wmesh(k) = rhs(k)
     tinew(k) = tseafrz
     rhs  (k) = tinew(k) - tiold(k) ! needed for posterio diagnostic

! FD debug
!  write(*,'(i4,20(1x,e10.3))') iteration,wmesh(0:lice_top)
!  write(*,'(i4,20(1x,e10.3))') iteration,tinew(0:lice_top)


!---------------------------------------------------------------------
! Calculation of ice and snow mass changes due to melting
!---------------------------------------------------------------------

! top
   k=lice_top
   slope = tinew(k)-tme(k)

    if (slope >0.d0) then    ! condifion on top melting  (first time)
      melt_ts=.true.
      tinew(k) = tme(k)
    endif


!---------------------------------------------------------------------
! reset top melting to 'off' if velocity has the wrong sign
!---------------------------------------------------------------------

if ( melt_ts ) then
  if ( lice_top==ns ) then
     if ( wmesh(lice_top) < -snow_precip ) then
       melt_ts = .false.
       wmesh(lice_top)=-sublimation_speed-snow_precip
     endif
  else
     if ( wmesh(lice_top) < 0.d0 ) then
       melt_ts = .false.
       wmesh(lice_top)=-sublimation_speed
     endif
  endif
endif

!---------------------------------------------------------------------
! need to get rid of snow if wmesh too large
!---------------------------------------------------------------------

   if ( wmesh(ns)*dtice > hsold + 1d-18 ) then
      sume2 = 0.d0
      do k=ni+1,ns
        sume2 = sume2 + em(k) * re(k) * heold(k)
      enddo
      energy_snow_melt = -sume2
      fwf = fwf + hsold / dtice / rhowat * rhosno
      dhs = - hsold
      lice_top=ni ! do not treat the snow layer
      melt_ts = .false.
      thin_snow_active = .true.
! FD debug
!write(*,*) 'wmesh(ns) < hsold', energy_snow_melt, hsold, dhs, wmesh(ns)*dtice
      wmesh(ns)=hsold/dtice
      sumrad=0.d0
      do k=1,lice_top
         sumrad = sumrad + rad(k)
      enddo
   endif


!---------------------------------------------------------------------
! compute vertical velocity in sigma coordinate
! and new thickness for ice and snow
!---------------------------------------------------------------------
! vertical regridding due to the sublimation/condensation at the top of the ice/snow
! vertical regridding due to the ice formation at the bottom of the ice

! apply a linear scaling between bottom and top ice growth/decay
   fac=0.d0
   do k=1,ni-1
     fac = fac + dzi(k)
     wmesh(k)=wmesh(0 ) * (1.d0-fac) + wmesh(ni) * fac
   enddo
   fac=0.d0
   do k=ni+1,ns-1
     fac = fac + dzi(k)
     wmesh(k)=wmesh(ni) * (1.d0-fac) + wmesh(ns) * fac
   enddo

! find final thickness after regridding
   dhi = ( wmesh(0 ) - wmesh(ni ) ) * dtice
   dhs = ( 0.d0      - wmesh(ns ) ) * dtice

   hsnew = max(hsold+dhs,0.d0)
   hinew = hiold+dhi
! FD debug
!write(*,*) 'snow',hsnew,dhs,wmesh(ns)*dtice

   he(   1:ni) = hinew*dzi(   1:ni)
   he(ni+1:ns) = hsnew*dzi(ni+1:ns)

!---------------------------------------------------------------------
! reset ice/snow temperature if needed

  do k=lice_top,1,-1
    tinew(k  ) = min(tinew(k  ),tme(k))
    tinew(k-1) = min(tinew(k-1),tme(k))
  enddo
  
! FD debug
!  write(*,'(i4,10(1x,e13.6))') iteration,dhi,dhs,dhi/hiold

  maxvel=0.d0
  do k=0,lice_top
     maxvel=max(maxvel,abs(tinew(k)-timid(k)))
  enddo
!write(*,*) 'maxvel',maxvel

  if (maxvel.lt.1d-12) exit

  timid  = tinew
  wold   = wmesh

 enddo ! convergence on ice temperature and wmesh
!---------------------------------------------------------------------
!---------------------------------------------------------------------

2000 continue

!---------------------------------------------------------------------
! freshwater and heat flux on melting to the ocean
!---------------------------------------------------------------------
   dq    = 0.0d0
   hmelt = 0.0d0

  if ( melt_ts ) then
k=lice_top
    dhi   =   wmesh(k  ) * re(k) / rhowat
    hmelt = hmelt + dhi
    dq    = dq + dhi *cheat*(sstz-tme(k))     ! (W/m2) heat required to warm or
                                              ! cool hmelt of melt water to sstz
! change to ocean
   fwf  = fwf  + hmelt   ! positive flux for input to ocean
   heat_melt_to_ocean = heat_melt_to_ocean  - dq            ! change in heat flux to ocean
  endif


!---------------------------------------------------------------------
! recompute the fluxes for energy check

! sensible heat flux with ocean
!  bshf = bshf + cheat*Ch*u_sea*theta_ther*rhs(0)

! surface flux
  topflux = Fnet - dzf * rhs(lice_top) - energy_snow_melt / dtice

if (energy_check) then
!print the top and bottom heat flux before accreation/deposition processes
! the two numbers should be equal in steady state.
 if (debug) then
  write(*,*) 'ice heat flux',bshf,topflux,flux2,energy_snow_melt/dtice
  write(*,*) 'hi,hs',hinew,hsnew,dhs/dtice,-wmesh(ni),wmesh(0)
 endif
! plus surface adjustment for disappearing thin snow
     ftot2 = dtice * ( topflux - bshf + flux2  + sumrad ) + energy_snow_melt
     if ( thin_snow_active ) ftot2 = ftot2 + dtice * Fprec  ! need to account for heat due to snow precip
     sume2=0.d0
 if ( .not.thin_snow_active ) then
  do k=1,ns
    tt1 = 0.5d0*(tinew(k-1)+tinew(k))
    qm(k) = func_qmelt(tme(k),tt1)
    em(k) = func_el(tme(k),tt1)
    sume2 = sume2 + he(k) * em(k) * re(k)
  enddo
 else ! for thin snow, the energy is assumed unchanged and uniform in the snow, ice is done the same way
  do k=1,ni
    tt1 = 0.5d0*(tinew(k-1)+tinew(k))
    qm(k) = func_qmelt(tme(k),tt1)
    em(k) = func_el(tme(k),tt1)
    sume2 = sume2 + he(k) * em(k) * re(k)
  enddo
  do k=ni+1,ns
    em(k) = func_el(tme(k),tinew(k))
    sume2 = sume2 + he(k) * em(k) * re(k)
  enddo
 endif
   fac=abs(ftot2-sume2+sume)
  write(*,'(i4,300(1x,e10.3))') iteration,wmesh(0:ns)
  write(*,'(i4,300(1x,e10.3))') iteration,tinew(0:ns)
  write(*,*) 'ftot fina',ftot2,sume2-sume,iteration,fac
   fac=abs(ftot2-sume2+sume)
   if (fac>1d-3) then 
       OPEN(1,file='restart.dat')
       WRITE(1,*) nlice,nlsno,ith_cond
       WRITE(1,*) tiold(0:nlice+nlsno),seasal
       WRITE(1,*) si(1:nlice)
       WRITE(1,*) hiold,hsold
       WRITE(1,*) snowfall,dwnlw,tair,qair,uair,swrad,oceflx,pres
       WRITE(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
       CLOSE(1)
    stop 'too large error'
   endif
endif

!---------------------------------------------------------------------
! reset top temperature for thin snow for uniformity
 if ( thin_snow_active) then
  do k=lice_top+1,ns
    tinew(k)=tinew(k-1)
  enddo
 endif

 if (hi .lt. hicut) then
   hinew = hicut
   tinew = tocn
   he(   1:ni) = hinew*dzi(   1:ni)
   he(ni+1:ns) = hsnew*dzi(ni+1:ns)
endif
!---------------------------------------------------------------------
! preparation for ocean thermodynamics

   bshf = bshf +  heat_melt_to_ocean

   hi=hinew
   hs=hsnew
   
   sim=0.D0
   do k=1,lice_top
      sim=sali(k)*he(k)
   enddo
   if (hi>0.d0) sim=sim/hi
   
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

 endif ! end of condition for hi > hslim

!---------------------------------------------------------------------
! prepare for next iteration:
!---------------------------------------------------------------------

    tiold = tinew

! convert back to LIM3
!     do k=1,ni
!       tt1 = 0.5d0*(tinew(k-1)+tinew(k))
!       ti(ni-k+1) = tt1 + tzero + 1d-2
!     enddo
     do k=ni+1,ns
       tt1 = 0.5d0*(tinew(k-1)+tinew(k))
       ts(ns-k+1) = tt1 + tzero + 1d-2
     enddo
     tbo = tinew(0 ) + tzero + 1d-2
     tsu = tinew(ns) + tzero + 1d-2
! FD keep the full FE solution for plotting reasons
     do k=1,ni
      ti(ni-k+1)=tinew(k-1) + tzero + 1d-2
     enddo
! call back for LIM3 of thickness (needed for radiative transfer)
     do k=1,ni
       dzi(ni-k+1) = dziold(k) * hi
     enddo
     do k=ni+1,ns
       dzs(ns-k+1) = dziold(k) * hs
     enddo
!write(*,*) 'ice heat flux',bshf

!because of potential change in steric height due to the presence of the ice, 
!we assume that the ice is not impacting the volume of water, hence
!the flux of volume between the ice and the ocean is set to zero,
!only the salt flux is impacted as a pseudo-volume flux
! (then, full precipitation should fall in water)

! check that this is actually the case in forcing.F90

end subroutine ice_thermo_FE


   !==============================================================================
   ! retrieve mean layer frome energy temperature
   !==============================================================================
   subroutine invert_energy_all(maxlay,em,ti,si)
   implicit none
   integer maxlay
   double precision em(maxlay),ti(0:maxlay),si(maxlay)
   integer k
   double precision ztmelts,zaaa,zbbb,zccc,zdiscrim,tt
      !-------------------
      ! Ice temperatures
      !-------------------
         DO k = 1, maxlay
             !Energy of melting q(S,T) [J.m-3]
             !Ice layer melt temperature
             ztmelts    =  -0.054d0*si(k)
             !Conversion q(S,T) -> T (second order equation)
             zaaa       =  cp_ice
             zbbb       =  ( cp_wat - cp_ice ) * ztmelts - em(k) - mlfus
             zccc       =  mlfus * ztmelts
             zdiscrim   =  SQRT( MAX(zbbb*zbbb - 4.d0*zaaa*zccc,0.d0) )
             tt = ( - zbbb - zdiscrim ) / ( 2.d0 *zaaa )
             tt = MIN( ztmelts, MAX(-80d0, tt ) )
             ti(k) = 2.d0 * tt - ti(k-1)
! FD debug   write(*,*) k,ti(k),ti(k-1),tt,em(k)
         END DO

   end subroutine invert_energy_all
   !==============================================================================


!-----------------------------------------------------------------
subroutine cal_energy(maxlay,ti,si,em)
implicit none
!arguments
  integer maxlay
  double precision, dimension(0:maxlay) :: ti
  double precision, dimension(maxlay)   :: em, si
!locals
  integer k
  double precision tme, tt1, tt2, qm
! specific heat capicity multiplied by mean sea density 
  double precision,save :: mlfus = 0.334d6        ! J/kg latent heat of fusion

     do k=1,maxlay
       tt1 = 0.5d0*(ti(k-1)+ti(k))
       tt2 = min(tt1,-1d-16)
       tme = Tfreeze1(si(k))
       qm  = cp_ice * ( tme - tt1 ) + mlfus * ( 1.d0 - tme/tt2 )
       em(k)= - qm + cp_wat * tme
     enddo

end subroutine cal_energy
!-----------------------------------------------------------------


!****************************************************************
!
subroutine degrade_to_two_layer_uniform(  &
                lice_top,&
                re,em,tme,heold,tiold,tinew,topflux,bshf,si,&
                fwf,hinew,hsnew,sstz,dtice,wmesh,theta_ther)
!
!****************************************************************
 implicit none
! arguments
 integer lice_top
 double precision, dimension(maxlay) :: re, em, tme, heold, si
 double precision, dimension(0:maxlay) :: tinew, tiold, wmesh
 double precision ftot,topflux,bshf,heat_melt_to_ocean, &
                  fwf,hinew,hsnew,dhs,hminice,sstz,dtice, &
                  theta_ther
! local
 double precision sume, dq, dhi, hmelt, sumh, sumr, sums, tt, efus, sum2, fac
 integer k

     wmesh=0.d0
     sume=0.d0
     sumh=0.d0
     sumr=0.d0
     sums=0.d0
     do k=1,ni
      sume = sume + heold(k) * em(k) * re(k)
      sumh = sumh + heold(k) * re(k)
      sumr = sumr + heold(k) * re(k)
      sums = sums + heold(k) * si(k) * re(k)
     enddo
     sumr = sumr / hi
     sums = sums / sumh
     theta_ther = 0.d0
! snow layer
     sum2 = sume
     do k=ni+1,ns
      sum2 = sum2 + heold(k) * em(k) * re(k)
     enddo

     ftot = dtice * ( topflux - bshf )

     if (ftot + sum2 > 0.d0 ) then ! total melt
       hmelt = sumh / rhowat
       fwf  = fwf  + hmelt / dtice  ! positive flux for input to ocean
       dq   = ftot - sume     ! surplus of energy to be given to ocean
!       heat_melt_to_ocean = heat_melt_to_ocean  + dq/dtice              ! change in heat flux to ocean
       tinew = min(tme(1)-0.01d0,sstz)
       hinew = 0.d0
       hsnew = 0.d0
       wmesh(0)=-hi/dtice
       wmesh(ns)=hs/dtice
     else
! partial ice melt
       efus = sume/sumh
       call invert_energy_one(efus,sums,tt)
       ftot = - bshf * dtice ! reverse sign for melting basal ice
       fac = max( ftot / (-sume), 0.d0 )
       dhi  = min(fac, 1.d0) * hi
       hmelt = dhi * sumr
       fwf  = fwf  + hmelt/dtice/rhowat  ! positive flux for input to ocean
       dq = max( ftot + sume, 0.d0 )
       hinew = hi - dhi
       tinew(0:ni) = tt
       wmesh(0)=-dhi/dtice
! partial snow case
! snow layer
       sume = 0.d0
       sumr = 0.d0
       do k=ni+1,ns
        sume = sume + heold(k) * em(k) * re(k)
        sumr = sumr + heold(k) * re(k)
       enddo
       efus = sume/sumr
       call invert_energy_one(efus,0.d0,tt)
       if (nlsno==1) then
         tinew(ni+1:ns) = 2.d0 * tt - tinew(ni) ! trick for nlnso=1
       else
         tinew(ni+1:ns) = tt
       endif
       sumr = sumr / hs
       ftot = topflux * dtice + dq ! add remaining energy after melting ice to energy available for melting snow
       fac = max ( ftot / (-sume), 0.d0 )
       dhs  = min(fac, 1.d0) * hs
       hmelt = dhs * sumr
       fwf  = fwf  + hmelt/dtice/rhowat  ! positive flux for input to ocean
       dq = max( ftot + sume, 0.d0 )
       wmesh(ns)=dhs/dtice
       hsnew = hs - dhs
!       heat_melt_to_ocean = heat_melt_to_ocean  + dq/dtice              ! change in heat flux to ocean
     endif

end subroutine degrade_to_two_layer_uniform



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      SUBROUTINE trisolverD(la,a,b,diag,rhs)
!--------------------------------------------------------------
! TRIAGONAL SYSTEM SOLVER
! for double float
!--------------------------------------------------------------
!
      IMPLICIT NONE
      integer la
      double precision A(la), B(la),DIAG(la),RHS(la)
      double precision tmp,tmpdiag
      INTEGER K

! first Gauss elimination
!
! At surface

          tmp      = 1.d0 / DIAG(1)
          a(1)     = A(1)   * tmp
          rhs(1)   = rhs(1) * tmp 
!
! at mid depth
!
          DO K = 2, LA-1
            tmpdiag  = DIAG(K) - A(K-1) * B(K)
            tmp      = 1.d0 / tmpdiag
            a(K)   = A(K) * tmp
            rhs(K) = (rhs(K) -rhs(K-1)*B(K) ) * tmp
          END DO
!
! at bottom, solution in rhs
!
          tmpdiag    = DIAG(LA) - A(LA-1) * B(LA)
          tmp        = 1.d0 / tmpdiag
          rhs(LA)    = (rhs(LA) -rhs(LA-1)*B(LA) ) * tmp

!--------------------------------------------------------------
! second and final Gauss elimination to surface
! solution in rhs

          DO K = LA-1,1,-1
            rhs(K)   = rhs(K) -rhs(K+1)*A(K)
          END DO

      RETURN
      END SUBROUTINE trisolverD

!---------------------------------------------------------------------


end module ice_thermodynamic_FE

