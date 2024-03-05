module ice_thermo_lim

  use var_thermo_vertical

implicit none

! FD debug
double precision sume1, sume2, dsume, sume3, sume4

      integer ipsnow, ipmelt

! restart
double precision, dimension(maxlay) :: ti0, ts0, dzi0, dzs0, si0
double precision hi0, hs0, tsu0, snowfall0

contains

      SUBROUTINE ice_thermo_diff(dtice,ln_write,numout)
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
        !!           number of layers in the ice/snow: nlice, nlsno
        !!           total ice/snow thickness : ht_i_b, ht_s_b
        !!
        !! ** External : 
        !!
        !! ** References :
        !!
        !! ** History :
        !!           (02-2003) Martin Vancoppenolle, Louvain-la-Neuve, Belgium
        !!
     implicit none
!arguments
      double precision dtice
      LOGICAL :: ln_write
      integer numout
      ! Local variables
      INTEGER    numeqmin, numeqmax, numeq
      double precision &
     &           ztcond_i(0:nlice), &
     &           zkappa_s(0:nlsno), &
     &           zkappa_i(0:nlice),ztstemp(0:nlsno),ztitemp(0:nlice), &
     &           zspeche_i(0:nlice),ztsold(0:nlsno),ztiold(0:nlice),  &
     &           zeta_s(0:nlsno),zeta_i(0:nlice),ztrid(2*maxlay+1,3), &
     &           zindterm(2*maxlay+1),zindtbis(2*maxlay+1), &
     &           zdiagbis(2*maxlay+1)

      double precision &
                 zrchu1, zrchu2, zqsat, zssdqw, ztmelt_i, zkimin, &
                 zg1s, ztsuold, ztsutemp, zg1, zeps, zerrit, &
                 zerrmax, zf, zbeta, zdifcase, zref
      double precision dzf
      double precision es
      integer layer
      integer nconv, nconv_max
      double precision sk, tk
      double precision tmelts
      double precision hslim

      ! Local constants	
      zeps      =  1.0d-20
      zg1s      =  2.d0
      zg1       =  2.d0
      zbeta     =  0.117d0   ! factor for thermal conductivity
      zerrmax   =  1.0d-11   ! max error at the surface
hslim=0.0005d0

      IF ( ln_write ) THEN
         WRITE(numout,*) ' ** ice_thermo_lim : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~ '
      ENDIF

! corrections by Martin
      ! new lines
!      zerrmax   = 1.0e-4
      nconv_max = 50

      ! Switches 
      ! ipsnow equals 1 if snow is present and 0 if absent
! FD pb with sign function      ipsnow   = int(1.d0-max(0.d0,sign(1.d0,-hs)))
      if ( hs > hslim ) then 
        ipsnow = 1
      else
        ipsnow = 0
      endif
      ! ipmelt equals 1 if surface is melting and 0 if not
      ipmelt   = int(max(0.d0,sign(1.d0,tsu-tp0)))

      tmelts = -fracsal*si(1) + tp0
      tmelt   = dble(ipsnow)*tmelt_sno+ (1.d0-dble(ipsnow))*tmelts

      ! Oceanic  heat flux and precipitations
      heatflx_bot = oceflx 
!
!------------------------------------------------------------------------------
!  1) Thermal conductivity at the ice interfaces
!------------------------------------------------------------------------------ 
!
      ! Pringle et al., JGR 2007 formula
      ! 2.11 + 0.09 S/T - 0.011.T

      ! thermal conductivity in the snow
      ztcond_i(0) = func_ki(si(1),ti(1)-tp0)

      DO layer = 1, nlice-1
         sk = 0.5d0 * ( si(layer) + si(layer+1) )
         tk = 0.5d0 * ( ti(layer) + ti(layer+1) - 2.d0 * tp0 )
         ztcond_i(layer) = func_ki(sk,tk)
      END DO

      ztcond_i(nlice) = func_ki(si(nlice),tbo-tp0)

! FD compute initial energy
      sume1 = 0.d0
      DO layer = 1, nlice
         tmelts = -fracsal*si(layer)
         sume1 = sume1 + rhoice * func_El(tmelts,ti(layer) - tp0) * dzi(layer)
      END DO
      DO layer = 1, nlsno
         tmelts = 0.d0
         sume1 = sume1 + rhosno * func_El(tmelts,ts(layer) - tp0) * dzs(layer) * ipsnow
      END DO

! FD debug, restart variable
      ti0(1:nlice) = ti(1:nlice)
      ts0(1:nlsno) = ts(1:nlsno)
      tsu0=tsu
      hi0=hi
      hs0=hs
      si0(1:nlice)=si(1:nlice)
      dzi0(1:nlice) = dzi(1:nlice)
      dzs0(1:nlsno) = dzs(1:nlsno)
      snowfall0 = snowfall

!
!------------------------------------------------------------------------------
!  3) kappa factors                                                             
!------------------------------------------------------------------------------ 
!
      ! snow
      zkappa_s(0)         = cond_sno/max( zeps, dzs(1) )
      do layer = 1, nlsno-1
         zkappa_s(layer)  = 2.d0*cond_sno/ &
     &   max(zeps, dzs(layer) + dzs(layer+1) )
      end do
      zkappa_s(nlsno)    = cond_sno/max(zeps, dzs(nlsno) )

      ! ice
      zkappa_i(0) = ztcond_i(0)/max(zeps, dzi(1) )
      do layer = 1, nlice-1
         zkappa_i(layer)  = 2.d0*ztcond_i(layer)/ &
     &   max(zeps, dzi(layer) + dzi(layer+1) ) 
      end do
      zkappa_i(nlice) = ztcond_i(nlice) / &
     &                          MAX(zeps, dzi(nlice) )

      ! interface
      zkappa_s(nlsno)   = 2.d0*cond_sno*ztcond_i(0)/max(zeps, &
     &             (ztcond_i(0)*dzs(nlsno) + &
     &             cond_sno*dzi(1)))
      zkappa_i(0)        = zkappa_s(nlsno)*dble(ipsnow)  &
     &                     + zkappa_i(0)*(1.d0-dble(ipsnow))

!
!------------------------------------------------------------------------------|
!  4) iterative procedure begins                                               |
!------------------------------------------------------------------------------|
!

        !------------------------------
        ! keeping old values in memory
        !------------------------------

        ztsuold =  tsu
        tsu = min(tsu,tmelt-0.00001d0)
        DO layer = 1, nlsno
           ztsold(layer)    =  ts(layer)
        END DO
        DO layer = 1, nlice
           ztiold(layer)    = ti(layer)
        END DO

      nconv     =  0

        !------------------------------
        ! Beginning of the loop
        !------------------------------
      
      zerrit = 10000.d0

! FD Martin correction      DO WHILE ( zerrit .GT. zerrmax )
      DO WHILE ( ( zerrit .GT. zerrmax ) .AND. ( nconv < nconv_max ) )

        nconv   =  nconv+1

        ztsutemp = tsu
        DO layer = 1, nlsno
           ztstemp(layer) = ts(layer)
        END DO
        DO layer = 1, nlice
           ztitemp(layer) = ti(layer)
        END DO
!
!------------------------------------------------------------------------------|
!  5) specific heat in the ice                                                 |
!------------------------------------------------------------------------------|
!
      DO layer = 1, nlice
         tmelts = -fracsal*si(layer)
         zspeche_i(layer) = func_cp( tmelts, ti(layer)-tp0, &
                                         ztiold(layer)-tp0 )
      END DO
!
!------------------------------------------------------------------------------|
!  6) eta factors                                                              |
!------------------------------------------------------------------------------|
!
      DO layer = 1, nlsno
         zeta_s(layer) = dtice / max(rhosno*cp_ice*dzs(layer),zeps)
      END DO

      DO layer = 1, nlice
        zeta_i(layer) = dtice / max(rhoice*dzi(layer) * &
     &                  zspeche_i(layer)             &
     &                  ,zeps)
      END DO

!
!------------------------------------------------------------------------------|
!  7) surface flux computation                                                 |
!------------------------------------------------------------------------------|
!
      tsu = min(tsu,tmelt)
      ipmelt   = int(max(0.d0,sign(1.d0,tsu-tp0)))
      ! pressure of water vapor saturation (Pa)
      es         =  611.d0*10.d0**(9.5d0*(tsu-273.16d0)/(tsu-7.66d0))
      ! net longwave radiative flux
      netlw = emi*(dwnlw - stefa*tsu*tsu*tsu*tsu)
      ! sensible and latent heat flux
      CALL flx(hi,tsu,tair,qair,fsens, &
     &         flat,q0,zrchu1,zrchu2,uair,zref)
      fsens   =  -fsens
      flat   =  MIN( -flat , 0.d0 ) ! always negative, as precip 
                                           ! energy already added

      ! intermediate variable
      zssdqw     =  q0*q0*pres/ &
     &              (0.622d0*es)*log(10.d0)*9.5d0* &
     &              ((273.16d0-7.66d0)/(tsu-7.66d0)**2.d0)
      ! derivative of the surface atmospheric net flux
          dzf    =  -4.d0*emi*stefa*tsu*tsu*tsu &
     &              -(zrchu1+zrchu2*zssdqw)
      ! surface atmospheric net flux
          zf     =  fac_transmi * swrad + netlw + fsens + flat
! FD debug
!write(*,*) 'atmo flux during nonlinear conv',zf

!
!------------------------------------------------------------------------------|
!  8) tridiagonal system terms                                                 |
!------------------------------------------------------------------------------|
!
      ! layer denotes the number of the layer in the snow or in the ice
      ! numeq denotes the reference number of the equation in the tridiagonal
      ! system, terms of tridiagonal system are indexed as following :
      ! 1 is subdiagonal term, 2 is diagonal and 3 is superdiagonal one

      ! ice interior terms (top equation has the same form as the others)
      DO numeq = 1, maxlay
         ztrid(numeq,1) = 0.d0
         ztrid(numeq,2) = 0.d0
         ztrid(numeq,3) = 0.d0
         zindterm(numeq) = 0.d0
         zindtbis(numeq) = 0.d0
         zdiagbis(numeq) = 0.d0
      END DO
      DO numeq = nlsno + 2, nlsno + nlice 
           layer = numeq - nlsno - 1
           ztrid(numeq,1)   =  - zeta_i(layer)*zkappa_i(layer-1) 
           ztrid(numeq,2)   =  1.d0 + zeta_i(layer)*(zkappa_i(layer-1) + &
     &                       zkappa_i(layer))
           ztrid(numeq,3)   =  - zeta_i(layer)*zkappa_i(layer)
           zindterm(numeq)  =  ztiold(layer) + zeta_i(layer)* &
     &                      swradab_i(layer)
      END DO
      ! ice bottom terms
      numeq   =  nlsno + nlice + 1
      ztrid(numeq,1)  =  - zeta_i(nlice)*zkappa_i(nlice-1)   
      ztrid(numeq,2)  =  1.d0 + zeta_i(nlice)*( zkappa_i(nlice)*zg1 &
     &                   + zkappa_i(nlice-1) )
      ztrid(numeq,3)  =  0.d0
      zindterm(numeq) =  ztiold(nlice) + zeta_i(nlice)* &
     &                   ( swradab_i(nlice) &
     &                   + zkappa_i(nlice)*zg1 &
     &                   *tbo )

      IF ( hs > hslim ) THEN
!
!------------------------------------------------------------------------------|
!  snow-covered cells                                                          |
!------------------------------------------------------------------------------|
!
      ! snow interior terms (bottom equation has the same form as the others)
      do numeq = 3, nlsno + 1
           layer =  numeq - 1
           ztrid(numeq,1)   =  - zeta_s(layer)*zkappa_s(layer-1)
           ztrid(numeq,2)   =  1.d0 + zeta_s(layer)*( zkappa_s(layer-1) + &
     &                         zkappa_s(layer) )
           ztrid(numeq,3)   =  - zeta_s(layer)*zkappa_s(layer)
           zindterm(numeq)  =  ztsold(layer) + zeta_s(layer)* &
     &                         swradab_s(layer)
      end do
      
      ! case of only one layer in the ice (ice equation is altered)
      if (nlice.eq.1) then
         ztrid(nlsno+2,3)    =  0.d0
         zindterm(nlsno+2)   =  zindterm(nlsno+2) + zkappa_i(1) * tbo 
      endif

      IF (tsu.LT.tmelt) THEN
!
!------------------------------------------------------------------------------|
!  case 1 : no surface melting - snow present                                  |
!------------------------------------------------------------------------------|
!
      zdifcase    =  1.d0
      numeqmin    =  1
      numeqmax    =  nlice + nlsno + 1

      ! surface equation
      ztrid(1,1) = 0.d0
      ztrid(1,2) = dzf - zg1s*zkappa_s(0)
      ztrid(1,3) = zg1s*zkappa_s(0)
      zindterm(1) = dzf*tsu - zf

      ! first layer of snow equation
      ztrid(2,1)  =  - zkappa_s(0)*zg1s*zeta_s(1)
      ztrid(2,2)  =  1.d0 + zeta_s(1)*(zkappa_s(1) + zkappa_s(0)*zg1s)
      ztrid(2,3)  =  - zeta_s(1)* zkappa_s(1)
      zindterm(2) =  ztsold(1) + zeta_s(1)*swradab_s(1)

      else
!
!------------------------------------------------------------------------------|
!  case 2 : surface is melting - snow present                                  |
!------------------------------------------------------------------------------|
!
      zdifcase    =  2.d0
      numeqmin    =  2
      numeqmax    =  nlice + nlsno + 1

      ! first layer of snow equation
      ztrid(2,1)  =  0.d0
      ztrid(2,2)  =  1.d0 + zeta_s(1)*(zkappa_s(1) + zkappa_s(0)*zg1s)
      ztrid(2,3)  =  - zeta_s(1)*zkappa_s(1)
      zindterm(2) =  ztsold(1) + zeta_s(1)*(swradab_s(1) + zkappa_s(0)* &
     &               zg1s*tsu)

      ENDIF
      ELSE
!
!------------------------------------------------------------------------------|
!  cells without snow                                                          |
!------------------------------------------------------------------------------|
!
      IF (tsu.lt.tmelt) THEN
!
!------------------------------------------------------------------------------|
!  case 3 : no surface melting - no snow                                       |
!------------------------------------------------------------------------------|
!
      zdifcase    =  3.d0
      numeqmin    =  nlsno + 1
      numeqmax    =  nlice + nlsno + 1

      ! surface equation	
      ztrid(numeqmin,1)   =  0.d0
      ztrid(numeqmin,2)   =  dzf - zkappa_i(0)*zg1    
      ztrid(numeqmin,3)   =  zkappa_i(0)*zg1
      zindterm(numeqmin)  =  dzf*tsu - zf

      ! first layer of ice equation
      ztrid(numeqmin+1,1) =  - zkappa_i(0)*zg1*zeta_i(1)
      ztrid(numeqmin+1,2) =  1.d0 + zeta_i(1)*(zkappa_i(1) + zkappa_i(0)* &
     &                       zg1)  
      ztrid(numeqmin+1,3) =  - zeta_i(1)*zkappa_i(1)  
      zindterm(numeqmin+1)=  ztiold(1) + zeta_i(1)* swradab_i(1)

      ! case of only one layer in the ice (surface & ice equations are altered)
      if (nlice.eq.1) then
      ztrid(numeqmin,1)    =  0.d0
      ztrid(numeqmin,2)    =  dzf - zkappa_i(0)*2.d0
      ztrid(numeqmin,3)    =  zkappa_i(0)*2.d0

      ztrid(numeqmin+1,1)  =  -zkappa_i(0)*2.d0*zeta_i(1)
      ztrid(numeqmin+1,2)  =  1.d0 + zeta_i(1)*(zkappa_i(0)*2.d0 + &
     &                        zkappa_i(1))
      ztrid(numeqmin+1,3)  =  0.d0

      zindterm(numeqmin+1) =  ztiold(1) + zeta_i(1)* &
     &                      ( swradab_i(1) + zkappa_i(1)*tbo )
      endif

      else
!
!------------------------------------------------------------------------------|
!  case 4 : surface is melting - no snow                                       |
!------------------------------------------------------------------------------|
!
      zdifcase    =  4.d0
      numeqmin    =  nlsno + 2
      numeqmax    =  nlice + nlsno + 1

      ! first layer of ice equation
      ztrid(numeqmin,1) =  0.d0
      ztrid(numeqmin,2) =  1.d0 + zeta_i(1)*(zkappa_i(1) + zkappa_i(0)* &
     &                       zg1)  
      ztrid(numeqmin,3) =  - zeta_i(1)* zkappa_i(1)
      zindterm(numeqmin)  =  ztiold(1) + zeta_i(1)*(swradab_i(1) + &
     &                       zkappa_i(0)*zg1*tsu )

      ! case of only one layer in the ice (surface & ice equations are altered)
      if (nlice.eq.1) then
         ztrid(numeqmin,1)  =  0.d0
         ztrid(numeqmin,2)  =  1.d0 + zeta_i(1)*(zkappa_i(0)*2.d0 + &
     &                         zkappa_i(1))
         ztrid(numeqmin,3)  =  0.d0
         zindterm(numeqmin) =  ztiold(1) + zeta_i(1)*           &
     &                         (swradab_i(1)                      &
     &                         + zkappa_i(1)*tbo)               &
     &                         + tsu*zeta_i(1)*zkappa_i(0)*2.d0
      endif

      endif
      endif
!
!------------------------------------------------------------------------------|
!  9) tridiagonal system solving                                               |
!------------------------------------------------------------------------------|
!
      ! solving the tridiagonal system with Gauss elimination
      ! Thomas algorithm, from Computational fluid Dynamics, J.D. ANDERSON, 
      ! McGraw-Hill 1984.	
      ! computing temporary results
      zindtbis(numeqmin) =  zindterm(numeqmin)
      zdiagbis(numeqmin) =  ztrid(numeqmin,2)

      do numeq = numeqmin+1, numeqmax
         zdiagbis(numeq)  =  ztrid(numeq,2) - ztrid(numeq,1)* &
     &                       ztrid(numeq-1,3)/zdiagbis(numeq-1)
         zindtbis(numeq)  =  zindterm(numeq) - ztrid(numeq,1)* &
     &                       zindtbis(numeq-1)/zdiagbis(numeq-1)
      end do

      ! ice temperatures
      ti(nlice)    =  zindtbis(numeqmax)/zdiagbis(numeqmax)
      do numeq = nlice + nlsno + 1, nlsno + 2, -1
         layer    =  numeq - nlsno - 1
         ti(layer)  =  (zindtbis(numeq) - ztrid(numeq,3)* &
     &                       ti(layer+1))/zdiagbis(numeq)
      end do

      ! snow temperatures      
      if ( hs > hslim ) then
      ts(nlsno)  =  (zindtbis(nlsno+1) - ztrid(nlsno+1,3)* &
     &                     ti(1))/zdiagbis(nlsno+1)
          if (nlsno.gt.1) then 
          do numeq = nlsno, 2, -1
             layer    =  numeq - 1
             ts(layer)  =  (zindtbis(numeq) - ztrid(numeq,3)* &
     &                         ts(layer+1))/zdiagbis(numeq)
          end do
          endif
      endif

      ! surface temperature
      if (tsu.lt.tmelt) then
         tsu =  ( zindtbis(numeqmin) - ztrid(numeqmin,3)* &
     &                  ( dble(ipsnow)*ts(1) +             &
     &                  (1.d0-dble(ipsnow))*ti(1) ) ) /     &
     &                  zdiagbis(numeqmin)
      endif
!
!--------------------------------------------------------------------------
!   10) Has the scheme converged ?, end of the iterative procedure        |
!--------------------------------------------------------------------------
!
      ! we verify that nothing has started to melt
      tsu          =  min(tsu,tmelt)
      DO layer  =  1, nlsno
         ts(layer)  =  min(ts(layer),tmelt)
      END DO
      DO layer  =  1, nlice
         ztmelt_i         =  min(-fracsal*si(layer) + tp0, 273.149999999d0) 
         ti(layer)  =  min(ti(layer),ztmelt_i)
      END DO
! FD debug
!write(*,*) 'diff ti/ts',tsu-tp0,ts(1:nlsno)-tp0,ti(1:nlice)-tp0

      ! zerrit is a residual which has to be under zerrmax
      zerrit   =  ABS(tsu-ztsutemp)            
      DO layer = 1, nlsno
         zerrit  =  MAX(zerrit,ABS(ts(layer) - ztstemp(layer)))
      END DO
      DO layer = 1, nlice
         zerrit  =  max(zerrit,abs(ti(layer) - ztitemp(layer)))
      END DO

!
!--------------------------------------------------------------------------
!   11) Heat conduction fluxes                                            |
!--------------------------------------------------------------------------
!

      ! surface conduction flux
      fcsu    =  - dble(ipsnow)*zkappa_s(0)*zg1s*(ts(1) - tsu)  &
     &            - (1.d0-dble(ipsnow))*zkappa_i(0)*zg1*         &
     &            (ti(1) - tsu)

      ! bottom conduction flux
      fcbo  =  - zkappa_i(nlice)* &
     &            ( zg1*(tbo - ti(nlice)) )

! FD mechanism for refreezing the surface
zf=zf+dzf*(tsu-ztsutemp)
if (zf-fcsu.lt.0.d0.and.tsu==tmelt) tsu=tsu-1d-3
! FD debug
!write(*,*) 'diff fcsu',fcsu,zf,zf-fcsu

      END DO ! FD end of nonlinear convergence

      ! internal conduction fluxes : snow
      !--upper snow value
      hfc_s(0) = - ipsnow* &
     &             zkappa_s(0) * zg1s * ( ts(1) - tsu )
      !--basal snow value
      hfc_s(1) = - ipsnow* &
     &             zkappa_s(1) * ( ti(1) - ts(1) )

      ! internal conduction fluxes : ice
      !--upper layer
      hfc_i(0) =  - ipsnow * &   ! interface flux if there is snow
     &         ( zkappa_i(0)  * ( ti(1) - ts(nlsno ) ) ) &
     &         - ( 1.0 - ipsnow ) * ( zkappa_i(0) * &
     &           zg1 * ( ti(1) - tsu ) ) ! upper flux if no
      !--internal ice layers
      DO layer = 1, nlice - 1
         hfc_i(layer) = - zkappa_i(layer) * ( ti(layer+1) - &
     &                      ti(layer) )
      END DO
      !--under the basal ice layer
      hfc_i(nlice) = fcbo

      ! case of only one layer in the ice
      IF (nlice.EQ.1) THEN
         fcsu = -dble(ipsnow)*(zkappa_s(0)*(zg1s*(ts(1)- &
     &                         tsu))) &
     &           -(1.d0-dble(ipsnow))*zkappa_i(0)*zg1  *(ti(1) - tsu)
         fcbo = -zkappa_i(nlice)*zg1* (tbo-ti(nlice))
      ENDIF
!
!--------------------------------------------------------------------------
!   12) Update atmospheric heat fluxes and energy of melting              |
!--------------------------------------------------------------------------
!
      ! pressure of water vapor saturation (Pa)
      es         =  611.d0*10.d0**(9.5d0*(tsu-273.16d0)/(tsu-7.66d0))
      ! net longwave radiative flux
      netlw = emi*(dwnlw-stefa*tsu*tsu*tsu*tsu)
      ! sensible and latent heat flux
      CALL flx(hi,tsu,tair,qair,fsens, &
     &         flat,q0,zrchu1,zrchu2,uair,zref)

      fsens   =  -fsens
      flat   =  MIN( -flat , 0.d0 ) ! always negative, as precip 
                                           ! energy already added
! FD debug compute energy
      sume3 = 0.d0
      DO layer = 1, nlice
         tmelts = -fracsal*si(layer)
         sume3 = sume3 + rhoice * func_El(tmelts,ti(layer) - tp0) * dzi(layer)
      END DO
      DO layer = 1, nlsno
         tmelts = 0.d0
         sume3 = sume3 + rhosno * func_El(tmelts,ts(layer) - tp0) * dzs(layer) * ipsnow
      END DO
! FD debug energy
      dsume = fcsu-fcbo
      do layer=1,nlice
        dsume = dsume + swradab_i(layer)
      enddo
      do layer=1,nlsno
        dsume = dsume + swradab_s(layer) * ipsnow ! FD
      enddo
! FD debug write(*,*) 'debug energy after diff',sume3,sume3-sume1,dtice*dsume

      do layer=1,nlice
        dsume = dsume + swradab_i(layer)
      enddo
      do layer=1,nlsno
        dsume = dsume + swradab_s(layer) * ipsnow ! FD
      enddo
! FD debug energy
      dsume = netlw + flat + fsens + fac_transmi * swrad + oceflx
      do layer=1,nlice
        dsume = dsume + swradab_i(layer)
      enddo
      do layer=1,nlsno
        dsume = dsume + swradab_s(layer) * ipsnow ! FD
      enddo

      ! ice energy of melting
      CALL ice_energy
      
      IF ( ln_write ) THEN
         WRITE(numout,*) ' nconv : ', nconv
         WRITE(numout,*) ' zerrit : ', zerrit
         WRITE(numout,*) ' t_su_b: ', tsu 
         WRITE(numout,*) ' t_s_b : ', ( ts(layer), &
     &                      layer = 1, nlsno )
         WRITE(numout,*) ' t_i_b : ', ( ti(layer), &
     &                      layer = 1, nlice )
         WRITE(numout,*) ' t_bo_b : ', tbo
         WRITE(numout,*)
      ENDIF

!
!------------------------------------------------------------------------------
! End of ice_th_diff
      END SUBROUTINE ice_thermo_diff

      SUBROUTINE ice_energy
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_energy *** 
      !!                 
      !! ** Purpose :   Computes sea ice energy of melting q_i (J.m-3)
      !!
      !! ** Method  :   Formula (Bitz and Lipscomb, 1999)
      !!
      !! history : Martin Vancoppenolle, May 2007
      !!-------------------------------------------------------------------

      double precision    :: &      !: goes to trash
     &   ztmelts                    !: sea ice freezing point in degC

      INTEGER              ::     &
     &   jk                         !: vertical index

      !!-------------------------------------------------------------------

      ! Sea ice energy of melting
      DO jk = 1, nlice
            ztmelts      =   - fracsal * si(jk)
            qi(jk) = rhoice * func_qm ( ztmelts, ti(jk)-tp0 )
      END DO !jk

      ! Snow energy of melting
      DO jk = 1, nlsno
            ztmelts      =   0.d0
            qs(jk) = rhosno * func_qm ( ztmelts, ts(jk)-tp0 )
      END DO !jk

!------------------------------------------------------------------------------      
      END SUBROUTINE ice_energy

      SUBROUTINE ice_thermo_dh(dtice,ln_write,numout)
       !!------------------------------------------------------------------
       !!                ***         ROUTINE ice_th_dh         ***
       !! ** Purpose :
       !!           This routine determines variations of ice and snow thicknesses.
       !! ** Method  :
       !!           Ice/Snow surface melting arises from imbalance in surface fluxes
       !!           Bottom accretion/ablation arises from flux budget
       !!           Snow thickness can increase by precipitation and decrease by 
       !!              sublimation
       !!           If snow load excesses Archmiede limit, snow-ice is formed by
       !!              the flooding of sea-water in the snow
       !! ** Steps  
       !!           1) Compute available flux of heat for surface ablation
       !!           2) Compute snow and sea ice enthalpies
       !!           3) Surface ablation and sublimation
       !!           4) Bottom accretion/ablation
       !!           5) Case of Total ablation
       !!           6) Snow ice formation
       !!
       !! ** Inputs / Outputs
       !!
       !! ** External
       !!
       !! ** References : Bitz and Lipscomb, JGR 99
       !!                 Fichefet T. and M. Maqueda 1997, 
       !!                 J. Geophys. Res., 102(C6), 12609-12646   
       !!                 Vancoppenolle, Fichefet and Bitz, GRL 2005
       !!                 Vancoppenolle et al., OM08
       !!
       !! ** History  : 
       !!   original code    01-04 (LIM)
       !!   original routine
       !!               (05-2003) M. Vancoppenolle, Louvain-La-Neuve, Belgium
       !!               (05-2008) BIO-LIM
       !!
       !!------------------------------------------------------------------
       !! * Arguments
       !!------------------------------------------------------------------
      double precision dtice
      LOGICAL ln_write
      integer numout

      ! Local Variables
      double precision z_f_surf
      double precision zdeltah(maxlay)
      double precision dh_s_prec, dh_s_melt, dhs_subl
      double precision zqt_s_ini, zqt_s_fin, zdqt_s
      double precision zqt_i_ini, zqt_i_fin, zdqt_i
      double precision s_i_max
      double precision cons_err, zeps
      double precision zqfont_su, zzf_surf
      double precision dhi_subl
      double precision hsold, hsnew
      double precision old_hi
      integer num_iter_max, iter
      double precision tmelts
      integer layer
      double precision zfracs, zdhcf, zds, zhnnew, zhn, zhgnewzswi1
      double precision zswi1, zswi2, zswi12, zqprec, zqfont_bo
      double precision zzf_base, zoldsinew, zhgnew, zgrr

      zqt_s_ini = 0.d0
      zqt_s_fin = 0.d0
      zdqt_s    = 0.d0
      zqt_i_ini = 0.d0
      zqt_i_fin = 0.d0
      zdqt_i    = 0.d0
      s_i_max   = 15.d0

      ! Local Constants
      zeps = 1.0d-20

      IF ( ln_write ) THEN
         WRITE(numout,*) ' ** ice_th_dh : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~ '
         WRITE(numout,*) 
      ENDIF
!
!------------------------------------------------------------------------------|
!  1) Calculate available heat for surface ablation         
!------------------------------------------------------------------------------|
!
      z_f_surf  = netlw + flat + fsens - fcsu + fac_transmi * swrad
      z_f_surf  = MAX(0.d0, z_f_surf)
! FD      z_f_surf  = z_f_surf * MAX(0.d0, SIGN(1.d0,tsu-tmelt))
      if (tsu < tmelt ) z_f_surf  = 0.d0

      IF ( ln_write ) THEN
         WRITE(numout,*) ' Available heat for surface ablation ... '
         WRITE(numout,*) 
         WRITE(numout,*) ' z_f_surf : ', z_f_surf
         WRITE(numout,*) ' fratsb   : ', netlw
         WRITE(numout,*) ' fleb     : ', flat
         WRITE(numout,*) ' fcsb     : ', fsens
         WRITE(numout,*) ' fc_su    : ', fcsu
         WRITE(numout,*) ' ab       : ', fac_transmi 
         WRITE(numout,*) ' fsolgb   : ', swrad
         WRITE(numout,*)
         WRITE(numout,*) ' ht_i_b    : ', hi
         WRITE(numout,*) ' ht_s_b    : ', hs
         WRITE(numout,*) ' t_su_b   : ', tsu
         WRITE(numout,*)
      ENDIF

!
!------------------------------------------------------------------------------|
!  2) Snowfall and surface melt                                                |
!------------------------------------------------------------------------------|
!
      ! total snow heat content for conservation
! FD !!! this assumes ns=1 !!!
! FD      zqt_s_ini = qs(1) * hs
      zqt_s_ini = 0.d0
      DO layer = 1, nlsno !in case of melting of more than 1 layer
         zqt_s_ini = zqt_s_ini + qs(layer) * dzs(layer) * ipsnow ! FD
      ENDDO

      IF ( ln_write ) THEN
         WRITE(numout,*) ' Surface ablation and sublimation ... '
         WRITE(numout,*) 
         WRITE(numout,*) ' zqt_s_ini : ', zqt_s_ini / dtice
         WRITE(numout,*) ' ht_s_b    : ', hs
         WRITE(numout,*) ' q_s_b(1)  : ', qs(1)
         WRITE(numout,*)
      ENDIF

      !----------
      ! Snowfall (m/s) ! FD
      !----------
      dh_s_prec    =  snowfall * dtice
      dh_s_melt    =  0.d0
      zqprec       =  rhosno * ( cp_ice * ( tp0 - tair ) + mlfus )
! FD check for tair < tp0
      if (tair > tp0 .or. z_f_surf > 0.d0 ) then
          snowfall = 0.d0
         dh_s_prec = 0.d0
      endif
      ! Conservation update
      zqt_s_ini    =  zqt_s_ini + zqprec * snowfall * dtice
      fprecip      =  - zqprec * snowfall * ipsnow ! FD
! FD debug energy
      dsume = dsume + fprecip

      IF ( ln_write ) THEN
         WRITE(numout,*) ' snow falls! '
         WRITE(numout,*) ' dh_s_prec : ', dh_s_prec
         WRITE(numout,*) ' flux of h : ', -fprecip
! FD         WRITE(numout,*) ' zqt_s_ini : ', zqt_s_ini / dtice
         WRITE(numout,*)
      ENDIF

      !-----------
      ! Snow melt
      !-----------
      ! Melt of fallen snow
      zqfont_su        =  z_f_surf * dtice 
      IF ( ln_write ) WRITE(numout,*) ' snow melts! '
      IF ( ln_write ) WRITE(numout,*) ' zqfont_su : ', zqfont_su / dtice 

! FD not sure that you can treat snow this way
! FD      zdeltah(1)       =  MIN( 0.d0 , - zqfont_su / MAX( zqprec , zeps ) )
! FD      zqfont_su        =  MAX( 0.d0 , - dh_s_prec  * ipsnow - zdeltah(1) ) * zqprec ! FD
! FD      zdeltah(1)       =  MAX( - dh_s_prec, zdeltah(1) )
! FD      dh_s_melt        =  dh_s_melt + zdeltah(1) * ipsnow ! FD

      ! Melt of snow
      DO layer = 1, nlsno !in case of melting of more than 1 layer
         zdeltah(layer) =  - zqfont_su / MAX(zeps, qs(layer) )
         zqfont_su      = MAX( 0.d0, - dzs(layer) * ipsnow - zdeltah(layer) ) * qs(layer)
         zdeltah(layer) = MAX( zdeltah(layer), - dzs(layer) )
         dh_s_melt      =  dh_s_melt + zdeltah(layer) * ipsnow !resulting melt of snow    !FD
      END DO
      IF ( ln_write ) WRITE(numout,*) ' dh_s_melt : ', dh_s_melt 
      dhs     =  dh_s_melt + dh_s_prec

      ! old and new snow thicknesses
      hsold   =  hs
      hsnew   =  hs + dhs

      ! if snow is still present zhn = 1, else zhn = 0
      zhn     =  1.d0 - MAX( 0.d0 , SIGN( 1.d0 , - hsnew ) )
      hs      =  MAX( 0.d0 , hsnew )

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Conservation test for snow
! FD this assumes nlsno=1 !!!!
! FD      zqt_s_fin = qs(1) * hs
! FD      zqt_s_fin = 0.d0
! FD pb here even if assumed multi-layer because dzs has been updated yet
! FD      DO layer = 1, nlsno !in case of melting of more than 1 layer
! FD         zqt_s_fin = zqt_s_fin + qs(layer) * dzs(layer)
! FD      ENDDO
! FD      zdqt_s = zqt_s_fin - zqt_s_ini

! FDIF ( ln_write ) THEN
! FD      WRITE(numout,*)
! FD      WRITE(numout,*) ' Conservation in snow... '
! FD      WRITE(numout,*) ' dh_s_melt : ', dh_s_melt
! FD      WRITE(numout,*) ' dh_s_prec : ', dh_s_prec
! FD      WRITE(numout,*) ' ht_s_b    : ', hs
! FD      WRITE(numout,*) ' zqt_s_ini : ', zqt_s_ini / dtice
! FD      WRITE(numout,*) ' zqt_s_fin : ', zqt_s_fin / dtice
! FD      WRITE(numout,*) ' zdqt_s    : ', zdqt_s / dtice
! FD      WRITE(numout,*) ' z_f_surf  : ', - z_f_surf
! FD      IF ( zqt_s_fin.GT.0.0 ) THEN
! FD         cons_err = ABS(zdqt_s / dtice  + z_f_surf )
! FD      ELSE! FD
! FD         cons_err = ABS(zqt_s_ini / dtice + zdqt_s / dtice )
! FD      ENDIF
! FD      WRITE(numout,*) ' Cons error, snow : ', cons_err
! FD      WRITE(numout,*)
! FDENDIF
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !------------------
      ! Ice surface melt
      !------------------
      IF ( ln_write ) WRITE(numout,*) ' ice melts!  '

      zzf_surf = zqfont_su / dtice
      zdqt_i       = 0.d0
      dhi_surf     = 0.d0
      DO layer = 1, nlice
         zdeltah(layer) =  - zqfont_su / qi(layer)
         zqfont_su      =  MAX( 0.d0 , - dzi(layer) - zdeltah(layer) ) *  qi(layer)
         zdeltah(layer) =  MAX( zdeltah(layer) , - dzi(layer) )
         dhi_surf       =  dhi_surf + zdeltah(layer) !resulting melt of ice
         zdqt_i         =  zdqt_i  + zdeltah(layer) * qi(layer) / dtice
      END DO

      cons_err = ABS( zzf_surf + zdqt_i )

      IF ( ln_write ) THEN
         WRITE(numout,*) ' Conservation in sea ice, surface '
         WRITE(numout,*) ' dh_i_surf: ', dhi_surf
         WRITE(numout,*) ' ht_i_b   : ', hi
         WRITE(numout,*) ' zzf_surf : ', zzf_surf
         WRITE(numout,*) ' zdqt_i   : ', zdqt_i
         WRITE(numout,*) ' Cons error, ice : ', cons_err
         WRITE(numout,*)
      ENDIF
      
!
!------------------------------------------------------------------------------|
!  3) Sublimation at the surface                                               |
!------------------------------------------------------------------------------|
!
      !------------------
      ! Snow sublimation
      !------------------

      !------------------
      ! Ice sublimation
      !------------------

      !-------------------
      ! Snow condensation
      !-------------------

      ! 4.3) Snow/ice sublimation
      !
      ! If fleb is negative, snow condensates at the surface.
      !
      dhs_subl =  - snosub*flat/(rhosno*mlsub)*dtice
      dhs      =  dhs + dhs_subl
      zdhcf    =  hs + dhs_subl 
      hs       =  MAX(0.d0, zdhcf)
      dhs      =  hs - hsold
      dhi_subl =  - MAX(0.d0,-zdhcf)*rhosno/rhoice

      dhi_surf =  dhi_surf + dhi_subl

      hsnew    =  hs

      IF (hs.le.0.0) THEN
         dhs =  MAX ( 0.d0 , dhs )
      ENDIF

      IF ( ln_write ) THEN
         WRITE(numout,*) ' Snow sublimation ... '
         WRITE(numout,*) ' '
         WRITE(numout,*) ' snosub : ', snosub
         WRITE(numout,*) ' dhs_subl : ', dhs_subl
         WRITE(numout,*) ' dhi_subl : ', dhi_subl
      ENDIF

!
!------------------------------------------------------------------------------|
!  4) Basal growth and melt                                                    |
!------------------------------------------------------------------------------|
!
      IF ( ln_write ) THEN
         WRITE(numout,*) ' Basal growth and melt ... '
         WRITE(numout,*) 
         WRITE(numout,*) ' fbbqb     : ', oceflx
         WRITE(numout,*) ' fc_bo     : ', fcbo
         WRITE(numout,*)
      ENDIF

      ! formation / melting of ice at the base is determined by the balance of
      ! the conductive heat flux in the ice (fc_bo_i), and the heat fluxes 
      ! from the ocean (fbbqb). Conductive heat flux is positive 
      ! downwards and fbbq is positive to the ice, i.e., upwards.
      ! Melt/formation rate is modulated by the enthalpy of the bottom ice layer. 
      ! accretion of ice at the base

      !--------------
      ! Basal growth
      !--------------
      ! Initial value (tested 1D, can be anything between 1 and 20)
! FD test      si_acc_new   = 10.d0
      si_acc_new   = 4.d0
      num_iter_max = 10

      ! the growth rate (dh_i_bott) is function of the new ice
      ! heat content (q_i_b(nlay_i+1)). q_i_b depends on the new ice
      ! salinity (snewice). snewice depends on dh_i_bott
      ! it converges quickly, so, no problem

      ! Iterative procedure
      IF ( ( fcbo+oceflx ) .LT. 0.d0 ) THEN

         IF ( ln_write ) WRITE(numout,*) ' Energy available for basal growth : ',fcbo+oceflx

         DO iter = 1, num_iter_max
               tmelts             =   - fracsal * si_acc_new
               ! New ice heat content (Bitz and Lipscomb, 1999)
               qi(nlice+1) = rhoice * func_qm ( tmelts, tbo-tp0 )
               dhi_bot    =  - dtice * ( fcbo + oceflx ) / qi(nlice+1) ! Bottom growth rate = - F*dt / q
               !
               ! New ice salinity ( Cox and Weeks, JGR, 1988 )
               !
               ! zswi2  (1) if dh_i_bott/rdt .GT. 3.6e-7
               ! zswi12 (1) if dh_i_bott/rdt .LT. 3.6e-7 and .GT. 2.0e-8
               ! zswi1  (1) if dh_i_bott/rdt .LT. 2.0e-8
               !
               zgrr   = MIN( 1.0d-3, MAX ( dhi_bot / dtice , zeps ) )
               zswi2  = MAX( 0.d0 , SIGN( 1.d0 , zgrr - 3.6d-7 ) )
               zswi12 = MAX( 0.d0 , SIGN( 1.d0 , zgrr - 2.0d-8 ) ) * &
     &                  ( 1.d0 - zswi2 )
               zswi1  = 1.d0 - zswi2 * zswi12
               zfracs = zswi1  * 0.12d0 +   &
     &                  zswi12 * ( 0.8925d0 + 0.0568d0 * &
     &                  LOG( 100.d0 * zgrr ) ) +  &
     &                  zswi2  * 0.26d0 /   &
     &                ( 0.26d0 + 0.74d0 * EXP ( - 724300.d0 * zgrr ) )
! FD too strong assumption               zfracs = 1. ! FD Martin shortcuts the above calculation
               zds     = zfracs * seasal - si_acc_new
               si_acc_new = zfracs * seasal

               ! the brine volume in the skeletal layer is equal to f
               eskel  = zfracs

               ! salt flux due to initial salt entrapment
               fsalt = seasal * ( 1.d0 - zfracs ) * dhi_bot / dtice * &
     &                rhoice / 1000.d0

         END DO ! iter
      ENDIF ! fc_bo_i

      si_acc_new = 4.d0 ! FD force the bottom salinity
!      IF ( ln_write ) THEN
!         WRITE(numout,*) ' eskel : ', eskel
!         WRITE(numout,*)
!      ENDIF

      ! Final values
      IF ( (fcbo+oceflx) .LT. 0.d0 ) THEN
      ! New ice salinity must not exceed 15 psu
         zoldsinew   = si_acc_new
         si_acc_new  = MIN( si_acc_new , s_i_max )
         ! Metling point in K
         tmelts             =   - fracsal * si_acc_new
         ! New ice heat content (Bitz and Lipscomb, 1999)
         qi(nlice+1) = rhoice * func_qm ( tmelts, tbo-tp0 )
         dhi_bot    =  - dtice * ( fcbo + oceflx ) / qi(nlice+1)
         IF ( ln_write ) WRITE(numout,*) ' dh_i_bott : ', dhi_bot

      ENDIF 

      !-----------------
      ! Basal melt
      !-----------------
      IF ( ( fcbo + oceflx ) .GE. 0.d0 ) THEN

         IF ( ln_write ) WRITE(numout,*) ' Energy available for basal melt   : ',fcbo + oceflx

         zqfont_bo   = dtice * ( fcbo + oceflx )
         zzf_base    = zqfont_bo / dtice
         zdqt_i      = 0.d0

         IF ( ln_write ) WRITE(numout,*) ' zqfont_bo : ', zqfont_bo

         dhi_bot     =  0.d0
         DO layer = nlice, 1, -1
            zdeltah(layer) =  - zqfont_bo / qi(layer)
            zqfont_bo      =  MAX ( 0.d0 , - dzi(layer) - zdeltah(layer) ) *  qi(layer)
            dhi_bot        =  dhi_bot + zdeltah(layer)
            zdqt_i         =  zdqt_i + zdeltah(layer) * qi(layer) / dtice
         END DO

         IF ( ln_write ) WRITE(numout,*) ' dh_i_bott : ', dhi_bot

         cons_err = ABS( zzf_base + zdqt_i )

         IF ( ln_write ) THEN
            WRITE(numout,*) ' Conservation in sea ice, base '
            WRITE(numout,*) ' dh_i_bott: ', dhi_bot
            WRITE(numout,*) ' ht_i_b   : ', hi
            WRITE(numout,*) ' zzf_base : ', zzf_base
            WRITE(numout,*) ' zdqt_i   : ', zdqt_i
!           WRITE(numout,*) ' Conservation error  ice surface : ', 
!    &                      cons_err
!           WRITE(numout,*)
         ENDIF

      ENDIF

      ! It can be than an internal temperature is greater than melt point
      ! then, see lim3 for correction

      ! new ice thickness
      zhgnew         = hi + dhi_surf + dhi_bot
      old_hi = hi

      hi = zhgnew

      !--------------------------------------
      ! Meltwater flow due to surface melt
      !--------------------------------------
      massmelt = ( - rhoice * MIN ( dhi_surf  , 0.d0 ) &
     &             - rhosno * MIN ( dh_s_melt , 0.d0 ) )  
      IF ( ln_write ) THEN
         WRITE(numout,*) ' massmelt : ', massmelt
         WRITE(numout,*)
      ENDIF

!
!------------------------------------------------------------------------------|
!  5) Formation of snow-ice                                                    |
!------------------------------------------------------------------------------|
!
      ! When snow load excesses Archimede's limit, snow-ice interface goes down
      ! under sea-level, flooding of seawater transforms snow into ice
      ! dh_snowice is positive for the ice

      dh_sni = MAX( 0.d0 , ( rhosno * hs + (rhoice - rhowat ) * hi) / ( rhosno + rhowat - rhoice ) )
        
      zhgnew         = MAX( zhgnew , zhgnew + dh_sni )
      zhnnew         = MIN( hs , hs - dh_sni )

      hs  = zhnnew
      hi  = zhgnew

      IF ( ln_write ) THEN
         WRITE(numout,*) ' dh_snowice : ', dh_sni
         WRITE(numout,*)
         WRITE(numout,*) ' At the end of the routine ... '
         WRITE(numout,*) ' ht_s_b : ', hs
         WRITE(numout,*) ' ht_i_b : ', hi
      ENDIF

!------------------------------------------------------------------------------|
! Fin de la subroutine ice_th_dh
      END SUBROUTINE ice_thermo_dh


      SUBROUTINE ice_thermo_remap(dtice,ln_write,numout)

!!-----------------------------------------------------------------------------
!! ** Purpose :
!!              This routine redistributes heat content and salt mass
!!              on the new grid
!!
!!
!! ** Method  : Linear redistribution
!!           
!! ** Steps   : switches, snow, ice heat content, ice salt content
!!
!! ** Arguments
!!
!! ** Inputs / Outputs
!!
!! ** External
!!
!! ** References
!!
!! ** History  : (05-2003) Martin Vancoppenolle, LIM1D
!!               (05-2008) Martin Vancoppenolle, LLN, revision BIOLIM
!!               (09-2009) Martin Vancoppenolle, LLN, restructuration BIOLIM
!!-----------------------------------------------------------------------------
      ! arguments
      double precision dtice
      LOGICAL ln_write
      integer numout

      ! Local Variables
      double precision, DIMENSION (maxlay) :: &
     &   zsh_i0              ,&   !: old ice salt content (ppt.m-2)
     &   zsh_i1                   !: new ice salt content (ppt.m-2)

      double precision, DIMENSION (maxlay) :: &
     &   zqh_i0              ,&   !: old ice heat content (J.m-2)
     &   zqh_i1                   !: new ice heat content (J.m-2)

      double precision, DIMENSION (maxlay) :: &
     &   zqh_s0              ,&   !: old snow heat content (J.m-2)
     &   zqh_s1                   !: new snow heat content (J.m-2)

      double precision, DIMENSION( maxlay ) :: &
     &  zthick1                  !: thickness of physical layers

      double precision, DIMENSION (maxlay+2) :: &
     &   zthick0                 !: thickness of old layers
! FD additions
      LOGICAL ln_con
      double precision zqt_s_ini,zqt_i_ini,zqt_s_fin,zqt_i_fin,zqt_ini,zqt_fin
      double precision zfmelt
      integer snind, icsuswi, icboind, snswi, icboswi, snicswi, icsuind
      integer snicind
      double precision deltah, zqsnow, zdeltah, zhsnow
      double precision s_i_snic, zsh_i_new, z_ms_i_ini
      double precision zeps, zerror
      integer layer, limsum,  ntop1, nbot1, layer_a, ntop0, nbot0
      double precision aaa,bbb,ccc,discrim
      double precision zfgrml, zfsigr, zqh_i_sni, zqh_i_new, zswitch
      double precision zsh_i_sni, z_ms_i_fin
      double precision tmelts

      ! Local Constants
      zeps   = 1.0d-20

      ln_con   = .TRUE.

! FD this should be move outside...
      si(:) = sinew(:) ! update salinity

      IF (ln_write) THEN
         WRITE(numout,*) ' ** ice_phy_remap : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) 
      ENDIF
!
!------------------------------------------------------------------------------|
!  CONSERVATION CHECKS                                                         |
!------------------------------------------------------------------------------|
!

      z_ms_i_ini = 0.d0
      DO layer = 1, nlice
         z_ms_i_ini = z_ms_i_ini + si(layer) * dzi(layer)
      END DO
      IF ( ln_write ) WRITE(numout,*) ' z_ms_i_ini : ', z_ms_i_ini

      zfmelt = ( si(1) * dhi_surf + si(nlice) * MIN(dhi_bot,0.d0) ) / dtice

      zqt_s_ini = 0.d0 ! total heat content in snow (initial value)
      zqt_i_ini = 0.d0 ! total heat content in ice  (initial value)
      zqt_s_fin = 0.d0 ! total heat content in snow (final value)
      zqt_i_fin = 0.d0 ! total heat content in ice  (final value)
      zqt_ini   = 0.d0 ! total heat content snow + ice (initial value)
      zqt_fin   = 0.d0 ! total heat content snow + ice (final value)

!
!------------------------------------------------------------------------------|
!  SWITCHES                                                                    |
!------------------------------------------------------------------------------|
!
      ! snow surface behaviour : calcul de snind-snswi
      ! snind : index tel que sa valeur est 0 si la neige s'accrète
      !                1 si la 1ere couche est entamée par la fonte
      !                2 si la 2eme ....
      !                                  etc ...
      snind    =  0  
      deltah   =  0.d0
      DO layer = 1, nlsno
         IF ( - dhs .GT. deltah ) THEN
            snind  =  layer 
         ENDIF
         deltah = deltah + dzs(layer)
      END DO
      snind = snind * ipsnow ! FD

      ! snswi : switch which value equals 1 if snow melts
      !                                 0 if not
      snswi    =  MAX(0,INT(-dhs/MAX(zeps,ABS(dhs)))) * ipsnow ! FD

      ! ice surface behaviour : computation of icsuind-icsuswi
      ! icsuind : index which values equals 0 if there is nothing
      !                                 1 if first layer has started to melt
      !                                 2 if first and second layer have started ...
      !                                 etc ...
      icsuind    =  0  
      deltah   =  0.d0
      DO layer = 1, nlice
         IF  ( -dhi_surf.GT.deltah) THEN
            icsuind  =  layer 
         ENDIF
         deltah = deltah + dzi(layer)
      END DO

      ! icsuswi : switch which value equals 1 if ice melts at the surface
      !                                   0 if not
      icsuswi    =  max(0,int(-dhi_surf/max(zeps,abs(dhi_surf))))

      ! ice bottom behaviour : computation of icboind-icboswi
      ! icboind : index which values equals 0 if accretion is on the way
      !                                     1 if last layer has started to melt
      !                               2 if penultiem layer has started ...
      !                               etc ...
      icboind    =  0  
      deltah   =  0.d0
      DO layer = nlice, 1, -1
         IF (-dhi_bot .GT. deltah) THEN
            icboind  =  nlice+1-layer 
         ENDIF
         deltah = deltah + dzi(layer)
      END DO

      ! icboswi : switch which value equals 1 if accretion of ice is on the way
      !                                         0 if ablation is on the way
      icboswi    =  max(0,int(dhi_bot/max(zeps,abs(dhi_bot))))

      ! snow-ice formation : calcul de snicind-snicswi
      ! snicind : index which values equals 0 if no snow-ice forms
      !                1 if last layer of snow has started to melt
      !                          2 if penultiem layer ...
      snicind    =  0
      deltah   =  0.d0
      DO layer = nlsno, 1, -1
         IF ( dh_sni .GT. deltah ) THEN
            snicind  =  nlsno+1-layer
         ENDIF
         deltah = deltah + dzs(layer)
      END DO

      ! snicswi : switch which value equals 1 if snow-ice forms
      !                                     0 if not
      snicswi    =  max(0,int(dh_sni/max(zeps,abs(dh_sni))))

      IF ( ln_write ) THEN
         WRITE(numout,*) ' - Switches ... '
         WRITE(numout,*) '  snind   : ', snind
         WRITE(numout,*) '  snswi   : ', snswi
         WRITE(numout,*) '  icsuind : ', icsuind
         WRITE(numout,*) '  icsuswi : ', icsuswi
         WRITE(numout,*) '  icboind : ', icboind
         WRITE(numout,*) '  icboswi : ', icboswi
         WRITE(numout,*) '  snicind : ', snicind
         WRITE(numout,*) '  snicswi : ', snicswi
         WRITE(numout,*) 
      ENDIF

!------------------------------------------------------------------------------!

      !------------------------------!
      !                              !
      !     SNOW redistribution      ! 
      !                              !
      !------------------------------!

  if (ipsnow == 1 .and. hs > 1d-12 ) then
!
!------------------------------------------------------------------------------|
!      S-1) 'Old' snow cotes
!------------------------------------------------------------------------------|
!
      !
      ! by 'old', I mean that layers coming from accretion are included, and 
      ! that surface layers which were partly melted are reduced ...

      ! indexes of the vectors
      ntop0    =  1
      nbot0    =  nlsno+1- snind+(1-snicind)*snicswi

      ! cotes of the top of the layers
      zs(0)   =  0.d0
      DO layer = 1, nbot0-1
         limsum    =  snswi*(layer+snind-1) + (1-snswi)*(layer-1)
         limsum    =  MIN( limsum , nlsno )
         zs(layer) =  dhs
         DO layer_a = 1, limsum
            zs(layer) = zs(layer) + dzs(layer_a)
         END DO
      END DO

      zs(nbot0) =  dhs - snicswi*dh_sni
      DO layer = 1, nlsno
        zs(nbot0) = zs(nbot0) + dzs(layer)
      END DO
      zs(1)     =  dhs*(1-snswi) + snswi*zs(1)

      ! thicknesses
      DO layer = ntop0, nbot0
         zthick0(layer)  =  zs(layer) - zs(layer-1)
      END DO
!
!------------------------------------------------------------------------------|
!      S-2) 'Old' enthalpies
!------------------------------------------------------------------------------|
!
      ! enthalpies
      zqh_s0(1) =  rhosno * ( cp_ice * ( tp0 - (1 - snswi) * tair - & ! fallen snow
     &             snswi * ts(1) ) + mlfus )* zthick0(1)
      zqt_s_ini = zqt_s_ini + zqh_s0(1)

      DO layer = 2, nbot0 ! remaining snow
         limsum = (1-snswi)*(layer-1) + snswi*(layer+snind-1)
         limsum = MIN( limsum, nlsno )
         zqh_s0(layer) = rhosno*( cp_ice*(tp0 - ts(limsum)) + mlfus) * zthick0(layer)
         zswitch   = 1.d0 - MAX (0.d0, SIGN ( 1.d0, zeps - hs ) )
         zqt_s_ini = zqt_s_ini + ( 1 - snswi ) * zqh_s0(layer) * zswitch
      END DO
      zqt_ini = zqt_s_ini

      IF (ln_write) THEN
         WRITE(numout,*) ' - Snow redistribution ... '
         WRITE(numout,*) ' ntop0 : ', ntop0
         WRITE(numout,*) ' nbot0 : ', nbot0
         WRITE(numout,*) ' zs    : ', ( zs(layer), layer = 0, nbot0)
         WRITE(numout,*) ' zqh_s0: ', ( zqh_s0(layer), layer = 1, nbot0)
         WRITE(numout,*) ' q_s_b1: ', qs(1)
         WRITE(numout,*)
      ENDIF

      !-----------------------------------------
      ! snow enthalpy lost in sea ice formation
      !-----------------------------------------
      zqsnow     =  0.d0
      zdeltah    =  0.d0

      DO layer =  nlsno, 1, -1
         zhsnow    =  MAX(0.d0, dh_sni-zdeltah)
         zqsnow    =  zqsnow + qs(layer) * zhsnow
      END DO

!
!------------------------------------------------------------------------------|
!      S-3) New snow grid
!------------------------------------------------------------------------------|
!
      ! indexes of the vectors
      ntop1    =  1
      nbot1    =  nlsno

      ! cotes and thicknesses
         DO layer = 1, nlsno
            dzs(layer) = hs / dble(nlsno)
         END DO
         zs(1)  = dzs(1) / 2.d0
         DO layer = 2, nlsno
            zs(layer)  = zs(layer-1) + ( dzs(layer-1) + dzs(layer) ) / 2.d0
         END DO

         IF ( ln_write ) THEN
            WRITE(numout,*) ' dzs : ', ( dzs(layer), layer = 1, nlsno )
            WRITE(numout,*) ' zs  : ', (  zs(layer), layer = 1, nlsno )
         ENDIF
!
!------------------------------------------------------------------------------|
!      S-4) New snow enthalpy
!------------------------------------------------------------------------------|
!
      CALL ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , &
     &                    zthick0, dzs, zqh_s0 , zqh_s1 )

      zqt_s_fin = 0.d0 ! total heat content
!
!------------------------------------------------------------------------------|
!      S-5) Recover snow temperature
!------------------------------------------------------------------------------|
!
      DO layer = 1, nlsno
         ts(layer)  =  ipsnow * ( tp0 - ( zqh_s1(layer) / rhosno/ &
     &                       MAX( dzs(layer) , zeps ) - mlfus )  &
     &                       / cp_ice ) +                        &
     &                       ( 1.d0 - ipsnow ) * tp0
         zqt_s_fin = zqt_s_fin + zqh_s1(layer)
         qs(layer) = ipsnow * zqh_s1(layer) / MAX( dzs(layer) , zeps )
      END DO
      zqt_fin = zqt_s_fin

      IF (ln_write) THEN
         WRITE(numout,*) ' - Heat conservation, ice_th_remap , snow '
         WRITE(numout,*) ' zqt_s_fin : ', zqt_s_fin
         WRITE(numout,*) ' zqt_s_ini : ', zqt_s_ini
         WRITE(numout,*) ' TS1. zqt_s_fin - zqt_s_ini : ', zqt_s_fin - zqt_s_ini
         WRITE(numout,*)
      ENDIF

  else
!
!------------------------------------------------------------------------------|
!      S-3) New snow grid
!------------------------------------------------------------------------------|
!
      ! indexes of the vectors
      ntop1    =  1
      nbot1    =  nlsno

      ! cotes and thicknesses
         DO layer = 1, nlsno
            dzs(layer) = hs / dble(nlsno)
         END DO
         zs(1)  = dzs(1) / 2.d0
         DO layer = 2, nlsno
            zs(layer)  = zs(layer-1) + ( dzs(layer-1) + dzs(layer) ) / 2.d0
         END DO
         DO layer = 1, nlsno
            ts(layer) = tsu
         ENDDO

         IF ( ln_write ) THEN
            WRITE(numout,*) ' dzs : ', ( dzs(layer), layer = 1, nlsno )
            WRITE(numout,*) ' zs  : ', (  zs(layer), layer = 1, nlsno )
         ENDIF

  endif ! ipsnow == 1 FD
!------------------------------------------------------------------------------!
 
      !------------------------------!
      !                              !
      !     ICE redistribution       ! 
      !                              !
      !------------------------------!

!------------------------------------------------------------------------------!

!
!------------------------------------------------------------------------------!
!  I-1) Old ice grid                                                           !
!------------------------------------------------------------------------------!
!
      ! indexes of the vectors
      ntop0    =  1
      nbot0    =  MAX(1, MIN( nlice + 1 - icboind + (1-icsuind)*icsuswi &
     &                      + snicswi, nlice + 2 ) )

      ! cotes of the top of the layers
      zi(0)   =  0.d0
      DO layer = 1, nbot0-1
         limsum    =  ((icsuswi*(icsuind+layer-1) + (1-icsuswi)*layer))* &
     &               (1-snicswi) + (layer-1)*snicswi
         zi(layer)=  icsuswi*dhi_surf + snicswi*dh_sni
         DO layer_a = 1, limsum
            zi(layer) = zi(layer) + dzi(layer_a)
         END DO
      END DO

      zi(nbot0) =  icsuswi*dhi_surf + snicswi*dh_sni + dhi_bot
      DO layer = 1, nlice
         zi(nbot0) = zi(nbot0) + dzi(layer)
      END DO
      zi(1)     =  snicswi*dh_sni + (1-snicswi)*zi(1)

      ! thicknesses
      DO layer = ntop0, nbot0
         zthick0(layer)  =  zi(layer) - zi(layer-1)
      END DO

      IF (ln_write) THEN
         WRITE(numout,*) ' - Ice redistribution ... '
         WRITE(numout,*) ' ntop0 : ', ntop0
         WRITE(numout,*) ' nbot0 : ', nbot0
         WRITE(numout,*) ' zi    : ', ( zi(layer), layer = 0, nbot0)
         WRITE(numout,*) ' zthick0: ', ( zthick0(layer), layer = ntop0,nbot0)
      ENDIF
!
!------------------------------------------------------------------------------|
!  I-2) Old ice heat content                                                   |
!------------------------------------------------------------------------------|
!
      !-------------------------
      ! sources of heat content
      !-------------------------
      !- new ice
      tmelts = - fracsal * si_acc_new
      zqh_i_new =  rhoice * func_qm ( tmelts, tbo-tp0 ) * zthick0(nbot0)

      !- snow ice
      zqh_i_sni =  rhowat * cp_wat * ( tp0 - tbo ) * dh_sni * &
     &            ( rhoice - rhosno ) / rhoice * snicswi  ! generally positive
      fcons_sni =  - zqh_i_sni / dtice
      zqh_i_sni =  zqsnow + zqh_i_sni

      !----------------------
      ! sea ice heat content
      !----------------------
      zqh_i0(:) = 0.d0
      !- internal ice layers
      DO layer = ntop0, nbot0
          limsum = MAX( 1 , MIN( snicswi*(layer-1)      &
     &           + icsuswi * ( layer - 1 +  icsuind )   &
     &           + ( 1 - icsuswi ) * ( 1 - snicswi ) * layer, nlice ) )
          zqh_i0(layer) = qi(limsum) * zthick0(layer)
      END DO
        
      !- boundary values 
      zqh_i0(1)     =  snicswi * zqh_i_sni + ( 1 - snicswi ) * zqh_i0(1)
      zqh_i0(nbot0) =  ( 1 - icboswi ) * zqh_i0(nbot0) + icboswi * zqh_i_new

      IF (ln_write) THEN
         WRITE(numout,*) ' zqh_i_new: ', zqh_i_new
         WRITE(numout,*) ' zqh_i_sni: ', zqh_i_sni
         WRITE(numout,*) ' zqh_i0: ', ( zqh_i0(layer), layer = ntop0,nbot0)
      ENDIF

      DO layer = ntop0, nbot0
         zqt_i_ini = zqt_i_ini + zqh_i0(layer) 
      END DO

      zqt_ini = zqt_ini + zqt_i_ini ! includes zqsnic and zqsnow
! FD debug
! write(*,*) 'debug energy avant remap',-zqt_ini,-zqt_ini-sume1,dsume*dtice

!
!------------------------------------------------------------------------------|
!  I-3) Old ice salt content                                                   |
!------------------------------------------------------------------------------|
!
      ! basal new ice salt content
! FD      zsh_i_new = eskel * seasal * MAX ( dhi_bot, 0.d0 )
      zsh_i_new = si_acc_new * MAX ( dhi_bot, 0.d0 ) ! FD force bottom salinity
      ! snow ice salt content
      s_i_snic  =  ( rhoice - rhosno ) / rhoice * seasal * frac_sni ! snow ice salinity
      zsh_i_sni = s_i_snic * dh_sni

      IF ( ln_write ) THEN
         WRITE(numout,*) ' Salt sources : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~ '
         WRITE(numout,*) ' zsh_i_new : ', zsh_i_new
         WRITE(numout,*) ' zsh_i_sni : ', zsh_i_sni
      ENDIF

      zsh_i0(:) = 0.d0

      DO layer = ntop0, nbot0
          limsum =  snicswi*(layer-1) + icsuswi*(layer-1+icsuind) &
     &            + (1-icsuswi)*(1-snicswi)*layer
          limsum =  MAX(1,MIN(snicswi*(layer-1) + icsuswi*(layer-1 &
     &           +  icsuind)                                       &
     &           + (1-icsuswi)*(1-snicswi)*layer, nlice))
          zsh_i0(layer) = si(limsum)*zthick0(layer)
      END DO

      zsh_i0(nbot0) =  (1-icboswi)* zsh_i0(nbot0) + icboswi * zsh_i_new
      zsh_i0(1)     =  snicswi * zsh_i_sni + ( 1 - snicswi ) * zsh_i0(1)
     
      z_ms_i_ini = 0.d0
      DO layer = ntop0, nbot0
         z_ms_i_ini = z_ms_i_ini + zsh_i0(layer)
      END DO

      IF ( ln_write ) THEN
         WRITE(numout,*) ' frac_sni:', frac_sni
         WRITE(numout,*) ' seasal  : ', seasal
         WRITE(numout,*) ' s_i_snic : ', s_i_snic
         WRITE(numout,*) ' zsh_i0   : ', zsh_i0(ntop0:nbot0)
         WRITE(numout,*)
      ENDIF
!
!------------------------------------------------------------------------------|
!  I-4) New ice profile                                                        |
!------------------------------------------------------------------------------|
!
      ! indexes of the vectors
      ntop1    =  1 
      nbot1    =  nlice

      ! cotes and thicknesses
         DO layer = 1, nlice
            dzi(layer) = hi / dble(nlice)
         END DO
         zi(1)  = dzi(1) / 2.d0
         DO layer = 2, nlice
            zi(layer)  = zi(layer-1) + ( dzi(layer-1) + dzi(layer) ) / 2.d0
         END DO
         IF ( ln_write ) THEN
            WRITE(numout,*) ' dzi : ', ( dzi(layer), layer = 1, nlice )
            WRITE(numout,*) ' zi  : ', (  zi(layer), layer = 1, nlice )
         ENDIF
!
!------------------------------------------------------------------------------|
!  I-5) Redistribute ice heat content                                          |
!------------------------------------------------------------------------------|
!
      ! Remapping ice enthalpy
      CALL ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , &
     &                    zthick0, dzi, zqh_i0 , zqh_i1 )

      IF (ln_write) THEN
         WRITE(numout,*) ' zqh_i1 : ',( zqh_i1(layer), layer = ntop1,nbot1)
         WRITE(numout,*)
      ENDIF

      DO layer = ntop1, nbot1
         zqt_i_fin = zqt_i_fin + zqh_i1(layer) 
      END DO

      zqt_fin = zqt_fin + zqt_i_fin ! includes zqsnic and zqsnow
!
!------------------------------------------------------------------------------|
!  I-6) Redistribute ice salt content                                          |
!------------------------------------------------------------------------------|
!
      CALL ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , &
     &                    zthick0, dzi, zsh_i0 , zsh_i1 )

!
!------------------------------------------------------------------------------|
!  I-7) Heat conservation test                                                 |
!------------------------------------------------------------------------------|
!
      IF ( ln_write ) THEN
         WRITE(numout,*) ' - Heat conservation, ice_th_remap , sea ice... '
         WRITE(numout,*) ' zqt_i_fin : ', zqt_i_fin
         WRITE(numout,*) ' zqt_i_ini : ', zqt_i_ini
         WRITE(numout,*) ' TI1. zqt_i_fin - zqt_i_ini : ',zqt_i_fin - zqt_i_ini
         WRITE(numout,*)
         WRITE(numout,*) ' - Heat conservation, total ... '
         WRITE(numout,*) ' zqt_ini         : ', zqt_ini
         WRITE(numout,*) ' zqt_fin         : ', zqt_fin
         WRITE(numout,*) ' zqt_fin-zqt_ini : ', zqt_fin-zqt_ini
      ENDIF
!
!------------------------------------------------------------------------------|
!  I-8) Recover energy of melting, salinity and temperature                    |
!------------------------------------------------------------------------------|
!
      DO layer = 1, nlice
         qi(layer) = zqh_i1(layer) / MAX( dzi(layer) , zeps )
      END DO

      DO layer = 1, nlice
         si(layer) = zsh_i1(layer) / MAX( dzi(layer) , zeps )
      END DO

      DO layer = 1, nlice
         tmelts = -fracsal*si(layer) + tp0
         aaa = cp_ice
         bbb = (cp_wat-cp_ice)*(tmelts-tp0) + qi(layer)/ rhoice - mlfus 
         ccc = mlfus * (tmelts-tp0)
         discrim = SQRT( bbb*bbb - 4.d0*aaa*ccc )
         ti(layer) = tp0 + (- bbb - discrim) / ( 2.d0*aaa )
      END DO
!
!------------------------------------------------------------------------------|
!  I-9) Salt conservation test                                                 |
!------------------------------------------------------------------------------|
!
      z_ms_i_fin = 0.d0
      DO layer = 1, nlice
         z_ms_i_fin = z_ms_i_fin + si(layer) * dzi(layer)
      END DO
      IF ( ln_write ) WRITE(numout,*) ' z_ms_i_fin : ', z_ms_i_fin

      ! Flux due to bottom formatiophy
      zfgrml = zsh_i_new / dtice
      zfsigr = zsh_i_sni / dtice

      zfgrml = zfgrml + zfmelt
      IF ( ln_write ) THEN
         WRITE(numout,*) ' zfgrml : ', zfgrml
         WRITE(numout,*) ' zfsigr : ', zfsigr
      ENDIF

      zerror = 1.0d-9
      CALL ice_sal_conserv(1,1,'ice_phy_remap: ',zerror,     &
     &                           z_ms_i_ini,z_ms_i_fin,  &
     &                           zfgrml, zfsigr, dtice)
!      si(:)=sinew(:) ! FD another trick to ensure that salinity remains constant

! FD compute final energy
      sume2 = 0.d0
      DO layer = 1, nlice
         tmelts = -fracsal*si(layer)
         sume2 = sume2 + rhoice * func_El(tmelts,ti(layer) - tp0) * dzi(layer)
      END DO
      DO layer = 1, nlsno
         tmelts = 0.d0
         sume2 = sume2 + rhosno * func_El(tmelts,ts(layer) - tp0) * dzs(layer) * ipsnow ! FD
      END DO
! FD debug energy
      write(*,*) 'debug energy',dsume * dtice, sume2 - sume1, sume2
      if ( abs(dsume*dtice-sume2+sume1) .gt. 1d4 .or. ts(1) .gt. tp0 ) then
          call flush(numout)
          call flush(6)
       OPEN(1,file='restart.dat')
       WRITE(1,*) nlice,nlsno,ith_cond
       WRITE(1,*) ti0(1:nlice),ts0(1:nlsno),tsu0,tbo
       WRITE(1,*) si0(1:nlice),seasal
       WRITE(1,*) hi0,hs0
       WRITE(1,*) dzi0(1:nlice),dzs0(1:nlsno)
       WRITE(1,*) snowfall0,dwnlw,tsu0,tair,qair,uair,swrad,oceflx,pres
       WRITE(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
       CLOSE(1)
          stop 'energy cons pb'
      endif

!-------------------------------------------------------------------------------
! Fin de la routine ice_th_remap
      END SUBROUTINE ice_thermo_remap

end module ice_thermo_lim
