      SUBROUTINE ice_rad(nlay_s,nlay_i,kideb,kiut,numofday)

!-----------------------------------------------------------------------------!
!     This routine computes absorption / transmission of radiation through
!     the snow ice system
!     (c) Martouf, UCL-ASTR, June 2007. Nadal Federer 1 set partout
!-----------------------------------------------------------------------------!

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'
      INCLUDE 'bio.com'

!-----------------------------------------------------------------------------!

! Local variables
      REAL(8), DIMENSION(maxnlay) ::  
     &   zkappa_alg , !: extinction radiation coefficient for algae 
     &   zkappa_det   !: extinction radiation coefficient for detritus 

!==============================================================================!
 
      WRITE(numout,*) ' ** ice_rad : '
      WRITE(numout,*) ' ~~~~~~~~~~~~ '

      DO ji = kideb, kiut
!
!------------------------------------------------------------------------------!
!  1) Radiation transmitted below the surface                                  !
!------------------------------------------------------------------------------!
!
!     WRITE(numout,*)
!     WRITE(numout,*) ' Radiation below the surface '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
      !-------------
      ! Snow switch
      !-------------
      ! Is there snow or not ?
      isnow   = INT( 1.0 - MAX( 0.0 , SIGN( 1.0 , - ht_s_b(ji) ) ) )

      !----------------------------
      ! Surface transmissivity i_0 
      !----------------------------
      ! ab = 1 - i_0 
      ! Grenfell and Maykut (1977)
      ab(ji) = 1.0 - ( 1.0 - isnow ) * rad_inot_i - isnow * rad_inot_s
!     ab(ji) = 1.0 - ( 1.0 - isnow ) * 0.30 - isnow * 0.15
!     ab(ji) = 1.0 - ( 1.0 - isnow ) * 0.30 - isnow * 0.
      ! Bitz and Lipscomb (1999)
!     ab(ji) = 1.0 - 0.30 / ( 10.0 * ht_s_b(ji) + 1.0 ) 

      !-----------------------------------------------------
      ! Solar radiation transmitted below the surface layer
      !-----------------------------------------------------
      ftrice      =  fsolgb(ji) * ( 1.0 - ab(ji) )
!
!------------------------------------------------------------------------------!
!  2) Extinction coefficients
!------------------------------------------------------------------------------!
!
!     WRITE(numout,*)
!     WRITE(numout,*) ' Extinction coefficients '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~'
!     WRITE(numout,*)

      zmurad = 0.656 ! angle factor

      zchla = 0.0        ! temporary value of chla in mg/m-3

      ! Chlorophyll a interpolation
! FD      CALL ice_bio_interp_bio2phy(kideb,kiut,nlay_i,.FALSE.)
! FD
! FD      WRITE(numout,*) ' chla_i : ', ( chla_i(layer), layer = 1, nlay_i )
! FD      WRITE(numout,*)
! FD modification
      chla_i(:)=0.0;

      DO layer = 1, nlay_i
         zchla = chla_i(layer)
         zkappa_alg(layer) = zchla * astar_alg / zmurad
         zkappa_det(layer) = zkappa_alg(layer) * fdet_alg
      END DO

!     WRITE(numout,*) ' zkappa_alg : ', ( zkappa_alg(layer), 
!    &                layer = 1, nlay_i )
!     WRITE(numout,*) ' zkappa_det : ', ( zkappa_det(layer), 
!    &                layer = 1, nlay_i )
!
!------------------------------------------------------------------------------!
!  3) Radiation transmitted through / absorbed by snow                         !
!------------------------------------------------------------------------------!
!
!     WRITE(numout,*)
!     WRITE(numout,*) ' Radiation absorption in snow '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
!     WRITE(numout,*)

      radtr_s(0) =  ftrice
      zzs = 0.0
      DO layer = 1, nlay_s
         zzs = zzs + deltaz_s_phy(layer)
         radtr_s(layer) = radtr_s(0) * exp( - rad_kappa_s*( MAX( 0.0 ,
     &                    zzs ) ) )
         radab_s(layer) = radtr_s(layer-1) - radtr_s(layer)
      END DO
!     WRITE(numout,*) ' radtr_s : ', 
!    &                ( radtr_s(layer) , layer = 0, nlay_s )
!     WRITE(numout,*) ' radab_s : ', 
!    &                ( radab_s(layer) , layer = 1, nlay_s )
!
!------------------------------------------------------------------------------!
!  4) Radiation transmitted through / absorbed by ice                          !
!------------------------------------------------------------------------------!
!
      WRITE(numout,*)
      WRITE(numout,*) ' Radiation absorption in ice '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*)
      ! transmitted at the upper ice layer
      radtr_i(0) = isnow * radtr_s(nlay_s)  + ( 1. - isnow ) * ftrice

      DO layer = 1, nlay_i
         ! extinction coefficient
         zraext_i   = rad_kappa_i + zkappa_alg(layer) 
     &              + zkappa_det(layer)

         ! layer thickness
         zdeltaz_i = deltaz_i_phy(layer)
         zdummy    = radtr_i(layer-1) * 
     &               EXP ( - zraext_i * deltaz_i_phy(layer) ) *
     &               deltaz_i_phy(layer)

         ! physicallay absorbed
         radab_phy_i(layer) = ( rad_kappa_i + zkappa_det(layer) )
     &                      * zdummy 
         radab_alg_i(layer) = zkappa_alg(layer) * zdummy 
         radtr_i(layer)     = radtr_i(layer-1) - radab_phy_i(layer)
     &                                         - radab_alg_i(layer)
         ! par
         par(layer)         = ( radtr_i(layer-1) + radtr_i(layer) ) 
     &                      / ( 2. * par_fsw )    ! average par in the layer
      END DO

      ! radiation sent to the ocean
      ftroce = radtr_i(nlay_i)
       
      WRITE(numout,*) ' i0     : ', 1.0-ab(ji)
      WRITE(numout,*) ' ftrice : ', ftrice
      WRITE(numout,*) ' ftroce : ', ftroce

!     WRITE(numout,*) ' radtr_i : ', 
!    &                ( radtr_i(layer) , layer = 0, nlay_i )
!     WRITE(numout,*) ' radab_phy_i : ', 
!    &                ( radab_phy_i(layer) , layer = 1, nlay_i )
!     WRITE(numout,*) ' radab_alg_i : ', 
!    &                ( radab_alg_i(layer) , layer = 1, nlay_i )

!     !----------------------
!     ! Additional radiation
!     !----------------------
!     zconc = 0.9
!     zkappaw = 0.5
!     zi0w   = 0.80
!     zincsw = fsolgb(ji) / ( 1.0 - albgb(ji) ) * zi0w

!     DO layer = 1, nlay_i
!        zdeltaz_i = z_i(layer) - z_i(layer-1)
!        zdummy    = zincsw * EXP ( - zkappaw * zdeltaz_i ) *
!    &               zdeltaz_i
!        zradal = ( 1.0 - zconc ) / zconc * zdummy
!        zindcsw = zindcsw - zradal

!        zdummy  = rkappa_i_phy + zkappa_det(layer) + zkappa_alg(layer)
!        radab_phy_i(layer) = radab_phy_i(layer) + zradal*
!    &                      ( rkappa_i_phy + zkappa_det(layer) ) /
!    &                        zdummy
!        radab_alg_i(layer) = radab_alg_i(layer) + zradal*
!    &                        zkappa_alg(layer) / zdummy

!     END DO

      !--------------------
      ! Conservation check
      !--------------------
! FD !!! here implicit assumption that nlay_s=1
! FD      sumrad = radab_s(1)
      sumrad = 0.d0
      DO layer = 1, nlay_s
         sumrad = sumrad + radab_s(layer)
      END DO
      DO layer = 1, nlay_i
         sumrad = sumrad + radab_phy_i(layer) + radab_alg_i(layer)
      END DO
      WRITE(numout,*) ' Conservation check '
      WRITE(numout,*) ' ftrice - ftroce : ', ftrice-ftroce
      WRITE(numout,*) ' sumrad          : ', sumrad

      END DO !ji

      WRITE(numout,*)
!     WRITE(numout,*) ' End of ice_rad '
!     WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!==============================================================================!
! end of the subroutine

      END SUBROUTINE
