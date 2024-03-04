      SUBROUTINE ice_bio_interp_bio2phy(kideb,kiut,nlay_i,ln_write)

! This routine interpolates chlorophyll a
! from the biological grid to the physical grid
! (c) Martin Vancoppenolle, May 2007
 
      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'
      INCLUDE 'bio.com'

      INTEGER :: 
     &  ji          , ! : index for space
     &  jk          , ! : index for ice layers
     &  jn          , ! : index for tracers
     &  layer1      , ! : relayering index
     &  layer2        ! : relayering index

      REAL(8), DIMENSION( 0:nlay_bio ) ::
     &  z0

      REAL(8), DIMENSION( 0:maxnlay ) ::
     &  z1

      REAL(8), DIMENSION( nlay_bio ) ::
     &  zqc         , ! : scalar content on the physical grid (input)
     &  zthick0       ! : thickness of biological layers

      REAL(8), DIMENSION( maxnlay ) ::
     &  zq1         , ! : scalar content on the biological grid (output)
     &  zthick1       ! : thickness of physical layers

      REAL(8), DIMENSION( maxnlay , nlay_bio ) ::
     &  zweight       ! : relayering matrix
! FD additions
      LOGICAL ln_write

!=============================================================================!

      IF ( ln_write ) THEN
         WRITE(numout,*) ' *** ice_bio_interp_bio2phy : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
      ENDIF
!
!-----------------------------------------------------------------------------!
! 1) Grids
!-----------------------------------------------------------------------------!
!
! compute the coordinates of the interfaces of the layers
      
      DO ji = kideb, kiut
      !---------------
      ! Physical grid
      !---------------
      z1(0) = 0.0
      DO layer = 1, nlay_i
         z1(layer) = ht_i_b(ji) / nlay_i * layer
         zthick1(layer) = z1(layer) - z1(layer-1)
      END DO ! layer

!     !+++++
!     WRITE(numout,*) ' z1      : ', ( z1(layer1) , 
!    &                layer1 = 0, nlay_i )
!     WRITE(numout,*) ' zthick1 : ', ( zthick1(layer1) ,
!    &                layer1 = 1, nlay_i )
!     !+++++

      !-----------------
      ! Biological grid
      !-----------------
      z0(0) = 0.0
      DO layer = 1, nlay_bio
         z0(layer) = z0(layer-1) + deltaz_i_bio(layer)
         zthick0(layer) = z0(layer) - z0(layer-1)
      END DO ! layer

!     !+++++
!     WRITE(numout,*) ' z0      : ', ( z0(layer1) , 
!    &                layer1 = 1, nlay_bio )
!     WRITE(numout,*) ' zthick0 : ', ( zthick0(layer1) , 
!    &                layer1 = 1, nlay_bio )
!     !+++++
!
!-----------------------------------------------------------------------------!
! 2) Scalar contents
!-----------------------------------------------------------------------------!
!
      DO layer = 1, nlay_bio
         zqc(layer) = chla_i_bio(layer) * zthick0(layer)
      END DO ! layer

!     !+++++
!     WRITE(numout,*) ' chla_i_bio :', ( chla_i_bio(layer) , 
!    &                layer = 1, nlay_bio )
!     WRITE(numout,*) ' zqc     : ', ( zqc(layer) , 
!    &                layer = 1, nlay_bio )
!     !+++++
!
!-----------------------------------------------------------------------------!
! 3) Weights 
!-----------------------------------------------------------------------------!
!
      ! weights of old layers on new ones
      DO layer1 = 1, nlay_i
         DO layer0 = 1, nlay_bio
            zweight(layer1,layer0) = MAX ( 0.0 , ( MIN ( z0(layer0) ,
     &      z1(layer1) ) - MAX ( z0 (layer0-1) , z1(layer1-1) ) ) / 
     &      zthick0(layer0) )
!           WRITE(numout,*) ' zweight : ', layer1, layer0, 
!    &                      zweight(layer1, layer0)
         END DO
      END DO
!
!-----------------------------------------------------------------------------!
! 4) Interpolation
!-----------------------------------------------------------------------------!
!
      !---------------
      ! Chlorophyll a
      !---------------
      DO layer1 = 1, nlay_i
         zq1(layer1) = 0.0
         DO layer0 = 1, nlay_bio
            zq1(layer1) = zq1(layer1) + zweight(layer1,layer0) *
     &                    zqc(layer0)
         END DO
      END DO
!     !+++++
!     WRITE(numout,*) ' Chlorophyll a : '
!     WRITE(numout,*) ' zq1     : ', ( zq1(layer1) , 
!    &                layer1 = 1, nlay_i )
!     !+++++

      DO layer1 = 1, nlay_i
         chla_i(layer1) = zq1(layer1) / zthick1(layer1)
      END DO

      END DO ! ji

!     WRITE(numout,*) ' chla_i  : ', ( chla_i(layer1) , 
!    &                layer1 = 1, nlay_i )

!     WRITE(numout,*) ' end '
!     WRITE(numout,*) ' ------------------------------------ '
!=============================================================================!
!-- End of ice_bio_interp_bio2phy --
 
      END
