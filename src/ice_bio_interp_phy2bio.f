      SUBROUTINE ice_bio_interp_phy2bio(kideb,kiut,nlay_i,ln_write)

! This routine interpolates salinity, temperature, brine salinity, brine volume
! on the biological grid
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

      REAL(8), DIMENSION( 0:maxnlay ) ::
     &  z0

      REAL(8), DIMENSION( 0:nlay_bio ) ::
     &  z1

      REAL(8), DIMENSION( maxnlay ) ::
     &  zqs         , ! : scalar content on the physical grid (input)
     &  zqt         , ! : scalar content on the physical grid (input)
     &  zqr         , ! : scalar content on the physical grid (input)
     &  zqpar       , ! : scalar content on the physical grid (input)
     &  zthick0       ! : thickness of biological layers

      REAL(8), DIMENSION( nlay_bio ) ::
     &  zq1         , ! : scalar content on the biological grid (output)
     &  zthick1       ! : thickness of physical layers

      REAL(8), DIMENSION( nlay_bio , maxnlay ) ::
     &  zweight       ! : relayering matrix

      REAL(8) ::    
     &  zaaa        , ! : dummyfactors for the computation of t_i_bio
     &  zbbb        ,
     &  zccc        ,
     &  zdiscrim    ,
     &  zsum0         ! : conservation test variable
     &  zsum1         ! : conservation test variable
! FD additions
      LOGICAL ln_write

!=============================================================================!

      IF ( ln_write ) THEN
         WRITE(numout,*) ' *** ice_bio_interp_phy2bio : '
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
      z0(0) = 0.0
      DO layer = 1, nlay_i
         z0(layer) = ht_i_b(ji) / nlay_i * layer
         zthick0(layer) = z0(layer) - z0(layer-1)
      END DO ! layer

!     !+++++
!     WRITE(numout,*) ' z0      : ', ( z0(layer0) , 
!    &                layer0 = 0, nlay_i )
!     WRITE(numout,*) ' zthick0 : ', ( zthick0(layer0) , 
!    &                layer0 = 1, nlay_i )
!     !+++++

      !-----------------
      ! Biological grid
      !-----------------
      z1(0) = 0.0
      DO layer = 1, nlay_bio
         z1(layer) = z1(layer-1) + deltaz_i_bio(layer)
         zthick1(layer) = z1(layer) - z1(layer-1)
      END DO ! layer

!     !+++++
!     WRITE(numout,*) ' z1      : ', ( z1(layer1) , 
!    &                layer1 = 1, nlay_bio )
!     WRITE(numout,*) ' zthick1 : ', ( zthick1(layer1) , 
!    &                layer1 = 1, nlay_bio )
!     !+++++
!
!-----------------------------------------------------------------------------!
! 2) Scalar contents
!-----------------------------------------------------------------------------!
!
      DO layer = 1, nlay_i
         zqs(layer)   = s_i_b(ji,layer) * zthick0(layer)
         zqt(layer)   = q_i_b(ji,layer) * zthick0(layer)
         zqr(layer)   = radab_alg_i(layer) ! no thickness ( for unit reasons )
         zqpar(layer) = par(layer) ! average par over the layer no thickness ( for unit reasons )
      END DO ! layer

!     !+++++
!     WRITE(numout,*) ' s_i_b   : ', ( s_i_b(ji,layer1) , 
!    &                layer1 = 1, nlay_i )
!     WRITE(numout,*) ' q_i_b   : ', ( q_i_b(ji,layer1) , 
!    &                layer1 = 1, nlay_i )
!     WRITE(numout,*) ' t_i_b   : ', ( t_i_b(ji,layer1) , 
!    &                layer1 = 1, nlay_i )
!     WRITE(numout,*) ' zqs     : ', ( zqs(layer1) , 
!    &                layer1 = 1, nlay_i )
!     WRITE(numout,*) ' zqt     : ', ( zqt(layer1) , 
!    &                layer1 = 1, nlay_i )
!     WRITE(numout,*) ' zqr     : ', ( zqr(layer1) , 
!    &                layer1 = 1, nlay_i )
!     !+++++
!
!-----------------------------------------------------------------------------!
! 3) Weights 
!-----------------------------------------------------------------------------!
!
      ! weights of old layers on new ones
      DO layer1 = 1, nlay_bio
         DO layer0 = 1, nlay_i
            zweight(layer1,layer0) = MAX ( 0.0 , ( MIN ( z0(layer0) ,
     &      z1(layer1) ) - MAX ( z0 (layer0-1) , z1(layer1-1) ) ) / 
     &      zthick0(layer0) )
         END DO
      END DO
!
!-----------------------------------------------------------------------------!
! 4) Interpolation
!-----------------------------------------------------------------------------!
!
      !--------------
      ! Ice salinity
      !--------------
      DO layer1 = 1, nlay_bio
         zq1(layer1) = 0.0
         DO layer0 = 1, nlay_i
            zq1(layer1) = zq1(layer1) + zweight(layer1,layer0) *
     &                    zqs(layer0)
         END DO
      END DO

!     !+++++
!     WRITE(numout,*) ' Salinity '
!     WRITE(numout,*) ' zq1     : ', ( zq1(layer1) , 
!    &                layer1 = 1, nlay_bio )
!     !+++++

      DO layer1 = 1, nlay_bio
         s_i_bio(layer1) = zq1(layer1) / zthick1(layer1)
      END DO

      IF ( ln_write ) THEN
         WRITE(numout,*) ' s_i_bio : ', ( s_i_bio(layer1) , 
     &                   layer1 = 1, nlay_bio )
      ENDIF

      !---------------------------------
      ! Biologically absorbed radiation
      !---------------------------------
      DO layer1 = 1, nlay_bio
         zq1(layer1) = 0.0
         DO layer0 = 1, nlay_i
            zq1(layer1) = zq1(layer1) + zweight(layer1,layer0) *
     &                    zqr(layer0)
         END DO
      END DO
      ! Conservation test
      zsum0 = 0.0
      DO layer0 = 1, nlay_i
         zsum0 = zsum0 + zqr(layer0)
      END DO

      zsum1 = 0.0
      DO layer1 = 1, nlay_bio
         zsum1 = zsum1 + zq1(layer1)
      END DO

      IF ( ABS( zsum0 - zsum1 ).GT. 1.0e-10 ) THEN
         WRITE(numout,*) ' ALERTE, Radiation not conserved '
         WRITE(numout,*) ' Step 1 '
         WRITE(numout,*) ' zsum0 : ', zsum0
         WRITE(numout,*) ' zsum1 : ', zsum1
      ENDIF

      ! Correction to avoid numerical diffusion
      ! if there is no chla, no absorption
      zdummy = 0.0
      DO layer = 1, nlay_bio
         IF ( chla_i_bio(layer) .LT. 1.0e-11 ) THEN
            zdummy = zdummy + zq1(layer)
            zq1(layer) = 0.0
         ENDIF
      END DO
      zq1(nlay_bio) = zq1(nlay_bio) + zdummy

      zsum1 = 0.0
      DO layer1 = 1, nlay_bio
         zsum1 = zsum1 + zq1(layer1)
      END DO

      IF ( ABS( zsum0 - zsum1 ).GT. 1.0e-10 ) THEN
         WRITE(numout,*) ' ALERTE, Radiation not conserved '
         WRITE(numout,*) ' Step 2 '
         WRITE(numout,*) ' zsum0 : ', zsum0
         WRITE(numout,*) ' zsum1 : ', zsum1
      ENDIF
         
      DO layer1 = 1, nlay_bio
         radab_alg_i_bio(layer1) = zq1(layer1) ! no thickness redis
      END DO

      !+++++
      WRITE(numout,*) ' radab_alg_i_bio : ', ( radab_alg_i_bio(layer1) ,
     &                layer1 = 1, nlay_bio )
      !+++++

      ! Photosynthetically usable radiation
      ! convert Wm-2 into mmol quanta m-2 s-1
      ! conversion describes how many quanta are contained in Fsw
      DO layer = 1, nlay_bio
         pur_bio(layer) = radab_alg_i_bio(layer) * 2.147
      END DO

      !+++++
      IF ( ln_write ) THEN
         WRITE(numout,*) ' pur_bio : ', ( pur_bio(layer1) ,
     &                   layer1 = 1, nlay_bio )
      ENDIF
      !+++++

      !-----
      ! PAR
      !-----
      ! probably not very well interpolated
      DO layer1 = 1, nlay_bio
         zq1(layer1) = 0.0
         DO layer0 = 1, nlay_i
            zq1(layer1) = zq1(layer1) + zweight(layer1,layer0) *
     &                    zqpar(layer0)
         END DO
      END DO
      DO layer1 = 1, nlay_bio
         par_bio(layer1) = zq1(layer1) ! no thickness redis
      END DO

      !--------------
      ! Heat content
      !--------------
      DO layer1 = 1, nlay_bio
         zq1(layer1) = 0.0
         DO layer0 = 1, nlay_i
            zq1(layer1) = zq1(layer1) + zweight(layer1,layer0) *
     &                    zqt(layer0)
         END DO
      END DO

!     !+++++
!     WRITE(numout,*) ' Heat content '
!     WRITE(numout,*) ' zq1     : ', ( zq1(layer1) , 
!    &                layer1 = 1, nlay_bio )
!     !+++++

      ! energy of melting
      DO layer1 = 1, nlay_bio
         zq1(layer1) = zq1(layer1) / zthick1(layer1)
      END DO

      ! invert energy of melting to get temperature back
      DO layer1 = 1, nlay_bio
         tmelts = - tmut * s_i_bio(layer1) + tpw
         zaaa = cpg
         zbbb = (cpw-cpg)*(tmelts-tpw) + zq1(layer1) / rhog
     &        - xlgm 
         zccc = xlgm * (tmelts-tpw)
         zdiscrim = SQRT( zbbb*zbbb - 4.0*zaaa*zccc )
         t_i_bio(layer1) = tpw + ( - zbbb - zdiscrim ) / (2.0*zaaa)
      END DO

      IF ( ln_write ) THEN
         WRITE(numout,*) ' t_i_bio : ', ( t_i_bio(layer1) , 
     &                   layer1 = 1, nlay_bio )
      ENDIF

      !--------------
      ! Brine volume
      !--------------
      DO layer1 = 1, nlay_bio
         e_i_bio(layer1) = - tmut * s_i_bio(layer1) /
     &                     ( t_i_bio(layer1) - tpw )
      END DO ! layer1

      IF ( ln_write ) THEN
         WRITE(numout,*) ' e_i_bio : ', ( e_i_bio(layer1) , 
     &                   layer1 = 1, nlay_bio )
      ENDIF

      END DO ! ji

!=============================================================================!
!-- End of ice_bio_interp_phy2bio --
 
      END
!
!=============================================================================!
!=============================================================================!
!

      SUBROUTINE ice_bio_interp_diffus(kideb,kiut,nlay_i,ln_write)

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
!    &  jn          , ! : index for tracers
     &  layer1      , ! : relayering index
     &  layer2        ! : relayering index

      REAL(8), DIMENSION( 0:maxnlay ) :: ! lower interface of the layer
     &  zz_phy

      REAL(8), DIMENSION( 0:nlay_bio ) :: ! lower interface of the layer
     &  zz_bio
! FD additions
      LOGICAL ln_write

!     REAL(8), DIMENSION( 0:nlay_bio ) ::
!    &  z1

!     REAL(8), DIMENSION( maxnlay ) ::
!    &  zqs         , ! : scalar content on the physical grid (input)
!    &  zqt         , ! : scalar content on the physical grid (input)
!    &  zqr         , ! : scalar content on the physical grid (input)
!    &  zqpar       , ! : scalar content on the physical grid (input)
!    &  zthick0       ! : thickness of biological layers

!     REAL(8), DIMENSION( nlay_bio ) ::
!    &  zq1         , ! : scalar content on the biological grid (output)

!=============================================================================!

      IF ( ln_write ) THEN
         WRITE(numout,*) ' *** ice_bio_interp_diffus : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
      ENDIF
!
!-----------------------------------------------------------------------------!
! 1) Grids
!-----------------------------------------------------------------------------!
!
! compute the coordinates of the interfaces of the layers
!
      zz_phy(0) = 0.
      DO layer = 1, nlay_i
         zz_phy(layer) = z_i_phy(layer) + deltaz_i_phy(layer) / 2.
      END DO

      zz_bio(0) = 0.
      DO layer = 1, nlay_bio
         zz_bio(layer) = z_i_bio(layer) + deltaz_i_bio(layer) / 2.
      END DO

      IF ( ln_write ) THEN
         WRITE(numout,*) ' zz_phy : ', ( zz_phy(layer), 
     &                                   layer = 0, nlay_i )
         WRITE(numout,*) ' zz_bio : ', ( zz_bio(layer), 
     &                                   layer = 0, nlay_bio )
      ENDIF

      DO layer_bio = 1, nlay_bio - 1
         zdist_max = 999.9
         zdist = zdist_max
         WRITE(numout,*) ' '
         WRITE(numout,*) ' layer_bio : ', layer_bio
         DO layer_phy = 1, nlay_i
            zdist = MIN ( zdist, zz_bio(layer_bio) - zz_phy(layer_phy) )
            IF ( ( zdist .GE. 0.0 ) .AND. ( zdist .LT. zdist_max ) ) 
     &      THEN
               index_mem = layer_phy
            ENDIF
            WRITE(numout,*) ' layer_phy : ', layer_phy
            WRITE(numout,*) ' zdist : ', zdist
            WRITE(numout,*) ' index_mem ', index_mem
         END DO ! layer_phy
         index_mem = MAX ( MIN( index_mem, nlay_i ) , 1 ) ! prevent absurd values sometimes reached in path cases
         zdummy1 = ( diff_br(index_mem+1) - diff_br(index_mem) ) /
     &             ( zz_phy(index_mem+1) - zz_phy(index_mem) )
         zdummy2 = zz_bio(layer_bio) - zz_phy(index_mem)

         diff_br_bio(layer_bio) = diff_br(index_mem) + zdummy1*zdummy2

      END DO ! layer_bio

      diff_br_bio(nlay_bio) = diff_br(nlay_i)

      IF ( ln_write ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' diff_br : ', ( diff_br(layer), layer = 1, 
     &                                    nlay_i ) 
         WRITE(numout,*) ' diff_br_bio : ', ( diff_br_bio(layer), 
     &                   layer = 1, nlay_bio ) 
      ENDIF
!      
!=============================================================================!
!-- End of ice_bio_interp_diff --
!
      END
