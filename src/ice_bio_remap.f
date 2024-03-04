      SUBROUTINE ice_bio_remap(nlay_s,nlay_i,kideb,kiut)
!-----------------------------------------------------------------------------!
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'thermo.com'
      include 'bio.com'

      REAL(8), DIMENSION( 0:maxnlay+2 ) ::
     &  z0

      REAL(8), DIMENSION( 0:maxnlay ) ::
     &  z1

      REAL(8), DIMENSION( maxnlay + 2 ) ::
     &  zq0         , ! : scalar content on the physical grid (input)
     &  zthick0       ! : thickness of biological layers

      REAL(8), DIMENSION( maxnlay ) ::
     &  zq1         , ! : scalar content on the biological grid (output)
     &  zthick1       ! : thickness of physical layers

      REAL(8), DIMENSION( nlay_bio , maxnlay ) ::
     &  zweight       ! : relayering matrix

!==============================================================================!

      WRITE(numout,*) 
      WRITE(numout,*) ' ice_bio_remap : '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~ '

      zeps   = 1.0e-20

!     cbu_i_bio(:,:) = cbun_i_bio(:,:)

      !--------------------
      ! Conservation check
      !--------------------
      CALL ice_bio_column(kideb,kiut,mt_i_bio_init,cbu_i_bio,
     &                    deltaz_i_bio, .FALSE.)

      DO ji = kideb, kiut
      WRITE(numout,*)
      WRITE(numout,*) ' Inputs     : '
      WRITE(numout,*) ' ~~~~~~       '
      WRITE(numout,*) ' dh_i_surf  : ', dh_i_surf(ji)
      WRITE(numout,*) ' dh_i_bott  : ', dh_i_bott(ji)
      WRITE(numout,*) ' dh_snowice : ', dh_snowice(ji)
!
!------------------------------------------------------------------------------|
!  1) Switches                                                                 |
!------------------------------------------------------------------------------|
!
      WRITE(numout,*)
      WRITE(numout,*) ' Switches : '
      WRITE(numout,*) ' ~~~~~~~~ '
      !----------------------------------
      ! Surface melt indicator : icsuind
      !----------------------------------
      ! ice surface behaviour : computation of icsuind-icsuswi
      ! icsuind : index which values equals 
      !           0 if there is nothing
      !           1 if first layer has started to melt
      !           2 if first and second layer have started ...
      !                             etc ...
      icsuind  =  0  
      deltah   =  0.0
      DO layer = 1, nlay_bio
         IF ( - dh_i_surf(ji) .GT. deltah ) THEN
            icsuind  =  layer 
         ENDIF
         deltah = deltah + deltaz_i_bio(layer)
      END DO

      !-------------------------------
      ! Surface melt switch : icsuswi
      !-------------------------------
      ! icsuswi : switch which value equals 1 if ice melts at the surface
      !                                     0 if not
      icsuswi    =  MAX( 0 , INT( - dh_i_surf(ji) 
     &              / MAX ( zeps , ABS( dh_i_surf(ji) ) ) ) )

      !--------------------------------
      ! Ice bottom indicator : icboind
      !--------------------------------
      ! icboind : index which values equals 0 if accretion is on the way
      !                                     1 if last layer has started to melt
      !                                     2 if penultiem layer has started ...
      !                                     etc ...
      icboind  =  0  
      deltah   =  0.0
      DO layer = nlay_bio, 1, -1
         IF ( - dh_i_bott(ji) .GT. deltah ) THEN
            icboind  =  nlay_bio + 1 - layer 
         ENDIF
         deltah = deltah + deltaz_i_bio(layer)
      END DO

      !-----------------------------
      ! Ice bottom switch : icboswi
      !-----------------------------
      ! icboswi : switch which value equals 1 if accretion of ice is on the way
      !                                     0 if ablation is on the way
      icboswi    =  MAX( 0 , INT( dh_i_bott(ji)
     &              / MAX( zeps , ABS( dh_i_bott(ji) ) ) ) )

      !------------------------------
      ! Snow ice indicator : snicind
      !------------------------------
      ! snow-ice formation : calcul de snicind-snicswi
      ! snicind : index which values equals 0 if no snow-ice forms
      !                1 if last layer of snow has started to melt
      !                2 if penultiem layer ...
      snicind  =  0
      deltah   =  0.0
      DO layer = nlay_s, 1, -1
         IF ( dh_snowice(ji) .GT. deltah) THEN
            snicind  =  nlay_s+1-layer
         ENDIF
         deltah = deltah + h_s(layer)
      END DO

      !---------------------------
      ! Snow ice switch : snicswi
      !---------------------------
      ! snicswi : switch which value equals 1 if snow-ice forms
      !                                     0 if not
      snicswi    =  MAX( 0 , INT( dh_snowice(ji) /
     &              MAX( zeps , ABS( dh_snowice(ji) ) ) ) )

      WRITE(numout,*) ' icsuind    : ', icsuind
      WRITE(numout,*) ' icsuswi    : ', icsuswi
      WRITE(numout,*) ' icboind    : ', icboind
      WRITE(numout,*) ' icboswi    : ', icboswi
      WRITE(numout,*) ' snicind    : ', snicind
      WRITE(numout,*) ' snicswi    : ', snicswi
      WRITE(numout,*) 
!
!------------------------------------------------------------------------------|
! 2) Old grid   
!------------------------------------------------------------------------------|
!
      WRITE(numout,*) ' Old grid : '
      WRITE(numout,*) ' ~~~~~~~~~~ '
      ! indexes of the vectors
      ntop0    =  1
      nbot0    =  nlay_bio + 1 - icboind + (1 - icsuind) * icsuswi +
     &            snicswi
      nbot0    =  MAX(1, MIN( nlay_bio + 1 - icboind + 
     &            (1-icsuind)*icsuswi
     &          + snicswi, nlay_bio + 2 ) )

      WRITE(numout,*) ' ntop0     : ', ntop0
      WRITE(numout,*) ' nbot0     : ', nbot0

      ! Cotes of the top of the layers
      z0(0)   =  0.0
      DO layer = 1, nbot0-1
        limsum    =  ( (icsuswi * (icsuind+layer-1) + 
     &                 (1-icsuswi)*layer) ) *
     &                 (1-snicswi) + (layer-1)*snicswi
         z0(layer) =  icsuswi*dh_i_surf(ji) + snicswi*dh_snowice(ji)
         DO layer_a = 1, limsum
            z0(layer) = z0(layer) + deltaz_i_bio(layer_a)
         END DO
      END DO

      z0(nbot0)  =  icsuswi*dh_i_surf(ji) + snicswi*dh_snowice(ji) + 
     &              dh_i_bott(ji)

      DO layer = 1, nlay_bio
         z0(nbot0) = z0(nbot0) + deltaz_i_bio(layer)
      END DO
       
      z0(1)      =  snicswi * dh_snowice(ji) + (1 - snicswi) * z0(1)

      DO layer = ntop0, nbot0
         zthick0(layer)  =  z0(layer) - z0(layer-1)
      END DO

      WRITE(numout,*) ' Old grid, ntop0 : ', ntop0, ' nbot0 : ', nbot0
      WRITE(numout,*) ' z0      : ', ( z0(layer), layer = 0, nbot0 )
      WRITE(numout,*) ' zthick0 : ', ( zthick0(layer) ,
     &                               layer = 1, nbot0 )
      WRITE(numout,*)
!
!------------------------------------------------------------------------------|
! 2) New grid   
!------------------------------------------------------------------------------|
!
      WRITE(numout,*) ' New grid : '
      WRITE(numout,*) ' ~~~~~~~~~~ '
      ! indexes of the vectors
      ntop1    =  1 
      nbot1    =  nlay_bio
      ! cotes and thicknesses
      CALL ice_bio_grid(kideb,kiut,.TRUE.) ! compute the biological grid

      ! cotes of the layer interfaces
      z1(0) = 0.0
      DO layer = 1, nlay_bio
         z1(layer) = z1(layer-1) + deltaz_i_bio(layer)
         zthick1(layer) = z1(layer) - z1(layer-1)
      END DO ! layer

      WRITE(numout,*) ' New grid, ntop1 : ', ntop1, ' nbot1 : ', nbot1
      WRITE(numout,*) ' z1      : ', ( z1(layer), layer = 0, nbot1 )
      WRITE(numout,*) ' zthick1 : ', ( zthick1(layer) ,
     &                               layer = 1, nbot1 )
!
!------------------------------------------------------------------------------|
! 3) Weights
!------------------------------------------------------------------------------|
!
      WRITE(numout,*)
      WRITE(numout,*) ' Weights  : '
      WRITE(numout,*) ' ~~~~~~~~~~ '
      DO layer1 = 1, nlay_bio
         DO layer0 = 1, nbot0
            zweight(layer1,layer0) = MAX ( 0.0 , ( MIN ( z0(layer0) ,
     &      z1(layer1) ) - MAX ( z0 (layer0-1) , z1(layer1-1) ) ) / 
     &      zthick0(layer0) )
!           WRITE(numout,*) ' layer0, layer1 : ', layer0, layer1
!           WRITE(numout,*) ' zweight        : ', zweight(layer1,layer0)
         END DO
      END DO
      WRITE(numout,*)
!
!------------------------------------------------------------------------------|
! 4) Tracer sources
!------------------------------------------------------------------------------|
!
      WRITE(numout,*) ' Tracer sources : '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~~ '
      
      DO jn = 1, ntra_bio

         IF ( flag_active (jn) ) THEN

         c_s_bio(jn)   = 0.0 !now no concentration in the snow
!        e_skel = 0.1        !sensitivity run
         ch_bo_bio(jn) = e_skel * c_skel_bio(jn) * MAX( dh_i_bott(ji) 
     &                 , 0.0 ) * stickfac
         ch_si_bio(jn) = ( ( rhog - rhon ) / rhog * c_w_bio(jn) 
     &                 * frtr_si_bio * stickfac
     &                   + ( rhon / rhog ) * c_s_bio(jn) ) * 
     &                   MAX( dh_snowice(ji) , 0.0 )
         IF ( i_flux .EQ. 3 ) THEN
            ch_bo_bio(jn) = cbu_i_bio(jn,nlay_bio)
     &                    * MAX( dh_i_bott(ji) , 0.0 )
            ch_si_bio(jn) = cbu_i_bio(jn,1)
     &                    * MAX( dh_snowice(ji) , 0.0 )
         ENDIF

         IF ( .NOT. ln_trbo ) ch_bo_bio = cbu_i_bio(jn,nlay_bio) * 
     &                          MAX ( dh_i_bott(ji) , 0.0 )
         IF ( .NOT. ln_trsi ) ch_si_bio = cbu_i_bio(jn,1) * 
     &                          MAX ( dh_snowice(ji) , 0.0 )

!        fcbp(jn) = ( c_skel_bio(jn) - e_skel * c_skel_bio(jn) * 
!    &                stick_fac ) * dh_i_bott(ji) / ddtb 
!        fcsi(jn) = ( c_w_bio(jn) * ( rhog - rhon ) / rhog * 
!    &              dh_snowice(ji) - ch_si_bio(jn) )
!    &              / ddtb
         fcbp(jn) = ch_bo_bio(jn) / ddtb
         fcsi(jn) = ch_si_bio(jn) / ddtb

         WRITE(numout,*) ' Tracer    : ', biotr_i_nam(jn)
!        WRITE(numout,*) ' e_skel    : ', e_skel
!        WRITE(numout,*) ' c_w_bio   : ', c_w_bio(jn)
!        WRITE(numout,*) ' c_skel_bio: ', c_skel_bio(jn)
         WRITE(numout,*) ' ch_bo_bio : ', ch_bo_bio(jn)
         WRITE(numout,*) ' ch_si_bio : ', ch_si_bio(jn)
!        WRITE(numout,*) ' frtr_si_bio ', frtr_si_bio
!        WRITE(numout,*) ' c_s_bio   : ', c_s_bio (jn)
!        WRITE(numout,*) ' c_ws_i    : ', 
!    &                   ( rhog - rhon ) / rhog * c_w_bio(jn)

         ENDIF

      END DO !jn
!
!------------------------------------------------------------------------------|
! 5) Tracer contents
!------------------------------------------------------------------------------|
!
      WRITE(numout,*)

      DO jn = 1, ntra_bio

         IF ( flag_active (jn) ) THEN

         WRITE(numout,*) ' ---------------------------- '
         WRITE(numout,*) ' Tracer : ', biotr_i_nam(jn)
         WRITE(numout,*) ' Tracer content : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~ '

         !--------------
         ! Inner layers
         !--------------
         mt_i_bio_init(jn) = 0.0

         DO layer = ntop0, nbot0
!           limsum = snicswi * ( layer - 1 ) + icsuswi * ( layer - 1 + 
!    &               icsuind ) + ( 1 - icsuswi ) * ( 1 - snicswi ) * 
!    &               layer
!           limsum = ((icsuswi*(icsuind+layer-1) + (1-icsuswi)*layer))*
!    &               (1-snicswi) + (layer-1)*snicswi
            limsum =  MAX(1,MIN(snicswi*(layer-1) + icsuswi*(layer-1
     &           +  icsuind)
     &           + (1-icsuswi)*(1-snicswi)*layer, nlay_bio))
            zq0(layer) = zthick0(layer) * cbu_i_bio(jn,limsum)

!           WRITE(numout,*) ' limsum : ', limsum
!           IF ( limsum .GE. 1 ) 
!    &         WRITE(numout,*) ' cbu_i_bio: ', cbu_i_bio(jn,limsum)

!           WRITE(numout,*) ' zq0    : ', zq0(layer)
!           WRITE(numout,*) ' zthick0: ', zthick0(layer)

         END DO
         !--------------
         ! Bottom layer
         !--------------
         ! bottom layer if ice forms at the bottom
         zq0(nbot0) = ( 1 - icboswi ) * zq0(nbot0) + icboswi *
     &                ch_bo_bio(jn)
         !-------------
         ! Upper layer
         !-------------
         ! first layer if snow ice forms
         zq0(1)     =  snicswi * ch_si_bio(jn) 
     &              + ( 1 - snicswi ) * zq0(1)

!        WRITE(numout,*) ' snicswi   : ', snicswi
!        WRITE(numout,*) ' c_si_bio  : ', ch_si_bio(jn) / zthick0(1)

         DO layer = ntop0, nbot0
            mt_i_bio_init(jn) = mt_i_bio_init(jn) +
     &                          zq0(layer)
         END DO

         WRITE(numout,*) ' zq0       : ', ( zq0(layer), 
     &                                    layer = 1, nbot0 )
!
!------------------------------------------------------------------------------|
! 6) Relayering
!------------------------------------------------------------------------------|
!
         WRITE(numout,*)
         WRITE(numout,*) ' Relayering : '
         WRITE(numout,*) ' ~~~~~~~~~~~~ '

         DO layer1 = 1, nlay_bio
            zq1(layer1) = 0.0
            DO layer0 = 1, nbot0
               zq1(layer1) = zq1(layer1) + zweight(layer1,layer0) *
     &                       zq0(layer0)
            END DO
         END DO

         WRITE(numout,*) ' zq1       : ', ( zq1(layer), 
     &                                    layer = 1, nlay_bio )
!
!------------------------------------------------------------------------------|
! 7) Recompute bulk / brine concentration from content
!------------------------------------------------------------------------------|
!
         WRITE(numout,*) ' Bulk concentration : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~ '

         !------------------------
         ! Bulk ice concentration
         !------------------------
         DO layer = 1, nlay_bio
            cbu_i_bio(jn,layer) = zq1(layer) / deltaz_i_bio(layer)
         END DO
         WRITE(numout,*) ' Tracer : ', biotr_i_nam(jn)
         WRITE(numout,*) ' cbu_i_bio   : ', ( cbu_i_bio(jn,layer), 
     &                                    layer = 1, nlay_bio )
      
         WRITE(numout,*) ' ---------------------------- '

         ENDIF

      END DO !jn

      !---------------
      ! Interpolation 
      !---------------
      CALL ice_bio_interp_phy2bio(kideb,kiut,nlay_i,.FALSE.) 
                                          ! interpolation of physical variables
                                          ! on the biological grid
      !---------------------
      ! Brine concentration
      !---------------------
      DO jn = 1, ntra_bio
!        IF ( flag_active (jn) .AND. biotr_i_typ(jn) .EQ. 'nut' ) THEN
         IF ( flag_active (jn) .AND. flag_diff(jn) ) THEN
!           DO layer = 1, nlay_bio
!              c_i_bio(jn,layer) = cbu_i_bio(jn,layer) / e_i_bio(layer)
!           END DO ! layer
            WRITE(numout,*) ' Tracer : ', biotr_i_nam(jn)
            WRITE(numout,*) ' c_i_bio : ', ( c_i_bio(jn,layer), 
     &                                       layer = 1, nlay_bio )
         ENDIF ! flags
      END DO !jn

      END DO  !ji

!
!-----------------------------------------------------------------------
! 8) Conservation check
!-----------------------------------------------------------------------
!
      !-------------------
      ! Conservation test
      !-------------------
      CALL ice_bio_column(kideb,kiut,mt_i_bio_final,cbu_i_bio,
     &                    deltaz_i_bio, .FALSE.)
      DO jn = 1, ntra_bio
         IF ( flag_active(jn) ) THEN
            WRITE(numout,*) ' mt_i_bio_final : ', mt_i_bio_final(jn)
         ENDIF
      END DO

      DO jn = 1, ntra_bio

         IF ( flag_active(jn) ) THEN

         ! Flux due to bottom formation
         f_bogr(jn) =  ch_bo_bio(jn) / ddtb
         f_sigr(jn) =  ch_si_bio(jn) / ddtb
         f_bogr(jn) =  0.0
         f_sigr(jn) =  0.0

         ! Flux due to snow ice formation

         ! Bottom flux ( positive upwards )
!        f_bo_tra(jn) =
!        f_bo_tra(jn) = zswitchw * ( - e_i_bio( nlay_bio )
!    &                * diff_bio * 2.0
!    &                / deltaz_i_bio(nlay_bio) * ( c_i_bio(jn,nlay_bio)
!    &                  - c_skel_bio(jn) ) )
!    &                + zswitchs * ( - qsummer * c_i_bio(jn,nlay_bio) )
!    &                / ddtb
!    &                * e_i_bio(nlay_bio)
!        f_su_tra(jn) = zswitchw * 0.0
!    &                + zswitchs * ( qsummer * c_s_bio(jn) )
!    &                / ddtb
         !+++++
         WRITE(numout,*) ' f_bogr : ', f_bogr(jn)
         WRITE(numout,*) ' f_sigr : ', f_sigr(jn)
         !+++++

         ENDIF

      END DO

      zerror = 1.0e-9
      CALL ice_bio_conserv(kideb,kiut,'ice_bio_remap: ',zerror,
     &                           mt_i_bio_init,mt_i_bio_final,
     &                           f_bogr, f_sigr, ddtb)
      CALL ice_bio_column(kideb,kiut,ct_i_bio,cbu_i_bio,
     &                    deltaz_i_bio, .FALSE.)

!==============================================================================|
! Fin de la routine ice_bio_remap
      WRITE(numout,*)

      END
