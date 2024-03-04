      SUBROUTINE ice_phy_remap(nlay_s,nlay_i,kideb,kiut)

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

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'
      INCLUDE 'bio.com'

      ! Local Variables
      LOGICAL ln_write

      REAL(8), DIMENSION (maxnlay) :: 
     &   zsh_i0              ,   !: old ice salt content (ppt.m-2)
     &   zsh_i1                  !: new ice salt content (ppt.m-2)

      REAL(8), DIMENSION (maxnlay) :: 
     &   zqh_i0              ,   !: old ice heat content (J.m-2)
     &   zqh_i1                  !: new ice heat content (J.m-2)

      REAL(8), DIMENSION (maxnlay) :: 
     &   zqh_s0              ,   !: old snow heat content (J.m-2)
     &   zqh_s1                  !: new snow heat content (J.m-2)

      REAL(8), DIMENSION( maxnlay ) ::
     &  zthick1                  !: thickness of physical layers

      REAL(8), DIMENSION (maxnlay+2) ::
     &   zthick0                 !: thickness of old layers
! FD additions
      LOGICAL ln_con

      ! Local Constants
      zeps   = 1.0e-20

      ln_write = .TRUE.
      ln_con   = .TRUE.

      s_i_b(1,:) = sn_i_b(:) ! update salinity

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

      ji = 1
      IF ( ln_con ) THEN
         CALL ice_sal_column(kideb,kiut,z_ms_i_ini,s_i_b(1,1:nlay_i),
     &                       deltaz_i_phy, nlay_i, .FALSE. )
         ! salt content loss due to ice melt
         zfmelt = ( s_i_b(ji,1) * dh_i_surf(ji) + s_i_b(ji,nlay_i) * 
     &              MIN(dh_i_bott(ji),0.0) ) / ddtb
      ENDIF

      zqt_s_ini = 0.0 ! total heat content in snow (initial value)
      zqt_i_ini = 0.0 ! total heat content in ice  (initial value)
      zqt_s_fin = 0.0 ! total heat content in snow (final value)
      zqt_i_fin = 0.0 ! total heat content in ice  (final value)
      zqt_ini   = 0.0 ! total heat content snow + ice (initial value)
      zqt_fin   = 0.0 ! total heat content snow + ice (final value)

      DO 10 ji = kideb, kiut
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
      deltah   =  0.0
      DO layer = 1, nlay_s
         IF ( - dh_s_tot(ji) .GT. deltah ) THEN
            snind  =  layer 
         ENDIF
         deltah = deltah + deltaz_s_phy(layer)
      END DO

      ! snswi : switch which value equals 1 if snow melts
      !                                 0 if not
      snswi    =  MAX(0,INT(-dh_s_tot(ji)/MAX(zeps,ABS(dh_s_tot(ji)))))

      ! ice surface behaviour : computation of icsuind-icsuswi
      ! icsuind : index which values equals 0 if there is nothing
      !                                 1 if first layer has started to melt
      !                                 2 if first and second layer have started ...
      !                                 etc ...
      icsuind    =  0  
      deltah   =  0.0
      DO layer = 1, nlay_i
         IF  ( -dh_i_surf(ji).GT.deltah) THEN
            icsuind  =  layer 
         ENDIF
         deltah = deltah + deltaz_i_phy(layer)
      END DO

      ! icsuswi : switch which value equals 1 if ice melts at the surface
      !                                   0 if not
      icsuswi    =  max(0,int(-dh_i_surf(ji)
     &              /max(zeps,abs(dh_i_surf(ji)))))

      ! ice bottom behaviour : computation of icboind-icboswi
      ! icboind : index which values equals 0 if accretion is on the way
      !                                     1 if last layer has started to melt
      !                               2 if penultiem layer has started ...
      !                               etc ...
      icboind    =  0  
      deltah   =  0.0
      DO layer = nlay_i, 1, -1
         IF (-dh_i_bott(ji) .GT. deltah) THEN
            icboind  =  nlay_i+1-layer 
         ENDIF
         deltah = deltah + deltaz_i_phy(layer)
      END DO

      ! icboswi : switch which value equals 1 if accretion of ice is on the way
      !                                         0 if ablation is on the way
      icboswi    =  max(0,int(dh_i_bott(ji)
     &              /max(zeps,abs(dh_i_bott(ji)))))

      ! snow-ice formation : calcul de snicind-snicswi
      ! snicind : index which values equals 0 if no snow-ice forms
      !                1 if last layer of snow has started to melt
      !                          2 if penultiem layer ...
      snicind    =  0
      deltah   =  0.0
      DO layer = nlay_s, 1, -1
         IF ( dh_snowice(ji) .GT. deltah ) THEN
            snicind  =  nlay_s+1-layer
         ENDIF
         deltah = deltah + deltaz_s_phy(layer)
      END DO

      ! snicswi : switch which value equals 1 if snow-ice forms
      !                                     0 if not
      snicswi    =  max(0,int(dh_snowice(ji)/
     &              max(zeps,abs(dh_snowice(ji)))))

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
      nbot0    =  nlay_s+1-snind+(1-snicind)*snicswi

      ! cotes of the top of the layers
      zm0(0)   =  0.0
      DO layer = 1, nbot0-1
         limsum    =  snswi*(layer+snind-1) + (1-snswi)*(layer-1)
         limsum    =  MIN( limsum , nlay_s )
         zm0(layer) =  dh_s_tot(ji)
         DO layer_a = 1, limsum
            zm0(layer) = zm0(layer) + deltaz_s_phy(layer_a)
         END DO
      END DO

      zm0(nbot0) =  dh_s_tot(ji) - snicswi*dh_snowice(ji)
      DO layer = 1, nlay_s
        zm0(nbot0) = zm0(nbot0) + deltaz_s_phy(layer)
      END DO
      zm0(1)     =  dh_s_tot(ji)*(1-snswi) + snswi*zm0(1)

      ! thicknesses
      DO layer = ntop0, nbot0
         zthick0(layer)  =  zm0(layer) - zm0(layer-1)
      END DO
!
!------------------------------------------------------------------------------|
!      S-2) 'Old' enthalpies
!------------------------------------------------------------------------------|
!
      ! enthalpies
      zqh_s0(1) =  rhon * ( cpg * ( tpw - (1 - snswi) * tabqb(ji) - ! fallen snow
     &             snswi * t_s_b(ji,1) ) + xlgn )* zthick0(1) 
      zqt_s_ini = zqt_s_ini + zqh_s0(1)

      DO layer = 2, nbot0 ! remaining snow
         limsum = (1-snswi)*(layer-1) + snswi*(layer+snind-1)
         limsum = MIN( limsum, nlay_s )
         zqh_s0(layer) = rhon*( cpg*(tpw - t_s_b(ji,limsum)) + xlgn) 
     &            *zthick0(layer)
         zswitch   = 1.0 - MAX (0.0, SIGN ( 1.0, zeps - ht_s_b(ji) ) )
         zqt_s_ini = zqt_s_ini + ( 1 - snswi ) * zqh_s0(layer) * zswitch
         WRITE(numout,*) 'limsum : ', limsum
      END DO
      zqt_ini = zqt_s_ini

      IF (ln_write) THEN
         WRITE(numout,*) ' - Snow redistribution ... '
         WRITE(numout,*) ' ntop0 : ', ntop0
         WRITE(numout,*) ' nbot0 : ', nbot0
         WRITE(numout,*) ' zm0   : ', ( zm0(layer), layer = 0, nbot0)
         WRITE(numout,*) ' zqh_s0: ', ( zqh_s0(layer), layer = 1, nbot0)
         WRITE(numout,*) ' q_s_b1: ', q_s_b(ji,1)
         WRITE(numout,*)
      ENDIF

      !-----------------------------------------
      ! snow enthalpy lost in sea ice formation
      !-----------------------------------------
      zqsnow     =  0.0
      zdeltah    =  0.0

      DO layer =  nlay_s, 1, -1
         zhsnow    =  MAX(0.0,dh_snowice(ji)-zdeltah)
         zqsnow    =  zqsnow + q_s_b(ji,layer) * zhsnow
      END DO

!
!------------------------------------------------------------------------------|
!      S-3) New snow grid
!------------------------------------------------------------------------------|
!
      ! indexes of the vectors
      ntop1    =  1
      nbot1    =  nlay_s

      ! cotes and thicknesses
      CALL ice_phy_grid( kideb , kiut , nlay_s , ht_s_b(ji), .FALSE., 
     &                   "sno" )
!
!------------------------------------------------------------------------------|
!      S-4) New snow enthalpy
!------------------------------------------------------------------------------|
!
      CALL ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , 
     &                    zthick0, deltaz_s_phy, zqh_s0 , zqh_s1 )

      zqt_s_fin = 0.0 ! total heat content
!
!------------------------------------------------------------------------------|
!      S-5) Recover snow temperature
!------------------------------------------------------------------------------|
!
      DO layer = 1, nlay_s
         isnow   = INT( 1.0 - MAX( 0.0 , SIGN( 1.0 , - ht_s_b(ji) ) ) )
         t_s_b(ji,layer)  =  isnow * ( tpw - ( zqh_s1(layer) / rhon/ 
     &                       MAX( deltaz_s_phy(layer) , zeps ) - xlgn )
     &                       / cpg ) +
     &                       ( 1.0 - isnow ) * tpw
         zqt_s_fin = zqt_s_fin + zqh_s1(layer)
         q_s_b(ji,layer) = isnow * zqh_s1(layer) / 
     &                     MAX( deltaz_s_phy(layer) , zeps )
      END DO
      zqt_fin = zqt_s_fin

      IF (ln_write) THEN
         WRITE(numout,*) ' - Heat conservation, ice_th_remap , snow '
         WRITE(numout,*) ' zqt_s_fin : ', zqt_s_fin
         WRITE(numout,*) ' zqt_s_ini : ', zqt_s_ini
         WRITE(numout,*) ' TS1. zqt_s_fin - zqt_s_ini : ', 
     &                   zqt_s_fin - zqt_s_ini
         WRITE(numout,*)
      ENDIF

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
      nbot0    =  MAX(1, MIN( nlay_i + 1 - icboind + (1-icsuind)*icsuswi
     &                      + snicswi, nlay_i + 2 ) )

      ! cotes of the top of the layers
      zm0(0)   =  0.0
      DO layer = 1, nbot0-1
         limsum    =  ((icsuswi*(icsuind+layer-1) + (1-icsuswi)*layer))*
     &               (1-snicswi) + (layer-1)*snicswi
         zm0(layer)=  icsuswi*dh_i_surf(ji) + snicswi*dh_snowice(ji)
         DO layer_a = 1, limsum
            zm0(layer) = zm0(layer) + deltaz_i_phy(layer_a)
         END DO
      END DO

      zm0(nbot0) =  icsuswi*dh_i_surf(ji) + snicswi*dh_snowice(ji) + 
     &              dh_i_bott(ji)
      DO layer = 1, nlay_i
         zm0(nbot0) = zm0(nbot0) + deltaz_i_phy(layer)
      END DO
      zm0(1)     =  snicswi*dh_snowice(ji) + (1-snicswi)*zm0(1)

      ! thicknesses
      DO layer = ntop0, nbot0
         zthick0(layer)  =  zm0(layer) - zm0(layer-1)
      END DO

      IF (ln_write) THEN
         WRITE(numout,*) ' - Ice redistribution ... '
         WRITE(numout,*) ' ntop0 : ', ntop0
         WRITE(numout,*) ' nbot0 : ', nbot0
         WRITE(numout,*) ' zm0   : ', ( zm0(layer), layer = 0, nbot0)
         WRITE(numout,*) ' zthick0: ', ( zthick0(layer), layer = ntop0, 
     &                                  nbot0)
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
      tmelts = - tmut * s_i_new + tpw
      zqh_i_new =  rhog * ( cpg * ( tmelts - t_bo_b(ji) ) 
     &             + xlgm * ( 1.0 - ( tmelts - tpw ) /
     &             ( t_bo_b(ji) - tpw ) ) - cpw * ( tmelts - tpw ) )
     &             * zthick0(nbot0)

      !- snow ice
      zqh_i_sni =  rho0 * cpw * ( tpw - t_bo_b(ji) ) * dh_snowice(ji) * 
     &            ( rhog - rhon ) / rhog * snicswi  ! generally positive
      fsnic  =  - zqh_i_sni / ddtb
      zqh_i_sni =  zqsnow + zqh_i_sni

      !----------------------
      ! sea ice heat content
      !----------------------
      zqh_i0(:) = 0.0
      !- internal ice layers
      DO layer = ntop0, nbot0
          limsum = MAX( 1 , MIN( snicswi*(layer-1) 
     &           + icsuswi * ( layer - 1 +  icsuind )
     &           + ( 1 - icsuswi ) * ( 1 - snicswi ) * layer, nlay_i ) )
          zqh_i0(layer) = q_i_b(ji,limsum) * zthick0(layer)
          WRITE(numout,*) ' limsum : ', limsum, ' layer : ', layer
      END DO
        
      !- boundary values 
      zqh_i0(1)     =  snicswi * zqh_i_sni + ( 1 - snicswi ) * zqh_i0(1)
      zqh_i0(nbot0) =  ( 1 - icboswi ) * zqh_i0(nbot0) + 
     &                 icboswi * zqh_i_new

      IF (ln_write) THEN
         WRITE(numout,*) ' zqh_i_new: ', zqh_i_new
         WRITE(numout,*) ' zqh_i_sni: ', zqh_i_sni
         WRITE(numout,*) ' zqh_i0: ', ( zqh_i0(layer), layer = ntop0, 
     &                                  nbot0)
      ENDIF

      DO layer = ntop0, nbot0
         zqt_i_ini = zqt_i_ini + zqh_i0(layer) 
      END DO

      zqt_ini = zqt_ini + zqt_i_ini ! includes zqsnic and zqsnow

!
!------------------------------------------------------------------------------|
!  I-3) Old ice salt content                                                   |
!------------------------------------------------------------------------------|
!
      ! basal new ice salt content
      zsh_i_new = e_skel * oce_sal * MAX ( dh_i_bott(ji), 0.0 )
      ! snow ice salt content
      s_i_snic  =  ( rhog - rhon ) / rhog * oce_sal * 
     &              frtr_si_phy ! snow ice salinity
      zsh_i_sni = s_i_snic * dh_snowice(ji)

      IF ( ln_write ) THEN
         WRITE(numout,*) ' Salt sources : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~ '
         WRITE(numout,*) ' zsh_i_new : ', zsh_i_new
         WRITE(numout,*) ' zsh_i_sni : ', zsh_i_sni
      ENDIF

      zsh_i0(:) = 0.0

      DO layer = ntop0, nbot0
          limsum =  snicswi*(layer-1) + icsuswi*(layer-1+icsuind)
     &            + (1-icsuswi)*(1-snicswi)*layer
          limsum =  MAX(1,MIN(snicswi*(layer-1) + icsuswi*(layer-1
     &           +  icsuind)
     &           + (1-icsuswi)*(1-snicswi)*layer, nlay_i))
          zsh_i0(layer) = s_i_b(ji,limsum)*zthick0(layer)
      END DO

      zsh_i0(nbot0) =  (1-icboswi)* zsh_i0(nbot0) + icboswi * zsh_i_new
      zsh_i0(1)     =  snicswi * zsh_i_sni + ( 1 - snicswi ) * zsh_i0(1)
     
      z_ms_i_ini = 0.0
      DO layer = ntop0, nbot0
         z_ms_i_ini = z_ms_i_ini + zsh_i0(layer)
      END DO

      IF ( ln_write ) THEN
         WRITE(numout,*) ' frtr_si_phy:', frtr_si_phy
         WRITE(numout,*) ' oce_sal  : ', oce_sal
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
      nbot1    =  nlay_i

      ! cotes and thicknesses
      CALL ice_phy_grid( kideb , kiut , nlay_i , ht_i_b(ji), .FALSE., 
     &                   "ice" )

! FD      IF (ln_write) THEN
! FD         WRITE(numout,*) ' z_i : ', ( z_i(layer), layer = ntop1, nbot1 )
! FD         WRITE(numout,*) ' h_i : ', ( h_i(layer), layer = ntop1, nbot1 )
! FD      ENDIF
!
!------------------------------------------------------------------------------|
!  I-5) Redistribute ice heat content                                          |
!------------------------------------------------------------------------------|
!
      ! Remapping ice enthalpy
      WRITE(numout,*) ' deltaz_i_phy : ', ( deltaz_i_phy(layer), 
     &                                      layer = 1, nlay_i )
      CALL ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , 
     &                    zthick0, deltaz_i_phy, zqh_i0 , zqh_i1 )

      IF (ln_write) THEN
         WRITE(numout,*) ' zqh_i1 : ',( zqh_i1(layer), layer = ntop1,
     &                                  nbot1)
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
      CALL ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , 
     &                    zthick0, deltaz_i_phy, zsh_i0 , zsh_i1 )

!
!------------------------------------------------------------------------------|
!  I-7) Heat conservation test                                                 |
!------------------------------------------------------------------------------|
!
      IF ( ln_write ) THEN
         WRITE(numout,*) ' - Heat conservation, ice_th_remap , 
     &   sea ice... '
         WRITE(numout,*) ' zqt_i_fin : ', zqt_i_fin
         WRITE(numout,*) ' zqt_i_ini : ', zqt_i_ini
         WRITE(numout,*) ' TI1. zqt_i_fin - zqt_i_ini : ', 
     &                   zqt_i_fin - zqt_i_ini
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
      DO layer = 1, nlay_i
         q_i_b(ji,layer) = zqh_i1(layer) / 
     &                     MAX( deltaz_i_phy(layer) , zeps )
      END DO

      DO layer = 1, nlay_i
         s_i_b(ji,layer) = zsh_i1(layer) / MAX( deltaz_i_phy(layer) , 
     &                     zeps )
      END DO

      DO layer = 1, nlay_i
         tmelts = -tmut*s_i_b(ji,layer) + tpw
         aaa = cpg
         bbb = (cpw-cpg)*(tmelts-tpw) + q_i_b(ji,layer)/ rhog
     &       - xlgm 
         ccc = xlgm * (tmelts-tpw)
         discrim = SQRT( bbb*bbb - 4.0*aaa*ccc )
         t_i_b(ji,layer) = tpw + (- bbb - discrim) / ( 2.0*aaa )
      END DO
!
!------------------------------------------------------------------------------|
!  I-9) Salt conservation test                                                 |
!------------------------------------------------------------------------------|
!
      CALL ice_sal_column(kideb,kiut,z_ms_i_fin,s_i_b(ji,1:nlay_i),
     &                    deltaz_i_phy, nlay_i, .FALSE.)
      IF ( ln_write ) WRITE(numout,*) ' z_ms_i_fin : ',
     &                z_ms_i_fin

      ! Flux due to bottom formatiophy
      zfgrml = zsh_i_new / ddtb
      zfsigr = zsh_i_sni / ddtb

      zfgrml = zfgrml + zfmelt
      IF ( ln_write ) THEN
         WRITE(numout,*) ' zfgrml : ', zfgrml
         WRITE(numout,*) ' zfsigr : ', zfsigr
      ENDIF

      zerror = 1.0e-9
      CALL ice_sal_conserv(kideb,kiut,'ice_phy_remap: ',zerror,
     &                           z_ms_i_ini,z_ms_i_fin,
     &                           zfgrml, zfsigr, ddtb)

 10   CONTINUE

      RETURN

!-------------------------------------------------------------------------------
! Fin de la routine ice_th_remap
      END
