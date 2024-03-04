      SUBROUTINE ice_sal_diff(nlay_i,kideb,kiut)

      !!------------------------------------------------------------------
      !!                ***         ROUTINE ice_sal_diff      ***
      !!
      !! ** Purpose :
      !!        This routine computes new salinities in the ice
      !!
      !! ** Method  : Vertical salinity profile computation 
      !!              Resolves brine transport equation
      !!           
      !! ** Steps
      !!
      !! ** Arguments
      !!
      !! ** Inputs / Outputs
      !!
      !! ** External
      !!
      !! ** References : Vancop. et al., 2008
      !!
      !! ** History  : 
      !!    (06-2003) Martin Vancop. LIM1D
      !!    (06-2008) Martin Vancop. BIO-LIM
      !!    (09-2008) Martin Vancop. Explicit gravity drainage
      !!
      !!------------------------------------------------------------------

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'

      REAL(8), DIMENSION(nlay_i) ::  
     &   z_ms_i             ,    !: mass of salt times thickness 
     &   z_sbr_i                 !: brine salinity

      REAL(8), DIMENSION(nlay_i) ::  !: dummy factors for tracer equation
     &   za                 ,    !: winter
     &   zb                 ,   
     &   ze                 ,    !: summer
     &   zf                 ,   
     &   zind               ,    !: independent term in the tridiag system
     &   zindw              ,    !: independent term in the tridiag system
     &   zinds              ,    !: independent term in the tridiag system
     &   zindtbis           ,    !:
     &   zdiagbis                !:

      REAL(8), DIMENSION(nlay_i,3) :: !: dummy factors for tracer equation
     &   ztrid              ,    !: tridiagonal matrix
     &   ztridw             ,    !: tridiagonal matrix, winter
     &   ztrids                  !: tridiagonal matrix, summer

      REAL(8) ::  
     &   zdummy1            ,    !: dummy factors
     &   zdummy2            ,    !: 
     &   zdummy3            ,    !: 
     &   zswitch_open       ,    !: switch for brine network open or not
     &   zswitchw           ,    !: switch for winter drainage 
     &   zswitchs           ,    !: switch for summer drainage
     &   zeps      = 1.0e-20     !: numerical limit

      ! Rayleigh number computation
      REAL(8) ::
     &   ze_i_min           ,    !: minimum brine volume
     &   zc                 ,    !: temporary scalar for sea ice specific heat
     &   zk                 ,    !: temporary scalar for sea ice thermal conductivity
     &   zalphara                !: multiplicator for diffusivity

      REAL(8), DIMENSION(nlay_i) ::  
     &   zsigma             ,    !: brine salinity at layer interfaces
     &   zperm              ,    !: permeability
     &   zpermin            ,    !: minimum permeability
     &   zrhodiff           ,    !: density difference
     &   zlevel             ,    !: height of the water column
     &   zthdiff                 !: thermal diffusivity

      INTEGER ::  
     &   layer2             ,    !: layer loop index
     &   indtr                   !: index of tridiagonal system

      CHARACTER(len=4)      ::  
     &   bc = 'conc'             !: Boundary condition 'conc' or 'flux'

      REAL(8) ::  
     &   z_ms_i_ini         ,    !: initial mass of salt
     &   z_ms_i_fin         ,    !: final mass of salt
     &   z_fs_b             ,    !: basal flux of salt
     &   z_fs_su            ,    !: surface flux of salt
     &   z_dms_i                 !: mass variation       
! FD additions
      LOGICAL ln_write
      LOGICAL ln_con
      LOGICAL ln_sal
      LOGICAL ln_be
      LOGICAL ln_gd
      LOGICAL ln_fl


      ln_write = .TRUE.  ! write outputs
      ln_con   = .TRUE.  ! conservation check
      ln_sal   = .TRUE.  ! compute salinity variations or not
      ln_be    = .TRUE.  ! compute brine expulsion or not
      ln_gd    = .TRUE.  ! compute gravity drainage or not
      ln_fl    = .TRUE.  ! compute flushing or not

      IF ( ln_write ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' ** ice_sal_diff : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) ' ln_sal = ', ln_sal
      ENDIF
      WRITE(numout,*) " nlay_i : ", nlay_i

      IF ( ln_sal ) THEN
!
!------------------------------------------------------------------------------|
! 1) Initialization
!------------------------------------------------------------------------------|
!
      IF ( ln_write ) THEN
         WRITE(numout,*) ' - Initialization ... '
      ENDIF

      DO 10 ji = kideb, kiut

      ! brine diffusivity
      diff_br(:) = 0.0

      !---------------------------
      ! Brine volume and salinity
      !---------------------------
      DO layer = 1, nlay_i
         e_i_b(layer) = - tmut * s_i_b(ji,layer) / ( t_i_b(ji,layer) 
     &                - tpw )
         z_sbr_i(layer) = s_i_b(ji,layer) / e_i_b(layer)
      END DO

      !----------
      ! Switches
      !----------
      ! summer switch = 1 if Tsu ge tpw and min brine volume superior than e_tres
      zswitchs = MAX( 0.0, SIGN ( 1.0 , t_su_b(ji) - tpw ) ) ! 0 si hiver 1 si ete
      zbvmin   = 1.0
      DO layer = 1, nlay_i
         zbvmin = MIN( e_i_b(layer) , zbvmin ) ! minimum brine volume
      END DO
      IF ( zbvmin .LT. e_tres ) zswitchs = 0.0

      ! winter switch
      zswitchw = 1.0 - zswitchs

      !------------------
      ! Percolating flux
      !------------------
      ! Percolating flow ( rho dh * beta * switch / rhow )
      qsummer = ( - rhog * MIN ( dh_i_surf(ji) , 0.0 ) 
     &            - rhon * MIN ( dh_s_tot(ji)  , 0.0 ) )  
      qsummer = qsummer * flu_beta * zswitchs / 1000.0

      !--------------------
      ! Conservation check
      !--------------------
      IF ( ln_con ) THEN 
         CALL ice_sal_column( kideb , kiut , z_ms_i_ini , 
     &                    s_i_b(1,1:nlay_i),
     &                    deltaz_i_phy, nlay_i, .FALSE. )
      ENDIF ! ln_con

      IF ( ln_write ) THEN
         WRITE(numout,*) ' nlay_i      : ', nlay_i
         WRITE(numout,*) ' kideb       : ', kideb
         WRITE(numout,*) ' kiut        : ', kiut 
         WRITE(numout,*) ' '
         WRITE(numout,*) ' deltaz_i_phy : ', ( deltaz_i_phy(layer),  
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' z_i_phy      : ', ( z_i_phy(layer),  
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' s_i_b       : ', ( s_i_b    (ji,layer), 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' t_i_b       : ', ( t_i_b    (ji,layer), 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' e_i_b       : ', ( e_i_b      (layer), 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' z_sbr_i     : ', ( z_sbr_i    (layer), 
     &                   layer = 1, nlay_i )
         WRITE(numout,*)
         WRITE(numout,*) ' zswitchs : ', zswitchs
         WRITE(numout,*) ' zswitchw : ', zswitchw
         WRITE(numout,*) 
      ENDIF ! ln_write

 10   CONTINUE

!
!------------------------------------------------------------------------------|
! 2) Rayleigh-number-based diffusivity
!------------------------------------------------------------------------------|
!
! Diffusivity is a function of the local Rayleigh number
! see Notz and Worster, JGR 2008
! Diffusivity, layer represents the interface'
! between layer and layer-1 '
!
      IF ( zswitchw .GE. 1.0 ) THEN

      IF ( ln_write ) THEN
         WRITE(numout,*) ' - Rayleigh-number based diffusivity ... '
         WRITE(numout,*) ' '
      ENDIF

      DO 20 ji = kideb, kiut 

      !-----------------------------------------
      ! Brine salinity between layer interfaces
      !-----------------------------------------
      DO layer = 1, nlay_i - 1
         zdummy1 = t_i_b(ji,layer) + deltaz_i_phy(layer) / 2. * 
     &             ( t_i_b(ji,layer+1) - t_i_b(ji,layer) ) /
     &             ( z_i_phy(layer+1) - z_i_phy(layer) ) - tpw
         zsigma(layer) = - zdummy1 / tmut
      END DO
      zsigma(nlay_i) = - ( t_i_b(ji,nlay_i) - tpw ) / tmut

      !--------------------
      ! Density difference
      !--------------------
      DO layer = 1, nlay_i  
         zrhodiff(layer) = - beta_ocs * ( oce_sal - zsigma(layer) )
      END DO

      !------------------------------------------
      ! Minimum permeability under current level
      !------------------------------------------
      DO layer = 1, nlay_i
         ze_i_min = 99999.0 
         DO layer2 = layer, nlay_i
            ze_i_min = MIN( ze_i_min , e_i_b(layer2) )
            zpermin(layer) = 1.0e-17 * ( ( 1000. * ze_i_min )**3.1 )
         END DO
      END DO ! layer

      !------------------------------------------------
      ! length of the water column under current level
      !------------------------------------------------
      DO layer = nlay_i, 1, -1
         zlevel(layer) = 0.0
         DO layer2 = layer, nlay_i
            zlevel(layer) = zlevel(layer) + deltaz_i_phy(layer2)
         END DO
      END DO
      zlevel(nlay_i) = deltaz_i_phy(nlay_i) / 2.0

      !---------------------
      ! Thermal diffusivity
      !---------------------
      zkimin = 0.1
      DO layer = 1, nlay_i - 1
         zdummy1 = t_i_b(ji,layer) + deltaz_i_phy(layer) / 2. * 
     &             ( t_i_b(ji,layer+1) - t_i_b(ji,layer) ) /
     &             ( z_i_phy(layer+1) - z_i_phy(layer) ) - tpw
         zdummy2 = s_i_b(ji,layer) + deltaz_i_phy(layer) / 2. * 
     &             ( s_i_b(ji,layer+1) - s_i_b(ji,layer) ) /
     &             ( z_i_phy(layer+1) - z_i_phy(layer) )
         zc = cpg + xlgm * tmut * zdummy2 / 
     &        MAX( zdummy1 * zdummy1 , zeps )
         zk = xkg + betak1 * zdummy2 / 
     &        MIN( -zeps , zdummy1 ) - betak2 * zdummy1
         zk = MAX( zk, zkimin )
         zthdiff(layer) = zk / ( rhog * zc )
      END DO

      zc = cpg + xlgm * tmut * s_i_b(ji,nlay_i) / 
     &     MAX( ( t_i_b(ji,nlay_i) - tpw ) * ( t_i_b(ji,nlay_i) - tpw ),
     &     zeps )
      zk = xkg + betak1 * s_i_b(ji,nlay_i) / 
     &     MIN( -zeps , t_i_b(ji,nlay_i) - tpw ) 
     &   - betak2 * ( t_i_b(ji,nlay_i) - tpw )
      zk = MAX( zk, zkimin )
      zthdiff(nlay_i) = zk / ( rhog * zc )

      !-----------------
      ! Rayleigh number
      !-----------------
      DO layer = 1, nlay_i
         rayleigh(layer) = gpes * MAX(zrhodiff(layer),0.0) * 
     &                     zpermin(layer) * zlevel(layer) / 
     &                     ( zthdiff(layer) * visc_br )
      END DO

      !-------------------
      ! Brine Diffusivity
      !-------------------
      DO layer = 1, nlay_i
         zalphara = ( TANH( ra_smooth * ( rayleigh(layer) - ra_c ) ) 
     &            + 1 ) / 2.0
!        diff_br(layer) = ( 1.0 - zalphara ) * d_br_mol + 
!    &                    zalphara * ( d_br_tur - d_br_mol )
         diff_br(layer) = ( 1.0 - zalphara ) * d_br_mol + 
     &                    zalphara * ( d_br_tur )
      END DO

      IF ( ln_write ) THEN
         WRITE(numout,*) ' zsigma   : ', ( zsigma(layer),  
     &                   layer = 1, nlay_i)
         WRITE(numout,*) ' zrhodiff : ', ( zrhodiff(layer),
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' zpermin : ', ( zpermin(layer), 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' zthdiff  : ', ( zthdiff(layer), 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' zlevel  : ', ( zlevel(layer), 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' rayleigh      : ', ( rayleigh(layer), 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' diff_br  : ', ( diff_br(layer), 
     &                   layer = 1, nlay_i )
      ENDIF

 20   CONTINUE

      ENDIF ! zswitchw

!
!------------------------------------------------------------------------------|
! 3) Compute dummy factors for tracer diffusion equation
!------------------------------------------------------------------------------|

      IF ( ln_write ) THEN
         WRITE(numout,*) ' - Compute dummy factors for tracer diffusion'
         WRITE(numout,*) ' '
      ENDIF

      DO 30 ji = kideb, kiut 
         
      !----------------
      ! Winter factors
      !----------------
      ! za factors
      zdummy1 = ddtb 
      DO layer = 1, nlay_i
         za(layer) = zdummy1 / ( deltaz_i_phy(layer) * e_i_b(layer) )
      END DO

      ! zb factors
      DO layer = 1, nlay_i - 1
         ! interpolate brine volume at the interface between layers
         zdummy1 = ( e_i_b(layer + 1 ) - e_i_b(layer) ) /  
     &             ( z_i_phy(layer + 1) - z_i_phy(layer) )
         zdummy2 = deltaz_i_phy(layer) / 2.0
         zdummy3 = e_i_b(layer) + zdummy1 * zdummy2
         zswitch_open = 0.0
         ! compute zswitch_open which equals 1 if the brine network is open
         IF ( zdummy3 .GE. e_tres ) zswitch_open = 1.0
         zb(layer) = zdummy3 * zswitch_open * diff_br(layer) / 
     &               ( z_i_phy(layer + 1) - z_i_phy(layer) )
      END DO

      zswitch_open = 0.0
      IF ( e_i_b(nlay_i) .GE. e_tres ) zswitch_open = 1.0

      ! Fixed boundary condition (imposed cc.)
      IF ( bc .EQ. 'conc' )  
     &   zb(nlay_i) = 2. * e_i_b(nlay_i) * zswitch_open / 
     &                  deltaz_i_phy(nlay_i) * diff_br(nlay_i)

      IF ( ln_write ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' -Winter factors '
         WRITE(numout,*) ' zswitch_open : ', zswitch_open
         WRITE(numout,*) ' za       : ', ( za (layer),  
     &                   layer = 1, nlay_i)
         WRITE(numout,*) ' zb       : ', ( zb (layer),  
     &                   layer = 1, nlay_i)
      ENDIF

      !----------------------
      ! Summer factors
      !----------------------
      ! ze factors
      DO layer = 1, nlay_i
         ze(layer) = qsummer * zswitchs / 
     &               ( e_i_b(layer) * deltaz_i_phy(layer) )
      END DO ! layer

! FD sure we'll do!
! FD      ! zf factors
! FD      ! could remove those, they are totally useless!!! ;-)
! FD      DO layer = 1, nlay_i - 1
! FD         zf(layer) = 1./2. * deltaz_i_phy(layer) / 
! FD     &               ( z_i_phy(layer+1)  - z_i_phy(layer) )
! FD      END DO ! layer

      IF ( ln_write ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' -Summer factors '
         WRITE(numout,*) ' ze : ', ( ze(layer), layer = 1, nlay_i ) 
! FD         WRITE(numout,*) ' zf : ', ( zf(layer), layer = 1, nlay_i ) 
      ENDIF
      WRITE(numout,*) ' truc 1 '

 30   CONTINUE
!
!-----------------------------------------------------------------------
! 4) Tridiagonal system terms for tracer diffusion equation, winter
!-----------------------------------------------------------------------
!
      DO 40 ji = kideb, kiut 

      !----------------
      ! first equation
      !----------------
      ztridw(1,1) = 0.0
      ztridw(1,2) = 1.0 + za(1) * zb(1)
      ztridw(1,3) = - za(1) * zb(1)
      zindw(1)    = z_sbr_i(1)
      WRITE(numout,*) ' truc 2 '

      !-----------------
      ! inner equations
      !-----------------
      DO layer = 2, nlay_i - 1
         ztridw(layer,1) = - za(layer) * zb(layer-1)
         ztridw(layer,2) = 1.0 + za(layer) * ( zb(layer-1) + 
     &                                         zb(layer) )
         ztridw(layer,3) = - za(layer) * zb(layer)
         zindw(layer)    = z_sbr_i(layer)
      END DO
      WRITE(numout,*) ' truc 3 '

      !----------------
      ! last equation
      !----------------
      WRITE(numout,*) " nlay_i : ", nlay_i
      ztridw(nlay_i,1) = - za(nlay_i) * zb(nlay_i-1)
      WRITE(numout,*) ' truc 4 '
      ztridw(nlay_i,2) = 1.0 + ( za(nlay_i) * ( zb(nlay_i-1) +
     &                   zb(nlay_i) ) )
      ztridw(nlay_i,3) = 0.
      WRITE(numout,*) ' truc 5 '
      zindw(nlay_i)    = z_sbr_i(nlay_i) + za(nlay_i) * zb(nlay_i) 
     &                 * oce_sal
      WRITE(numout,*) ' truc 6 '

      IF ( ln_write ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' -Tridiag terms, winter ... '
         WRITE(numout,*)
         DO layer = 1, nlay_i
            WRITE(numout,*) ' layer : ', layer
            WRITE(numout,*) ' ztridw   : ', ztridw(layer,1), 
     &                      ztridw(layer,2), ztridw(layer,3)
            WRITE(numout,*) ' zindw    : ', zindw(layer)
         END DO
      ENDIF

 40   CONTINUE
!
!-----------------------------------------------------------------------
! 5) Tridiagonal system terms for tracer diffusion equation, summer
!-----------------------------------------------------------------------
!
      DO 50 ji = kideb, kiut 

      DO layer = 1, nlay_i
         ztrids(layer,1) = - ze(layer)
         ztrids(layer,2) = 1.0 + ze(layer)
         ztrids(layer,3) = 0.0
         zinds(layer) = z_sbr_i(layer)
      END DO
      ztrids(1,1) = 0.0

      IF ( ln_write ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' -Tridiag terms, summer ... '
         WRITE(numout,*)
         DO layer = 1, nlay_i
            WRITE(numout,*) ' layer : ', layer
            WRITE(numout,*) ' ztrids   : ', ztrids(layer,1), 
     &                ztrids(layer,2), ztrids(layer,3)
            WRITE(numout,*) ' zinds     : ',zinds(layer)
         END DO
      ENDIF

 50   CONTINUE
      
!
!-----------------------------------------------------------------------
! 6) Partitionning tridiag system between summer and winter
!-----------------------------------------------------------------------
!
      DO 60 ji = kideb, kiut 

      DO indtr = 1, 3
         DO layer = 1, nlay_i
            ztrid(layer,indtr) = zswitchw * ztridw(layer,indtr) + 
     &                           zswitchs * ztrids(layer,indtr)
         END DO ! layer
      END DO ! indtr

      DO layer = 1, nlay_i
         zind(layer) = zswitchw * zindw(layer) + 
     &                 zswitchs * zinds(layer)
      END DO ! layer

      IF ( ln_write ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' -Tridiag terms...'
         WRITE(numout,*)
         WRITE(numout,*) ' zswitchw : ', zswitchw
         WRITE(numout,*) ' zswitchs : ', zswitchs
         DO layer = 1, nlay_i
            WRITE(numout,*) ' layer    : ', layer
            WRITE(numout,*) ' ztrid    : ', ztrid(layer,1), 
     &                        ztrid(layer,2), ztrid(layer,3)
            WRITE(numout,*) ' zind     : ', zind(layer)
         END DO ! layer
      ENDIF

 60   CONTINUE
!
!-----------------------------------------------------------------------
! 7) Solving the tridiagonal system
!-----------------------------------------------------------------------
!
      DO 70 ji = kideb, kiut 

      ! The tridiagonal system is solved with Gauss elimination
      ! Thomas algorithm, from Computational fluid Dynamics, J.D. ANDERSON, 
      ! McGraw-Hill 1984.	

      zindtbis(1) =  zind(1)
      zdiagbis(1) =  ztrid(1,2)

      DO layer = 2, nlay_i
         zdiagbis(layer)  =  ztrid(layer,2) - ztrid(layer,1) *
     &                       ztrid(layer-1,3) / zdiagbis(layer-1)
         zindtbis(layer)  =  zind(layer) - ztrid(layer,1) *
     &                       zindtbis(layer-1) / zdiagbis(layer-1)
      END DO

      ! Recover brine salinity 
      z_sbr_i(nlay_i)     =  zindtbis(nlay_i) / zdiagbis(nlay_i) 
      DO layer = nlay_i - 1 , 1 , -1
         z_sbr_i(layer)   =  ( zindtbis(layer) - ztrid(layer,3)*
     &                       z_sbr_i(layer+1)) / zdiagbis(layer)
      END DO
      ! Recover ice salinity
      DO layer = 1, nlay_i
         sn_i_b(layer)  =  z_sbr_i(layer) * e_i_b(layer)
!        s_i_b(ji,layer)  =  z_sbr_i(layer) * e_i_b(layer)
      END DO

      IF ( ln_write ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' -Solving the tridiagonal system ... '
         WRITE(numout,*)
         WRITE(numout,*) ' zdiagbis: ', ( zdiagbis(layer) , 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' zindtbis: ', ( zdiagbis(layer) , 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' z_sbr_i : ', ( z_sbr_i(layer) , 
     &                   layer = 1, nlay_i )
         WRITE(numout,*) ' sn_i_b   : ', ( sn_i_b(layer) , 
     &                   layer = 1, nlay_i )
      ENDIF

 70   CONTINUE
!
!-----------------------------------------------------------------------
! 8) Mass of salt conserved ?
!-----------------------------------------------------------------------
!
      zerror = 1.0d-8

      DO 80 ji = kideb, kiut 

      ! Final mass of salt
      CALL ice_sal_column( kideb , kiut , z_ms_i_fin , 
     &                    sn_i_b(1:nlay_i),
     &                    deltaz_i_phy, nlay_i, .FALSE. )

      ! Bottom flux ( positive upwards )
      zswitch_open = 0.0
      IF ( e_i_b(nlay_i) .GE. e_tres ) zswitch_open = 1.0
      zfb      = zswitchw * ( - e_i_b( nlay_i ) ! had a minus before
     &                * diff_br(nlay_i) * 2.0
     &                / deltaz_i_phy(nlay_i) * ( z_sbr_i(nlay_i)
     &                  - oce_sal ) ) * zswitch_open
     &                + zswitchs * ( - qsummer * z_sbr_i(nlay_i) )
     &                / ddtb

      fsb = - zfb * rhog / 1000. ! ice-ocean salt flux
      WRITE(numout,*) ' fsb : ', fsb

      ! Surface flux of salt
      zfsu     = zswitchw * 0.0

      ! conservation check
      CALL ice_sal_conserv(kideb,kiut,'ice_sal_diff : ',zerror,
     &                           z_ms_i_ini,z_ms_i_fin,
     &                           zfb   , zfsu  , ddtb)

 80   CONTINUE

      ENDIF ! ln_sal
!
!------------------------------------------------------------------------------|
! End of la sous-routine
      WRITE(numout,*)
      END SUBROUTINE 
