      SUBROUTINE ice_sal_diff_CW(nlay_i,kideb,kiut)

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
     &   zdiagbis           ,    !:
     &   zflux                   !: flux of tracer under layer i

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
     &   zthdiff            ,    !: thermal diffusivity
     &   zgrad_t                 !: temperature gradient in the ice

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
         WRITE(numout,*) ' ** ice_sal_diff_CW : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) ' Cox and weeks based gravity drainage '
         WRITE(numout,*) ' ln_sal = ', ln_sal
      ENDIF

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

      ! gravity drainage parameter (Cox and Weeks 88)
      zeta            = 20.

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
! 2) Gravity drainage as from Cox and Weeks (1988)
!------------------------------------------------------------------------------|
!
      IF ( zswitchw .EQ. 1.0 ) THEN

      DO 20 ji = kideb, kiut

         !----------------------
         ! temperature gradient
         !----------------------
         z_h_lay    = ht_i_b(ji) / nlay_i
         zgrad_t(1) = 2. * ( t_i_b(ji,1) - ( ht_s_b(ji)*t_i_b(ji,1) +
     &                z_h_lay*t_s_b(ji,1) ) / 
     &                ( z_h_lay + ht_s_b(ji) ) ) / z_h_lay

         DO layer = 2, nlay_i - 1
            zgrad_t(layer) = ( t_i_b(ji,layer+1) - t_i_b(ji,layer-1) ) /
     &                       z_h_lay
         END DO
         zgrad_t(nlay_i) = - 2.*( t_i_b(ji,nlay_i) - t_bo_b(ji) ) /
     &                       z_h_lay

         IF ( ln_write ) WRITE(numout,*) ' zgrad_t : ', 
     &                   ( zgrad_t(layer), layer = 1, nlay_i )

         igrd = 1 ! switch for gravity drainage ( 1 if yes )
         zdummysb = - FLOAT(igrd) * rhog / 1000. * ht_i_b(ji) ! salt flux [ kg NaCl.m-2.s-1 ]
         zdsdt = 0.0
         DO layer = nlay_i, 1, -1
            IF ( e_i_b(layer) .LE. e_tres) igrd = 0
            IF ( zgrad_t(layer) .LE. 0 )   igrd = 0 ! temperature gradient must
                                                    ! be directed downwards
            zds_grd = MIN ( 0.0, delta_cw * ( 1.0 - zeta * e_i_b(layer)
     &              * zgrad_t(layer) * ddtb * igrd ) )
            sn_i_b(layer) = s_i_b(ji,layer) + zds_grd
            zdsdt = zdsdt + FLOAT(igrd) * zds_grd / ddtb / FLOAT(nlay_i)
         END DO
         fsb = zdummysb * zdsdt

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
      IF ( zswitchw .EQ. 1. ) THEN
         za(:) = 0.0
         zb(:) = 0.0
      ENDIF

      !----------------------
      ! Summer factors
      !----------------------
      ! ze factors
      DO layer = 1, nlay_i
         ze(layer) = qsummer * zswitchs / 
     &               ( e_i_b(layer) * deltaz_i_phy(layer) )
      END DO ! layer

      ! zf factors
      ! could remove those, they are totally useless!!! ;-)
      DO layer = 1, nlay_i - 1
         zf(layer) = 1./2. * deltaz_i_phy(layer) / 
     &               ( z_i_phy(layer+1)  - z_i_phy(layer) )
      END DO ! layer

      IF ( ln_write ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' -Summer factors '
         WRITE(numout,*) ' ze : ', ( ze(layer), layer = 1, nlay_i ) 
         WRITE(numout,*) ' zf : ', ( zf(layer), layer = 1, nlay_i ) 
      ENDIF

 30   CONTINUE
!
!-----------------------------------------------------------------------
! 4) Tridiagonal system terms for tracer diffusion equation, winter
!-----------------------------------------------------------------------
!
      DO 40 ji = kideb, kiut 

         ztridw(:,:) = 0.

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
            ztrid(layer,indtr) = zswitchs * ztrids(layer,indtr)
         END DO ! layer
      END DO ! indtr

      DO layer = 1, nlay_i
         zind(layer) = zswitchs * zinds(layer)
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
      IF ( zswitchs .EQ. 1.0 ) THEN
      DO layer = 1, nlay_i
         sn_i_b(layer)  =  z_sbr_i(layer) * e_i_b(layer)
!        s_i_b(ji,layer)  =  z_sbr_i(layer) * e_i_b(layer)
      END DO
      ENDIF

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

      zflux(nlay_i) = zfb

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
