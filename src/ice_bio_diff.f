      SUBROUTINE ice_bio_diff(kideb,kiut,nlay_i)

!------------------------------------------------------------------------------!
!                                 ice_bio_diff
!
!    Transport and diffusion of tracers
!    (c) Martin Vancoppenolle, May 2007
!        1.1 Rayleigh number based diffusivity, Oct 2008
!------------------------------------------------------------------------------!
 
      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'
      INCLUDE 'bio.com'

      INTEGER :: 
     &   ji                 ,    ! : horizontal space index
     &   jn                      ! : horizontal space index jn

      REAL(8), DIMENSION( maxnlay ) ::  !: dummy factors for tracer equation
     &   za                 ,    !: winter
     &   zb                 ,   
     &   ze                 ,    !: summer
     &   zind               ,    !: independent term in the tridiag system
     &   zindw              ,    !: independent term in the tridiag system
     &   zinds              ,    !: independent term in the tridiag system
     &   zindtbis           ,    !:
     &   zdiagbis

      REAL(8), DIMENSION(20,3) ::!: dummy factors for tracer equation
     &   ztrid              ,    !: tridiagonal matrix
     &   ztridw             ,    !: tridiagonal matrix, winter
     &   ztrids                  !: tridiagonal matrix, summer

      REAL(8) ::  
     &   zdummy1            ,    !: dummy factors
     &   zdummy2            ,    !: 
     &   zdummy3            ,    !: 
     &   zswitch_open       ,    !: switch for brine network open or not
     &   zswitchw           ,    !: switch for winter drainage 
     &   zswitchs                !: switch for summer drainage

      INTEGER ::  
     &   indtr              ,    !: index of tridiagonal system
     &   iter                    !: time step  index

      CHARACTER(len=4)      ::  
         !conc works, with time step of 3600s and diff of 1.0e-8
         ! flux makes problems
     &   bc = 'conc'             !: Boundary condition 'conc' or 'flux'
! FD additions
      LOGICAL ln_write_bio
      LOGICAL ln_con_bio

      ln_write_bio = .TRUE.
      ln_con_bio   = .TRUE.


!=======================================================================

      WRITE(numout,*) 
      WRITE(numout,*) ' ** ice_bio_diff : '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*)

      DO ji = kideb, kiut
!
!-----------------------------------------------------------------------
! 1) Initialization
!-----------------------------------------------------------------------
!
      IF ( ln_write_bio ) THEN
         WRITE(numout,*) ' Initialization '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~ '
         WRITE(numout,*)
         WRITE(numout,*) ' diff_br_bio : ', ( diff_br_bio(layer), 
     &                   layer = 1, nlay_bio )
      ENDIF

!     beta_sal = flu_beta ! remove that
      
      !---------------
      ! Interpolation 
      !---------------
      CALL ice_bio_interp_phy2bio(kideb,kiut,nlay_i,.FALSE.) 
                             ! interpolation of physical variables
                             ! on the biological grid
                             ! mass of salt, heat content, brine volume, Rb, PAR, PUR

      CALL ice_bio_interp_diffus(kideb,kiut,nlay_i,.TRUE.) 

      !--------------------------------
      ! Brine concentration of tracers
      !--------------------------------
      ! Brine conservation is diffused
      DO jn = 1, ntra_bio
         IF ( flag_diff(jn) .AND. flag_active(jn) ) THEN
            DO jk = 1, nlay_bio
               c_i_bio(jn,jk) = cbu_i_bio(jn,jk) / e_i_bio(jk)
            END DO
         ENDIF
         IF ( flag_adsorb(jn) .AND. flag_active(jn) ) THEN
            DO jk = 1, nlay_bio
               c_i_bio(jn,jk) = 0.
            END DO
         ENDIF
      END DO

      !--------------------
      ! Conservation check
      !--------------------

      CALL ice_bio_column(kideb,kiut,mt_i_bio_init,cbu_i_bio,
     &                    deltaz_i_bio, .FALSE.)

      IF ( ln_write_bio ) THEN
         DO jn = 1, ntra_bio
            WRITE(numout,*) ' mt_i_bio_init : ', mt_i_bio_init(jn)
         END DO
      ENDIF

      ! layer by layer
      DO jn = 1, ntra_bio
         DO jk = 1, nlay_bio
            m_i_bio_init(jn,jk) = cbu_i_bio(jn,jk)*deltaz_i_bio(jk)
         END DO
      END DO

      !----------
      ! Switches
      !----------
      ! summer switch
      zswitchs = MAX( 0.0, SIGN ( 1.0 , t_su_b(ji) - tpw ) ) ! 0 si hiver 1 si ete
      
      zbvmin   = 1.0
      DO layer = 1, nlay_bio
         zbvmin = MIN( e_i_b(layer) , zbvmin ) ! minimum brine volume
      END DO
      IF ( zbvmin .LT. e_tres ) zswitchs = 0.0

      ! winter switch
      zswitchw = 1.0 - zswitchs

      IF ( ln_write_bio ) THEN
         WRITE(numout,*) ' zswitchs : ', zswitchs
         WRITE(numout,*) ' zswitchw : ', zswitchw
         WRITE(numout,*) 
      ENDIF
!
!-----------------------------------------------------------------------
! 2) Compute dummy factors for tracer diffusion equation
!-----------------------------------------------------------------------
!
      DO jn = 1, ntra_bio 

      IF ( ( flag_diff(jn) .OR. flag_adsorb(jn) ) 
     &   .AND. flag_active(jn) ) THEN

      IF ( ln_write_bio ) THEN
         WRITE(numout,*) ' --------------------------------- '
         WRITE(numout,*) '  Diffusion for ', biotr_i_nam(jn)
         WRITE(numout,*) ' --------------------------------- '

         WRITE(numout,*) 
         WRITE(numout,*) ' Dummy factors '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~ '
         WRITE(numout,*)
      ENDIF

      !----------------------
      ! Winter factors
      !----------------------
      ! za factors
!     zdummy1 = ddtb * diff_bio ! CHANGE
      zdummy1 = ddtb            ! CHANGE
      DO layer = 1, nlay_bio
         za(layer) = zdummy1 / ( deltaz_i_bio(layer) * e_i_bio(layer) )
      END DO

      ! zb factors
      DO layer = 1, nlay_bio - 1
         ! interpolate brine volume at the interface between layers
         zdummy1 = ( e_i_bio(layer + 1) - e_i_bio(layer) ) /  
     &             ( z_i_bio(layer + 1) - z_i_bio(layer) )
         zdummy2 = deltaz_i_bio(layer) / 2.0
         zdummy3 = e_i_bio(layer) + zdummy1 * zdummy2
         zswitch_open = 0.0
         ! compute zswitch_open which equals 1 if the brine network is open
         IF ( zdummy3 .GE. e_tres ) zswitch_open = 1.0
!        zb(layer) = zdummy3 * zswitch_open /  ! CHANGE
!    &               ( z_i_bio(layer+1) - z_i_bio(layer) ) * 
         zb(layer) = zdummy3 * zswitch_open /  ! CHANGE
     &               ( z_i_bio(layer+1) - z_i_bio(layer) ) * 
     &               diff_br_bio(layer)
      END DO

      zswitch_open = 0.0
      IF ( e_i_bio(nlay_bio) .GE. e_tres ) zswitch_open = 1.0

      ! Cw fixed boundary condition (imposed cc.)
      IF ( bc .EQ. 'conc' )  
!    &   zb(nlay_bio) = 2. * e_i_bio(nlay_bio) * zswitch_open / 
!    &                  deltaz_i_bio(nlay_bio)
     &   zb(nlay_bio) = 2. * e_i_bio(nlay_bio) * zswitch_open / 
     &                  deltaz_i_bio(nlay_bio) * diff_br_bio(nlay_bio)

      !----------------------
      ! Summer factors
      !----------------------

      !------------------
      ! Percolating flux
      !------------------
!     ! Percolating flow ( rho dh * beta * switch / rhow )
!     qsummer = ( - rhog * MIN ( dh_i_surf(ji) , 0.0 ) 
!    &            - rhon * MIN ( dh_s_tot(ji)  , 0.0 ) )  
!     qsummer = qsummer * beta_sal * zswitchs / 1000.0

      ! ze factors
      DO layer = 1, nlay_bio
         ze(layer) = qsummer * zswitchs / 
     &             ( e_i_bio(layer) * deltaz_i_bio(layer) )
      END DO ! layer

      IF ( ln_write_bio ) THEN
         WRITE(numout,*) ' Winter factors '
         WRITE(numout,*) ' za       : ', ( za (layer),  
     &                   layer = 1, nlay_bio)
         WRITE(numout,*) ' zb       : ', ( zb (layer),  
     &                   layer = 1, nlay_bio)
         WRITE(numout,*)
         WRITE(numout,*) ' Summer factors '
         WRITE(numout,*) ' zswitchs : ', zswitchs
         WRITE(numout,*) ' qsummer  : ', qsummer 
         WRITE(numout,*) ' ze : ', ( ze(layer), layer = 1, nlay_bio ) 
      ENDIF
!
!-----------------------------------------------------------------------
! 3) Tridiagonal system terms for tracer diffusion equation, winter
!-----------------------------------------------------------------------
!
      !----------------
      ! first equation
      !----------------
      ztridw(1,1) = 0.0
      ztridw(1,2) = 1.0 + za(1) * zb(1)
      ztridw(1,3) = - za(1) * zb(1)
      zindw(1)    = c_i_bio(jn,1)

      !-----------------
      ! inner equations
      !-----------------
      DO layer = 2, nlay_bio - 1
         ztridw(layer,1) = - za(layer) * zb(layer-1)
         ztridw(layer,2) = 1.0 + za(layer) * ( zb(layer-1) + 
     &                                         zb(layer) )
         ztridw(layer,3) = - za(layer) * zb(layer)
         zindw(layer)    = c_i_bio(jn,layer)

      END DO

      !----------------
      ! last equation
      !----------------
      ztridw(nlay_bio,1) = - za(nlay_bio) * zb(nlay_bio-1)
      ztridw(nlay_bio,2) = 1.0 + ( za(nlay_bio) * ( zb(nlay_bio-1) + 
     &                     zb(nlay_bio) ) )
      ztridw(nlay_bio,3) = 0.
      zindw(nlay_bio)    = c_i_bio(jn,nlay_bio) + 
     &                     za(nlay_bio) * zb(nlay_bio) * c_skel_bio(jn)

      IF ( i_flux .EQ. 2 ) THEN
          ztridw(nlay_bio,2) = 1.0 + ( za(nlay_bio) * zb(nlay_bio-1) )
          zindw(nlay_bio)    = c_i_bio(jn,nlay_bio)
      ENDIF

      IF ( ln_write_bio ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' Tridiag terms, Winter : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*)
         DO layer = 1, nlay_bio
            WRITE(numout,*) ' layer : ', layer
            WRITE(numout,*) ' ztridw   : ', ztridw(layer,1), 
     &                      ztridw(layer,2),
     &                      ztridw(layer,3)
            WRITE(numout,*) ' zindw     : ',zindw(layer)
         END DO
      ENDIF
!
!-----------------------------------------------------------------------
! 4) Tridiagonal system terms for tracer diffusion equation, winter
!-----------------------------------------------------------------------
!
      DO layer = 1, nlay_bio 
         ztrids(layer,1) = - ze(layer)
         ztrids(layer,2) = 1.0 + ze(layer)
         ztrids(layer,3) = 0.0
         zinds(layer) = c_i_bio(jn,layer)
      END DO
      ztrids(1,1) = 0.0

      IF ( ln_write_bio ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' Tridiag terms, summer : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*)
         DO layer = 1, nlay_bio
            WRITE(numout,*) ' layer : ', layer
            WRITE(numout,*) ' ztrids   : ', ztrids(layer,1), 
     &                      ztrids(layer,2),
     &                      ztrids(layer,3)
            WRITE(numout,*) ' zinds     : ',zinds(layer)
         END DO
      ENDIF
!
!-----------------------------------------------------------------------
! 5) Partitionning tridiag system between summer and winter
!-----------------------------------------------------------------------
!
      DO indtr = 1, 3
         DO layer = 1, nlay_bio
            ztrid(layer,indtr) = zswitchw * ztridw(layer,indtr) + 
     &                           zswitchs * ztrids(layer,indtr)
         END DO ! layer
      END DO ! indtr

      DO layer = 1, nlay_bio
         zind(layer) = zswitchw * zindw(layer) + 
     &                 zswitchs * zinds(layer)
      END DO ! layer

      IF ( ln_write_bio ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' Tridiag terms : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~ '
         WRITE(numout,*)
         DO layer = 1, nlay_bio
            WRITE(numout,*) ' layer    : ', layer
            WRITE(numout,*) ' ztrid    : ', ztrid(layer,1), 
     &                        ztrid(layer,2), ztrid(layer,3)
            WRITE(numout,*) ' zind     : ', zind(layer)
         END DO ! layer
      ENDIF
!
!-----------------------------------------------------------------------
! 6) Solving the tridiagonal system
!-----------------------------------------------------------------------
!
      ! The tridiagonal system is solved with Gauss elimination
      ! Thomas algorithm, from Computational fluid Dynamics, J.D. ANDERSON, 
      ! McGraw-Hill 1984.	
      zindtbis(1) =  zind(1)
      zdiagbis(1) =  ztrid(1,2)
      DO layer = 2, nlay_bio
         zdiagbis(layer)  =  ztrid(layer,2) - ztrid(layer,1) *
     &                       ztrid(layer-1,3) / zdiagbis(layer-1)
         zindtbis(layer)  =  zind(layer) - ztrid(layer,1) *
     &                       zindtbis(layer-1) / zdiagbis(layer-1)
      END DO

      !-----------------------------
      ! Tracer brine concentrations
      !-----------------------------
      c_i_bio(jn,nlay_bio) =  zindtbis(nlay_bio) / zdiagbis(nlay_bio) 
      DO layer = nlay_bio - 1 , 1 , -1
         c_i_bio(jn,layer)  =  (zindtbis(layer) - ztrid(layer,3)*
     &                       c_i_bio(jn,layer+1)) / zdiagbis(layer)
      END DO

      IF ( ln_write_bio ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' Resolution  '
         WRITE(numout,*) ' ~~~~~~~~~~~ '
         WRITE(numout,*)
         DO layer = 1, nlay_bio
            WRITE(numout,*) ' layer    : ', layer
            WRITE(numout,*) ' zdiagbis : ', zdiagbis(layer)
            WRITE(numout,*) ' zindtbis : ', zindtbis(layer)
         END DO
      ENDIF

      ENDIF ! biotr_i_typ .EQ. 'nut' or 'org'

      END DO ! jn
!
!-----------------------------------------------------------------------
! 7) Recover bulk tracer concentrations
!-----------------------------------------------------------------------
!
      !--------------------------------
      ! Tracer bulk ice concentrations
      !--------------------------------
      DO jn = 1, ntra_bio

         IF ( flag_diff(jn) .AND. flag_active(jn) ) THEN

            DO layer = 1, nlay_bio
!           cbun_i_bio(jn,layer) = c_i_bio(jn,layer) * e_i_bio(layer)
               cbu_i_bio(jn,layer) = c_i_bio(jn,layer) * e_i_bio(layer)
            END DO

         ENDIF ! flag_diff

         IF ( flag_adsorb(jn) .AND. flag_active(jn) ) THEN
            DO layer = 1, nlay_bio
               cbu_i_bio(jn,layer) = cbu_i_bio(jn,layer) 
     &                            + c_i_bio(jn,layer) * e_i_bio(layer)
            END DO
         ENDIF ! flag_adsorb

      END DO ! jn

      IF ( ln_write_bio ) THEN
         DO jn = 1, ntra_bio
         IF (      ( flag_diff(jn) .OR. flag_adsorb(jn) ) 
     &        .AND.         flag_active(jn)               ) THEN
  
            WRITE(numout,*) 
            WRITE(numout,*) ' Tracer concentrations '
            WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~ '
            WRITE(numout,*)
            WRITE(numout,*) ' Tracer : ', biotr_i_nam(jn)
            WRITE(numout,*) ' c_i_bio   : ', ( c_i_bio(jn,layer), 
     &                      layer = 1, nlay_bio )
            WRITE(numout,*) ' cbun_i_bio : ', ( cbun_i_bio(jn,layer), 
     &                      layer = 1, nlay_bio )
            WRITE(numout,*)
         ENDIF ! flag_diff
         END DO ! jn
      ENDIF ! ln_write_bio
!
!-----------------------------------------------------------------------
! 8) Conservation check
!-----------------------------------------------------------------------
!
      IF ( ln_con_bio ) THEN

      IF ( ln_write_bio ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Conservation check : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*)
      ENDIF ! ln_write_bio

      CALL ice_bio_column(kideb,kiut,mt_i_bio_final,cbu_i_bio,
     &                    deltaz_i_bio, .FALSE.)

      zswitch_open = 0.0
      IF ( e_i_bio(nlay_bio) .GE. e_tres ) zswitch_open = 1.0

      zerror = 1.0d-8

      DO jn = 1, ntra_bio

         IF ( flag_diff(jn) .AND. flag_active(jn) ) THEN

         IF ( ln_write_bio ) THEN
            WRITE(numout,*) ' mt_i_bio_final : ', mt_i_bio_final(jn)
         ENDIF

         ! Bottom flux ( positive upwards )
!        f_bo_tra(jn) = zswitchw * ( - e_i_bio( nlay_bio ) ! had a minus before
!    &                * diff_bio * 2.0
!    &                / deltaz_i_bio(nlay_bio) * ( c_i_bio(jn,nlay_bio)
!    &                  - c_skel_bio(jn) ) ) * zswitch_open
!    &                + zswitchs * ( - qsummer * c_i_bio(jn,nlay_bio) )
!    &                / ddtb
      ! info dans ice_sal_diff
!     zswitch_open = 0.0
!     IF ( e_i_b(nlay_i) .GE. e_tres ) zswitch_open = 1.0
!     zfb      = zswitchw * ( - e_i_b( nlay_i ) ! had a minus before
!    &                * diff_br(nlay_i) * 2.0
!    &                / deltaz_i_phy(nlay_i) * ( z_sbr_i(nlay_i)
!    &                  - oce_sal ) ) * zswitch_open
!    &                + zswitchs * ( - qsummer * z_sbr_i(nlay_i) )
!    &                / ddtb
!    !
! traduction en bio
      zswitch_open = 0.0
      IF ( e_i_bio(nlay_bio) .GE. e_tres ) zswitch_open = 1.0
      zfb      = zswitchw * ( - e_i_bio( nlay_bio ) ! had a minus before
     &                * diff_br_bio(nlay_bio) * 2.0
     &                / deltaz_i_bio(nlay_bio) * ( c_i_bio(jn,nlay_bio)
     &                  - c_skel_bio(jn) ) ) * zswitch_open
     &                + zswitchs * ( - qsummer * c_i_bio(jn,nlay_bio) )
     &                / ddtb
      f_bo_tra(jn) = zfb
      fcb(jn) = - zfb ! ice-ocean tracer flux
      IF ( i_flux .EQ. 2 ) fcb(jn) = 0.

      ! surface flux
      f_su_tra(jn) = zswitchw * 0.0
     &             + zswitchs * ( qsummer * c_s_bio(jn) )
     &             / ddtb

      ! fcb_max = idealized flux if c_i_bio(nlay_bio) .EQ. 0.0 all the time
      fcb_max(jn)= - zswitchw * ( - e_i_bio( nlay_bio ) ! had a minus before
     &                * diff_br_bio(nlay_bio) * 2.0
!    &                / deltaz_i_bio(nlay_bio) * ( c_i_bio(jn,nlay_bio)
     &                / deltaz_i_bio(nlay_bio) * ( 0.0                 
     &                  - c_skel_bio(jn) ) ) * zswitch_open
     &                - zswitchs * ( - qsummer * c_i_bio(jn,nlay_bio) )
     &                / ddtb
!     fcb_max(jn)= zswitchw * ( - e_i_bio( nlay_bio ) ! had a minus before
!    &                * diff_br_bio(nlay_bio) * 2.0
!    &                / deltaz_i_bio(nlay_bio) * ( 0.0
!    &                  - c_skel_bio(jn) ) ) * zswitch_open / ddtb


         IF ( ln_write_bio ) THEN
            WRITE(numout,*) ' f_bo_tra : ', f_bo_tra(jn)
            WRITE(numout,*) ' f_su_tra : ', f_su_tra(jn)
            WRITE(numout,*) ' zswitchs : ', zswitchs
            WRITE(numout,*) ' qsummer  : ', qsummer 
            WRITE(numout,*) ' c_i_bio  : ', c_i_bio(jn,nlay_bio)
            WRITE(numout,*) ' mt_i_bio_init  : ', mt_i_bio_init(jn)
            WRITE(numout,*) ' mt_i_bio_final : ', mt_i_bio_final(jn)

         ENDIF ! ln_write_bio

         ENDIF ! flag_diff

      END DO ! jn

      CALL ice_bio_conserv(kideb,kiut,'ice_bio_diff : ',zerror,
     &                           mt_i_bio_init,mt_i_bio_final,
     &                           f_bo_tra, f_su_tra, ddtb)

      ENDIF ! ln_con_bio
!
!-----------------------------------------------------------------------
! 8) Layer by layer conservation check
!-----------------------------------------------------------------------
!

!     WRITE(numout,*)
!     WRITE(numout,*) ' Conservation check, layer by layer : '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
!     WRITE(numout,*)

!     DO jn = 1, ntra_bio

!     IF ( flag_diff(jn) .AND. flag_active(jn) ) THEN

!     !--------------------
!     ! inner layer fluxes
!     !--------------------
!     WRITE(numout,*) ' layer by layer conservation check '
!     fdiff(jn,0) = 0.0 
!     DO layer = 1, nlay_bio - 1
!        ! interpolate brine volume at the interface between layers
!        zdummy1 = ( e_i_bio(layer + 1 ) - e_i_bio(layer) ) /  
!    &             ( z_i_bio(layer + 1) - z_i_bio(layer) )
!        zdummy2 = deltaz_i_bio(layer) / 2.0
!        zdummy3 = e_i_bio(layer) + zdummy1 * zdummy2
!        zswitch_open = 0.0
!        ! compute zswitch_open which equals 1 if the brine network is open
!        IF ( zdummy3 .GE. e_tres ) zswitch_open = 1.0
!        fdiff(jn,layer) = - zswitch_open * zdummy3 * diff_bio * 
!    &             ( c_i_bio(jn,layer+1) - c_i_bio(jn,layer) )
!        WRITE(numout,*) ' layer        : ', layer
!        WRITE(numout,*) ' fdiff(layer) : ', fdiff(jn,layer)
!        WRITE(numout,*) ' zswitch_open : ', zswitch_open
!        WRITE(numout,*) ' zdummy3      : ', zdummy3
!        WRITE(numout,*) ' diff_bio     : ', diff_bio
!     END DO

!     DO layer = 1, nlay_bio - 1
!        WRITE(numout,*) ' fdiff(layer)    : ', fdiff(jn,layer)
!     END DO

!     !------------------
!     ! lower layer flux
!     !------------------

!     zswitch_open = 0.0
!     IF ( e_i_bio(nlay_bio) .GE. e_tres ) zswitch_open = 1.0
!     fdiff(jn,nlay_bio) = - 2.0 * e_i_bio(nlay_bio) * diff_bio *
!    &               ( c_skel_bio(jn) - c_i_bio(jn,nlay_bio) )
!    &               / deltaz_i_bio(nlay_bio) * zswitch_open
!     WRITE(numout,*) ' c_skel_bio : ', c_skel_bio(jn)
!     WRITE(numout,*) ' c_i_bio(N) : ', c_i_bio(jn,nlay_bio)
!     WRITE(numout,*) ' deltaz_i_bio : ', deltaz_i_bio(nlay_bio)
!     WRITE(numout,*) ' zswitch_open : ', zswitch_open

!     DO layer = 1, nlay_bio
!        WRITE(numout,*) ' fdiff(layer) : ', fdiff(jn,layer)
!     END DO
!     !------------------
!     ! tracer content
!     !------------------
!     ! layer by layer
!     DO layer = 1, nlay_bio
!        m_i_bio_final(jn,layer) = cbu_i_bio(jn,layer)*
!    &                             deltaz_i_bio(layer)
!     END DO

!     DO layer = 1, nlay_bio
!        WRITE(numout,*) ' fdiff(layer) : ', fdiff(jn,layer)
!     END DO

!     DO layer = 1, nlay_bio
!        IF ( zswitchw .GT. 0.9 ) THEN
!           WRITE(numout,*)
!           WRITE(numout,*) 'layer : ', layer
!           zdiff = m_i_bio_final(jn,layer) - m_i_bio_init(jn,layer)
!           zfluxdt =  ( - fdiff(jn,layer) + fdiff(jn,layer-1) ) * ddtb
!           WRITE(numout,*) ' m_i_bio_init : ', m_i_bio_init(jn,layer) 
!           WRITE(numout,*) ' m_i_bio_final: ', m_i_bio_final(jn,layer)
!           WRITE(numout,*) ' difference   : ', zdiff
!           WRITE(numout,*) ' fdiff(jk)    : ', fdiff(jn,layer)
!           WRITE(numout,*) ' fdiff(jk-1)  : ', fdiff(jn,layer-1)
!           WRITE(numout,*) ' zfluxdt      : ', zfluxdt
!           WRITE(numout,*) ' error        : ', zdiff-zfluxdt
!        ENDIF
!     END DO

!     ENDIF ! IF biotr_i_typ EQ 'nut' or 'org'

!     END DO ! jn

      END DO ! ji

      IF ( ln_write_bio ) THEN

         WRITE(numout,*)
         WRITE(numout,*) ' *** After diffusion of tracers *** '
         WRITE(numout,*) '     model output '

         DO jn = 1, ntra_bio
            IF ( flag_active(jn) ) THEN
               WRITE(numout,*) ' biotr_i_nam : ', biotr_i_nam(jn)
               WRITE(numout,*) ' cbun_i_bio : ', ( cbun_i_bio(jn, jk), 
     &                         jk = 1, nlay_bio )
            ENDIF ! flag_active
         END DO ! jn
         WRITE(numout,*)

      ENDIF ! ln_write_bio

      WRITE(numout,*)
      WRITE(numout,*) ' End of ice_bio_diff '
      WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!=============================================================================!
!-- End of ice_bio_diff --
 
      END
