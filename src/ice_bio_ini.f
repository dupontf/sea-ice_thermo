      SUBROUTINE ice_bio_ini(kideb,kiut,nlay_i)

! This routine reads namelist and 
! initializes biogeochemistry routines
! (c) Martin Vancoppenolle, May 2007
!
! version: 2.21
 
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
     &  numbio = 500  ! : reference number for bio.param

!     dimension ain(imax,jmax),zinda(imax,jmax),ifvt(imax,jmax)
!
      WRITE(numout,*) ' ** ice_bio_ini : '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~~ '
      WRITE(numout,*) 
!
!-----------------------------------------------------------------------------!
! 1) Reads namelist
!-----------------------------------------------------------------------------!
!
      WRITE(numout,*) ' Biological parameters ... '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*)

      OPEN( unit = numbio , file='bio.param', status='old' )
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      ! Process switches for biogeochemistry
      WRITE(numout,*), " 0a"
      READ(numbio,*) c_mod
      READ(numbio,*)
      WRITE(numout,*), " 0a"
      READ(numbio,*) i_flux
      WRITE(numout,*), " 0c "
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*) ln_trbo     ! activate basal entrapment for tracers
      WRITE(numout,*), " 1 "
      READ(numbio,*)
      READ(numbio,*) ln_trsi     ! activate snow-ice entrapment for tracers
      WRITE(numout,*), " 2 "
      READ(numbio,*)
      READ(numbio,*) ln_trdiff   ! activate tracer diffusion
      WRITE(numout,*), " 3 "
      READ(numbio,*)
      READ(numbio,*) ln_trremp   ! activate tracer remapping
      WRITE(numout,*), " 4 "
      WRITE(numout,*) 'ln_trremp : ', ln_trremp
      READ(numbio,*)
      READ(numbio,*) ln_lim_dsi  ! activate DSi limitation
      WRITE(numout,*), " 5 "
      WRITE(numout,*) 'ln_lim_dsi: ', ln_lim_dsi
      READ(numbio,*)
      READ(numbio,*) ln_lim_no3  ! activate NO3 limitation
      WRITE(numout,*), " 5b"
      READ(numbio,*)
      READ(numbio,*) ln_lim_lig  ! activate light limitation
      WRITE(numout,*), " 6 "
      READ(numbio,*)
      READ(numbio,*) ln_lim_tem  ! activate temperature inhibition
      WRITE(numout,*), " 7 "
      READ(numbio,*)
      READ(numbio,*) ln_lim_sal  ! activate brine salinity inhibition
      WRITE(numout,*), " 8 "
      READ(numbio,*)
      READ(numbio,*) ln_lys      ! activate lysis
      WRITE(numout,*), " 9 "
      READ(numbio,*)
      READ(numbio,*) ln_rem      ! remineralization or not
      WRITE(numout,*), " 10 "
      READ(numbio,*)
      READ(numbio,*) ln_syn      ! diatom synthesis or not
      WRITE(numout,*), " 11 "
      READ(numbio,*)
      READ(numbio,*) ln_dormant  ! diatom dormancy or not
      WRITE(numout,*), " 12 "
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      ! initial values
      READ(numbio,*) cdsi_ini    ! initial DSi bulk concentration (mmol m-3)
      WRITE(numout,*), " 13 "
      READ(numbio,*)
      READ(numbio,*) cdfe_ini    ! initial DFe bulk concentration (micrmol m-3)
      WRITE(numout,*), " 14 "
      READ(numbio,*)
      READ(numbio,*) cdar_ini    ! initial DAR bulk concentration (mmol C m-3)
      WRITE(numout,*), " 15 "
      READ(numbio,*)
      READ(numbio,*) cdaf_ini    ! initial DAF bulk concentration (mmol C m-3)
      WRITE(numout,*), " 16 "
      READ(numbio,*)
      READ(numbio,*) ctoc_ini    ! initial TOC bulk concentration (mmol C m-3)
      WRITE(numout,*), " 17 "
      READ(numbio,*)
      READ(numbio,*) cpfe_ini    ! initial PFe bulk concentration (micrmol m-3)
      WRITE(numout,*), " 18 "
      READ(numbio,*)
      READ(numbio,*) cno3_ini    ! initial NO3 bulk concentration (micrmol m-3)
      WRITE(numout,*), " 18b "
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      ! biophysical parameters
      READ(numbio,*) frtr_si_bio ! fractionation coeff in snow ice
      WRITE(numout,*), " 19 "
      READ(numbio,*)
      READ(numbio,*) stickfac    ! sticking constant
      WRITE(numout,*), " 20 "
      READ(numbio,*)
      READ(numbio,*) par_fsw     ! conversion factor between par and fsw (muE/J)
      WRITE(numout,*), " 21 "
      READ(numbio,*)
      READ(numbio,*) diff_da     ! diffusivity for diatoms (SL/BA only)
      WRITE(numout,*), " 21b "
      READ(numbio,*)
      READ(numbio,*) diff_nut    ! diffusivity for nuts (SL/BA only)
      WRITE(numout,*), " 21c "
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      ! physiological parameters
      READ(numbio,*) pp_presc    ! prescribed primary production (mmol C m-3 s-1)
      WRITE(numout,*), " 22 "
      READ(numbio,*)
      READ(numbio,*) chla_c      ! chlorophyll to carbon ratio (g/mol)
      WRITE(numout,*), " 22b"
      READ(numbio,*)
      READ(numbio,*) si_c        ! silicate to carbon ratio (mol/mol)
      WRITE(numout,*), " 23 "
      READ(numbio,*)
      READ(numbio,*) no3_c       ! nitrate to carbon ratio (mol/mol)
      WRITE(numout,*), " 23b "
      READ(numbio,*)
      READ(numbio,*) mumax_bio   ! maximal specific growth rate (s-1)
      WRITE(numout,*), " 24 "
      READ(numbio,*)
      READ(numbio,*) klys_bio    ! autolysis rate (hours-1)
      WRITE(numout,*), " 25 "
      READ(numbio,*)
      READ(numbio,*) kmin_bio    ! fraction of diatom loss that is remineralized in Si
      WRITE(numout,*), " 25b"
      READ(numbio,*)
      READ(numbio,*) khs_si_bio  ! half saturatio ct for Si uptake (mmol m-3)
      WRITE(numout,*), " 26 "
      READ(numbio,*)
      READ(numbio,*) khs_n_bio   ! half saturatio ct for NO3 uptake (mmol m-3)
      WRITE(numout,*), " 26b "
      READ(numbio,*)
      READ(numbio,*) ek_bio      ! light adaptation (micrmol quanta m-2 s-1)
      WRITE(numout,*), " 27 "
      READ(numbio,*)
      READ(numbio,*) astar_alg   ! specific att coeff ( m-1 / ( mg chla m-3 ) )
      WRITE(numout,*), " 28 "
      READ(numbio,*)
      READ(numbio,*) fdet_alg    ! fraction of attenuation due to detrital mater compared to algal attenuation (Arrigo)
      WRITE(numout,*), " 28b "
      READ(numbio,*)
      READ(numbio,*) rg_bio      ! prescribed primary production ( K-1 )
      WRITE(numout,*), " 29 "
      READ(numbio,*)
      READ(numbio,*) lim_sal_swi ! Type of salinity limitation 1=Arrigo, 2=4th order gaussian, 3=sine, 4=gaussian, 5=gaussian with 3rd order arg
      READ(numbio,*)
      WRITE(numout,*), " 29b"
      READ(numbio,*) lim_sal_wid ! Width of salinity limitation (case 4)
      READ(numbio,*)
      WRITE(numout,*), " 29c"
      READ(numbio,*) lim_sal_smax! Salinity at which limitation equals 1 (case 4)
      WRITE(numout,*), " 29d"
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      READ(numbio,*)
      ! parameters that are not yet used but might be for future simco version
      READ(numbio,*) kmax_bio    ! maximal specific photo (hours-1)
      WRITE(numout,*), " 31 "
      READ(numbio,*)
      READ(numbio,*) maint_bio   ! maintenance rate (hours-1)
      WRITE(numout,*), " 32 "
      READ(numbio,*)
      READ(numbio,*) khs_sr_bio  ! half saturation ct for SR synthesis (-)
      WRITE(numout,*), " 33 "
      READ(numbio,*)
      READ(numbio,*) ecg_bio     ! energetic cost growth (-)
      WRITE(numout,*), " 34 "
      READ(numbio,*)
      READ(numbio,*) ke_bio      ! EPS excrestion constant (hours-1)
      WRITE(numout,*), " 35 "
      READ(numbio,*)
      READ(numbio,*) khs_fe_bio  ! half saturation ct for Fe uptake (micrmol m-3)
      WRITE(numout,*), " 36 "
      READ(numbio,*)
      READ(numbio,*) kdis_bio    ! frustule dissolution rate (hours-1)
      WRITE(numout,*), " 37 "
      READ(numbio,*)
      READ(numbio,*) kads_bio    ! iron adsorption rate (hours-1)
      WRITE(numout,*), " 38 "

      CLOSE( numbio )

      ! Writes in the output file
      WRITE(numout,*) ' c_mod     : ', c_mod
      WRITE(numout,*) ' i_flux    : ', i_flux
      WRITE(numout,*) ' ln_trbo   : ', ln_trbo
      WRITE(numout,*) ' ln_trsi   : ', ln_trsi
      WRITE(numout,*) ' ln_trdiff : ', ln_trdiff
      WRITE(numout,*) ' ln_trremp : ', ln_trremp
      WRITE(numout,*) ' ln_lim_dsi: ', ln_lim_dsi
      WRITE(numout,*) ' ln_lim_no3: ', ln_lim_no3
      WRITE(numout,*) ' ln_lim_lig: ', ln_lim_lig
      WRITE(numout,*) ' ln_lim_tem: ', ln_lim_tem
      WRITE(numout,*) ' ln_lim_sal: ', ln_lim_sal
      WRITE(numout,*) ' ln_lys    : ', ln_lys
      WRITE(numout,*) ' ln_rem    : ', ln_rem
      WRITE(numout,*) ' ln_syn    : ', ln_syn
      WRITE(numout,*) ' ln_dormant: ', ln_dormant
      WRITE(numout,*)
      WRITE(numout,*) ' cdsi_ini  : ', cdsi_ini
      WRITE(numout,*) ' cdfe_ini  : ', cdfe_ini
      WRITE(numout,*) ' cdar_ini  : ', cdar_ini
      WRITE(numout,*) ' cdaf_ini  : ', cdaf_ini
      WRITE(numout,*) ' ctoc_ini  : ', ctoc_ini
      WRITE(numout,*) ' cpfe_ini  : ', cpfe_ini
      WRITE(numout,*) ' cno3_ini  : ', cno3_ini
      WRITE(numout,*)
      WRITE(numout,*) ' frtr_si_bio:', frtr_si_bio
      WRITE(numout,*) ' stickfac  : ', stickfac
      WRITE(numout,*) ' par_fsw   : ', par_fsw
      WRITE(numout,*) ' diff_da   : ', diff_da
      WRITE(numout,*) ' diff_nut  : ', diff_nut
      WRITE(numout,*)
      WRITE(numout,*) ' pp_presc  : ', pp_presc
      WRITE(numout,*) ' chla_c    : ', chla_c
      WRITE(numout,*) ' si_c      : ', si_c
      WRITE(numout,*) ' no3_c     : ', no3_c
      WRITE(numout,*) ' mumax_bio : ', mumax_bio   
      WRITE(numout,*) ' klys_bio  : ', klys_bio
      WRITE(numout,*) ' khs_si_bio: ', khs_si_bio
      WRITE(numout,*) ' ek_bio    : ', ek_bio
      WRITE(numout,*) ' astar_alg : ', astar_alg
      WRITE(numout,*) ' fdet_alg  : ', fdet_alg
      WRITE(numout,*) ' rg_bio    : ', rg_bio
      WRITE(numout,*) ' lim_sal_swi  : ', lim_sal_swi
      WRITE(numout,*) ' lim_sal_wid  : ', lim_sal_wid
      WRITE(numout,*) ' lim_sal_smax : ', lim_sal_smax
      WRITE(numout,*)
      WRITE(numout,*) ' kmin_bio  : ', kmin_bio 
      WRITE(numout,*) ' kmax_bio  : ', kmax_bio
      WRITE(numout,*) ' maint_bio : ', maint_bio
      WRITE(numout,*) ' khs_sr_bio: ', khs_sr_bio
      WRITE(numout,*) ' ecg_bio   : ', ecg_bio
      WRITE(numout,*) ' ke_bio    : ', ke_bio
      WRITE(numout,*) ' khs_fe_bio: ', khs_fe_bio
      WRITE(numout,*) ' kdis_bio  : ', kdis_bio 
      WRITE(numout,*) ' kads_bio  : ', kads_bio 
      WRITE(numout,*) ' fe_c      : ', fe_c
      WRITE(numout,*) ' fe_c_eps  : ', fe_c_eps
      WRITE(numout,*) 
!
!-----------------------------------------------------------------------------!
! 2) Converts units
!-----------------------------------------------------------------------------!
!
      WRITE(numout,*) ' Converting units of params... '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*)
      ! moves from hours-1 to seconds-1
!     mumax_bio  = mumax_bio  / 3600.0
!     kmax_bio   = kmax_bio   / 3600.0
!     maint_bio  = maint_bio  / 3600.0
!     klys_bio   = klys_bio   / 3600.0
!     ke_bio     = ke_bio     / 3600.0
!     kmin_bio   = kmin_bio   / 3600.0
!     kdis_bio   = kdis_bio   / 3600.0
!     kads_bio   = kads_bio   / 3600.0
!
!-----------------------------------------------------------------------------!
! 3) Create grid and interpolate physical variables
!-----------------------------------------------------------------------------!
!
      WRITE(numout,*) ' Creates grids and do interpolation... '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*)

      CALL ice_bio_grid(kideb,kiut,.TRUE.)   ! compute the biological grid
      CALL ice_bio_interp_phy2bio(kideb,kiut,nlay_i, .FALSE.) ! Interpolation of physical
                                             ! variables on the bio grid
!
!-----------------------------------------------------------------------------!
! 4) Initialize bulk concentrations
!-----------------------------------------------------------------------------!
!
      WRITE(numout,*) ' Initialize tracers... '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*)

      !------------------------
      ! DSi : Dissolved silica
      !------------------------
      jn = 1
      flag_active(jn) = .true.  ! Activated
      flag_adsorb(jn) = .false. ! Not adsorbed
      ! Units = ( mmole Si m-3 ) or ( µmol Si l-1)
      biotr_i_nam(jn) = 'DSi' ! Name of the tracer
      biotr_i_typ(jn) = 'nut' ! Type = Nutrient
      biotr_i_uni(jn) = 'mmol_m_3' ! Units
      IF ( flag_active(jn) ) THEN
         DO layer = 1, nlay_bio
            cbu_i_bio(jn,layer) = cdsi_ini
         END DO
         IF ( c_mod .NE. 'ML' ) cbu_i_bio(jn,1:nlay_bio-1) = 0.
      ELSE
         cbu_i_bio(jn,:) = 0.0
      ENDIF
      flag_diff(jn) = .true. ! DSi is diffused

!     !----------------------
!     ! DFe : Dissolved iron
!     !----------------------
      jn = 2
      ! Units = ( nmole Fe m-3 ) or ( micrmol Fe l-1 )
      flag_active(jn) = .false. ! Not activated
      flag_adsorb(jn) = .false. ! Not adsorbed
      biotr_i_nam(jn) = 'DFe'   ! Name of the tracer
      biotr_i_typ(jn) = 'nut'   ! Type = Nutrient
      biotr_i_uni(jn) = 'mumol_m3' ! Units
      IF ( flag_active(jn) ) THEN
         DO layer = 1, nlay_bio
            cbu_i_bio(jn,layer) = cdfe_ini
         END DO
         IF ( c_mod .NE. 'ML' ) cbu_i_bio(jn,1:nlay_bio-1) = 0.
      ELSE
         cbu_i_bio(jn,:) = 0.0
      ENDIF
      flag_diff(jn) = .false.

      !------------------------------------
      ! DAR : Reserve products for diatoms
      !------------------------------------
      jn = 3
      ! Units = ( mmole C m-3 )
      flag_active(jn) = .false. ! Not activated
      flag_adsorb(jn) = .false. ! Not adsorbed
      biotr_i_nam(jn) = 'DSR'   ! Name of the tracer
      biotr_i_typ(jn) = 'alg'   ! Type = Nutrient
      biotr_i_uni(jn) = 'mmol_m_3' ! Units
      DO layer = 1, nlay_bio
         cbu_i_bio(jn,layer) = 0.0
      END DO
      layer = nlay_bio
         IF ( flag_active(jn) ) THEN
            cbu_i_bio(jn,layer) = cdar_ini
         ELSE
            cbu_i_bio(jn,layer) = 0.0
         ENDIF
      flag_diff(jn) = .false.

      !--------------------------------------
      ! DAF : Funtional products for diatoms
      !--------------------------------------
      jn = 4
      ! Units = ( mmole C m-3 )
      flag_active(jn) = .true.     ! Activated ?
      flag_adsorb(jn) = .false.    ! Adsorbed ?
      flag_diff(jn)   = .false.    ! Diffused ?
      biotr_i_nam(jn) = 'DAF'      ! Name of the tracer
      biotr_i_typ(jn) = 'alg'      ! Type = algal
      biotr_i_uni(jn) = 'mmol_m_3' ! Units
      DO layer = 1, nlay_bio
         cbu_i_bio(jn,layer) = 0.0
      END DO

      IF ( flag_active(jn) ) THEN
         DO layer = 1, nlay_bio
            cbu_i_bio(jn,layer) = cdaf_ini
         END DO
         IF ( c_mod .NE. 'ML' ) cbu_i_bio(jn,1:nlay_bio-1) = 0.
      ELSE
         cbu_i_bio(jn,:) = 0.0 
      ENDIF

      !----------------------------
      ! TOC : Total organic carbon
      !----------------------------
      jn = 5
      ! Units = ( mmole C m-3 )
      flag_active(jn) = .false. ! Not activated
      flag_adsorb(jn) = .false. ! Not adsorbed
      biotr_i_nam(jn) = 'TOC' ! Name of the tracer
      biotr_i_typ(jn) = 'org' ! Type = Nutrient
      biotr_i_uni(jn) = 'mmol_m_3' ! Units
      DO layer = 1, nlay_bio
         IF ( flag_active(jn) ) THEN
            cbu_i_bio(jn,layer) = ctoc_ini
         ELSE
            cbu_i_bio(jn,layer) = 0.0
         ENDIF
      END DO
      flag_diff(jn) = .false.

      !------------------------
      ! PFe : Particulate iron
      !------------------------
      jn = 6
      ! Units = ( micrmole C m-3 ) or ( nmol Fe l-1 )
      flag_active(jn) = .false. ! Not activated
      flag_adsorb(jn) = .false. ! Not adsorbed
      biotr_i_nam(jn) = 'PFe' ! Name of the tracer
      biotr_i_typ(jn) = 'org' ! Type = Nutrient
      biotr_i_uni(jn) = 'mumol_m3' ! Units
      DO layer = 1, nlay_bio
         IF ( flag_active(jn) ) THEN
            cbu_i_bio(jn,layer) = cpfe_ini
         ELSE
            cbu_i_bio(jn,layer) = 0.0
         ENDIF
      END DO
      flag_diff(jn) = .false.

      !--------------------------
      ! NO3 : dissolved nitrates
      !--------------------------
      jn = 7
      ! Units = ( mmol N m-3 )
      flag_active(jn) = .true.  ! Not activated
      flag_adsorb(jn) = .false. ! Not adsorbed
      flag_diff(jn) = .true.
      biotr_i_nam(jn) = 'NO3' ! Name of the tracer
      biotr_i_typ(jn) = 'nut' ! Type = Nutrient
      biotr_i_uni(jn) = 'mmol_m3' ! Units
      IF ( flag_active(jn) ) THEN
         DO layer = 1, nlay_bio
            cbu_i_bio(jn,layer) = cno3_ini
         END DO
         IF ( c_mod .NE. 'ML' ) cbu_i_bio(jn,1:nlay_bio-1) = 0.
      ELSE
         cbu_i_bio(jn,:) = 0.0
      ENDIF
!
!-----------------------------------------------------------------------------!
! 5) Brine concentrations
!-----------------------------------------------------------------------------!
!
      WRITE(numout,*) ' Compute brine concentrations... '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*)

      DO jn = 1, ntra_bio

         ! Dissolved tracers
!        IF ( ( biotr_i_nam(jn) .EQ. 'DSi' ) .OR.
!    &        ( biotr_i_nam(jn) .EQ. 'DFe' ) ) THEN
!        IF ( flag_active(jn) .AND. flag_diff(jn) ) THEN
!           DO jk = 1, nlay_bio
!              c_i_bio(jn,jk)    = cbu_i_bio(jn,jk) / e_i_bio(jk) 
!           END DO
!        ENDIF

         ! Particulate tracers
!        IF ( ( biotr_i_nam(jn) .EQ. 'PFe' ) .OR. 
!    &        ( biotr_i_nam(jn) .EQ. 'DAR' ) .OR. 
!    &        ( biotr_i_nam(jn) .EQ. 'DAF' ) ) THEN
!        IF ( flag_active(jn) .AND. ( .NOT. flag_diff(jn) ) ) THEN
!           DO jk = 1, nlay_bio
!              c_i_bio(jn,jk) = 0.0
!           END DO
!        ENDIF

         ! TOC

      END DO
!
!-----------------------------------------------------------------------------!
! 6) Seawater, skeletal layer concentrations and chlorophyll a
!-----------------------------------------------------------------------------!
!
      WRITE(numout,*) ' Seawater / skeltal concentrations, chla... '
!     WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*)

      !----------------------------------------------
      ! Compute seawater and skeletal concentrations
      !----------------------------------------------
      ! Initialization follows V1 specifications
      
      c_w_bio(1) = 51.25 ! DSi (Zemmelink et al 08) SIMBA (60)
!     c_w_bio(1) = 5.000 ! sensitivity run
      c_w_bio(2) = 1.97  ! DFe 
      c_w_bio(3) = 0.03  ! DSR
      c_w_bio(4) = 0.60  ! DAF (Zemmelink et al 08)
      c_w_bio(4) = 0.05  ! GRL runs
!     c_w_bio(4) = 0.10  ! sensitivity run
      c_w_bio(5) = 54.05 ! TOC
      c_w_bio(6) = 0.33  ! PFe
      c_w_bio(7) = 26.9  ! NO3 (Zemmelink et al 08)

      DO jn = 1, ntra_bio
         IF ( flag_active(jn) ) THEN
            c_skel_bio(jn) = c_w_bio(jn)
         ELSE
            c_skel_bio(jn) = 0.0
         ENDIF
      END DO

      DO jn = 1, ntra_bio
         IF ( flag_active(jn) ) THEN
            WRITE(numout,*) ' Tracer     : ', biotr_i_nam(jn)
            WRITE(numout,*) ' ~~~~~~~      '
            WRITE(numout,*) ' Type       : ', biotr_i_typ(jn)
            WRITE(numout,*) ' cbu_i_bio  : ', ( cbu_i_bio(jn,jk), 
     &                   jk = 1, nlay_bio )
            WRITE(numout,*) ' c_w_bio    : ', c_w_bio(jn)
            WRITE(numout,*) ' c_skel_bio : ', c_skel_bio(jn)
            WRITE(numout,*)
         ENDIF
      END DO ! jn

      !---------------
      ! Chlorophyll a
      !---------------
      ! Units, mg m-3 (micro g.l-1)
      DO layer = 1, nlay_bio
         chla_i_bio(layer) = cbu_i_bio(4,layer) / 2.0
      END DO

      WRITE(numout,*) ' Chla '
      WRITE(numout,*) ' ~~~~ '
      WRITE(numout,*) ' chla_i_bio : ', 
     &                ( chla_i_bio(layer), layer = 1, nlay_bio )

      WRITE(numout,*)
      WRITE(numout,*) '    *** Initial values *** '
      WRITE(numout,*) '    model output '
      DO jn = 1, ntra_bio
         IF ( flag_active(jn) ) THEN
           WRITE(numout,*) ' biotr_i_nam : ', biotr_i_nam(jn)
           WRITE(numout,*) ' cbu_i_bio : ', ( cbu_i_bio(jn, jk), jk = 1,
     &                                        nlay_bio )
         ENDIF
      END DO
      WRITE(numout,*)

      WRITE(numout,*)
      WRITE(numout,*) ' End of ice_bio_ini '
      WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!
!-----------------------------------------------------------------------------!
!-- End of ice_bio_ini --
      RETURN
 
      END
