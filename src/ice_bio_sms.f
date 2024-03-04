      SUBROUTINE ice_bio_sms(nlay_i,kideb,kiut)
!
!-----------------------------------------------------------------------------!
!                           *** ice_bio_sms ***
!     Biological sources minus sinks
!     (c) Martin Vancoppenolle, UCL-ASTR, June 2007
!     equations from Christiane and Sylvie
!     plus valuable input from Hugues and Anne
!     version : source_2.21
!-----------------------------------------------------------------------------!
!
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'thermo.com'
      include 'bio.com'

      INTEGER :: zcase
      DIMENSION zalpha(6) ! coefficients for salinity limitation (Arrigo and Sullivan, 92)
      INTEGER, DIMENSION(3) :: zindex
! FD additions
      LOGICAL ln_write_bio

      DATA zalpha /1.1e-2,3.012e-2,1.0342e-3,-4.6033e-5, 
     &            4.926e-7,-1.659e-9/
      ln_write_bio = .TRUE.
      zeps = 1.0e-14

      !==============================================================================!

      bio_do_fsw =  5.     ! FSW  dormancy threshold
      bio_do_chl = 0.1     ! chla dormancy threshold

      !==============================================================================!

      WRITE(numout,*) 
      WRITE(numout,*) ' ** ice_bio_sms : '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~~ '
      ji = 1

      DO jn = 1, ntra_bio
         IF ( flag_active(jn) ) THEN
            WRITE(numout,*) ' Tracer :   ', biotr_i_nam(jn)
            WRITE(numout,*) ' cbu_i_bio: ', ( cbu_i_bio(jn,layer), 
     &                                      layer = 1, nlay_i )
         ENDIF
      END DO

      zcase = 0

      IF ( flag_active(1) ) zcase = 1
      IF ( flag_active(4) .OR. flag_active(7) ) zcase = 2
      IF ( c_mod .NE. 'ML' ) zcase = 3
      
      WRITE(numout,*) ' flag_active(1) : ', flag_active(1)
      WRITE(numout,*) ' flag_active(4) : ', flag_active(4)
      WRITE(numout,*) ' flag_active(7) : ', flag_active(7)
      WRITE(numout,*) ' zcase          : ', zcase

      SELECT CASE ( zcase )
!
!------------------------------------------------------------------------------!
! 1) Only DSi
!------------------------------------------------------------------------------!
!
      CASE(1)

      WRITE(numout,*)
      WRITE(numout,*) ' 1 : Only DSi  '
      WRITE(numout,*) ' ~~~~~~~~~~~~~ '
      WRITE(numout,*) ' Prescribed biology: PP = ', pp_presc, 
     &                ' mmol C/m3/s '

      ! flag for solar radiation
      zflag = 1.0
      IF ( fsolgb(ji) .LE. 0.0 ) zflag = 0.0 

      DO layer = 1, nlay_bio
         zsf = 0.0
         IF ( layer .GE. 8 ) zsf = pp_presc
         zd_dsi = - zsf * ddtb * zflag
         cbu_i_bio(1,layer) = MAX(cbu_i_bio(1,layer) + zd_dsi,0.0)
      END DO ! layer
!
!------------------------------------------------------------------------------!
! 2) Interactive component
!------------------------------------------------------------------------------!
!
      CASE(2)

      WRITE(numout,*)
      WRITE(numout,*) ' 2 : DSi + NO3 + DAF '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*) ' si_c      : ', si_c
      WRITE(numout,*) ' no3_c     : ', no3_c
      WRITE(numout,*) ' mumax_bio : ', mumax_bio   
      WRITE(numout,*) ' klys_bio  : ', klys_bio
      WRITE(numout,*) ' khs_si_bio: ', khs_si_bio
      WRITE(numout,*) ' ek_bio    : ', ek_bio
      WRITE(numout,*) ' astar_alg : ', astar_alg
      WRITE(numout,*) ' rg_bio    : ', rg_bio
      WRITE(numout,*) ' lim_sal_swi :', lim_sal_swi
      WRITE(numout,*) ' lim_sal_wid :', lim_sal_wid
      WRITE(numout,*) ' lim_sal_smax:', lim_sal_smax

      ! syn_bio synthesis
      ! lys_bio lysis
      ! rem_bio remineralisation

      lim_dsi(:) = 1.; lim_lig(:) = 1.; lim_tem(:) = 1.; lim_sal(:) = 1.
      lim_no3(:) = 1.; lys_bio(:) = 0.; rem_bio(:) = 0.; syn_bio(:) = 0.

      DO layer = 1, nlay_bio
         !--------------
         ! Diatom lysis
         !--------------
         IF ( ln_lys )
     &      lys_bio(layer) = klys_bio * cbu_i_bio(4,layer)

         !------------------
         ! Remineralization
         !------------------
         IF ( ln_rem )
     &      rem_bio(layer) = kmin_bio * lys_bio(layer)

      END DO
      !-----------------------
      ! Limitation-inhibition
      !----------------------
      DO layer = 1, nlay_bio
         ! DSi limitation
         IF ( ln_lim_dsi ) ! local DSi concentration counts
     &      lim_dsi(layer) = c_i_bio(1,layer) 
     &                     / ( khs_si_bio + c_i_bio(1,layer) )
         ! NO3 limitation
         IF ( ln_lim_no3 ) ! local NO3 concentration counts
     &      lim_no3(layer) = c_i_bio(7,layer) 
     &                     / ( khs_n_bio + c_i_bio(7,layer) )
         ! light limitation
         IF ( ln_lim_lig )
!    &      lim_lig(layer) = 1. - EXP ( - par_bio(layer) / ek_bio )
         ! Platt and Jasby 76 is more accurate
     &      lim_lig(layer) = par_bio(layer) / 
     &                       SQRT( ek_bio**2 + par_bio(layer)**2 )

         ! temperature inhibition
         IF ( ln_lim_tem )
     &      lim_tem(layer) = EXP( rg_bio * ( t_i_bio(layer) - 273.15 ) )
      END DO

      ! Salinity inhibition
      IF ( ln_lim_sal ) THEN

         DO layer = 1, nlay_bio
            zbrine_sal = - ( t_i_bio(layer) - 273.15 ) / tmut

            ! 1 --- Arrigo and Sullivan (1992) ! ok
            IF ( lim_sal_swi .EQ. 1 ) THEN
               lim_sal(layer) =    zalpha(1) + zalpha(2) * zbrine_sal
     &         + zalpha(3) * zbrine_sal**2.0 + zalpha(4) * zbrine_sal**3.0
     &         + zalpha(5) * zbrine_sal**4.0 + zalpha(6) * zbrine_sal**5.0
            ENDIF

            ! 2 --- 4th order gaussian
            IF ( lim_sal_swi .EQ. 2 ) THEN
               zdummmy_arg = ( ( zbrine_sal - lim_sal_smax ) 
     &                         / lim_sal_wid )**4 
               lim_sal(layer) = EXP (- zdummy_arg )
            ENDIF

            ! 3 --- Sine
            IF ( lim_sal_swi .EQ. 3 ) THEN
               za_sin = 1.
               zpi = ACOS(-1.)
               zc_sin = LOG(2.) / ( LOG(100.) - LOG( lim_sal_smax) )
               zb_sin = zpi / ( 2. * lim_sal_smax**zc_sin )
               lim_sal(layer) = za_sin * 
     &                          ASIN ( zb_sin * zbrine_sal**c_sin )
            ENDIF

            ! 4 --- Gaussian ( Martouf and Grozny, 2010 )
            IF ( lim_sal_swi .EQ. 4 ) THEN
               WRITE(numout,*) ' 4th switch for salinity limitation '
               WRITE(numout,*) ' layer no     : ', layer
               WRITE(numout,*) ' lim_sal_smax : ', lim_sal_smax
               WRITE(numout,*) ' lim_sal_wid  : ', lim_sal_smax
               zdummy_arg = ( ( zbrine_sal - lim_sal_smax ) / 
     &                         lim_sal_wid )**2.
               WRITE(numout,*) ' zdummy_arg   : ', zdummy_arg
               WRITE(numout,*) ' zbrine_sal   : ', zbrine_sal
               lim_sal(layer) = EXP( - zdummy_arg )
               WRITE(numout,*) ' lim_sal      : ', lim_sal(layer)
            ENDIF

            ! 5 --- Gaussian with cubic argument
            IF ( lim_sal_swi .EQ. 5 ) THEN
               zUU = 1.9
               zTT = 100.
               zSS = lim_sal_smax
               za0 = -zUU
               zVV = ( zTT**3 - 2*zSS**3 ) /  ( 2.*zSS - zTT )
               zd0 = zUU / ( zSS * ( zVV + zSS**2 ) )
               zb0 = zd0 * zVV
               zdummy_arg = za0 + zb0 * zbrine_sal + zd0 * zbrine_sal**3
               lim_sal(layer) = EXP(-zdummy_arg*2)
            ENDIF

            ! Arrigo's regression is sometimes over 1 (near 34 psu). Hence we add a max min
            lim_sal(layer) = MAX( 0., MIN( lim_sal(layer), 1. ) )
            WRITE(numout,*) ' lim_sal      : ', lim_sal(layer)
         END DO

      ENDIF

      !------------------
      ! Diatom synthesis
      !------------------
      DO layer = 1, nlay_bio
         z_syn = mumax_bio * lim_tem(layer) * lim_sal(layer) * 
     &           MIN ( lim_dsi(layer), lim_no3(layer), lim_lig(layer) )

         IF ( ln_syn ) THEN
            zsyn1 = z_syn * cbu_i_bio(4,layer)                ! normal rate
            zsyn2 = zsyn1
            zsyn3 = zsyn1
            IF ( ln_lim_dsi ) ! local DSi concentration counts
     &         zsyn2 = cbu_i_bio(1,layer) / ddtb / si_c - zeps   ! if DSi depleted
            IF ( ln_lim_no3 ) ! local NO3 concentration counts
     &         zsyn3 = cbu_i_bio(7,layer) / ddtb / no3_c - zeps  ! if NO3 depleted
            syn_bio(layer) = MIN( zsyn1, zsyn2, zsyn3 ) 

            WRITE(numout,*) ' zsyn1 : ', zsyn1
            WRITE(numout,*) ' zsyn2 : ', zsyn2
            WRITE(numout,*) ' zsyn3 : ', zsyn3
         ENDIF

         ! Dormancy
         !----------
         ! see Lancelot et al., 2009) ! a priori it's a AND, maybe OR ?
         IF ( ln_dormant ) THEN
           IF ( ( fsolgb(ji)        .LT. bio_do_fsw ) .AND. 
     &          ( chla_i_bio(layer) .LT. bio_do_chl ) ) THEN
               syn_bio(layer) = 0.0
               rem_bio(layer) = 0.0
               lys_bio(layer) = 0.0
            ENDIF
         ENDIF

         IF ( ln_write_bio ) THEN
            WRITE(numout,*) ' Initial values. layer = ', layer
            WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
            WRITE(numout,*) ' DSi : ', cbu_i_bio(1,layer)
            WRITE(numout,*) ' DAF : ', cbu_i_bio(4,layer)
            WRITE(numout,*) ' NO3 : ', cbu_i_bio(7,layer)
            WRITE(numout,*)
            WRITE(numout,*) ' Terms '
            WRITE(numout,*) ' ~~~~~ '
            WRITE(numout,*) ' syn_bio   : ', syn_bio(layer)
            WRITE(numout,*) ' lys_bio   : ', lys_bio(layer)
            WRITE(numout,*) ' rem_bio   : ', rem_bio(layer)
            WRITE(numout,*) ' lim_dsi   : ', lim_dsi(layer)
            WRITE(numout,*) ' lim_no3   : ', lim_no3(layer)
            WRITE(numout,*) ' lim_lig   : ', lim_lig(layer)
            WRITE(numout,*) ' lim_tem   : ', lim_tem(layer)
            WRITE(numout,*) ' lim_sal   : ', lim_sal(layer)
            WRITE(numout,*)
            WRITE(numout,*) ' t_i_bio   : ', t_i_bio(layer)
            WRITE(numout,*) ' par_bio   : ', par_bio(layer)
            WRITE(numout,*)
         ENDIF

         ! Update variables
         !------------------
         zd_dsi1 = si_c * ( rem_bio(layer) - syn_bio(layer) ) * ddtb
         zd_dsi2 = - cbu_i_bio(1,layer) 
         zd_dsi = MAX ( zd_dsi1, zd_dsi2 )

         zd_no31 = no3_c * ( rem_bio(layer) - syn_bio(layer) ) * ddtb
         zd_no32 = - cbu_i_bio(7,layer) 
         zd_no3 = MAX ( zd_no31, zd_no32 )

         zd_daf1 = ( syn_bio(layer) - lys_bio(layer) ) * ddtb
         zd_daf2 = - cbu_i_bio(4,layer)
         zd_daf = MAX ( zd_daf1, zd_daf2 )

         cbu_i_bio(1,layer) = cbu_i_bio(1,layer) + zd_dsi
         cbu_i_bio(7,layer) = cbu_i_bio(7,layer) + zd_no3
         cbu_i_bio(4,layer) = cbu_i_bio(4,layer) + zd_daf

         IF ( ln_write_bio ) THEN
            WRITE(numout,*) ' tendency terms '
            WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
            WRITE(numout,*) ' zd_dsi : ', zd_dsi
            WRITE(numout,*) ' zd_no3 : ', zd_no3
            WRITE(numout,*) ' zd_daf : ', zd_daf
            WRITE(numout,*)
            WRITE(numout,*) ' Final values '
            WRITE(numout,*) ' ~~~~~~~~~~~~ '
            WRITE(numout,*) ' DSi : ', cbu_i_bio(1,layer)
            WRITE(numout,*) ' NO3 : ', cbu_i_bio(7,layer)
            WRITE(numout,*) ' DAF : ', cbu_i_bio(4,layer)
            WRITE(numout,*)
         ENDIF

      END DO ! layer
!
!------------------------------------------------------------------------------!
! 3) Skeletal layer or biologically-active layer
!------------------------------------------------------------------------------!
!
      CASE(3)

      WRITE(numout,*)
      WRITE(numout,*) ' 3 : One layer for a model of type ', c_mod
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
      WRITE(numout,*) ' si_c      : ', si_c
      WRITE(numout,*) ' no3_c     : ', no3_c
      WRITE(numout,*) ' mumax_bio : ', mumax_bio   
      WRITE(numout,*) ' klys_bio  : ', klys_bio
      WRITE(numout,*) ' khs_si_bio: ', khs_si_bio
      WRITE(numout,*) ' khs_n_bio : ', khs_n_bio
      WRITE(numout,*) ' ek_bio    : ', ek_bio
      WRITE(numout,*) ' astar_alg : ', astar_alg
      WRITE(numout,*) ' rg_bio    : ', rg_bio
      WRITE(numout,*) 

      IF ( c_mod .EQ. 'BA' ) CALL ice_bio_grid(kideb,kiut,.TRUE.) ! compute the biological grid for BAL case

!     CALL ice_bio_interp_phy2bio(kideb,kiut,nlay_i,.FALSE.) 
!     WRITE(numout,*) ' par       : ', ( par    (layer), 
!    &                layer = 1, nlay_i  )
!     WRITE(numout,*) ' par_bio   : ', ( par_bio(layer), 
!    &                layer = 1, nlay_bio  )

      lim_dsi(:) = 1.; lim_lig(:) = 1.; lim_tem(:) = 1.; lim_sal(:) = 1.
      lim_no3(:) = 1.; lys_bio(:) = 0.; rem_bio(:) = 0.; syn_bio(:) = 0.

      layer = nlay_bio ! only one layer

      !--------------
      ! Flux
      !--------------
      zindex(1) = 1
      zindex(2) = 4
      zindex(3) = 7
  
      DO i = 1, 3
         ii = zindex(i)
         ! diffusivity
         IF ( ii .NE. 4 ) THEN 
            zdiff = diff_nut ! for nutrients
         ELSE
            zdiff = diff_da  ! for diatoms
         ENDIF

         ! diffusive flux
         fcb(ii)  = 2. * zdiff * ( c_w_bio(ii) - cbu_i_bio(ii,layer) ) ! diffusive flux
     &        / h_bio 

         ! growth-melt flux
         IF ( i_flux .NE. 4 )  THEN  ! growth-melt
           fcbp(ii) = MAX( dh_i_bott(ji), 0. ) * c_w_bio(ii) / ddtb    ! growth
     &      + MIN( dh_i_bott(ji), 0. ) * cbu_i_bio(ii,layer) / ddtb    ! melt
         ELSE ! growth only
           fcbp(ii) = MAX( dh_i_bott(ji), 0. ) * c_w_bio(ii) / ddtb    ! growth
         ENDIF

         ! snow ice and max fluxes are irrelevant in SL / BAL case
         fcsi(ii)    = 0.
         fcb_max(ii) = 0.

      END DO

!     fcb(4)  = 2. * diff_da * ( c_w_bio(4) - cbu_i_bio(4,layer) )
!    &        / h_bio                                                ! diffusion
!     fcbp(4) = MAX( dh_i_bott(ji), 0. ) * c_w_bio(4) / ddtb         ! growth
!    &        + MIN( dh_i_bott(ji), 0. ) * cbu_i_bio(4,layer) / ddtb ! melt
!     fcsi(4)    = 0.
!     fcb_max(4) = 0.

!     fcb(1)  = 2. * diff_nut * ( c_w_bio(1) - cbu_i_bio(1,layer) ) 
!    &        / h_bio                                                ! diffusion (mmol/m2/s)
!     fcbp(1) = MAX( dh_i_bott(ji), 0. ) * c_w_bio(1) / ddtb         ! growth
!    &        + MIN( dh_i_bott(ji), 0. ) * cbu_i_bio(1,layer) / ddtb ! melt
!     fcsi(1)    = 0.
!     fcb_max(1) = 0.

!     fcb(7)  = 2. * diff_nut * ( c_w_bio(7) - cbu_i_bio(7,layer) ) 
!    &        / h_bio                                                ! diffusion
!     fcbp(7) = MAX( dh_i_bott(ji), 0. ) * c_w_bio(7) / ddtb         ! growth
!    &        + MIN( dh_i_bott(ji), 0. ) * cbu_i_bio(7,layer) / ddtb ! melt
!     fcsi(7)    = 0.
!     fcb_max(7) = 0.

      IF ( i_flux .EQ. 1 ) THEN ! growth-melt + diffusion
         zf_da  = fcb(4) + fcbp(4)
         zf_dsi = fcb(1) + fcbp(1)
         zf_no3 = fcb(7) + fcbp(7)
          WRITE(numout,*) ' fcb(1) : ', fcb(1), ' fcbp(1) : ', fcbp(1)
          WRITE(numout,*) ' fcb(7) : ', fcb(7), ' fcbp(7) : ', fcbp(7)
          WRITE(numout,*) ' c_w_bio(1) : ', c_w_bio(1)
          WRITE(numout,*) ' cbu_i_bio(1,layer) :', cbu_i_bio(1,layer)
          WRITE(numout,*) ' layer : ', layer
      ENDIF

      IF ( ( i_flux .EQ. 2 ) .OR. ( i_flux .EQ. 4 ) ) THEN ! growth-melt 
         zf_da  = fcbp(4)
         zf_dsi = fcbp(1)
         zf_no3 = fcbp(7)
         fcb(:) = 0.
      ENDIF

      IF ( i_flux .EQ. 3 ) THEN ! diffusion
         zf_da  = fcb(4)
         zf_dsi = fcb(1)
         zf_no3 = fcb(7)
         fcbp(:) = 0.
      ENDIF
      
      !--------------
      ! Diatom lysis
      !--------------
      IF ( ln_lys )
     &   lys_bio(layer) = klys_bio * cbu_i_bio(4,layer)

      !------------------
      ! Remineralization
      !------------------
      IF ( ln_rem )
     &   rem_bio(layer) = kmin_bio * lys_bio(layer)

      !-----------------------
      ! Limitation-inhibition
      !----------------------
      ! DSi limitation
      IF ( ln_lim_dsi ) ! local DSi concentration counts
     &   lim_dsi(layer) = cbu_i_bio(1,layer) 
     &                  / ( khs_si_bio + cbu_i_bio(1,layer) )
      ! NO3 limitation
      IF ( ln_lim_no3 ) ! local NO3 concentration counts
     &   lim_no3(layer) = cbu_i_bio(7,layer) 
     &                  / ( khs_n_bio + cbu_i_bio(7,layer) )
      ! light limitation
      IF ( ln_lim_lig )
     &   lim_lig(layer) = par(layer) / 
     &                    SQRT( ek_bio**2 + par(layer)**2 )
      WRITE(numout,*) ' layer      : ', layer

      ! temperature inhibition
      IF ( ln_lim_tem )
     &   lim_tem(layer) = EXP( rg_bio * ( t_i_bio(layer) - 273.15 ) )

      !------------------
      ! Diatom synthesis
      !------------------
      z_syn = mumax_bio * lim_tem(layer) * 
     &        MIN ( lim_dsi(layer), lim_no3(layer), lim_lig(layer) )

      IF ( ln_syn ) THEN
         zsyn1 = z_syn * cbu_i_bio(4,layer)                ! normal rate
         zsyn2 = zsyn1
         zsyn3 = zsyn1
         IF ( ln_lim_dsi ) ! local DSi concentration counts
     &      zsyn2 = cbu_i_bio(1,layer) / ddtb / si_c - zeps   ! if DSi depleted
         IF ( ln_lim_no3 ) ! local NO3 concentration counts
     &      zsyn3 = cbu_i_bio(7,layer) / ddtb / no3_c - zeps  ! if NO3 depleted
         syn_bio(layer) = MIN( zsyn1, zsyn2, zsyn3 ) 

         WRITE(numout,*) ' zsyn1 : ', zsyn1
         WRITE(numout,*) ' zsyn2 : ', zsyn2
         WRITE(numout,*) ' zsyn3 : ', zsyn3
      ENDIF

      IF ( ln_write_bio ) THEN
         WRITE(numout,*) ' Initial values. layer = ', layer
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) ' DSi : ', cbu_i_bio(1,layer)
         WRITE(numout,*) ' DAF : ', cbu_i_bio(4,layer)
         WRITE(numout,*) ' NO3 : ', cbu_i_bio(7,layer)
         WRITE(numout,*)
         WRITE(numout,*) ' Terms '
         WRITE(numout,*) ' ~~~~~ '
         WRITE(numout,*) ' syn_bio   : ', syn_bio(layer)
         WRITE(numout,*) ' lys_bio   : ', lys_bio(layer)
         WRITE(numout,*) ' rem_bio   : ', rem_bio(layer)
         WRITE(numout,*) ' lim_dsi   : ', lim_dsi(layer)
         WRITE(numout,*) ' lim_no3   : ', lim_no3(layer)
         WRITE(numout,*) ' lim_lig   : ', lim_lig(layer)
         WRITE(numout,*) ' lim_tem   : ', lim_tem(layer)
         WRITE(numout,*) ' lim_sal   : ', lim_sal(layer)
         WRITE(numout,*)
         WRITE(numout,*) ' t_i_bio   : ', t_i_bio(layer)
         WRITE(numout,*) ' par       : ', par(layer)
         WRITE(numout,*)
      ENDIF

      ! Update variables
      !------------------
      zd_dsi1 = ( si_c * ( rem_bio(layer) - syn_bio(layer) ) 
     &        + zf_dsi / h_bio ) * ddtb
      zd_dsi2 = - cbu_i_bio(1,layer) 
      zd_dsi = MAX ( zd_dsi1, zd_dsi2 )
      WRITE(numout,*) ' zf_dsi  : ', zf_dsi
      WRITE(numout,*) ' h_bio   : ', h_bio
      WRITE(numout,*) ' flux term : ', zf_dsi / h_bio * ddtb
      WRITE(numout,*) ' zd_dsi1 : ', zd_dsi1
      WRITE(numout,*) ' zd_dsi2 : ', zd_dsi2

      zd_no31 = ( no3_c * ( rem_bio(layer) - syn_bio(layer) )
     &        + zf_no3 / h_bio ) * ddtb
      zd_no32 = - cbu_i_bio(7,layer) 
      zd_no3 = MAX ( zd_no31, zd_no32 )

      zd_daf1 = ( syn_bio(layer) - lys_bio(layer) 
     &        + zf_da / h_bio ) * ddtb
      zd_daf2 = - cbu_i_bio(4,layer)
      zd_daf = MAX ( zd_daf1, zd_daf2 )

      cbu_i_bio(1,layer) = cbu_i_bio(1,layer) + zd_dsi
      cbu_i_bio(7,layer) = cbu_i_bio(7,layer) + zd_no3 
      cbu_i_bio(4,layer) = cbu_i_bio(4,layer) + zd_daf 

      IF ( ln_write_bio ) THEN
         WRITE(numout,*) ' tendency terms '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
         WRITE(numout,*) ' zd_dsi : ', zd_dsi
         WRITE(numout,*) ' zd_no3 : ', zd_no3
         WRITE(numout,*) ' zd_daf : ', zd_daf
         WRITE(numout,*) ' zf_dsi : ', zf_dsi
         WRITE(numout,*) ' zf_no3 : ', zf_no3
         WRITE(numout,*) ' zf_da  : ', zf_da 
         WRITE(numout,*)
         WRITE(numout,*) ' Final values '
         WRITE(numout,*) ' ~~~~~~~~~~~~ '
         WRITE(numout,*) ' DSi : ', cbu_i_bio(1,layer)
         WRITE(numout,*) ' NO3 : ', cbu_i_bio(7,layer)
         WRITE(numout,*) ' DAF : ', cbu_i_bio(4,layer)
         WRITE(numout,*)
      ENDIF
!
!------------------------------------------------------------------------------!
! 4) End of the routine
!------------------------------------------------------------------------------!
!
      CASE DEFAULT
         WRITE(numout,*) ' BIG PROBLEM '
      END SELECT

      DO jn = 1, ntra_bio
         IF ( flag_active(jn) ) THEN
            WRITE(numout,*) ' Tracer :   ', biotr_i_nam(jn)
            WRITE(numout,*) ' cbu_i_bio: ', ( cbu_i_bio(jn,layer), 
     &                                      layer = 1, nlay_i )
         ENDIF
      END DO
      !+++++

      WRITE(numout,*)
      WRITE(numout,*) '    *** After biological sources and sinks *** '
      WRITE(numout,*) '    model output '
      DO jn = 1, ntra_bio
         IF ( flag_active(jn) ) THEN
           WRITE(numout,*) ' biotr_i_nam : ', biotr_i_nam(jn)
           WRITE(numout,*) ' cbu_i_bio : ', ( cbu_i_bio(jn, jk), jk = 1,
     &                                      nlay_bio )
         ENDIF
      END DO
      WRITE(numout,*)

      WRITE(numout,*)
      WRITE(numout,*) ' End of ice_bio_sms '
      WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
     

!==============================================================================|
! Fin de la routine ice_bio_sms
      WRITE(numout,*)

      END
