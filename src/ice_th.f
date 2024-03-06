      SUBROUTINE ice_th(ntrmax,nlay_i,nlay_s,numofday)

        !!------------------------------------------------------------------
        !!                ***         ROUTINE ice_th       ***
        !! ** Purpose :
        !!           Ice thermodynamics 
        !! ** Method  :
        !!    This routine calls the thermodynamic and biological routines
        !!
        !! ** Arguments :
        !!           ntrmax, nlay_i, nlay_s, numofday
        !!
        !! ** Inputs / Ouputs : (global commons)
        !!
        !! ** External : 
        !!
        !! ** References : Vancoppenolle et al., JGR 2007
        !!
        !! ** History :
        !!       (1) CLIO, Goosse and Fichefet, JGR, 1999.
        !!       (2) LIM-1D, Vancoppenolle et al., JGR, 2007.
        !!       (3) BIO-LIM, Martin Vancoppenolle, 2008
        !!
        !!------------------------------------------------------------------
        !! * Arguments

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'
      INCLUDE 'bio.com'

      REAL(8) :: ain(imax,jmax),zinda(imax,jmax),ifvt(imax,jmax)
      REAL(8) :: qdtcn(imax,jmax),qlbsbq(imax,jmax),zfwat(imax,jmax)
      REAL(8) :: zhgbqp(imax,jmax)
      REAL(8) :: fdtcn(imax,jmax)
      COMMON/lowat/ zfwat
      INTEGER :: ntrmax

      ! Energy conservation
      LOGICAL con_i

      zeps0 = 1.0e-16
      zeps1 = 1.0e-20
      zeps2 = 1.0e-04

      con_i = .true. ! conservation check in the ice or not
      jl = 1         ! category number (temporary)
! FD debug
      con_i = .false. ! conservation check in the ice or not

! FD      WRITE(numout,*) ' * ice_th : '
! FD      WRITE(numout,*) ' ~~~~~~~~~~ '
! FD      WRITE(numout,*) 
! FD      WRITE(numout,*) ' nlay_i : ', nlay_i
! FD      WRITE(numout,*) ' nlay_s : ', nlay_s
!
!-------------------------------------------------------------------------------
!  1. Convert vectors
!-------------------------------------------------------------------------------
!
      ! This step is a residual from the CLIO 3D code that 
      ! should be removed in principle
      nbpb = 1
      npb(nbpb) = 1
 
      CALL gather(nbpb,hnpbqb,hnpbq,npb)   ! snowfall
      CALL gather(nbpb,fsolgb,fsolg,npb)   ! net solar flux
      CALL gather(nbpb,fbbqb,fbbq,npb)     ! oceanic heat flux
      CALL gather(nbpb,ratbqb,ratbqg,npb)  ! downwelling longwave flux
      CALL gather(nbpb,psbqb,psbq,npb)     ! atm pressure
      CALL gather(nbpb,tabqb,tabq,npb)     ! air temperature
      CALL gather(nbpb,qabqb,qabq,npb)     ! air humidity
      CALL gather(nbpb,vabqb,vabq,npb)     ! wind velocity
      CALL gather(nbpb,fscbqb,fsbbq,npb)  ! sensible heat flux
      CALL gather(nbpb,fltbqb,ffltbq,npb)  ! latent heat flux

      CALL gather(nbpb,tdewb,tdew,npb)     ! relative air humidity (output)
      CALL gather(nbpb,albgb,albg,npb)     ! albedo (output)

      CALL gather(nbpb,t_su_b,t_su,npb)    ! surface temperature
      CALL gather(nbpb,t_bo_b,t_bo,npb)    ! bottom temperature
      CALL gather(nbpb,ht_i_b,ht_i,npb)    ! ice thickness
      CALL gather(nbpb,ht_s_b,ht_s,npb)    ! snow depth

      DO 259 k = 1, maxnlay
         CALL gather(nbpb,t_i_b(1,k),t_i(1,1,k),npb) ! ice temperature
         CALL gather(nbpb,t_s_b(1,k),t_s(1,1,k),npb) ! snow temperature
         CALL gather(nbpb,s_i_b(1,k),s_i(1,1,k),npb) ! ice salinity
 259  CONTINUE
!
!-------------------------------------------------------------------------------
!  2) Prevent high values for temperatures, ice / snow enthalpy
!-------------------------------------------------------------------------------
!
      DO 262 ji = 1, nbpb

         DO layer = 1, nlay_s
            t_s_b(ji,layer)  =  MIN( tpw , t_s_b(ji,layer) )
         END DO
         DO layer = 1, nlay_i
            tmelts           =  - tmut*s_i_b(ji,layer) + tpw
            t_i_b(ji,layer)  =  MIN(tmelts,t_i_b(ji,layer))
         END DO

         CALL ice_th_enmelt(1,nbpb, nlay_s, nlay_i) ! ice enthalpy


 262  CONTINUE

      ! Initialize total heat content
      IF ( con_i ) CALL ice_th_glohec( qt_i_in , qt_s_in ,              
     &                                 q_i_layer_in , 1 , nbpb , jl, 
     &                                 nlay_s , nlay_i )
!
!-------------------------------------------------------------------------------
!  3) Call of the thermodynamic subroutines                                    
!-------------------------------------------------------------------------------
!
! FD      IF ( numit .EQ. nstart ) CALL ice_bio_ini( 1 , nbpb , nlay_i )    ! bio initialization

      CALL ice_rad( nlay_s , nlay_i , 1 , nbpb , numofday )        ! radiative transfer

      CALL ice_th_diff( nlay_s , nlay_i , 1 , nbpb , numofday )    ! heat diffusion

! FD      IF ( con_i ) THEN                                            ! conservation test
! FD         CALL ice_th_glohec( qt_i_fin , qt_s_fin , q_i_layer_fin ,
! FD     &                       1 , nbpb , jl , nlay_s , nlay_i )
! FD         CALL ice_th_con_dif( 1 , nbpb , nlay_s , nlay_i , jl )
! FD      ENDIF
! FD
! FD      IF ( gravdr .EQ. 'CW' )
! FD     &   CALL ice_sal_diff_CW(nlay_i,1,nbpb)                       ! salt transport (Cox and Weeks)
! FD
! FD      IF ( gravdr .EQ. 'RA' )
! FD     &   CALL ice_sal_diff(nlay_i,1,nbpb)                          ! salt transport (Rayleigh-number based)
! FD      sn_i_b(:) = s_i_b(1,:)  ! new salinity in absence of salinity model
! FD
! FD      ! the new salinities and brine volume should be updated in phy_remap
! FD
! FD      IF ( ln_trdiff .AND. ( c_mod .EQ. 'ML' ) ) 
! FD     &   CALL ice_bio_diff( 1 , nbpb , nlay_i )                    ! bio transport
! FD
! FD      CALL ice_bio_sms(nlay_i,1,nbpb)                              ! bio sources minus sinks
! FD
! FD      CALL ice_th_dh(nlay_s,nlay_i,1,nbpb)                         ! ice growth and melt
! FD
! FD      CALL ice_phy_remap(nlay_s,nlay_i,1,nbpb)                     ! physical remapping of heat
                                                                   ! content and salinity
      IF ( con_i ) THEN                                            ! conservation test
         CALL ice_th_glohec( qt_i_fin , qt_s_fin , q_i_layer_fin ,
     &                       1 , nbpb , jl, nlay_s, nlay_i )
         CALL ice_th_con_dh(1,nbpb,nlay_s,nlay_i,jl)
      ENDIF

! FD      IF ( ln_trremp .AND. ( c_mod .EQ. 'ML' ) )  
! FD     &   CALL ice_bio_remap(nlay_s, nlay_i, 1, nbpb)               ! bio remapping
! FD
! FD      !---------------
! FD      ! Chlorophyll a
! FD      !---------------
! FD      ! Units, mg m-3 (micro g.l-1)
! FD      DO layer = 1, nlay_bio
! FD         chla_i_bio(layer) = cbu_i_bio(4,layer) * chla_c
! FD      END DO
! FD
! FD      WRITE(numout,*)
! FD      WRITE(numout,*) '    *** After tracer remapping *** '
! FD      WRITE(numout,*) '    model output '
! FD
! FD      DO jn = 1, ntra_bio
! FD         IF ( flag_active(jn) ) THEN
! FD           WRITE(numout,*) ' biotr_i_nam : ', biotr_i_nam(jn)
! FD           WRITE(numout,*) ' cbu_i_bio : ', ( cbu_i_bio(jn, jk), jk = 1,
! FD     &                                        nlay_bio )
! FD         ENDIF
! FD      END DO
! FD      WRITE(numout,*) ' chla_i_bio : ', ( chla_i_bio(layer), layer = 1, 
! FD     &                                    nlay_bio )
! FD      WRITE(numout,*)
!
!-------------------------------------------------------------------------------
!  4) Outputs 
!-------------------------------------------------------------------------------
!
      !---------------
      ! Netcdf Output
      !---------------
      CALL ice_output(nlay_i,nlay_s)
!
!-------------------------------------------------------------------------------
!  5. Convert vectors
!-------------------------------------------------------------------------------
!
      ! This step is a residual from the CLIO 3D code that 
      ! should be removed in principle
      CALL scater(nbpb,firg,npb,fratsb)
      CALL scater(nbpb,fcsg,npb,fcsb)
      CALL scater(nbpb,fleg,npb,fleb)

      DO 273 k = 1, maxnlay
         CALL scater(nbpb,t_i(1,1,k),npb,t_i_b(1,k),npb)
         CALL scater(nbpb,t_s(1,1,k),npb,t_s_b(1,k),npb)
         CALL scater(nbpb,s_i(1,1,k),npb,s_i_b(1,k),npb)
 273  CONTINUE

      CALL scater(nbpb,t_su,npb,t_su_b)
      CALL scater(nbpb,t_bo,npb,t_bo_b)
      CALL scater(nbpb,ht_s,npb,ht_s_b)
      CALL scater(nbpb,ht_i,npb,ht_i_b)

      RETURN
!
!------------------------------------------------------------------------------
!- end of ice_th
      END
