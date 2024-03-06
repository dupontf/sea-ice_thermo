      SUBROUTINE ice_output(nlay_i,nlay_s)

!=============================================================================!
! ice_output : Creates and write in the ice output file
! Vancoppenolle and Fettweise, UCL-ASTR, June 2007
!=============================================================================!

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'
      INCLUDE 'bio.com'

      CHARACTER(len=10) :: 
     &   filenc='ice.out.nc'
      CHARACTER(len=1) ::
     &   c1
      CHARACTER(len=4) ::
     &   c4

      REAL*4, DIMENSION(nend)     :: dimval
      REAL*4, DIMENSION(nlay_i)   :: dimval2
      REAL*4, DIMENSION(nlay_bio) :: dimval3
      REAL*4, DIMENSION(nlay_s)   :: dimval4
      REAL*4, DIMENSION(nlay_i)   :: dummy_vector

!=============================================================================!

!
!-----------------------------------------------------------------------------!
! 1) File definition
!-----------------------------------------------------------------------------!
!
      IF ( numit .EQ. nstart ) THEN

         WRITE(numout,*) ' ice_output : '
         WRITE(numout,*) ' ~~~~~~~~~~ : '

         WRITE(numout,*) ' Initialization of the NetCdf file : '
         !---------------------
         ! File initialization
         !---------------------
         CALL CF_INI_FILE(filenc,"Output NETCDF for S(L)IMCO")

         !---------------------
         ! Dimensions creation
         !---------------------
         ! time 
         !------
         CALL CF_CREATE_DIM("time", "days since 0000-01-01 00:00:0.0", 
     &      0, dimval)
         WRITE(numout,*) ' time dimension created '

         ! snow vertical coordinate (moving)
         !-----------------------------------
         DO i = 1, nlay_s
            dimval4(i) = REAL(i)
         ENDDO
         CALL CF_CREATE_DIM("z_s_p"   ,"m"       ,nlay_s,
     &   dimval4(1:nlay_s))
         WRITE(numout,*) ' z_s_p dimension created ',nlay_s

         ! physical vertical coordinate (moving)
         !---------------------------------------
         DO i = 1, nlay_i
            dimval2(i) = REAL(i)
         ENDDO
         CALL CF_CREATE_DIM("z_i_p"   ,"m"       ,nlay_i,
     &   dimval2(1:nlay_i))
         WRITE(numout,*) ' z_i_p dimension created ',nlay_i

!        ! biological vertical coordinate (moving)
!        !-----------------------------------------
! FD         DO i = 1, nlay_bio 
! FD            dimval3(i) = REAL(i)
! FD         ENDDO
! FD         CALL CF_CREATE_DIM("z_i_b"   ,"m"       ,nlay_bio,
! FD     &   dimval3(1:nlay_bio))
! FD         WRITE(numout,*) ' z_i_b dimension created ',nlay_bio
! FD         CALL FLUSH(numout)

         !--------------------
         ! Variables creation
         !--------------------
         ! Times
         !-------
         CALL CF_CREATE_VAR("ts_m","Time step in months",
     .                   "mon","time","-","-","-" )

         CALL CF_CREATE_VAR("ts_y","Time step in years",
     .                   "y","time","-","-","-" )

         ! Grids
         !-------
         CALL CF_CREATE_VAR("z_ip","Vertical cote, physics",
     .                   "m","time","z_i_p","-","-" )

! FD         CALL CF_CREATE_VAR("z_ib","Vertical cote, biology",
! FD     .                   "m","time","z_i_b","-","-" )

         ! Ice / snow thickness and variations
         !-------------------------------------
         CALL CF_CREATE_VAR("h_i","Ice thickness",
     .                   "m","time","-","-","-" )
         WRITE(numout,*) ' h_i dimension created '

         CALL CF_CREATE_VAR("dhib","Bottom accretion",
     .                   "m/s","time","-","-","-" )

         CALL CF_CREATE_VAR("dhisu","Surface ice melt",
     .                   "m/s","time","-","-","-" )

         CALL CF_CREATE_VAR("dhisi","Snow ice formation",
     .                   "m/s","time","-","-","-" )

         CALL CF_CREATE_VAR("dhs","Snow variations",
     .                   "m/s","time","-","-","-" )

         CALL CF_CREATE_VAR("h_s","Snow depth   ",
     .                   "m","time","-","-","-" )
         WRITE(numout,*) ' h_s dimension created '

         ! Inner ice variables  
         !---------------------
         CALL CF_CREATE_VAR("t_i","Ice Temperature",
     .                   "C","time","z_i_p","-","-" )
         WRITE(numout,*) ' t_i dimension created '

         CALL CF_CREATE_VAR("s_i","Ice Salinity",
     .                   "ppt","time","z_i_p","-","-" )
         WRITE(numout,*) ' s_i dimension created '

         CALL CF_CREATE_VAR("e_i","Brine volume",
     .                   "%","time","z_i_p","-","-" )
         WRITE(numout,*) ' e_i dimension created '

         CALL CF_CREATE_VAR("ra","Rayleigh number",
     .                   "-","time","z_i_p","-","-" )
         WRITE(numout,*) ' Ra dimension created '

         CALL CF_CREATE_VAR("ub","brine flushing velocity",
     .                   "m/s","time","z_i_p","-","-" )
         WRITE(numout,*) ' Ub dimension created '

         CALL CF_CREATE_VAR("diff","Salt diffusivity", 
     .                   "m2/s","time","z_i_p","-","-" )
         WRITE(numout,*) ' Diff dimension created '

         CALL CF_CREATE_VAR("t_s","Snow Temperature",
     .                   "C","time","z_s_p","-","-" )
         WRITE(numout,*) ' t_s dimension created '

         CALL CF_CREATE_VAR("t_su","Surface Temperature",
     .                   "C","time","-","-","-" )
         WRITE(numout,*) ' t_su dimension created '

         CALL CF_CREATE_VAR("fsb","Ice-ocean bd salt flux",
     .                   "kgNaCl/m2/s","time","-","-","-" )
         WRITE(numout,*) ' fsb  dimension created '

         CALL CF_CREATE_VAR("fsbp","Ice-ocean se salt flux",
     .                   "kgNaCl/m2/s","time","-","-","-" )
         WRITE(numout,*) ' fsbp dimension created '

         ! Forcing
         !---------
         CALL CF_CREATE_VAR("albe","Surface albedo",
     .                   "/","time","-","-","-" )
         CALL CF_CREATE_VAR("F_sw","Incoming solar",
     .                   "W.m-2","time","-","-","-" )
         CALL CF_CREATE_VAR("Fswn","Net solar",
     .                   "W.m-2","time","-","-","-" )
         CALL CF_CREATE_VAR("Flwd","Incoming longwave",
     .                   "W.m-2","time","-","-","-" )
         CALL CF_CREATE_VAR("Flwu","Emitted longwave",
     .                   "W.m-2","time","-","-","-" )
         CALL CF_CREATE_VAR("F_sh","Sensible heat flux",
     .                   "W.m-2","time","-","-","-" )
         CALL CF_CREATE_VAR("F_lh","Latent heat flux",
     .                   "W.m-2","time","-","-","-" )
         CALL CF_CREATE_VAR("Tair","Air temperature ",
     .                   "C","time","-","-","-" )
         CALL CF_CREATE_VAR("qair","Air humidity    ",
     .                   "kg/kg","time","-","-","-" )
         CALL CF_CREATE_VAR("wspd","Wind speed",
     .                   "m/s","time","-","-","-" )
         CALL CF_CREATE_VAR("qsfc","Surface humidity", 
     .                   "kg/kg","time","-","-","-" )
         CALL CF_CREATE_VAR("tdew","Relative air humidity", 
     .                   "%","time","-","-","-" )

         ! Radiative transfer 
         !--------------------
         CALL CF_CREATE_VAR("T_is","Transmitted rad flx through ice",
     .                   "W.m-2","time","-","-","-" )

         CALL CF_CREATE_VAR("T_oc","Transmitted rad flux to ocean",
     .                   "W.m-2","time","-","-","-" )

         CALL CF_CREATE_VAR("A_s","Absorbed rad flux by snow",
     .                   "W.m-2","time","z_s_p","-","-" )

         CALL CF_CREATE_VAR("A_ib","Absorbed rad in ice by algae",
     .                   "W.m-2","time","z_i_p","-","-" )

         CALL CF_CREATE_VAR("A_ip","Absorbed rad in ice physically",
     .                   "W.m-2","time","z_i_p","-","-" )

! FD         ! Biological variables
! FD         !----------------------
! FD         DO jn = 1, ntra_bio
! FD
! FD         IF ( flag_diff(jn) .AND. flag_active(jn) ) THEN
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         CALL CF_CREATE_VAR(biotr_i_nam(jn)//"d",
! FD     &        "Brine cc., "//biotr_i_nam(jn),
! FD     &        biotr_i_uni(jn),"time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Brine tracer conc created ', 
! FD     &   biotr_i_nam(jn)//"d"
! FD         ENDIF
! FD
! FD         IF ( flag_active(jn) ) THEN
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         CALL CF_CREATE_VAR(biotr_i_nam(jn)//"b",
! FD     &        "Bulk cc., "//biotr_i_nam(jn),
! FD     &        biotr_i_uni(jn),"time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Brine tracer conc created ',
! FD     &   biotr_i_nam(jn)//"b"
! FD         ENDIF
! FD
! FD         IF ( flag_active(jn) ) THEN
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         CALL CF_CREATE_VAR(biotr_i_nam(jn)//"t",
! FD     &        "Total content, "//biotr_i_nam(jn),
! FD     &        biotr_i_uni(jn)//".m","time","-","-","-" )
! FD         WRITE(numout,*) ' Total tracer content created ',
! FD     &   biotr_i_nam(jn)//"t"
! FD         ENDIF
! FD
! FD         IF ( flag_active(jn) ) THEN
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         CALL CF_CREATE_VAR("F"//biotr_i_nam(jn),
! FD     &        "Ice-ocean "//biotr_i_nam(jn)//" flux",
! FD     &        biotr_i_uni(jn)//".m/s","time","-","-","-")
! FD         WRITE(numout,*) ' Tracer flux dimension created ',
! FD     &   biotr_i_nam(jn)//"t"
! FD         ENDIF
! FD
! FD         IF ( flag_active(jn) ) THEN
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         CALL CF_CREATE_VAR("F"//biotr_i_nam(jn)//"b",
! FD     &        "Ice-ocean "//biotr_i_nam(jn)//" basal flux",
! FD     &        biotr_i_uni(jn)//".m/s","time","-","-","-")
! FD         WRITE(numout,*) ' Tracer flux dimension created ',
! FD     &   biotr_i_nam(jn)//"t"
! FD         ENDIF
! FD
! FD         IF ( flag_active(jn) ) THEN
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         CALL CF_CREATE_VAR("F"//biotr_i_nam(jn)//"si",
! FD     &        "Ice-ocean "//biotr_i_nam(jn)//" snow ice flux",
! FD     &        biotr_i_uni(jn)//".m/s","time","-","-","-")
! FD         WRITE(numout,*) ' Tracer flux dimension created ',
! FD     &   biotr_i_nam(jn)//"t"
! FD         ENDIF
! FD
! FD         IF ( flag_active(jn) ) THEN
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         CALL CF_CREATE_VAR("F"//biotr_i_nam(jn)//"bmax",
! FD     &        "Ice-ocean "//biotr_i_nam(jn)//" maximal flux",
! FD     &        biotr_i_uni(jn)//".m/s","time","-","-","-")
! FD         WRITE(numout,*) ' Tracer flux dimension created ',
! FD     &   biotr_i_nam(jn)//"t"
! FD         ENDIF
! FD
! FD         END DO ! jn
! FD
! FD         ! Chlorophyll a
! FD         CALL CF_CREATE_VAR("Chla_bio", "Chlorophyll a on bio grid",
! FD     &        "mmol.m-3","time","z_i_b","-","-")
! FD         WRITE(numout,*) ' Chla bio created '
! FD
! FD         CALL CF_CREATE_VAR("Chla_phy", "Chlorophyll a on phy grid",
! FD     &        "mmol.m-3","time","z_i_p","-","-" )
! FD         WRITE(numout,*) ' Chla phy created '
! FD
! FD         ! PAR - PUR
! FD         CALL CF_CREATE_VAR("PAR","Photosynth. available radiation",
! FD     &                      "\B5E.m-2.s-1","time","z_i_p","-","-" )
! FD         WRITE(numout,*) ' PAR created '
! FD
! FD         CALL CF_CREATE_VAR("PUR","Photosynth. used radiation",
! FD     &                      "\B5E.m-2.s-1","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' PUR created '
! FD
! FD         IF ( flag_active(3) .OR. flag_active(4) ) THEN
! FD         ! Sources and sinks
! FD         CALL CF_CREATE_VAR("syn","Diatom synthesis",
! FD     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Diatom growth created '
! FD
! FD         CALL CF_CREATE_VAR("lys","Diatom lysis",
! FD     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Diatom lysis created '
! FD
! FD         CALL CF_CREATE_VAR("rem","Remineralization",
! FD     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Remineralization created '
! FD
! FD         CALL CF_CREATE_VAR("lim_lig","Light limitation",
! FD     .                   "-","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Light limitation created '
! FD
! FD         CALL CF_CREATE_VAR("lim_dsi","DSi limitation",
! FD     .                   "-","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' DSi limitation created '
! FD
! FD         CALL CF_CREATE_VAR("lim_no3","NO3 limitation",
! FD     .                   "-","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' NO3 limitation created '
! FD
! FD         CALL CF_CREATE_VAR("lim_tem","Temp limitation",
! FD     .                   "-","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Temp limitation created '
! FD
! FD         CALL CF_CREATE_VAR("lim_sal","Sal limitation",
! FD     .                   "-","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Sal limitation created '
! FD
! FD         CALL CF_CREATE_VAR("phi","Photosynthesis",
! FD     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Photosynthesis created '
! FD
! FD         CALL CF_CREATE_VAR("resp","Respiration",
! FD     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Respiration created '
! FD
! FD         CALL CF_CREATE_VAR("excr","Excretion",
! FD     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Respiration created '
! FD
! FD         CALL CF_CREATE_VAR("l_sr","Reservoir lysis",
! FD     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
! FD         WRITE(numout,*) ' Reservoir lysis created '
! FD
! FD         ENDIF

         !---------------
         ! File creation
         !---------------
         CALL CF_CREATE_FILE(filenc)

      DO layer = 1, nlay_s
         dummy_vector(layer) =  REAL( layer ) - 0.5 
      END DO
      CALL CF_WRITE (filenc, 'z_s_p', 1, nlay_s, 1, 1, 
     &               dummy_vector )
      DO layer = 1, nlay_i
         dummy_vector(layer) =  REAL( layer ) - 0.5 
      END DO
      CALL CF_WRITE (filenc, 'z_i_p', 1, nlay_i, 1, 1, 
     &               dummy_vector )
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) =  REAL( layer ) - 0.5 
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'z_i_b', 1, nlay_bio, 1, 1, 
! FD     &               dummy_vector )
      ENDIF ! FD end of numit=nstart

!
!-----------------------------------------------------------------------------!
! 2) Write in the netcdf
!-----------------------------------------------------------------------------!
!
      IF (mod(numit-nstart+1,nfr_out).ne.0) RETURN ! only do below every nfr_out timesteps


      WRITE(numout,*) ' Write into the NetCdf file ',numit
      ji = 1
      !-----------
      ! Open file
      !-----------
      CALL CF_OPEN  (filenc,id)
      !-------------------
      ! Write in the file
      !-------------------
      ! Time
      !------
! FD modifications
      dimval(1)=REAL(numit)*real(ddtb)/86400.
      CALL CF_WRITE (filenc, 'time', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
      dimval(1)=REAL(numit)/86400.0*REAL(ddtb)/30.0
      CALL CF_WRITE (filenc, 'ts_m', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
      dimval(1)=REAL(numit)/86400.*REAL(ddtb)/360.0
      CALL CF_WRITE (filenc, 'ts_y', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)

      ! Grids
      !-------
      ! vertical cotes, physical
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL( ht_i_b(ji) ) / REAL( nlay_i ) * 
     &                       ( REAL( layer ) - 1./2. )
      END DO
      CALL CF_WRITE (filenc, 'z_ip', 
     & (numit-nstart+1)/nfr_out, nlay_i, 1, 1, 
     &               dummy_vector )

! FD      ! vertical cotes, biological
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL( z_i_bio(layer) )
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'z_ib',
! FD     & (numit-nstart+1)/nfr_out, nlay_bio, 1, 1, 
! FD     &               dummy_vector )

      ! Ice / snow thickness and variations
      !-------------------------------------
! FD debug
      dimval(1)=REAL(ht_i_b(ji))
      CALL CF_WRITE (filenc, 'h_i ', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(dh_i_bott(ji))
      CALL CF_WRITE (filenc, 'dhib', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(dh_i_surf(ji))
      CALL CF_WRITE (filenc, 'dhisu', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
! FD      WRITE(numout,*) ' ice output, dh_snowice1: ', dh_snowice(ji)
! FD debug
      dimval(1)=REAL(dh_snowice(ji))
      CALL CF_WRITE (filenc, 'dhisi', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
! FD      WRITE(numout,*) ' ice output, dh_snowice2: ', 
! FD     &               REAL(dh_snowice(ji) / ddtb * 86400.0 * 100.0)

! FD debug
      dimval(1)=REAL(ht_s_b(ji))
      CALL CF_WRITE (filenc, 'h_s ', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(dh_s_tot(ji) / ddtb * 86400.0 * 100.0)
      CALL CF_WRITE (filenc, 'dhs', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(t_su_b(ji) - 273.16)
      CALL CF_WRITE (filenc, 't_su', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(fsb)
      CALL CF_WRITE (filenc, 'fsb', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)
      dimval(1)=REAL(fsbp)
      CALL CF_WRITE (filenc, 'fsbp', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval)


      ! Forcing
      !---------
! FD modifications
      dimval(1)=REAL(emig*ratbqb(ji) )
      CALL CF_WRITE (filenc, 'Flwd', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      dimval(1)=- REAL(fratsb(ji) - emig*ratbqb(ji))
      CALL CF_WRITE (filenc, 'Flwu', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(fsolgb(ji) / ( 1.0 - albgb(ji) ) )
      CALL CF_WRITE (filenc, 'F_sw', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(ab(ji)*fsolgb(ji) )
      CALL CF_WRITE (filenc, 'Fswn', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(fcsb(ji))
      CALL CF_WRITE (filenc, 'F_sh', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(fleb(ji)) 
      CALL CF_WRITE (filenc, 'F_lh', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      WRITE(numout,*) ' albg  : ', albg(1,1)
      WRITE(numout,*) ' albgb : ', albgb(ji)
      dimval(1)=REAL(albgb(ji))
      CALL CF_WRITE (filenc, 'albe', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )

      dimval(1)=REAL(tabqb(ji) - 273.15 )
      CALL CF_WRITE (filenc, 'Tair', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(qabqb(ji))
      CALL CF_WRITE (filenc, 'qair', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(vabqb(ji))
      CALL CF_WRITE (filenc, 'wspd', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(qsfcb(ji))
      CALL CF_WRITE (filenc, 'qsfc', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(tdewb(ji))*100.0
      CALL CF_WRITE (filenc, 'tdew', (numit-nstart+1)/nfr_out, 1, 1, 1, 
     &               dimval )

      ! Inner ice variables  
      !---------------------
      ! temperatures
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL(t_i_b(ji,layer) - 273.16 )
      END DO
      CALL CF_WRITE (filenc, 't_i ', (numit-nstart+1)/nfr_out, nlay_i,
     & 1, 1, 
     &               dummy_vector )

      ! salinity
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL( s_i_b(ji,layer) )
      END DO
      CALL CF_WRITE (filenc, 's_i ', (numit-nstart+1)/nfr_out, nlay_i,
     & 1, 1, 
     &               dummy_vector )

      ! brine volume
      DO layer = 1, nlay_i
! FD         dummy_vector(layer) = REAL ( - tmut * s_i_b(ji,layer) / 
! FD     &                         ( t_i_b(ji,layer) - 273.16) ) * 100.0
! FD replace the analytical form based on T = - mu S
         dummy_vector(layer) = REAL ( brvolum(layer) ) * 100.0
      END DO
      CALL CF_WRITE (filenc, 'e_i ', (numit-nstart+1)/nfr_out, nlay_i,
     & 1, 1, 
     &               dummy_vector )

      ! Rayleigh number
      DO layer = 1, nlay_i
         dummy_vector(layer) = rayleigh(layer)
      END DO
      CALL CF_WRITE (filenc, 'ra', (numit-nstart+1)/nfr_out, nlay_i,
     & 1, 1, 
     &               dummy_vector )
     
      ! Brine flushing velocity
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL( brflush(layer) )
      END DO
      CALL CF_WRITE (filenc, 'ub', (numit-nstart+1)/nfr_out, nlay_i,
     & 1, 1, 
     &               dummy_vector )
     
      ! Salt diffusivity in brine
      DO layer = 1, nlay_i
         dummy_vector(layer) =diff_br(layer)
      END DO
      CALL CF_WRITE (filenc, 'diff', (numit-nstart+1)/nfr_out, nlay_i,
     & 1, 1, 
     &               dummy_vector )

      ! snow temperature
      DO layer = 1, nlay_s
         dummy_vector(layer) = REAL(t_s_b(ji,layer) - 273.16 )
      END DO
      CALL CF_WRITE (filenc, 't_s ', (numit-nstart+1)/nfr_out, nlay_s,
     & 1, 1, 
     &               dummy_vector )

      ! Radiative transfer
      !--------------------
      ! Transmitted radiation in the ice
! FD modifications
      dimval(1)=REAL(ftrice)
      CALL CF_WRITE (filenc, 'T_is', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval )
      dimval(1)=REAL(ftroce)
      CALL CF_WRITE (filenc, 'T_oc', (numit-nstart+1)/nfr_out, 
     &               1, 1, 1, dimval )
      ! Absorbed radiation in snow, physically
      DO layer = 1, nlay_s
         dummy_vector(layer) = REAL(radab_s(layer))
      END DO
      CALL CF_WRITE (filenc, 'A_s', (numit-nstart+1)/nfr_out, nlay_s,
     & 1, 1, 
     &               dummy_vector )
      ! Absorbed radiation in ice, physically
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL(radab_phy_i(layer))
      END DO
      CALL CF_WRITE (filenc, 'A_ip', (numit-nstart+1)/nfr_out, nlay_i,
     & 1, 1, 
     &               dummy_vector )
      ! Absorbed radiation in ice, biologically
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL(radab_alg_i(layer))
      END DO
      CALL CF_WRITE (filenc, 'A_ib', (numit-nstart+1)/nfr_out, nlay_i,
     & 1, 1, 
     &               dummy_vector )

! FD      ! Biological variables  
! FD      !----------------------
! FD      ! tracer brine concentration
! FD      DO jn = 1, ntra_bio
! FD
! FD         IF ( flag_diff(jn) .AND. flag_active(jn) ) THEN
! FD            DO layer = 1, nlay_bio
! FD               dummy_vector(layer) = REAL( c_i_bio(jn,layer) )
! FD            END DO
! FD            WRITE(numout,*) ' Tracer number : ', jn
! FD            WRITE(numout,*) ' dummy_vector : ', 
! FD     &                      ( dummy_vector(layer), layer = 1, nlay_bio )
! FD            CALL CF_WRITE (filenc, biotr_i_nam(jn)//"d", 
! FD     &                     (numit-nstart+1)/nfr_out, nlay_bio, 1, 1, 
! FD     &                     dummy_vector )
! FD         ENDIF
! FD
! FD         ! tracer bulk ice concentration
! FD         IF ( flag_active(jn) ) THEN
! FD         DO layer = 1, nlay_bio
! FD            dummy_vector(layer) = REAL( cbu_i_bio(jn,layer) )
! FD         END DO
! FD
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         WRITE(numout,*) ' dummy_vector : ', 
! FD     &                   ( dummy_vector(layer), layer = 1, nlay_bio )
! FD
! FD         CALL CF_WRITE (filenc, biotr_i_nam(jn)//"b", 
! FD     &                  (numit-nstart+1)/nfr_out, nlay_bio, 1, 1, 
! FD     &                  dummy_vector )
! FD
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         WRITE(numout,*) ' dummy_vector : ', 
! FD     &                   REAL( ct_i_bio(jn) )
! FD
! FD! FD modifications
! FD      dimval(1)=REAL( ct_i_bio(jn) )
! FD         CALL CF_WRITE (filenc, biotr_i_nam(jn)//"t",
! FD     &                  (numit-nstart+1)/nfr_out, 1, 1, 1, 
! FD     &                  dimval )
! FD
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         WRITE(numout,*) ' dummy_vector : ', 
! FD     &                   REAL( fcb(jn) )
! FD
! FD! FD modifications
! FD      dimval(1)=REAL( fcb(jn) )
! FD         CALL CF_WRITE (filenc, "F"//biotr_i_nam(jn),
! FD     &                  (numit-nstart+1)/nfr_out, 1, 1, 1, 
! FD     &                  dimval )
! FD
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         WRITE(numout,*) ' dummy_vector : ', 
! FD     &                   REAL( fcbp(jn) )
! FD
! FD! FD modifications
! FD      dimval(1)=REAL( fcbp(jn) )
! FD         CALL CF_WRITE (filenc, "F"//biotr_i_nam(jn)//"b",
! FD     &                  (numit-nstart+1)/nfr_out, 1, 1, 1, 
! FD     &                  dimval )
! FD
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         WRITE(numout,*) ' dummy_vector : ', 
! FD     &                   REAL( fcsi(jn) )
! FD
! FD! FD modifications
! FD      dimval(1)=REAL( fcsi(jn) )
! FD         CALL CF_WRITE (filenc, "F"//biotr_i_nam(jn)//"si",
! FD     &                  (numit-nstart+1)/nfr_out, 1, 1, 1, 
! FD     &                  dimval )
! FD
! FD         WRITE(numout,*) ' Tracer number : ', jn
! FD         WRITE(numout,*) ' dummy_vector : ', 
! FD     &                   REAL( fcb_max(jn) )
! FD
! FD! FD modifications
! FD      dimval(1)=REAL( fcb_max(jn) )
! FD         CALL CF_WRITE (filenc, "F"//biotr_i_nam(jn)//"bmax",
! FD     &                  (numit-nstart+1)/nfr_out, 1, 1, 1, 
! FD     &                  dimval )
! FD
! FD         ENDIF
! FD      END DO ! jn
! FD
! FD      ! Chla, physical grid
! FD      DO layer = 1, nlay_i
! FD         dummy_vector(layer) = REAL(chla_i(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'Chla_phy', (numit-nstart+1)/nfr_out,
! FD     & nlay_i, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! Chla, biological grid
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(chla_i_bio(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'Chla_bio', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1,
! FD     &               dummy_vector )
! FD
! FD      ! PAR-PUR
! FD      DO layer = 1, nlay_i
! FD         dummy_vector(layer) = REAL(par(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'PAR', (numit-nstart+1)/nfr_out,
! FD     & nlay_i, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(pur_bio(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'PUR', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      !-------------------
! FD      ! Sources and sinks
! FD      !-------------------
! FD
! FD      IF ( flag_active(3) .OR. flag_active(4) ) THEN
! FD
! FD      ! diatom growth
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(syn_bio(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'syn', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! diatom lysis
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(lys_bio(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'lys', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! remineralization
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(rem_bio(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'rem', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! Light limitation
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(lim_lig(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'lim_lig', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! DSi limitation
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(lim_dsi(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'lim_dsi', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! NO3 limitation
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(lim_no3(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'lim_no3', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! Temperature limitation
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(lim_tem(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'lim_tem', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! Salinity limitation
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(lim_sal(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'lim_sal', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! Photosynthesis
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(phi_bio(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'phi', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! Respiration
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(resp_bio(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'resp', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! Excretion
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(excr_bio(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'excr', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ! Reservoir lysis
! FD      DO layer = 1, nlay_bio
! FD         dummy_vector(layer) = REAL(lsr_bio(layer))
! FD      END DO
! FD      CALL CF_WRITE (filenc, 'l_sr', (numit-nstart+1)/nfr_out,
! FD     & nlay_bio, 1, 1, 
! FD     &               dummy_vector )
! FD
! FD      ENDIF

      !------------
      ! Close file
      !------------
      CALL CF_CLOSE (filenc)
!
!-----------------------------------------------------------------------------!
!
! FD      WRITE(numout,*) 

      END SUBROUTINE
