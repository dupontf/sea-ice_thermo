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

      WRITE(numout,*) ' ice_output : '
      WRITE(numout,*) ' ~~~~~~~~~~ : '
!
!-----------------------------------------------------------------------------!
! 1) File definition
!-----------------------------------------------------------------------------!
!
      IF ( numit .EQ. nstart ) THEN

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
         DO i = nstart , nend
            dimval(i) = REAL(i) 
         ENDDO
      !  WRITE(numout,*) ' dimval : ', ( dimval(i), i = nstart, nend )
         IF ( ddtb .EQ. 3600 ) 
     &      CALL CF_CREATE_DIM("time"    ,"hours"     ,nitrun,
     &      dimval(nstart:nend))
         IF ( ddtb .EQ. 86400 ) 
     &      CALL CF_CREATE_DIM("time"    ,"days"     ,nitrun,
     &      dimval(nstart:nend))
         WRITE(numout,*) ' time dimension created '

         ! snow vertical coordinate (moving)
         !-----------------------------------
         DO i = 1, nlay_s
            dimval4(i) = REAL(i)
         ENDDO
         CALL CF_CREATE_DIM("z_s_p"   ,"m"       ,nlay_s,
     &   dimval4(1:nlay_s))
         WRITE(numout,*) ' z_s_p dimension created '

         ! physical vertical coordinate (moving)
         !---------------------------------------
         DO i = 1, nlay_i
            dimval2(i) = REAL(i)
         ENDDO
         CALL CF_CREATE_DIM("z_i_p"   ,"m"       ,nlay_i,
     &   dimval2(1:nlay_i))
         WRITE(numout,*) ' z_i_p dimension created '

!        ! biological vertical coordinate (moving)
!        !-----------------------------------------
         DO i = 1, nlay_bio 
            dimval3(i) = REAL(i)
         ENDDO
         CALL CF_CREATE_DIM("z_i_b"   ,"m"       ,nlay_bio,
     &   dimval3(1:nlay_bio))
         WRITE(numout,*) ' z_i_b dimension created '

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

         CALL CF_CREATE_VAR("z_ib","Vertical cote, biology",
     .                   "m","time","z_i_b","-","-" )

         ! Ice / snow thickness and variations
         !-------------------------------------
         CALL CF_CREATE_VAR("h_i","Ice thickness",
     .                   "m","time","-","-","-" )
         WRITE(numout,*) ' h_i dimension created '

         CALL CF_CREATE_VAR("dhib","Bottom accretion",
     .                   "cm/day","time","-","-","-" )

         CALL CF_CREATE_VAR("dhisu","Surface ice melt",
     .                   "cm/day","time","-","-","-" )

         CALL CF_CREATE_VAR("dhisi","Snow ice formation",
     .                   "cm/day","time","-","-","-" )

         CALL CF_CREATE_VAR("dhs","Snow variations",
     .                   "cm/day","time","-","-","-" )

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

         ! Biological variables
         !----------------------
         DO jn = 1, ntra_bio

         IF ( flag_diff(jn) .AND. flag_active(jn) ) THEN
         WRITE(numout,*) ' Tracer number : ', jn
         CALL CF_CREATE_VAR(biotr_i_nam(jn)//"d",
     &        "Brine cc., "//biotr_i_nam(jn),
     &        biotr_i_uni(jn),"time","z_i_b","-","-" )
         WRITE(numout,*) ' Brine tracer conc created ', 
     &   biotr_i_nam(jn)//"d"
         ENDIF

         IF ( flag_active(jn) ) THEN
         WRITE(numout,*) ' Tracer number : ', jn
         CALL CF_CREATE_VAR(biotr_i_nam(jn)//"b",
     &        "Bulk cc., "//biotr_i_nam(jn),
     &        biotr_i_uni(jn),"time","z_i_b","-","-" )
         WRITE(numout,*) ' Brine tracer conc created ',
     &   biotr_i_nam(jn)//"b"
         ENDIF

         IF ( flag_active(jn) ) THEN
         WRITE(numout,*) ' Tracer number : ', jn
         CALL CF_CREATE_VAR(biotr_i_nam(jn)//"t",
     &        "Total content, "//biotr_i_nam(jn),
     &        biotr_i_uni(jn)//".m","time","-","-","-" )
         WRITE(numout,*) ' Total tracer content created ',
     &   biotr_i_nam(jn)//"t"
         ENDIF

         IF ( flag_active(jn) ) THEN
         WRITE(numout,*) ' Tracer number : ', jn
         CALL CF_CREATE_VAR("F"//biotr_i_nam(jn),
     &        "Ice-ocean "//biotr_i_nam(jn)//" flux",
     &        biotr_i_uni(jn)//".m/s","time","-","-","-")
         WRITE(numout,*) ' Tracer flux dimension created ',
     &   biotr_i_nam(jn)//"t"
         ENDIF

         IF ( flag_active(jn) ) THEN
         WRITE(numout,*) ' Tracer number : ', jn
         CALL CF_CREATE_VAR("F"//biotr_i_nam(jn)//"b",
     &        "Ice-ocean "//biotr_i_nam(jn)//" basal flux",
     &        biotr_i_uni(jn)//".m/s","time","-","-","-")
         WRITE(numout,*) ' Tracer flux dimension created ',
     &   biotr_i_nam(jn)//"t"
         ENDIF

         IF ( flag_active(jn) ) THEN
         WRITE(numout,*) ' Tracer number : ', jn
         CALL CF_CREATE_VAR("F"//biotr_i_nam(jn)//"si",
     &        "Ice-ocean "//biotr_i_nam(jn)//" snow ice flux",
     &        biotr_i_uni(jn)//".m/s","time","-","-","-")
         WRITE(numout,*) ' Tracer flux dimension created ',
     &   biotr_i_nam(jn)//"t"
         ENDIF

         IF ( flag_active(jn) ) THEN
         WRITE(numout,*) ' Tracer number : ', jn
         CALL CF_CREATE_VAR("F"//biotr_i_nam(jn)//"bmax",
     &        "Ice-ocean "//biotr_i_nam(jn)//" maximal flux",
     &        biotr_i_uni(jn)//".m/s","time","-","-","-")
         WRITE(numout,*) ' Tracer flux dimension created ',
     &   biotr_i_nam(jn)//"t"
         ENDIF

         END DO ! jn

         ! Chlorophyll a
         CALL CF_CREATE_VAR("Chla_bio", "Chlorophyll a on bio grid",
     &        "mmol.m-3","time","z_i_b","-","-")
         WRITE(numout,*) ' Chla bio created '

         CALL CF_CREATE_VAR("Chla_phy", "Chlorophyll a on phy grid",
     &        "mmol.m-3","time","z_i_p","-","-" )
         WRITE(numout,*) ' Chla phy created '

         ! PAR - PUR
         CALL CF_CREATE_VAR("PAR","Photosynth. available radiation",
     &                      "µE.m-2.s-1","time","z_i_p","-","-" )
         WRITE(numout,*) ' PAR created '

         CALL CF_CREATE_VAR("PUR","Photosynth. used radiation",
     &                      "µE.m-2.s-1","time","z_i_b","-","-" )
         WRITE(numout,*) ' PUR created '

         IF ( flag_active(3) .OR. flag_active(4) ) THEN
         ! Sources and sinks
         CALL CF_CREATE_VAR("syn","Diatom synthesis",
     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
         WRITE(numout,*) ' Diatom growth created '

         CALL CF_CREATE_VAR("lys","Diatom lysis",
     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
         WRITE(numout,*) ' Diatom lysis created '

         CALL CF_CREATE_VAR("rem","Remineralization",
     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
         WRITE(numout,*) ' Remineralization created '

         CALL CF_CREATE_VAR("lim_lig","Light limitation",
     .                   "-","time","z_i_b","-","-" )
         WRITE(numout,*) ' Light limitation created '

         CALL CF_CREATE_VAR("lim_dsi","DSi limitation",
     .                   "-","time","z_i_b","-","-" )
         WRITE(numout,*) ' DSi limitation created '

         CALL CF_CREATE_VAR("lim_no3","NO3 limitation",
     .                   "-","time","z_i_b","-","-" )
         WRITE(numout,*) ' NO3 limitation created '

         CALL CF_CREATE_VAR("lim_tem","Temp limitation",
     .                   "-","time","z_i_b","-","-" )
         WRITE(numout,*) ' Temp limitation created '

         CALL CF_CREATE_VAR("lim_sal","Sal limitation",
     .                   "-","time","z_i_b","-","-" )
         WRITE(numout,*) ' Sal limitation created '

         CALL CF_CREATE_VAR("phi","Photosynthesis",
     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
         WRITE(numout,*) ' Photosynthesis created '

         CALL CF_CREATE_VAR("resp","Respiration",
     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
         WRITE(numout,*) ' Respiration created '

         CALL CF_CREATE_VAR("excr","Excretion",
     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
         WRITE(numout,*) ' Respiration created '

         CALL CF_CREATE_VAR("l_sr","Reservoir lysis",
     .                   "mmol C m-3 h-1","time","z_i_b","-","-" )
         WRITE(numout,*) ' Reservoir lysis created '

         ENDIF

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
      DO layer = 1, nlay_bio
         dummy_vector(layer) =  REAL( layer ) - 0.5 
      END DO
      CALL CF_WRITE (filenc, 'z_i_b', 1, nlay_bio, 1, 1, 
     &               dummy_vector )
      ENDIF ! FD end of numit=nstart

!
!-----------------------------------------------------------------------------!
! 2) Write in the netcdf
!-----------------------------------------------------------------------------!
!
      WRITE(numout,*) ' Write into the NetCdf file '
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
      dimval(1)=REAL(numit)
      CALL CF_WRITE (filenc, 'time', numit-nstart+1, 
     &               1, 1, 1, dimval)
      dimval(1)=REAL(numit)/86400.0*REAL(ddtb)/30.0
      CALL CF_WRITE (filenc, 'ts_m', numit-nstart+1, 
     &               1, 1, 1, dimval)
      dimval(1)=REAL(numit)/86400.*REAL(ddtb)/365.0
      CALL CF_WRITE (filenc, 'ts_y', numit-nstart+1, 
     &               1, 1, 1, dimval)

      ! Grids
      !-------
      ! vertical cotes, physical
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL( ht_i_b(ji) ) / REAL( nlay_i ) * 
     &                       ( REAL( layer ) - 1./2. )
      END DO
      CALL CF_WRITE (filenc, 'z_ip', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )

      ! vertical cotes, biological
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL( z_i_bio(layer) )
      END DO
      CALL CF_WRITE (filenc, 'z_ib', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! Ice / snow thickness and variations
      !-------------------------------------
! FD debug
      dimval(1)=REAL(ht_i_b(ji))
      CALL CF_WRITE (filenc, 'h_i ', numit-nstart+1, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(dh_i_bott(ji) / ddtb * 86400.0 * 100.0)
      CALL CF_WRITE (filenc, 'dhib', numit-nstart+1, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(dh_i_surf(ji) / ddtb * 86400.0 * 100.0)
      CALL CF_WRITE (filenc, 'dhisu', numit-nstart+1, 
     &               1, 1, 1, dimval)
      WRITE(numout,*) ' ice output, dh_snowice1: ', dh_snowice(ji)
! FD debug
      dimval(1)=REAL(dh_snowice(ji) / ddtb * 86400.0 * 100.0)
      CALL CF_WRITE (filenc, 'dhisi', numit-nstart+1, 
     &               1, 1, 1, dimval)
      WRITE(numout,*) ' ice output, dh_snowice2: ', 
     &               REAL(dh_snowice(ji) / ddtb * 86400.0 * 100.0)

! FD debug
      dimval(1)=REAL(ht_s_b(ji))
      CALL CF_WRITE (filenc, 'h_s ', numit-nstart+1, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(dh_s_tot(ji) / ddtb * 86400.0 * 100.0)
      CALL CF_WRITE (filenc, 'dhs', numit-nstart+1, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(t_su_b(ji) - 273.16)
      CALL CF_WRITE (filenc, 't_su', numit-nstart+1, 
     &               1, 1, 1, dimval)
! FD debug
      dimval(1)=REAL(fsb)
      CALL CF_WRITE (filenc, 'fsb', numit-nstart+1, 
     &               1, 1, 1, dimval)
      dimval(1)=REAL(fsbp)
      CALL CF_WRITE (filenc, 'fsbp', numit-nstart+1, 
     &               1, 1, 1, dimval)


      ! Forcing
      !---------
! FD modifications
      dimval(1)=REAL(emig*ratbqb(ji) )
      CALL CF_WRITE (filenc, 'Flwd', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      dimval(1)=- REAL(fratsb(ji) - emig*ratbqb(ji))
      CALL CF_WRITE (filenc, 'Flwu', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(fsolgb(ji) / ( 1.0 - albgb(ji) ) )
      CALL CF_WRITE (filenc, 'F_sw', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(ab(ji)*fsolgb(ji) )
      CALL CF_WRITE (filenc, 'Fswn', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(fcsb(ji))
      CALL CF_WRITE (filenc, 'F_sh', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(fleb(ji)) 
      CALL CF_WRITE (filenc, 'F_lh', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      WRITE(numout,*) ' albg  : ', albg(1,1)
      WRITE(numout,*) ' albgb : ', albgb(ji)
      dimval(1)=REAL(albgb(ji))
      CALL CF_WRITE (filenc, 'albe', numit-nstart+1, 1, 1, 1, 
     &               dimval )

      dimval(1)=REAL(tabqb(ji) - 273.15 )
      CALL CF_WRITE (filenc, 'Tair', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(qabqb(ji))
      CALL CF_WRITE (filenc, 'qair', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(vabqb(ji))
      CALL CF_WRITE (filenc, 'wspd', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(qsfcb(ji))
      CALL CF_WRITE (filenc, 'qsfc', numit-nstart+1, 1, 1, 1, 
     &               dimval )
      dimval(1)=REAL(tdewb(ji))*100.0
      CALL CF_WRITE (filenc, 'tdew', numit-nstart+1, 1, 1, 1, 
     &               dimval )

      ! Inner ice variables  
      !---------------------
      ! temperatures
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL(t_i_b(ji,layer) - 273.16 )
      END DO
      CALL CF_WRITE (filenc, 't_i ', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )

      ! salinity
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL( s_i_b(ji,layer) )
      END DO
      CALL CF_WRITE (filenc, 's_i ', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )

      ! brine volume
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL ( - tmut * s_i_b(ji,layer) / 
     &                         ( t_i_b(ji,layer) - 273.16) ) * 100.0
      END DO
      CALL CF_WRITE (filenc, 'e_i ', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )

      ! Rayleigh number
      DO layer = 1, nlay_i
         dummy_vector(layer) = rayleigh(layer)
      END DO
      CALL CF_WRITE (filenc, 'ra', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )
     
      ! Salt diffusivity in brine
      DO layer = 1, nlay_i
         dummy_vector(layer) =diff_br(layer)
      END DO
      CALL CF_WRITE (filenc, 'diff', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )

      ! snow temperature
      DO layer = 1, nlay_s
         dummy_vector(layer) = REAL(t_s_b(ji,layer) - 273.16 )
      END DO
      CALL CF_WRITE (filenc, 't_s ', numit-nstart+1, nlay_s, 1, 1, 
     &               dummy_vector )

      ! Radiative transfer
      !--------------------
      ! Transmitted radiation in the ice
! FD modifications
      dimval(1)=REAL(ftrice)
      CALL CF_WRITE (filenc, 'T_is', numit-nstart+1, 
     &               1, 1, 1, dimval )
      dimval(1)=REAL(ftroce)
      CALL CF_WRITE (filenc, 'T_oc', numit-nstart+1, 
     &               1, 1, 1, dimval )
      ! Absorbed radiation in snow, physically
      DO layer = 1, nlay_s
         dummy_vector(layer) = REAL(radab_s(layer))
      END DO
      CALL CF_WRITE (filenc, 'A_s', numit-nstart+1, nlay_s, 1, 1, 
     &               dummy_vector )
      ! Absorbed radiation in ice, physically
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL(radab_phy_i(layer))
      END DO
      CALL CF_WRITE (filenc, 'A_ip', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )
      ! Absorbed radiation in ice, biologically
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL(radab_alg_i(layer))
      END DO
      CALL CF_WRITE (filenc, 'A_ib', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )

      ! Biological variables  
      !----------------------
      ! tracer brine concentration
      DO jn = 1, ntra_bio

         IF ( flag_diff(jn) .AND. flag_active(jn) ) THEN
            DO layer = 1, nlay_bio
               dummy_vector(layer) = REAL( c_i_bio(jn,layer) )
            END DO
            WRITE(numout,*) ' Tracer number : ', jn
            WRITE(numout,*) ' dummy_vector : ', 
     &                      ( dummy_vector(layer), layer = 1, nlay_bio )
            CALL CF_WRITE (filenc, biotr_i_nam(jn)//"d", 
     &                     numit-nstart+1, nlay_bio, 1, 1, 
     &                     dummy_vector )
         ENDIF

         ! tracer bulk ice concentration
         IF ( flag_active(jn) ) THEN
         DO layer = 1, nlay_bio
            dummy_vector(layer) = REAL( cbu_i_bio(jn,layer) )
         END DO

         WRITE(numout,*) ' Tracer number : ', jn
         WRITE(numout,*) ' dummy_vector : ', 
     &                   ( dummy_vector(layer), layer = 1, nlay_bio )

         CALL CF_WRITE (filenc, biotr_i_nam(jn)//"b", 
     &                  numit-nstart+1, nlay_bio, 1, 1, 
     &                  dummy_vector )

         WRITE(numout,*) ' Tracer number : ', jn
         WRITE(numout,*) ' dummy_vector : ', 
     &                   REAL( ct_i_bio(jn) )

! FD modifications
      dimval(1)=REAL( ct_i_bio(jn) )
         CALL CF_WRITE (filenc, biotr_i_nam(jn)//"t",
     &                  numit-nstart+1, 1, 1, 1, 
     &                  dimval )

         WRITE(numout,*) ' Tracer number : ', jn
         WRITE(numout,*) ' dummy_vector : ', 
     &                   REAL( fcb(jn) )

! FD modifications
      dimval(1)=REAL( fcb(jn) )
         CALL CF_WRITE (filenc, "F"//biotr_i_nam(jn),
     &                  numit-nstart+1, 1, 1, 1, 
     &                  dimval )

         WRITE(numout,*) ' Tracer number : ', jn
         WRITE(numout,*) ' dummy_vector : ', 
     &                   REAL( fcbp(jn) )

! FD modifications
      dimval(1)=REAL( fcbp(jn) )
         CALL CF_WRITE (filenc, "F"//biotr_i_nam(jn)//"b",
     &                  numit-nstart+1, 1, 1, 1, 
     &                  dimval )

         WRITE(numout,*) ' Tracer number : ', jn
         WRITE(numout,*) ' dummy_vector : ', 
     &                   REAL( fcsi(jn) )

! FD modifications
      dimval(1)=REAL( fcsi(jn) )
         CALL CF_WRITE (filenc, "F"//biotr_i_nam(jn)//"si",
     &                  numit-nstart+1, 1, 1, 1, 
     &                  dimval )

         WRITE(numout,*) ' Tracer number : ', jn
         WRITE(numout,*) ' dummy_vector : ', 
     &                   REAL( fcb_max(jn) )

! FD modifications
      dimval(1)=REAL( fcb_max(jn) )
         CALL CF_WRITE (filenc, "F"//biotr_i_nam(jn)//"bmax",
     &                  numit-nstart+1, 1, 1, 1, 
     &                  dimval )

         ENDIF
      END DO ! jn

      ! Chla, physical grid
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL(chla_i(layer))
      END DO
      CALL CF_WRITE (filenc, 'Chla_phy', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )

      ! Chla, biological grid
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(chla_i_bio(layer))
      END DO
      CALL CF_WRITE (filenc, 'Chla_bio', numit-nstart+1, nlay_bio, 1, 1,
     &               dummy_vector )

      ! PAR-PUR
      DO layer = 1, nlay_i
         dummy_vector(layer) = REAL(par(layer))
      END DO
      CALL CF_WRITE (filenc, 'PAR', numit-nstart+1, nlay_i, 1, 1, 
     &               dummy_vector )

      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(pur_bio(layer))
      END DO
      CALL CF_WRITE (filenc, 'PUR', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      !-------------------
      ! Sources and sinks
      !-------------------

      IF ( flag_active(3) .OR. flag_active(4) ) THEN

      ! diatom growth
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(syn_bio(layer))
      END DO
      CALL CF_WRITE (filenc, 'syn', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! diatom lysis
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(lys_bio(layer))
      END DO
      CALL CF_WRITE (filenc, 'lys', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! remineralization
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(rem_bio(layer))
      END DO
      CALL CF_WRITE (filenc, 'rem', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! Light limitation
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(lim_lig(layer))
      END DO
      CALL CF_WRITE (filenc, 'lim_lig', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! DSi limitation
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(lim_dsi(layer))
      END DO
      CALL CF_WRITE (filenc, 'lim_dsi', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! NO3 limitation
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(lim_no3(layer))
      END DO
      CALL CF_WRITE (filenc, 'lim_no3', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! Temperature limitation
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(lim_tem(layer))
      END DO
      CALL CF_WRITE (filenc, 'lim_tem', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! Salinity limitation
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(lim_sal(layer))
      END DO
      CALL CF_WRITE (filenc, 'lim_sal', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! Photosynthesis
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(phi_bio(layer))
      END DO
      CALL CF_WRITE (filenc, 'phi', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! Respiration
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(resp_bio(layer))
      END DO
      CALL CF_WRITE (filenc, 'resp', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! Excretion
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(excr_bio(layer))
      END DO
      CALL CF_WRITE (filenc, 'excr', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ! Reservoir lysis
      DO layer = 1, nlay_bio
         dummy_vector(layer) = REAL(lsr_bio(layer))
      END DO
      CALL CF_WRITE (filenc, 'l_sr', numit-nstart+1, nlay_bio, 1, 1, 
     &               dummy_vector )

      ENDIF

      !------------
      ! Close file
      !------------
      CALL CF_CLOSE (filenc)
!
!-----------------------------------------------------------------------------!
!
      WRITE(numout,*) 

      END SUBROUTINE
