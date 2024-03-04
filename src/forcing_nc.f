      SUBROUTINE forcing_nc(xjour)

        !!------------------------------------------------------------------
        !!                ***       ROUTINE forcing      ***
        !! ** Purpose :
        !!      This routine computes the model forcing
        !!       forc_swi = 0 -> read
        !!                  1 -> computed
        !!                 99 -> prescribed
        !!
        !! ** Method  :
        !!
        !! ** Arguments :
        !!
        !! ** Inputs / Ouputs : (global commons)
        !! 
        !! ** External : 
        !!
        !! ** References : 
        !!
        !! ** History :
        !!
        !!------------------------------------------------------------------
        !! * Arguments

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'forcing.com'

      INTEGER :: 
     &  ji          , ! : index for space
     &  jk          , ! : index for ice layers
     &  jf          , ! : index for forcing field
     &  numforc= 600  ! : reference number for bio.param

      CHARACTER(len=10) :: 
     &   filenc='forcing.nc'

      REAL(4)        zforc(1), zforc2(2) ! forcing field dummy array
      DIMENSION ws(96),zmue(96),zalcnp(96) ! for solar flux formula

      DIMENSION budyko(19)
! FD addition
      LOGICAL ln_write_forc

      DATA budyko /1.00,0.98,0.95,0.92,0.89,0.86,0.83,0.80,0.78,0.75,
     &             0.72,0.69,0.67,0.64,0.61,0.58,0.56,0.53,0.50/
      ln_write_forc = .TRUE.

      WRITE(numout,*) ' * forcing_nc : '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~ '

      ! Control parameters number of steps in the integration of shortwave rad
      nintsr =  24   ! 24 ideally 
      zsolar = 1368. ! solar constant
!
!-----------------------------------------------------------------------------!
! 1) Names of the forcing variables
!-----------------------------------------------------------------------------!
!
      forc_nam(1) = 'fsw'
      forc_nam(2) = 'flw'
      forc_nam(3) = 'par'
      forc_nam(4) = 'tair'
      forc_nam(5) = 'pres'
      forc_nam(6) = 'qair'
      forc_nam(7) = 'wspd'
      forc_nam(8) = 'cld'
      forc_nam(9) = 'foce'
      forc_nam(10)= 'sfal'
      forc_nam(11)= 'albe'
!
!-----------------------------------------------------------------------------!
! 2) Reads namelist
!-----------------------------------------------------------------------------!
!
      IF ( nit .EQ. nstart ) THEN

      IF ( ln_write_forc ) THEN
         WRITE(numout,*) ' Forcing parameters ... '
         WRITE(numout,*)
         WRITE(numout,*) ' forc_nam    : ',( forc_nam(i), i = 1, n_forc)
      ENDIF

      OPEN( unit = numforc, file='forcing.param', status='old' )
      READ(numforc,*)
      READ(numforc,*)
      READ(numforc,*)
      READ(numforc,*)
      READ(numforc,*)
      READ(numforc,*) ts_forc    ! Forcing time step
      READ(numforc,*)
      READ(numforc,*) n0_forc    ! Number of the first time step in the forcing file
      READ(numforc,*)
      READ(numforc,*) n1_forc    ! Number of the last time step in the forcing file
      IF ( ln_write_forc ) WRITE(numout,*) '  ts_forc : ', ts_forc

      READ(numforc,*)
      READ(numforc,*)
      READ(numforc,*)
      READ(numforc,*)
      READ(numforc,*)
      READ(numforc,*)
      READ(numforc,*)

      IF ( ln_write_forc ) WRITE(numout,*) '  forc_swi, forc_uni, 
     & forc_val : '
      DO i = 1, n_forc
         READ(numforc,*) idum, dummy1, dummy2
         READ(numforc,*)
         forc_swi(i) = idum
         forc_uni(i) = dummy1
         forc_val(i) = dummy2
         IF ( ln_write_forc ) WRITE(numout,*) forc_nam(i), forc_swi(i), 
     &                                        forc_uni(i), forc_val(i)
      END DO

      CALL CF_OPEN  (filenc,id) ! open forcing file
      
      IF ( ln_write_forc ) WRITE(numout,*)

      ENDIF ! nit=nstart
!
!-----------------------------------------------------------------------------!
! 3) Read the variables
!-----------------------------------------------------------------------------!
!

      IF ( nit .EQ. nstart ) THEN
         n_fofr = INT( ts_forc / ddtb ) ! forcing frequency
         n_forc_min = n_fofr / 2 + nstart 
         n_forc_max = ( nitrun - INT( FLOAT(n_fofr) / 2. ) ) + nstart
         WRITE(numout,*) ' n_fofr : ', n_fofr
         WRITE(numout,*) ' nstart     : ', nstart
         WRITE(numout,*) ' n_forc_min : ', n_forc_min
         WRITE(numout,*) ' n_forc_max : ', n_forc_max
         WRITE(numout,*) ' nend       : ', nend
         i_forc_day = nstart / n_fofr - n0_forc
      ENDIF

      i_forc = MOD( ( nit - nstart ) - n_fofr / 2 , n_fofr) ! indicates if forcing is to be read
      WRITE(numout,*) ' nit - nstart : ', nit - nstart
      WRITE(numout,*) ' n_fofr       : ', n_fofr
      WRITE(numout,*) ' i_forc : ', i_forc

      !-------------------------------------
      ! Time steps at which values are READ
      !-------------------------------------
      IF ( ( nit .EQ. nstart ) .OR. ( i_forc .EQ. 0.0 ) ) THEN
         i_forc_count = 0
         i_forc_day   = i_forc_day + 1

         IF ( i_forc_day .GT. n1_forc ) i_forc_day = 1

         WRITE(numout,*) ' Forcing read at this time step, nit : ', nit
         WRITE(numout,*) ' Forcing step, i_forc_day : ', i_forc_day

         DO i = 1, n_forc
            IF ( forc_swi(i) .EQ. 0 ) THEN 
               IF ( nit .GT. nstart ) forc_arr_old(i) = forc_arr_new(i)
               IF ( nit .LT. n_forc_max ) THEN
                  CALL CF_READ1D ( filenc, forc_nam(i), 
     &                 i_forc_day , 1, zforc )
                  forc_arr_new(i) = DBLE(zforc(1)) * forc_uni(i)
               ENDIF
               IF ( nit .EQ. nstart ) forc_arr_old(i) = forc_arr_new(i)
               forc_coeff(i) = ( forc_arr_new(i) - forc_arr_old(i) ) /
     &                           FLOAT(n_fofr)
!              WRITE(numout,*) ' forc_coeff : ', forc_coeff(i)
!              WRITE(numout,*) ' forc_arr_old:', forc_arr_old(i)
!              WRITE(numout,*) ' forc_arr_new:', forc_arr_new(i)
            ENDIF
         END DO
       
      ENDIF
      !-----------------------------------------
      ! Time steps at which values are COMPUTED
      !-----------------------------------------
      i_forc_count = i_forc_count + 1
      DO i = 1, n_forc
         IF ( forc_swi(i) .EQ. 0 ) forc_arr(i) = forc_arr_old(i) 
     &                     + forc_coeff(i) * FLOAT(i_forc_count)
!        WRITE(numout,*) ' forc_arr   : ', forc_arr(i)
      END DO

      ! Recover arrays (dirty patch to remove in the end)
      IF ( forc_swi(1) .EQ. 0 ) zfsw = forc_arr(1)
      IF ( forc_swi(2) .EQ. 0 ) zflw = forc_arr(2)
      IF ( forc_swi(3) .EQ. 0 ) zpar = forc_arr(3)
      IF ( forc_swi(4) .EQ. 0 ) ztair= forc_arr(4)
      IF ( forc_swi(5) .EQ. 0 ) zpres= forc_arr(5)
      IF ( forc_swi(6) .EQ. 0 ) zqair= forc_arr(6)
      IF ( forc_swi(7) .EQ. 0 ) zwspd= forc_arr(7)
      IF ( forc_swi(8) .EQ. 0 ) zcld = forc_arr(8)
      IF ( forc_swi(9) .EQ. 0 ) zfoce= forc_arr(9)
      IF ( forc_swi(10).EQ. 0 ) zsfal= forc_arr(10)
      IF ( forc_swi(11).EQ. 0 ) THEN 
         zalbe = forc_arr(11)
         zfswn = ( 1. - zalbe ) * zfsw
      ENDIF
!     ! sensitivity test 
!     ztair = ztair + 2.
!     zsfal = zsfal + zsfal
!     zcld  = zcld + zcld / 10.
!
!-----------------------------------------------------------------------------!
! 4) Compute some of the forcing fields if required
!-----------------------------------------------------------------------------!
!
      !----------------------
      ! 4.1) FLW computation
      !----------------------
      i = 1 ; j = 1
      IF ( ( forc_swi(2) .GT. 0 ) .AND. ( forc_swi(2) .LT. 99 ) ) THEN
         ze = zqair / ( 1. - zqair ) * zpres / 100.     ! wvp in hPa
     &      / ( 0.622 - zqair / ( 1. - zqair ) )
         zta4 = ztair * ztair * ztair * ztair
      ENDIF

      ! Efimova (61) and Key et al (96)
      !---------------------------------
      IF ( forc_swi(2) == 1 ) THEN
         zflw = emig * stefan * zta4 * 
     &          ( 0.746 + 0.0066 * ze ) + ( 1 + 0.26 * zcld )
      ENDIF

      ! Berliand 
      !----------
      IF ( forc_swi(2) == 2 ) THEN
         alat    = asin(covrai(i,j))/radian
         clat    = (90.0-alat)/10.0
         indx    = 1+int(clat)
         zflw = stefan*zta4*
     &            (1.0-(0.39-0.05*sqrt(ze))*(1.0-budyko(indx)*
     &              zcld*zcld))
         zcorr_fac_lw = 0.1
         zflw = zflw + zflw * zcorr_fac_lw
      ENDIF ! forc_swi(2)

      !-------------------------
      ! 4.2) Albedo computation
      !-------------------------
      ! Shine 
      !-------
      IF ( forc_swi(11) == 1 ) THEN
         CALL shine(tfsn, tfsg, t_su(i,j), ht_i(i,j), ht_s(i,j), 
     &              alb_c, alb_o)
         zalbe =  ( ( 1. - zcld ) * alb_c + zcld * alb_o )
      ENDIF ! 

      ! Albedo from observed values , ISPOL 
      !-------------------------------------
      IF ( forc_swi(11) == 2 ) THEN
         IF ( ht_s(i,j) .LE. 0.0 ) THEN
            zalbe = 0.50 !bare ice albedo
         ELSE
            IF ( t_su(i,j) .GE. 273.15 ) THEN
               zalbe = 0.65 ! melting snow albedo
            ELSE
               zalbe = 0.80 ! dry snow albedo
            ENDIF
         ENDIF
         alb_c = zalbe
         alb_o = zalbe + 0.06
         zalbe =  ( ( 1. - zcld ) * alb_c + zcld * alb_o )
      ENDIF ! 

      ! net SW flux in case SW is read in a file
      !------------------------------------------
      IF ( ( forc_swi(11) .GT. 0 ) .AND. ( forc_swi(11) .LT. 99 ) .AND. 
     &     ( forc_swi(1) .EQ. 0 ) ) zfswn = zfsw * zalbe

      !----------------------
      ! 4.3) FSW computation
      !----------------------
      IF ( ( forc_swi(1) .GT. 0 ) .AND. ( forc_swi(1) .LT. 99 ) ) THEN
         zeps0 = 1.d-13
         dpi   = 2*pi
         indaet = 1
         ijour  = INT(xjour)
         dec    = pdecli(indaet,ijour) * radian
         sdec   = sin(dec)
         cdec   = cos(dec)
         DO j = js1 , js2
            DO i = is1(j) , is2(j)
            ! geometric factors
            !-------------------
            ! covrai = sinus of latitude
            slat   = covrai(i,j)
            zps    = slat*sdec
            zpc    = COS(ASIN(slat))*cdec
            zljour = ACOS(-SIGN(one,zps)
     &             * MIN(one,SIGN(one,zps)*(zps/zpc)))
            dws    = (2.0*zljour)/REAL(nintsr)
            zlmidi = ASIN(( zps +zpc ))/radian
            zalcnq = 0.0
            DO k = 1, nintsr
              ws(k)     = zljour-(REAL(k)-0.5)*dws
              zmue(k)   = MAX(zero,zps+zpc*COS(ws(k)))
              zalcnp(k) = 0.05/(1.1*zmue(k)**1.4+0.15)
              zalcnq    = zalcnp(k)
            END DO
            zalcnq = zalcnq/MAX(2.0*zljour,zeps0)
            zmudum = 0.4

            ! Irradiance 
            !------------
            ze = zqair / ( 1. - zqair ) * zpres      ! wvp in Pa
     &         / ( 0.622 - zqair / ( 1. - zqair ) )
            ztc = ztair - 273.15
            zesw = 611.*EXP(17.269*ztc/(ztc+237.3)) ! Pa
            ztdew = ze / zesw ! no units

            !-------
            ! Shine 
            !-------
            IF ( forc_swi(1) == 1 ) THEN
               !-----------------
               ! DAILY time step
               !-----------------
               IF ( ddtb .EQ. 86400.0 ) THEN
                  frsdrg = 0.0
                  frsdfg = 0.0
                  frsdro = 0.0
                  frsdfo = 0.0
                  opt_dept = 16.297 ! to put in the namelist
                  DO k = 1 , nintsr ! integrate over the whole day
                     frsdrg = frsdrg+dws*      ! clear 
     &                       (zsolar*zmue(k)*zmue(k)*(1.0-alb_c))/
     &                     (1.2*zmue(k)+(1.0+zmue(k))*ze*1.0e-05+0.0455)
                     frsdfg  = frsdfg+dws*     ! overcast
     &                     ( ( 53.5 + 1274.5*zmue(k) ) * SQRT(zmue(k) )
     &                     * ( 1.0 - 0.996*alb_o ) )  
     &                     / (1.0+0.139*(1.0-0.9435*alb_o)*opt_dept)
                  END DO
                  ! net solar heat flux (1-a)FSW
                  zfswn = ( ( 1.0 - zcld ) * frsdrg + zcld * frsdfg ) / 
     &                    dpi
               ENDIF ! ddtb

               !--------------------
               ! SUBDAILY time step
               !--------------------
               IF ( ddtb .LT. 86400.00 ) THEN
!                  WRITE(numout,*) ' Shine formula, subdaily time step '
!                  WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '

                   ! Solar angle
                   ijour  = INT(xjour)
                   zdecl  = pdecli(indaet,ijour) * radian
                   zsdec  = SIN(zdecl)
                   zcdec  = COS(zdecl)
                   zslat  = covrai(i,j) ! covrai=sinus of latitude
                   zps    = zslat*zsdec
                   zpc    = COS(ASIN(slat))*cdec
                   WRITE(numout,*) ' ijour : ', ijour
                   WRITE(numout,*) ' xjour : ', xjour
                   WRITE(numout,*) ' zdecl : ', zdecl
                   WRITE(numout,*) ' zsdec : ', zsdec
                   WRITE(numout,*) ' zcdec : ', zcdec
                   WRITE(numout,*) ' zslat : ', zslat
                   WRITE(numout,*) ' zps   : ', zps
                   WRITE(numout,*) ' zpc   : ', zpc
                   ihour = ( xjour - FLOAT(ijour) ) * 24.
                   WRITE(numout,*) ' ihour : ', ihour

                   ! Hour angle
                   zljour = ACOS(-SIGN(one,zps)
     &                    * MIN(one,SIGN(one,zps)*(zps/zpc)))
                   zhourang = - pi + 2.*pi / 24. * FLOAT(ihour)
                   WRITE(numout,*) ' zljour   : ', zljour
                   WRITE(numout,*) ' zhourang : ', zhourang
                   ! Cosine of solar angle
                   zcosz = 0.0
                   IF ( ( zhourang .GT. -zljour ) .AND. 
     &                  ( zhourang .LT. zljour ) ) 
     &                zcosz = MAX(0.0,zps+zpc*COS(zhourang))
                   zcosz2 = zcosz * zcosz
                   WRITE(numout,*) ' zcosz : ', zcosz
                   WRITE(numout,*) ' zcosz2: ', zcosz2
                   ! Irradiance
                   zqsr_clear = zsolar * zcosz2 
     &                        / ( 1.2 * zcosz 
     &                          + ( 1.0 + zcosz ) * ze * 1.0e-5 
     &                          + 0.0455 )
                   opt_dept = 16.297 ! to put in the namelist
                   zqsr_cloud = ( 53.5 + 1274.5 * zcosz ) * 
     &                          SQRT( zcosz ) / ( 1.0 + 0.139 * 
     &                          ( 1.0 - 0.9435 * zalbe ) * opt_dept )
                   zfsw = ( 1. - zcld ) * zqsr_clear 
     &                  + zcld          * zqsr_cloud
                   zfswn = ( 1. - zcld ) * ( 1. - alb_c ) * zqsr_clear
     &                   + zcld * ( 1. - alb_o ) * zqsr_cloud
                   WRITE(numout,*) ' zqsr_clear : ', zqsr_clear
                   WRITE(numout,*) ' zqsr_cloud : ', zqsr_cloud

               ENDIF ! ddtb

            ENDIF ! Shine

            !----------
            ! Zillmann
            !----------
            IF ( forc_swi(1) == 2 ) THEN
               IF ( ddtb .EQ. 86400.0 ) THEN
               frsdtg = 0.0
               frsdto = 0.0
               DO k = 1 , nintsr
                  albo   = (1.0-zcld)*zalcnp(k)+zcld*zalbe
                  frsdtg = frsdtg+dws*(1.0-zalbe)*
     &                     (zsolar*zmue(k)*zmue(k))/
     &                     ((zmue(k)+2.7)*ze*1.0e-05+1.085*zmue(k)+0.10)
               END DO
               ! rewrite this next line correctly
               zfswn      = 0.9*min(one,(1-.62*zcld+.0019*zcld))*
     &                      frsdtg/dpi
               ENDIF
            ENDIF ! 
            END DO
         END DO
      ENDIF ! forc_swi(1)

      IF ( forc_swi(3) == 1 ) THEN
         !-----------------
         ! 4.2 PAR formula
         !-----------------
         zpar = 0.43
      ENDIF
!
!-----------------------------------------------------------------------------!
! 5) Case of prescribed values
!-----------------------------------------------------------------------------!
!
      IF ( forc_swi(1) .EQ. 99 ) zfsw  = forc_val(1) * forc_uni(1)
      IF ( forc_swi(2) .EQ. 99 ) zflw  = forc_val(2) * forc_uni(2)
      IF ( forc_swi(3) .EQ. 99 ) zpar  = forc_val(3) * forc_uni(3)
      IF ( forc_swi(4) .EQ. 99 ) ztair = forc_val(4) * forc_uni(4)
      IF ( forc_swi(5) .EQ. 99 ) zpres = forc_val(5) * forc_uni(5)
      IF ( forc_swi(6) .EQ. 99 ) zqair = forc_val(6) * forc_uni(6)
      IF ( forc_swi(7) .EQ. 99 ) zwspd = forc_val(7) * forc_uni(7)
      IF ( forc_swi(8) .EQ. 99 ) zcld  = forc_val(8) * forc_uni(8)
      IF ( forc_swi(9) .EQ. 99 ) zfoce = forc_val(9) * forc_uni(9)
      IF ( forc_swi(10).EQ. 99 ) zsfal = forc_val(10) * forc_uni(10)
      IF ( forc_swi(11).EQ. 99 ) THEN
         zalbe = forc_val(11) * forc_uni(11)
         zfswn = ( 1. - zalbe ) * zfsw
      ENDIF
!
!-----------------------------------------------------------------------------!
! 6) Temporary dirty plug
!-----------------------------------------------------------------------------!
!
      i = 1
      j = 1

      fsolg(i,j)  = zfswn 
      ratbqg(i,j) = zflw
      ! par is not yet included, should be
      tabq(i,j)   = ztair
      psbq(i,j)   = zpres
      qabq(i,j)   = zqair
      vabq(i,j)   = zwspd
      cloud(i,j)  = zcld
      oce_flx     = zfoce
      hnpbq(i,j)  = zsfal
      albg(i,j)   = zalbe
      ze = zqair / ( 1. - zqair ) * zpres      ! wvp in Pa
     &   / ( 0.622 - zqair / ( 1. - zqair ) )
      ztc = ztair - 273.15
      zesw = 611.*EXP(17.269*ztc/(ztc+237.3)) ! Pa
      ztdew = ze / zesw ! no units
      tdew(i,j)   = ztdew

      IF ( ln_write_forc ) THEN
         WRITE(numout,*) ' Forcing fields ... '
         WRITE(numout,*)
         WRITE(numout,*) ' fsolg       : ', fsolg(i,j)
         WRITE(numout,*) ' ratbqg      : ', ratbqg(i,j)
         WRITE(numout,*) ' tabq        : ', tabq(i,j)
         WRITE(numout,*) ' psbq        : ', psbq(i,j)
         WRITE(numout,*) ' qabq        : ', qabq(i,j)
         WRITE(numout,*) ' vabq        : ', vabq(i,j)
         WRITE(numout,*) ' cloud       : ', cloud(i,j)
         WRITE(numout,*) ' oce_flx     : ', oce_flx
         WRITE(numout,*) ' hnpbq       : ', hnpbq(i,j)
         WRITE(numout,*) ' albg        : ', albg(i,j)
         WRITE(numout,*)
      ENDIF

!------------------------------------------------------------------------------
! end of forcing 
      RETURN
      END
