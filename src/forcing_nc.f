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
      real(8) :: zps

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
         WRITE(numout,*) ' oce_flx     : ', oce_flx
         WRITE(numout,*) ' hnpbq       : ', hnpbq(i,j)
         WRITE(numout,*) ' albg        : ', albg(i,j)
         WRITE(numout,*)
      ENDIF

!------------------------------------------------------------------------------
! end of forcing 
      RETURN
      END

      SUBROUTINE forcing_semtner(xjour)

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

      real, dimension(12) :: sw12=(/0.,0.,30.3,157.6,281.8,305.7,
     &216.5,143.3,58.9,0.6,0.,0./)
      real, dimension(12) :: lw12=(/165.6, 164.0, 164.0, 184.7, 240.4, 
     &               286.6, 304.1, 297.7, 262.7, 221.3, 178.3, 173.6/)
      real, dimension(12) :: sh12=(/18.8,12.1,11.5,4.6,-7.2,
     &-6.2,-4.8,-6.4,-2.7,0.2,8.9,12.6/)
      real, dimension(12) :: lh12=(/0.0,-0.3,-0.5,-1.4,-7.3,
     &-11.2,-10.2,-10.5,-6.2,-0.3,-0.2,-0.2/)
      real, dimension(12) :: al12=(/0.84,0.84,0.83,0.81,0.80,
     &0.78,0.64,0.69,0.84,0.84,0.84,0.84/)
      real, dimension(12) :: fo12=(/2.,2.,2.,2.,2.,2.,2.,2.,
     &2.,2.,2.,2./)
      real, dimension(12) :: pr12=(/1.,1.,1.,1.,1.,0.,0.,0.,
     &0.,1.,1.,1./)
      real, dimension(12) :: xmon=(/0.5,1.5,2.5,3.5,4.5,5.5,
     &                              6.5,7.5,8.5,9.5,10.5,11.5/)

      real xd
! FD addition
      LOGICAL ln_write_forc

! FD debug      ln_write_forc = .TRUE.
      ln_write_forc = .FALSE.

! FD      WRITE(numout,*) ' * forcing_nc : '
! FD      WRITE(numout,*) ' ~~~~~~~~~~~~~~ '

      xd=mod(xjour,360.)/30.
      call xcubic_per(xd, fsolg(1,1) ,  12, xmon, sw12, 12.)
      call xcubic_per(xd, ratbqg(1,1),  12, xmon, lw12, 12.)
      call xcubic_per(xd, fsbbq(1,1) ,  12, xmon, sh12, 12.)
      call xcubic_per(xd, ffltbq(1,1),  12, xmon, lh12, 12.)
      call xcubic_per(xd, albg(1,1)  ,  12, xmon, al12, 12.)
      call xcubic_per(xd, oce_flx    ,  12, xmon, fo12, 12.)
      call xcubic_per(xd, hnpbq(1,1) ,  12, xmon, pr12, 12.)
! FD corrections
! FD at time of first paper 1.05 factor correction on SW and LW
      fsolg(1,1) = MAX(fsolg(1,1)*1.03,0.) ! solar radiation
      fsolg(1,1) = ( 1. - albg(1,1) ) * fsolg(1,1) 
      hnpbq(1,1) = MAX(hnpbq(1,1) * 0.01*0.2,0.) ! precipitation in m/day
      oce_flx = oce_flx * 3.0 ! ocean flux
      ratbqg(1,1) = ratbqg(1,1)*1.03 ! long wave descending heat flux
      


      IF ( ln_write_forc ) THEN
         WRITE(numout,*) ' Forcing fields ... '
         WRITE(numout,*)
         WRITE(numout,*) ' fsol        : ', fsolg(1,1)
         WRITE(numout,*) ' lwh         : ', ratbqg(1,1)
         WRITE(numout,*) ' shf         : ', fsbbq(1,1)
         WRITE(numout,*) ' lhf         : ', ffltbq(1,1)
         WRITE(numout,*) ' oce_flx     : ', oce_flx
         WRITE(numout,*) ' hnpbq       : ', hnpbq(1,1)
         WRITE(numout,*) ' albg        : ', albg(1,1)
         WRITE(numout,*)
      ENDIF

!------------------------------------------------------------------------------
! end of forcing 
      RETURN
      END

      subroutine xcubic_per(x, fup, na, xa, fa, xperiodic)
      implicit none
      integer na
      real x, fup, fpup, fppup, xa(na), fa(na), xperiodic
! locals
      real capx,f00,f10,c0,c1,c2,c3,c4
      integer i,ix,ix1,ixp1,ixm1
      real, parameter :: pt5=0.5,one=1.0,two=2.0,six=6.0,ov6=one/six

            capx = (x-xa(1))/(xa(2)-xa(1))
            if (capx<0.0) capx=capx+xperiodic/(xa(2)-xa(1))
            ix   = int(capx) + 1
            capx = capx - REAL(ix-1)
            ix1  = mod(ix,na) + 1

            c3 = capx
            c1 = one - c3
            c0 = - ov6 * c1 * c3
            c2 = c0 * ( two - c3 )
            c4 = c0 * ( one + c3 )

            ixp1 = ix1
            ixm1 = mod(ix-2+na,na) + 1

            f00 = ( fa(ixp1) - fa(ix) )
     &           +( fa(ixm1) - fa(ix) )

            ixp1 = mod(ix1,na) + 1
            ixm1 = ix

            f10 = ( fa(ixp1) - fa(ix1) )
     &           +( fa(ixm1) - fa(ix1) )
!*
!*         final reconstruction
            fup   = c1 * fa(ix)  + c2 * f00 + c3 * fa(ix1)  + c4 * f10

      end

