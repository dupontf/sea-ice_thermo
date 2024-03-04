      PROGRAM biolim

      ! This is the BIOLIM code
      ! composed of the halo-thermodynamic model of Vancoppenolle et al. (JGR 2007)
      ! and of a newly developed biogeochemical component

      ! (c) Martin Vancoppenolle (UCL-ASTR, 2002-2008)

      INCLUDE 'type.com'
      INCLUDE 'const.com'
      INCLUDE 'para.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'dynami.com'
      INCLUDE 'forcing.com'
 
!------------------------------------------------------------------------------
!     1) CREATES THE OUTPUT FILES
!------------------------------------------------------------------------------


      OPEN(numout,file='run.out.txt',status='replace')        ! output text file
      OPEN(503   ,file='cons_bio.out.txt',status='replace')   ! output conservation file
      OPEN(85    ,file='forcing.out.txt',status='replace')
      OPEN(86    ,file='time.step', status='replace')

      WRITE(numout,*) ' ---------------------------------------------- '
      WRITE(numout,*) '         *** BIO-LIM sea ice model ***          '
      WRITE(numout,*) '  (c) UCL-ASTR, Louvain-La-Neuve Belgium        '
      WRITE(numout,*) '  Martin Vancoppenolle, 2008.                   '
      WRITE(numout,*) ' ---------------------------------------------- '
      WRITE(numout,*) 
      WRITE(numout,*) ' biolim :  '
      WRITE(numout,*) ' ~~~~~~~~~ '
      WRITE(numout,*) 
!
!------------------------------------------------------------------------------
!     2) INITIALIZING THE RUN    
!------------------------------------------------------------------------------
!
      nn99=0

      !-----------------------------------
      ! Physical constants and parameters
      !-----------------------------------
      CALL defcst(nn99)
      CALL defgrid  ! have to suppress this

      !-----------------------
      ! Do the initialization 
      !-----------------------
      CALL init_nc(n_s,n_i)

      ntrmax = int(dts(ks2)/ddtb)

 
      ! Permet de fixer le jour ou on veut commencer l''iteration
      ! premjour est fixer dans run.param
      premjour= nstart
! FD      tpstot = (nstart    - 1) * 86400
      tpstot = DBLE( nstart - 1 ) * ddtb

      ! nstart indique la premiere iteration qui doit etre effectuee
      ! numit indique le nombre d''iteration
      numit  = nstart    - 1 
!     nstart = numit + 1
      nlast  = numit + nitrun
      dtsec  = 86400.
      dtsd2 = 0.5*dts(ks2)

      IF ( ninfo .EQ. 0 ) ninfo = nlast + 1
      IF ( nsav .EQ. 0 )  nsav  = nlast + 1
 
      ! Le lstab est definit dans le run.param. Il definit
      ! la valeur de mixage qui dans le cas 1d sera fixee a 0

      mixage = MAX(MIN(lstab+lstab,2),-lstab)
! FD      numofday = nday1     - 1  
      numofday = INT( tpstot / 86400.d0)
      yeaday   = 365.0
      yrsec  = yeaday*86400.
!
!-----------------------------------------------------------------------
!     3) MAIN LOOP                                                     |
!-----------------------------------------------------------------------
!
      WRITE(numout,*)
      WRITE(numout,*) ' Loop begins, iteration : ', nstart
      WRITE(numout,*)
      WRITE(numout,*) '  nstart = ', nstart
      WRITE(numout,*) '  nend   = ', nend

      DO 300 nit= nstart, nend  

         !---------------------
         ! Time step variables
         !---------------------
         numit = numit + 1
         WRITE(numout,*) 
     &   '============================================================='
         WRITE(numout,*)
         WRITE(numout,*) '   *** Beginning the time step *** ', numit
         WRITE(numout,*)
         WRITE(numout,*) 
     &   '============================================================='
         WRITE(numout,*)
! FD         numofday = numofday + 1
         tpstot = tpstot + dtsd2
         numofday = INT( tpstot / 86400.d0 )
         ! The date corresponds to the middle of the time step.
         !   - iyear is the year, the first year is the year 1.
         !   - integer part of xjour is the day of the year(1-yeaday)
         !   - non-integer part of xjour is fraction of the day

         ! Number of days in the year
         ! and number of the year
         IF (numofday.eq.int(yeaday+1)) THEN
            nyear1   = nyear1   + 1
            yeaday = 365.0
            IF (mod(nyear1  ,4).eq.0) THEN
               yeaday = 366.0
            ENDIF
            yrsec  = yeaday*86400.
            numofday = 1
         ENDIF
            
         iyear  = 1+int(tpstot/yrsec)
! FD         xjour  = mod(tpstot,yrsec)/dtsec + 1.0
! FD         xjour  = real(numofday) + 0.5
! FD         ! xjour is not good
         xjour  = REAL(nit) * ddtb / 86400.

         y = (tpstot/yrsec)

         !       Sinon avance de 1jour
         d = (tpstot/yrsec)*yeaday + 1.0
         !	  d = (tpstot/yrsec)*yeaday
          
         year = int(y)
         !	   d = (y-year)*yeaday
         day = int(d)
         tpstot = tpstot + dtsd2
         ! ecriture du pas de temps tous les 30 jours 

         IF ((REAL(numit)/100.0).eq.(real(Int(real(numit)/100.0)))) THEN
           WRITE(86,*) numit
         ENDIF

         !---------
         ! Forcing
         !---------
         CALL forcing_nc(xjour)

         !--------------------
         ! Ice thermodynamics
         !--------------------
         WRITE(numout,*) ' n_i :', n_i
         WRITE(numout,*) ' n_s :', n_s
         CALL ice_th(ntrmax,n_i,n_s,numofday)
 
 300  CONTINUE ! End of the loop

!-----------------------------------------------------------------------
!     4) END OF THE RUN                                                |
!-----------------------------------------------------------------------
      ! close output files
      CLOSE(numout) 
      CLOSE(85)
      CLOSE(86)
      CLOSE(999) 
      CLOSE(24)
 
!-----------------------------------------------------------------------------!
! End of the main routine
!
      END
