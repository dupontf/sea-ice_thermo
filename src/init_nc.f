      SUBROUTINE init_nc(nlay_s,nlay_i)

      ! CLIO, UCL-ASTR, 2002 (H. Goosse)
      ! BIO-LIM1D, 2008 (M. Vancoppenolle)

      ! This routine initializes ice variables
      
!------------------------------------------------------------------------------!

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'dynami.com'
      INCLUDE 'moment.com'
      INCLUDE 'thermo.com'

      CHARACTER(len=10) :: 
     &   filenc='init.nc'

      REAL(4)        zini(1)  ! forcing field dummy array
      REAL(4)        zzi(maxnlay), zti(maxnlay), zsi(maxnlay) ! forcing field dummy array
      REAL(4)        zzs(maxnlay), zts(maxnlay) ! forcing dummy array for snow temperature profile
      INTEGER        nzs, nzi

      REAL(8), DIMENSION (maxnlay) :: 
     &   zsh_i0              ,   !: old ice salt content (ppt.m-2)
     &   zsh_i1                  !: new ice salt content (ppt.m-2)

      REAL(8), DIMENSION (maxnlay) :: 
     &   zqh_i0              ,   !: old ice heat content (J.m-2)
     &   zqh_i1                  !: new ice heat content (J.m-2)

      REAL(8), DIMENSION (maxnlay+2) ::
     &   zthick0                 !: thickness of old layers

      WRITE(numout,*) ' * init_nc : '
      WRITE(numout,*) ' ~~~~~~~~~~~ '
      WRITE(numout,*) 
      WRITE(numout,*) ' Variables values prescribed initially '

      WRITE(numout,*) ' nlay_s : ', nlay_s
      WRITE(numout,*) ' nlay_i : ', nlay_i

      zeps = 1.0e-10

!
!------------------------------------------------------------------------------|
!  1) Case of a NETCDF file
!------------------------------------------------------------------------------|
!
      i = 1
      j = 1
      !-----------
      ! Open file
      !-----------
      CALL CF_OPEN  (filenc,id) ! open forcing file

      !-----------
      ! Read data
      !-----------
! FD read snow
!      CALL CF_READDIM ( filenc, 'nlay_i', nzi )
!      CALL CF_READDIM ( filenc, 'nlay_s', nzs )

      zn = 0.0
      CALL CF_READ1D ( filenc, 'nlay_s', 1, 1, zini )
      zn = REAL(zini(1))
      nzs = NINT(zn)

      DO layer = 1, nzs
         CALL CF_READ1D ( filenc, 'z_s', layer, 1, zini)
         zzs(layer) = zini(1)
         CALL CF_READ1D ( filenc, 't_s', layer, 1, zini)
! FD         t_s(i,j,layer) = REAL(zini(1))
         zts(layer) = zini(1)
! FD debug
      write(*,*) 'zts',layer,zzs(layer),zts(layer)
      END DO

      zn = 0.0
      CALL CF_READ1D ( filenc, 'nlay_i', 1, 1, zini )
      zn = REAL(zini(1))
      nzi = NINT(zn)

      zj_d = 0.0
      CALL CF_READ1D ( filenc, 'j_d', 1, 1, zini )
      zj_d = REAL(zini(1))

      CALL CF_READ1D ( filenc, 'h_s', 1, 1, zini)
      ht_s(i,j) = REAL(zini(1))
      ! forced
      ht_s(i,j) = 0.05

      CALL CF_READ1D ( filenc, 'h_i', 1, 1, zini)
      ht_i(i,j) = REAL(zini(1))

      CALL CF_READ1D ( filenc, 't_su', 1, 1, zini)
      t_su(i,j) = REAL(zini(1))

      DO layer = 1, nzi
         CALL CF_READ1D ( filenc, 'z_i', layer, 1, zini)
         zzi(layer) = zini(1)
         CALL CF_READ1D ( filenc, 't_i', layer, 1, zini)
         zti(layer) = zini(1)
         CALL CF_READ1D ( filenc, 's_i', layer, 1, zini)
         zsi(layer) = zini(1)
      END DO

      t_bo(i,j)   =  273.15-tmut*oce_sal 

      !----------------------
      ! Interpolate profiles
      !----------------------

            !------
            ! Snow
            !------
            WRITE(numout,*) ' i = ', i, ' j = ', j
            CALL ice_phy_grid(1,1,nlay_s,ht_s(i,j), .FALSE., "sno" )   ! compute the physical grid

            ntop0  = 1
            nbot0  = nzs
            ntop1  = 1
            nbot1  = nlay_s

            !--- Temperature
            zm0(0) = 0.0         ! layer interfaces cotes
            DO layer = 1, nbot0-1
               zm0(layer) = REAL( zzs(layer) + zzs(layer+1) ) / 2.
            END DO
            zm0(nbot0) = ht_s(i,j)

            DO layer = 1, nbot0  !layer thickness and salt content
               zthick0(layer) = zm0(layer)-zm0(layer-1)
               zsh_i0(layer) = zthick0(layer)*REAL(zts(layer)) !volume of snow
            END DO

            CALL ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , 
     &                       zthick0, deltaz_s_phy, zsh_i0 , zsh_i1 )

            DO layer = 1, nlay_s
               t_s(i,j,layer) = zsh_i1(layer) /
     &                                     deltaz_s_phy(layer)
! FD debug
      write(*,*) 't_s',layer,t_s(i,j,layer)
            END DO

            !-----
            ! Ice
            !-----
            CALL ice_phy_grid(1,1,nlay_i,ht_i(i,j), .FALSE., "ice" )   ! compute the physical grid

            ! grid indexes for redistribution
            ntop0  = 1 
! FD            nbot0  = INT(zn)
            nbot0  = nzi
            ntop1  = 1
            nbot1  = nlay_i

            !--- Salinity
            zm0(0) = 0.0         ! layer interfaces cotes
            DO layer = 1, nbot0-1
               zm0(layer) = REAL( zzi(layer) + zzi(layer+1) ) / 2.
            END DO
            zm0(nbot0) = ht_i(i,j)

            DO layer = 1, nbot0  !layer thickness and salt content
               zthick0(layer) = zm0(layer)-zm0(layer-1)
               zsh_i0(layer) = rhog*zthick0(layer)*REAL(zsi(layer)) !mass of salt
            END DO

            CALL ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , 
     &                       zthick0, deltaz_i_phy, zsh_i0 , zsh_i1 )

            DO layer = 1, nlay_i
               s_i(i,j,layer) = zsh_i1(layer) / 
     &         (rhog*deltaz_i_phy(layer))
            END DO

            DO layer = 1, nbot0  ! heat content
               ztmelts    =   - tmut * REAL(zsi(layer)) + tpw 
               zq         = rhog * ( cpg * ( ztmelts -REAL(zti(layer)) )
     &                      + xlgm*( 1.0 - (ztmelts-tpw) / 
     &                        MIN( (REAL( zti(layer) ) -tpw),-zeps) )  
     &                      - cpw      * ( ztmelts-tpw  ) ) 
               zqh_i0(layer) = zq * zthick0(layer)
           END DO
           CALL ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , 
     &                       zthick0, deltaz_i_phy, zqh_i0 , zqh_i1 )

           DO layer = 1, nlay_i
              zq = zqh_i1(layer) / MAX( deltaz_i_phy(layer) , zeps )
              ztmelts = -tmut*s_i(i,j,layer) + tpw
              aaa = cpg
              bbb = (cpw-cpg)*(ztmelts-tpw) + zq / rhog - xlgm
              ccc = xlgm * (ztmelts-tpw)
              discrim = SQRT( bbb*bbb - 4.0*aaa*ccc )
              t_i(i,j,layer) = tpw + (- bbb - discrim) / ( 2.0*aaa )
           END DO
!
!------------------------------------------------------------------------------|
!  1.3) setting some fluxes to zero.                                           |
!------------------------------------------------------------------------------|

! check!!!
            i = 1
            j = 1
            firg(i,j)   = 0.0
            fcsg(i,j)   = 0.0
            fleg(i,j)   = 0.0
            fsbbq(i,j)  = 0.0
            qstobq(i,j) = 0.0

            WRITE(numout,*) ' Surface temperature t_su : ', t_su(i,j)
            WRITE(numout,*) ' Basal   temperature t_bo : ', t_bo(i,j)
            WRITE(numout,*) ' Snow temperatures   t_s  : ', 
     &                      ( t_s(i,j,jk), jk = 1, nlay_s )
            WRITE(numout,*) ' Ice temperatures    t_i  : ', 
     &                      ( t_i(i,j,jk), jk = 1, nlay_i )
            WRITE(numout,*) ' Ice salinities      s_i  : ', 
     &                      ( s_i(i,j,jk), jk = 1, nlay_i )
            WRITE(numout,*) ' Ice thickness       ht_i : ', ht_i(i,j)
            WRITE(numout,*) ' Snow depth          ht_s : ', ht_s(i,j)

!------------------------------------------------------------------------------|
!- end of init.f
      RETURN
      END
