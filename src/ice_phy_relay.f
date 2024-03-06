      SUBROUTINE ice_phy_relay( nbot0 , nbot1 , ntop0 , ntop1 , !zl0,zl1,
     &                          hl0, hl1, ql0, ql1 )

      !!------------------------------------------------------------------
      !!                ***         ROUTINE ice_phy_relay     ***
      !!
      !! ** Purpose :
      !! This routine redistributes the scalar quantity ql0 - distributed on layers of
      !! lower boundary limit zl0 - onto new scalars ql1 on a new grid - described by 
      !! zl1.
      !!       rl01(i,j) represents linear weights of old layers into the new ones.
      !!       index of top layer --> ntop0/1
      !!       index of bottom layer --> nbot0/1
      !!
      !! ** Method  : Relayering
      !!           
      !! ** Steps
      !!
      !! ** Arguments
      !!
      !! ** Inputs / Outputs
      !!
      !! ** External
      !!
      !! ** References : Vancop. et al., 2007
      !!
      !! ** History  : 
      !!    (12-2002) Martin Vancop. First test
      !!    (06-2003) Martin Vancop. LIM1D
      !!    (06-2008) Martin Vancop. BIO-LIM
      !!    (09-2008) Martin Vancop. Explicit gravity drainage
      !!
      !!------------------------------------------------------------------

      INCLUDE 'para.com'
      INCLUDE 'type.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'dynami.com'
      INCLUDE 'moment.com'
      INCLUDE 'thermo.com'
 
      REAL(8), DIMENSION ( maxnlay ) ::
     &  ql0         !: old scalar

      REAL(8), DIMENSION ( maxnlay + 2 ) ::
     &  hl0         !: old layer thickness 

      REAL(8), DIMENSION ( maxnlay ) ::
     &  ql1      ,  !: new scalar
     &  hl1         !: old layer thickness 

      REAL(8), DIMENSION ( 0:maxnlay ) ::
     &  zl0      ,  !: old layer interfaces
     &  zl1         !: new layer interfaces

      REAL(8) ::
     &  rl01        !: relayering matrix

      INTEGER ::
     &  layer1   ,  !: first layer index
     &  layer2      !: second layer index
! FD additions
      LOGICAL ln_write

      zlimit = 1.0e-10   !: limiting factor to avoid divisions per 0
      ln_write = .FALSE.
      ql1(:) = 0.0
!
!==============================================================================!

        zl0(0)     = 0.0
        zl0(ntop0) = hl0(ntop0)
        DO layer0 = ntop0+1, nbot0
           zl0(layer0) = zl0(layer0-1) + hl0(layer0)
        END DO

        zl1(0)     = 0.0
        zl1(ntop1) = hl1(ntop1)
        DO layer1 = ntop1+1, nbot1
           zl1(layer1) = zl1(layer1-1) + hl1(layer1)
        END DO

      IF ( ln_write ) THEN
         WRITE(numout,*) ' ** ice_phy_relay : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) 
         WRITE(numout,*) ' Input... '
         WRITE(numout,*) ' ql0 : ', ( ql0(jk), jk = 1, nbot0 )
         WRITE(numout,*) ' nbot0 : ', nbot0 
         WRITE(numout,*) ' nbot1 : ', nbot1 
         WRITE(numout,*) ' ntop0 : ', ntop0
         WRITE(numout,*) ' ntop1 : ', ntop1
         WRITE(numout,*) ' zl0   : ', ( zl0(jk), jk = ntop0, nbot0 )
         WRITE(numout,*) ' zl1   : ', ( zl1(jk), jk = ntop1, nbot1 )
         WRITE(numout,*) ' hl0   : ', ( hl0(layer0), 
     &                                layer0 = ntop0, nbot0 )
         WRITE(numout,*) ' hl1   : ', ( hl1(layer1), 
     &                                layer1 = ntop1, nbot1 )
      ENDIF
!
!------------------------------------------------------------------------------|
!  1) Relayering procedure                                                     |
!------------------------------------------------------------------------------|
!
        !-------------
        ! new scalars
        !-------------
        IF (ln_write) WRITE(numout,*) ' Redistribution of the scalar '
! FD this version works
!        DO layer1 = ntop1, nbot1
!           ql1(layer1) = 0.0
!           DO layer0 = ntop0, nbot0
!     
!               rl01 = MAX( 0.0 , ( MIN(zl0(layer0), 
!     &         zl1(layer1)) - MAX(zl0(layer0-1), zl1(layer1-1) ) ) / 
!     &         MAX( hl0(layer0) , zlimit ) ) 
!              ql1(layer1) = ql1(layer1) + rl01 * ql0(layer0)
!              IF (ln_write) WRITE(numout,*) ' ql1 : ', ql1(layer1)
!              IF (ln_write) WRITE(numout,*) ' ql0 : ', ql0(layer0)
!              IF (ln_write) WRITE(numout,*) ' rl01: ', rl01
!           END DO
!        END DO
         ql1(ntop1:nbot1) = 0.0
         layer0 = ntop0
         layer1 = ntop1
         do while (layer0 <= nbot0 .and. layer1 <= nbot1)
               rl01 = MAX( 0.0 , ( MIN(zl0(layer0), 
     &         zl1(layer1)) - MAX(zl0(layer0-1), zl1(layer1-1) ) ) / 
     &         MAX( hl0(layer0) , zlimit ) ) 
              ql1(layer1) = ql1(layer1) + rl01 * ql0(layer0)
              IF (ln_write) WRITE(numout,*) ' ql1 : ', ql1(layer1)
              IF (ln_write) WRITE(numout,*) ' ql0 : ', ql0(layer0)
              IF (ln_write) WRITE(numout,*) ' rl01: ', rl01
            if (zl0(layer0) > zl1(layer1)) then
               layer1 = layer1 + 1
            else
               layer0 = layer0 + 1
            endif
        END DO
        IF (ln_write) WRITE(numout,*) ' ql1 : ', 
     &                ( ql1(jk), jk = ntop1, nbot1 )

        RETURN
!------------------------------------------------------------------------------!
      END SUBROUTINE
