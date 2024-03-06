      SUBROUTINE ice_phy_grid(kideb,kiut,nlay,zht_i,ln_write,snoorice)

! This routine creates thermodynamical vertical grid
! (c) Martin Vancoppenolle, May 2007
 
      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'
      INCLUDE 'bio.com'

      INTEGER :: 
     &  ji          , ! : index for space
     &  jk            ! : index for ice layers

      REAL(8) :: zht_i
      INTEGER :: kideb,kiut,nlay
      CHARACTER(len=3) :: 
     &   snoorice     ! : message indicating whether we have to compute ice or snow
! FD additions
      LOGICAL ln_write

!=============================================================================!

      IF ( ln_write ) THEN
         ji = 1
         WRITE(numout,*) ' *** ice_phy_grid : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*)
         WRITE(numout,*) ' ht_i_b : ', zht_i
         WRITE(numout,*) ' nlay   : ', nlay   
         WRITE(numout,*) ' which? : ', snoorice
      ENDIF

      DO ji = kideb, kiut

      !---------------------------
      ! layer thickness
      !---------------------------
      IF ( snoorice .EQ. 'ice' ) THEN
         DO layer = 1, nlay
            deltaz_i_phy(layer) = zht_i / REAL(nlay,8)
         END DO
      ENDIF

      IF ( snoorice .EQ. 'sno' ) THEN
         DO layer = 1, nlay
            deltaz_s_phy(layer) = zht_i / REAL(nlay,8)
         END DO
      ENDIF

      !---------------------------
      ! layer cotes
      !---------------------------
      IF ( snoorice .EQ. 'ice' ) THEN
         z_i_phy(1)  = deltaz_i_phy(1) / 2.0
         DO layer = 2, nlay
            z_i_phy(layer)  = z_i_phy(layer-1) + ( deltaz_i_phy(layer-1)
     &                      + deltaz_i_phy(layer) ) / 2.0
         END DO

         IF ( ln_write ) THEN
            WRITE(numout,*) ' deltaz_i_phy : ', ( deltaz_i_phy(layer),  
     &                      layer = 1, nlay )
            WRITE(numout,*) ' z_i_phy      : ', ( z_i_phy(layer),  
     &                      layer = 1, nlay )
            WRITE(numout,*) ' end '
         ENDIF

      ENDIF

      IF ( snoorice .EQ. 'sno' ) THEN
         z_s_phy(1)  = deltaz_s_phy(1) / 2.0
         DO layer = 2, nlay
            z_s_phy(layer)  = z_s_phy(layer-1) + ( deltaz_s_phy(layer-1)
     &                      + deltaz_s_phy(layer) ) / 2.0
         END DO

         IF ( ln_write ) THEN
            WRITE(numout,*) ' deltaz_s_phy : ', ( deltaz_s_phy(layer),  
     &                      layer = 1, nlay )
            WRITE(numout,*) ' z_s_phy      : ', ( z_s_phy(layer),  
     &                      layer = 1, nlay )
            WRITE(numout,*) ' end '
         ENDIF

      ENDIF

      END DO ! ji

!=============================================================================!
!-- End of ice_phy_grid --
 
      END
