      SUBROUTINE ice_bio_grid(kideb,kiut,ln_write)

! This routine creates biogeochemical vertical grid
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
     &  jk          , ! : index for ice layers
     &  jn          , ! : index for tracers
     &  numbio = 500  ! : reference number for bio.param

      REAL(8) ::  
     &  zh_bio1 = 0.1 , ! : thickness of the upper bio layer
     &  zh_bio2 = 0.1   !: thickness of the lower bio layer
! FD additions
      LOGICAL ln_write

!=============================================================================!

      IF ( ln_write ) THEN
         WRITE(numout,*) ' *** ice_bio_grid : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*)
         WRITE(numout,*) ' c_mod : ', c_mod
      ENDIF

      DO ji = kideb, kiut

      !-----------------
      ! Layer thickness
      !-----------------

      ! Uniform multi-layer case
      IF ( c_mod .EQ. 'ML' ) THEN
         DO layer = 1, nlay_bio
            deltaz_i_bio(layer) = ht_i_b(ji) / nlay_bio
         END DO
      ENDIF

      ! Deformed case (uppermost and lowermost have fixed thickness)
      IF ( c_mod .EQ. 'DL' ) THEN
        deltaz_i_bio(1) = MIN( ht_i_b(ji) / REAL(nlay_bio) , zh_bio1 )
        deltaz_i_bio(nlay_bio) = MIN( ht_i_b(ji) / REAL(nlay_bio) ,
     &                           zh_bio2 )
        DO layer = 2, nlay_bio - 1
           deltaz_i_bio(layer) = ( ht_i_b(ji) - deltaz_i_bio(1) - 
     &                             deltaz_i_bio(nlay_bio) ) / 
     &                           FLOAT( nlay_bio - 2 )
        END DO
      ENDIF

      ! Skeletal layer case
      IF ( c_mod .EQ. 'SL' ) THEN
         h_bio = 0.02
         deltaz_i_bio(nlay_bio) = MIN( h_bio, ht_i_b(ji) / 2. )
         DO layer = 1, nlay_bio - 1
            deltaz_i_bio(layer) = ( ht_i_b(ji) - h_bio ) / 
     &                            ( nlay_bio - 1 )
         END DO
      ENDIF

      ! Biologically-active layer case
      IF ( c_mod .EQ. 'BA' ) THEN
         i_bal = 1
         
         DO layer = 2, n_i
            IF ( ( e_i_b(layer-1) .LT. flu_bvtr ) .AND. 
     &           ( e_i_b(layer)   .GT. flu_bvtr ) ) i_bal = layer
         END DO ! layer

         WRITE(numout,*) ' e_i_b   : ', ( e_i_b(layer), 
     &                     layer = 1, n_i )
         WRITE(numout,*) ' flu_bvtr: ', flu_bvtr
         WRITE(numout,*) ' i_bal   : ', i_bal
         WRITE(numout,*)

         ! compute the location of BAL
         zdh = ht_i_b(ji) / REAL(n_i)
         zde = ( e_i_b(i_bal) - e_i_b(i_bal-1) )
         zm  = zde / zdh
         zdz = ( flu_bvtr - e_i_b(i_bal - 1 ) ) / zm
         WRITE(numout,*) ' zdh : ', zdh
         WRITE(numout,*) ' zde : ', zde
         WRITE(numout,*) ' zm  : ', zm 
         WRITE(numout,*) ' zdz : ', zdz
         WRITE(numout,*)
         IF ( i_bal .GT. 1 ) zbal = z_i_phy(i_bal-1) + zdz
         IF ( i_bal .EQ. 1 ) zbal = 0. ! location of the 5% contour
         h_bio = MAX( ht_i_b(ji) - zbal - 0.02, 0.02 )

         WRITE(numout,*) ' ht_i_b  : ', ht_i_b(ji)
         WRITE(numout,*) ' zbal  = ', zbal 
         WRITE(numout,*) ' h_bio = ', h_bio
         WRITE(numout,*)


         ! vertical coordinates
         deltaz_i_bio(nlay_bio) = h_bio
         DO layer = 1, nlay_bio - 1
            deltaz_i_bio(layer) = ( ht_i_b(ji) - h_bio ) / 
     &                            ( nlay_bio - 1 )
         END DO ! layer
      ENDIF ! c_mod .EQ. 'BA'

      !-------------
      ! Layer cotes
      !-------------
      z_i_bio(1)  = deltaz_i_bio(1) / 2.0
      DO layer = 2, nlay_bio
         z_i_bio(layer)  = z_i_bio(layer-1) + ( deltaz_i_bio(layer-1) +
     &                     deltaz_i_bio(layer) ) / 2.0
      END DO

      END DO ! ji

      IF ( ln_write ) THEN
         WRITE(numout,*) ' deltaz_i_bio : ', ( deltaz_i_bio(layer),  
     &                   layer = 1, nlay_bio)
         WRITE(numout,*) ' z_i_bio  : ', ( z_i_bio(layer),  
     &                   layer = 1, nlay_bio)

         WRITE(numout,*) ' end '
      ENDIF

!=============================================================================!
!-- End of ice_bio_ini --
 
      END
