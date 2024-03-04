      SUBROUTINE ice_bio_column(kideb,kiut,zmt_i,zc_i,zdeltaz,ln_write)

! This routine sums the mass of tracers over the whole column
! (c) Martin Vancoppenolle, June 2007
 
      INCLUDE 'type.com'
      INCLUDE 'para.com'

      INTEGER :: 
     &  ji          , ! : index for space
     &  jk          , ! : index for ice layers
     &  jn            ! : index for tracers

      REAL(8), DIMENSION(ntra_bio) ::
     &   zmt_i
      REAL(8), DIMENSION(ntra_bio,nlay_bio) ::
     &   zc_i
      REAL(8), DIMENSION(nlay_bio) ::
     &   zdeltaz
! FD additions
      LOGICAL ln_write

!=============================================================================!

      IF ( ln_write ) THEN
         WRITE(numout,*) ' *** ice_bio_column : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~ '
      ENDIF

      DO jn = 1, ntra_bio
         zmt_i(jn) = 0.0
      END DO

      DO jn = 1, ntra_bio
         DO jk = 1, nlay_bio
            zmt_i(jn) = zmt_i(jn) + zc_i(jn,jk) *
     &      zdeltaz(jk)
         END DO
      END DO

      RETURN

!=============================================================================!
!-- End of ice_bio_column --
 
      END
