      SUBROUTINE ice_sal_column( kideb,kiut,zmt_i,zs_i,zdeltaz,nlay_i,
     &                           ln_write )

      !-----------------------------------------------------------------------!
      ! This routine sums the mass of salt over the whole sea ice column
      ! (c) Martin Vancoppenolle, September 2008
      !-----------------------------------------------------------------------!
      INCLUDE 'type.com'
      INCLUDE 'para.com'

      INTEGER :: 
     &  ji          , ! : index for space
     &  jk            ! : index for ice layers

      REAL(8) ::
     &   zmt_i        ! : total mass of salt
      REAL(8), DIMENSION(maxnlay) ::
     &   zs_i         ! : bulk salinity
      REAL(8), DIMENSION(maxnlay) ::
     &   zdeltaz      ! : thickness of the layers
! FD additions
      LOGICAL ln_write

      IF ( ln_write ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' *** ice_sal_column : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) ' kideb : ', kideb
         WRITE(numout,*) ' kiut  : ', kiut 
         WRITE(numout,*) ' zmt_i : ', zmt_i
         WRITE(numout,*) ' zs_i  : ', ( zs_i(jk), jk = 1, nlay_i )
         WRITE(numout,*) ' zdeltaz:', ( zdeltaz(jk), jk = 1, nlay_i )
      ENDIF

      zmt_i = 0.0
      DO jk = 1, nlay_i
         zmt_i = zmt_i + zs_i(jk) * zdeltaz(jk)
      END DO

      IF ( ln_write ) THEN
         WRITE(numout,*) ' zmt_i : ', zmt_i
      ENDIF


      !=============================================================================!
      !-- End of ice_sal_column --
      RETURN

      END
