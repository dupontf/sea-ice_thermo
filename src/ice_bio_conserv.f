      SUBROUTINE ice_bio_conserv(kideb,kiut,message,err,
     &                           zmb0,zmb1,zfb,zfsu,zdt)

! This routine tests conservation of the mass of tracers
! (c) Martin Vancoppenolle, June 2007
 
      INCLUDE 'type.com'
      INCLUDE 'para.com'

      CHARACTER(len=15) :: 
     &   message       ! : message indicating the name of the routine calling

      REAL(8) ::
     &   err

      REAL(8), DIMENSION(ntra_bio) ::
     &   zmb0 ,
     &   zmb1 ,
     &   zfb  ,
     &   zfsu

      REAL(8) ::
     &   zdt 
         
      INTEGER :: 
     &   ji          , ! : index for space
     &   jk          , ! : index for ice layers
     &   jn            ! : index for tracers

      REAL(8) ::
     &   zdm           ! : actual mass variation
     &   zdmf

!=============================================================================!

      WRITE(numout,*) ' ice_bio_conserv : '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~ '
      WRITE(numout,*) ' message     : ', message
      WRITE(numout,*) ' error max   : ', err
      WRITE(numout,*) ' kideb, kiut : ', kideb, kiut
      WRITE(numout,*) ' ddtb        : ', zdt

      DO jn = 1, ntra_bio
         WRITE(numout,*) ' jn : ', jn
         WRITE(numout,*) ' ntra_bio : ', ntra_bio
         WRITE(numout,*) ' mt_i_bio_init  : ', zmb0(jn)
         WRITE(numout,*) ' mt_i_bio_final : ', zmb1(jn)
!     END DO

!     DO jn = 1, ntra_bio
         zdm = ( zmb1(jn) - zmb0(jn) ) / zdt
         zdmf = zfb(jn) + zfsu(jn)

         WRITE(numout,*) ' Actual mass variation zdm       : ', zdm
         WRITE(numout,*) ' Mass variation from fluxes zdmf : ', zdmf
         WRITE(numout,*) ' Bio conserv error : ', ABS(zdm-zdmf)
         WRITE(503,*) zmb1(jn), ABS(zdm-zdmf)*zdt
         
         IF ( ABS ( zdm - zdmf ) .GT. err ) THEN
            WRITE(numout,*) ' Conservation error after ', message
            WRITE(numout,*) ' Error                           : ', 
     &                      ABS( zdm - zdmf )
            WRITE(numout,*) ' Actual mass variation zdm       : ', zdm
            WRITE(numout,*) ' Mass variation from fluxes zdmf : ', zdmf
            WRITE(numout,*)
            WRITE(numout,*) ' mt_i_bio_init   : ', zmb0(jn)
            WRITE(numout,*) ' mt_i_bio_final  : ', zmb1(jn)
!           WRITE(numout,*) ' c_i_bio   : ', ( c_i_bio(jn,layer), 
!    &                      layer = 1, nlay_bio )
!           WRITE(numout,*) ' cbu_i_bio : ', ( cbu_i_bio(jn,layer), 
!    &                      layer = 1, nlay_bio )

         ENDIF
      END DO

      RETURN

!=============================================================================!
!-- End of ice_bio_conserv --
 
      END
