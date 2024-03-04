      SUBROUTINE ice_sal_conserv(kideb,kiut,message,err,
     &                           zmb0,zmb1,zfb,zfsu,zdt)

!=============================================================================!

! This routine tests conservation of the mass of salt
! (c) Martin Vancoppenolle, September 2008
 
      INCLUDE 'type.com'
      INCLUDE 'para.com'

      CHARACTER(len=15) :: 
     &   message       ! : message indicating the name of the routine calling

      REAL(8) ::
     &   err

      REAL(8) ::
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

      WRITE(numout,*) ' ice_sal_conserv : '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~ '
      WRITE(numout,*) ' message     : ', message
      WRITE(numout,*) ' error       : ', err
      WRITE(numout,*) ' kideb, kiut : ', kideb, kiut
      WRITE(numout,*) ' ddtb        : ', zdt
      WRITE(numout,*) ' ms_i_init  : ', zmb0
      WRITE(numout,*) ' ms_i_final : ', zmb1

      zdm = ( zmb1 - zmb0 ) / zdt
      zdmf = zfb + zfsu

      WRITE(numout,*) ' Actual mass variation zdm       : ', zdm
      WRITE(numout,*) ' Mass variation from fluxes zdmf : ', zdmf
      WRITE(numout,*) ' Sal conserv error : ', ABS(zdm-zdmf)

      WRITE(503,*) zmb1, ABS(zdm-zdmf)*zdt
         
      IF ( ABS ( zdm - zdmf ) .GT. err ) THEN
         WRITE(numout,*) ' Conservation error after ', message
         WRITE(numout,*) ' Error                           : ', 
     &                   ABS( zdm - zdmf )
         WRITE(numout,*) ' Actual mass variation zdm       : ', zdm
         WRITE(numout,*) ' Mass variation from fluxes zdmf : ', zdmf
         WRITE(numout,*)
         WRITE(numout,*) ' ms_i_init       : ', zmb0
         WRITE(numout,*) ' ms_i_final      : ', zmb1
!           WRITE(numout,*) ' s_i_b     : ', ( s_i_b(ji,layer), 
!    &                      layer = 1, nlay_i )
      ENDIF
      RETURN

!=============================================================================!
!-- End of ice_sal_conserv --
      END
