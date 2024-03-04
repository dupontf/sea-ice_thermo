!
!		COMMONS FOR FORCING
!		====================
!
!
! forc_nam      : name of the forcing fields
! ts_forc       : time step of the forcing
! forc_swi      : switch for computation of the forcing
! forc_uni      : unit conversion factor for a given forcing field
! forc_val      : prescribed value of a given forcing field if any

      CHARACTER(len=8), DIMENSION(n_forc) ::         
     &   forc_nam

      INTEGER ::
     &   n_fofr,
     &   i_forc_day,
     &   i_forc_count,
     &   i_forc,
     &   n_forc_min,
     &   n_forc_max,
     &   n0_forc,
     &   n1_forc

      INTEGER, DIMENSION(n_forc) ::
     &   forc_swi

      REAL(8)                    ::
     &   ts_forc

      REAL(8), DIMENSION(n_forc) ::
     &   forc_uni,
     &   forc_val,
     &   forc_arr_old,
     &   forc_arr_new,
     &   forc_coeff,
     &   forc_arr

      COMMON /forcing/
     &   forc_nam,
     &   n_fofr,
     &   n0_forc,
     &   n1_forc,
     &   i_forc_day,
     &   i_forc_count,
     &   i_forc,
     &   n_forc_min,
     &   n_forc_max,
     &   forc_swi,
     &   ts_forc,
     &   forc_uni,
     &   forc_val, 
     &   forc_arr_old,
     &   forc_arr_new,
     &   forc_coeff,
     &   forc_arr

