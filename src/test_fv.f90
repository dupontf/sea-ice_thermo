program test

  use var_thermo_vertical
      use ice_thermodynamic_fv


implicit none

logical :: ln_write =.true.
integer k
real(8) :: dtice



       OPEN(1,file='restart.0.dat',status='old',action='read')
       READ(1,*) nlice,nlsno,ith_cond,dtice
call allocate_thermo_1d
!       READ(1,*) ti(1:nlice),ts(1:nlsno+1),tsu,tbo
       READ(1,*) ti(1:nlice),ts(1:nlsno),tsu,tbo
       READ(1,*) si(1:nlice)
       READ(1,*) hi,hs
       READ(1,*) snowfall,dwnlw,tsu,tair,qair,uair,swrad,oceflx,pres
       READ(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
       READ(1,*) fsens,flat
       CLOSE(1)

      call init_sigma_ice_FV
      call ice_thermo(dtice)

end
