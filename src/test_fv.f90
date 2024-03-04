program test

  use var_thermo_vertical
      use ice_thermodynamic_fv


implicit none

logical :: ln_write =.true.


       OPEN(1,file='restart.dat',status='old',action='read')
       READ(1,*) nlice,nlsno,ith_cond
       READ(1,*) ti(1:nlice),ts(1:nlsno+1),tsu,tbo
       READ(1,*) si(1:nlice)
       READ(1,*) hi,hs
       READ(1,*) snowfall,dwnlw,tsu,tair,qair,uair,swrad,oceflx
       READ(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
       CLOSE(1)

      call init_sigma_ice_FV
      call ice_thermo(3600.d0)

end
