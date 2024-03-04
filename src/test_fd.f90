program test

  use var_thermo_vertical
      use ice_thermo_lim


implicit none

logical :: ln_write =.true.


       OPEN(1,file='restart.dat',status='old',action='read')
       READ(1,*) nlice,nlsno
       READ(1,*) ti(1:nlice),ts(1:nlsno),tsu,tbo
       READ(1,*) si(1:nlice)
       READ(1,*) hi,hs
       READ(1,*) dzi(1:nlice),dzs(1:nlsno)
       READ(1,*) snowfall,dwnlw,tsu,tair,qair,uair,swrad,oceflx
       READ(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
       CLOSE(1)

      call ice_thermo_diff (3600.d0,ln_write,6)
      sinew=si ! default, no change
      call ice_thermo_dh   (3600.d0,ln_write,6)
      call ice_thermo_remap(3600.d0,ln_write,6)

end
