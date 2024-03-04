program test

  use var_thermo_vertical
  use ice_thermodynamic_FE


implicit none

integer k

       OPEN(1,file='restart.0.dat',status='old',action='read')
       READ(1,*) nlice,nlsno
       READ(1,*) tiold(0:nlice+nlsno)
       READ(1,*) si(1:nlice)
       READ(1,*) hi,hs
       READ(1,*) snowfall,dwnlw,tsu,tair,qair,uair,swrad,oceflx,seasal
       READ(1,*) fac_transmi,swradab_i(1:nlice),swradab_s(1:nlsno)
       CLOSE(1)

 ni=nlice
 ns=nlice+nlsno

 do k=0,ni-1
    zi(k) = dble(k)/dble(nlice)
 enddo
   zi(ni) = 1.d0
 do k=ni+1,ns-1
    zi(k) = 1.d0 + dble(k-ni)/dble(nlsno)
 enddo
   zi(ns) = 2.d0

!---------------------------------------------------------------------
! other quantities
!---------------------------------------------------------------------

 do k=1,ns
   dzi(k)=zi(k)-zi(k-1)
 enddo
  ziold =  zi
 dziold = dzi

!---------------------------------------------------------------------
!---------------------------------------------------------------------
 call ice_thermo_FE(3600.d0)

end
