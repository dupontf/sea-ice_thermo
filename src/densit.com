c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  bloc "densit.com" : 'include' in the routines etat and outave
c  modif : 25/02/99
 
      common / etatloc / gravit,
     &  cstrho(0:10), cfb1z(kmax), cfb1z4(kmax+1), bref(kmax),
     &  dztau(kmax), zwtau(1+kmax),
     &  rho0dz(kmax),
     &  cfm2up(kmax), cfm2dw(kmax), rhozdz(0:1,kmax),
     &  dzsdz(kmax,kmax)
c
c--end of file "densit.com"
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
