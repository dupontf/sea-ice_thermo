c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  fichier "const.com", incorpore par instruction 'include' dans les programmes
c   (et les routines des programmes) :
c      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN, TRSBATH.
c  modif : 06/02/98
 
c--blocs common :

      REAL(8) :: zero, one, epsil, cstmin, cstmax, untour, pi, radian,
     &  yeaday, cpo, omega, rterre, unsrt, gpes, rho0, svrdrp
     &  iyea_num
 
      common / cstfix /
     &  zero, one, epsil, cstmin, cstmax, untour, pi, radian,
     &  yeaday, cpo, omega, rterre, unsrt, gpes, rho0, svrdrp
     &  iyea_num
c--fin du fichier "const.com"
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
