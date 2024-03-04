      SUBROUTINE defgrid
 
        !!------------------------------------------------------------------
        !!                ***  ROUTINE defgrid      ***
        !! ** Purpose :
        !!           Defines dummy boundary constants
        !!           Residual from 3D version of CLIO
        !!           Should be removed
        !! ** Method  :
        !!           Definitions !!! 
        !!
        !! ** Arguments :
        !!
        !! ** Inputs / Ouputs : (global commons)
        !!
        !! ** External : 
        !!
        !! ** References : Vancoppenolle et al., JGR 2007
        !!
        !! ** History :
        !!       (1) CLIO, Goosse and Fichefet, JGR, 1999.
        !!       (2) LIM-1D, Vancoppenolle et al., JGR, 2007.
        !!       (3) BIO-LIM, Martin Vancoppenolle, 2008
        !!
        !!------------------------------------------------------------------
        !! * Arguments

        INCLUDE 'type.com'
        INCLUDE 'const.com'
        INCLUDE 'para.com'
        INCLUDE 'bloc.com'
        INCLUDE 'ice.com'
        INCLUDE 'reper.com'
        INCLUDE 'dynami.com'
 
!------------------------------------------------------------------------------
!  1) Local and global variables
!------------------------------------------------------------------------------

      WRITE(numout,*) ' * defgrid : '
      WRITE(numout,*) ' ~~~~~~~~~~~ '
      WRITE(numout,*) 
 
      do 11 j=1,jmax
       is1(j) = imax
       is2(j) = 1
 11   continue
 
      do 12 k=1,kmax
        z(k)  = 0.0
        unsdz(k)  = 0.0
 12   continue
      z(kmax+1) = 0.0
 
      do 20 j=1,jmax
       do 20 i=1,imax
        hu(i,j) = 0.0
 20   continue
 
      do 21 ns=1,nsmax
       do 21 j=1,jmax
        do 21 i=1,imax
         phiss(i,j,ns)  = 0.0
 21   continue
 
      do 30 k=1,kmax
       do 30 j=1,jmax
        do 30 i=1,imax
         q(i,j,k) = 0.0
          tms(i,j,k) = 1.0
          tmu(i,j,k) = 1.0
 30   continue
 
      do 31 k=1,kmax+1
       do 31 j=1,jmax
        do 31 i=1,imax
         w(i,j,k) = 0.0
 31   continue
 
c--Initialisation des champs scalaires et des forcings, de l'energie
c  cinetique turbulente
c  useful for ocesla
      do 32 ns=1,nsmax
       do 32 k=1,kmax
        do 32 j=1,jmax
         do 32 i=1,imax
          rappes(i,j,k,ns) = 0.0
          scalr(i,j,k,ns)=0.0
          phivs(i,j,k,ns) = 0.0
 32   continue
 
!------------------------------------------------------------------------------
!  2 ) Grid and resolution constants
!------------------------------------------------------------------------------
 
        !--maximum size of the basin
        ims1 = 1
        ims2 = 1
 
        js1 = 1
        js2 = 1
 
        ks1 = 1
        ks2 = kmax
 
        jcl1 = 1
        jcl2 = 0
 
      !--covrai and coriolis factor
      do 540 j=1,jmax
        do 540 i=1,imax
          covrai(i,j)=sin(zlatz*radian)
          fs2cor(i,j)=omega*covrai(i,j)
 540  continue
 
      RETURN
!------------------------------------------------------------------------------
! end of subroutine
      END
