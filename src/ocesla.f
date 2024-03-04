      subroutine ocesla
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  Programme simulant un ocean statique agissant uniquement comme un reservoir 
c  de chaleur. 
c  modif : 01/04/94
 
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'

c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|

c- 1) FLUX DE SURFACE
c--------------------
c- Prise en compte des Flux en surface pour les scalaires.
c- phiss <-deja.mult par deltaT/Dz, donc il ne faut multiplier phiss
c- par un coefficient
      do 30 ns=1,nsmax
         do 20 j=js1,js2
            do 10 i=is1(j),is2(j)
               scal(i,j,ks2,ns)=scal(i,j,ks2,ns) - phiss(i,j,ns)
C              scal(i,j,ks2,ns)=scal(i,j,ks2,ns)
C    &                         -dts(ks2)*unsdz(ks2)*phiss(i,j,ns)
 10         continue
 20      continue
 30   continue


c- 2) FLUX DE VOLUME
c--------------------
c- Prise en compte des flux de volume pour les scalaires.
c- phivs <-deja.mult par deltaT/Dz, donc il ne faut pas multiplier phivs
c- par un coefficient
      do 40 ns=1,nsmax 
	  do 50 j=1,jmax 
	    do 60 k=1,kmax
	       do 70 i=1,imax
	          scal(i,j,k,ns)= scal(i,j,k,ns) - phivs(i,j,k,ns)
 70            continue
 60        continue
 50     continue
 40   continue


c- 3) RAPPEL (VERSION EXPLICITE)
c--------------------------------
      do 80 ns=1,nsmax
       do 80 j=1,jmax
	   do 80 i=1,imax
	     do 80 k=1,kmax
              fflx = rappes(i,j,k,ns)*(scalr(i,j,k,ns)-scal(i,j,k,ns))
              scal(i,j,k,ns)= scal(i,j,k,ns) + fflx
 80   continue


      return
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- fin de la routine ocesla -
      end
