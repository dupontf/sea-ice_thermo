      subroutine gather(n,a,b,index)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c
      include 'type.com'
      REAL(8) :: a(n)
      INTEGER :: index(n)
      REAL(8) :: b(*)
c
      do 1 ji=1,n
      a(ji)=b(index(ji))
 1    continue
c
      return
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- Fin de la routine gather -
      end
