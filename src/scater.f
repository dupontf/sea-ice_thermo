      subroutine scater(n,a,index,b)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c
      include 'type.com'
      REAL(8) :: b(n)
      INTEGER :: index(n)
      REAL(8) :: a(*)
c
      do 1 ji=1,n
      a(index(ji))=b(ji)
 1    continue
c
      return
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- Fin de la routine scater -
      end
