c
c     ********************
      function pdecli(i,j)
c     ********************
c
c		THIS FUNCTION COMPUTES THE SOLAR DECLINATION FOR THE DAY
c		j (IN DECIMAL DEGREES).
c		j = DAY OF THE YEAR (j=1 ON JANUARY 1)
c		i = -1, 0, 1 FOR ODD, NORMAL AND LEAP YEARS,
c		    RESPECTIVELY
c
c---
c Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
c---
      include 'type.com'
      include 'const.com'
      xj     = j
      a0     = 0.39507671
      a1     = 22.85684301
      a2     = -0.38637317
      a3     = 0.15096535
      a4     = -0.00961411
      b1     = -4.29692073
      b2     = 0.05702074
      b3     = -0.09028607
      b4     = 0.00592797
      pi     = 3.1415927
      if (i) 1,2,3
 1    xj     = xj-0.5
      goto 4
 3    xj     = xj-1.
      goto 4
 2    xj     = j
 4    p      = pi*(2.0*xj-367.0)
      p      = p/yeaday
      p2     = p*2.0
      p3     = p*3.0
      p4     = p*4.0
      dcl1   = a1*cos(p)+a2*cos(p2)+a3*cos(p3)+a4*cos(p4)
      dcl2   = b1*sin(p)+b2*sin(p2)+b3*sin(p3)+b4*sin(p4)
      pdecli = a0+dcl1+dcl2
c
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Fin de la routine pdecli -
      end
