	subroutine flx(hi2,thsfc,tair,qair,fsens,flat,qsfc,
     &			zchu1,zchu2,uair2,zref)
c inputs :hi2(m),thsfc(K),qsfc(kg/kg),tair(K),qair(kg/kg),uair(m/s),zref(m)
c outputs:fsens (W/m^2), flat (W/m^2)
c--- Computes turbulent fluxes of sensible and latent heat (z=10m).
c--- Based on Andreas, E.L., 1987: A Theory for the Scalar Roughness 
c--- and the Scalar Transfer Coefficients over Snow and Sea Ice,
c--- Boundary Layer Meteorology, v.38, 159-184.
c 
c       implicit double precision (a-h,o-z)
        implicit real (a-h,o-z)        
        integer ii
c--- Constants:
creal aw,cp,cw2,epsi,epss,epsw,gd,grav,hilead,pci,one
c       real pcs,pi,qb2,qi,rd,rv,qs,qs0,rho2,rhoice,rhosf,rhosm,rhow
c       real sigma,t0, ai2, bi, ts2, ta2
c
creal alphae,alphah2,alrs,CEN,CHN2,cte,ctf
c       real kn,L2,lnztz0,lnzqz0,nu,psih2,psim2,psiq,uair
c       real qair,qsfc,Ri2,Rstar,sqrtCD,tair,thsfc,uair,usave,ustar
c       real Ri2,Rstar,sqrtCD,usave,ustar
        double precision tair,thsfc,uair2,hi2,qair,fsens,flat,qsfc
        double precision zchu1,zchu2,CH2,CE2,CD,CDN
c       real u_data,x2,z0,zeta,zq2,zref,zt
        dimension b0t(3),b1t(3),b2t(3), b0q(3),b1q(3),b2q(3)
        real lv,kn,L2,lnztz0,lnzqz0,nu
        
ccc	include 'const.cmn'
c
c--- Values of coefficients for polynomials, Table I, p.177
	data b0t/1.250, 0.149, 0.317/
	data b1t/0.000,-0.550,-0.565/
	data b2t/0.000, 0.000,-0.183/
	data b0q/1.610, 0.351, 0.396/
	data b1q/0.000,-0.628,-0.512/
	data b2q/0.000, 0.000,-0.180/
c--03/08/2001      
        uair=uair2
c--- Constants
	kn=0.4
	nu=1.461e-5
	alphah2=1.0
	alphae=1.0
	ai2   = 21.8746   
	bi    =-265.5
        grav  = 9.81            !  gravitational acceleration (m/s2)
        lv    = 2.501e+6        !  latent heat of vaporization (J/kg)
        rho2  = 1.275           !  density of dry air (kg/m3)
        cp    = 1005.           !  heat capacity of dry air (J/kg.K)
        zref  = 10.             !  (m)
c--- cst rajoutees le 02/08/2001---
        one   = 1.d0
        pi    = 4.0 * atan(one)
c--- get surface sphum
        ts2=thsfc-273.16
        qsfc   =  0.622*6.11/1013.*exp(min(ai2*ts2/(ts2-bi),10.))
c
         
	if(uair.le.0.) then
	  fsens=0.
	  flat=0.
          ctf=0.
          cte=0.
	  return

        else if (uair.lt.0.5) then
cTEA save value (will linearly interpolate down from 0.5 m/s)
          usave = uair
          uair = 0.5
        else
          usave = -1.
	endif
c
	Ri2=grav*zref*(tair-thsfc)/(tair*uair*uair)
c
c--- Thickness dependent roughness length (Guest & Davidson 1991 JGR)
	if(hi2.lt.0.01) then
	   z0=8.0e-4
	elseif(hi2.ge.0.01.and.hi2.le.0.10) then
	   z0=4.5e-4
	elseif(hi2.gt.0.10.and.hi2.le.0.30) then 
	   z0=2.4e-3
	elseif(hi2.gt.0.30.and.hi2.le.2.00) then
	   z0=1.3e-3
	elseif(hi2.gt.2.00) then
	   z0=2.0e-3
	endif
c
c--- Using reference height of 10 m, get u* from uair and z0 using (1)
c
	ustar = uair * kn / alog(zref/z0)
c
c--- Get roughness Reynolds number R* = u* z0 / v  (p.163)
c
	Rstar = ustar * z0 / nu
c
c--- Get ln(zT/z0) and ln(zQ/z0) from (53) and Table I  Andreas (p.177)
c
	if(Rstar.le.0.135) then
	   ii=1
	elseif(Rstar.gt.0.135.and.Rstar.lt.2.5) then
	   ii=2
	elseif(Rstar.ge.2.5) then
	   ii=3
	endif

	alrs=alog(Rstar)
	lnztz0 = b0t(ii) + b1t(ii)*alrs + b2t(ii)*alrs*alrs
	lnzqz0 = b0q(ii) + b1q(ii)*alrs + b2q(ii)*alrs*alrs
	zt = z0*exp(lnztz0)     ! Roughness length for temperature
	zq2 = z0*exp(lnzqz0)     ! Roughness length for q

c
c--- Get neutral drag coefficient CD from (10)  Andreas p. 162
c
	CDN = kn**2 / (alog(zref/z0))**2
c
c--- Get neutral CH2 and CE2 from (11), (12) using alphah2=alphae=1.0
c(p.162)
c
	sqrtCD = CDN**0.5
	CHN2 = alphah2*kn*sqrtCD / (kn/sqrtCD - lnztz0)
	CEN = alphae*kn*sqrtCD / (kn/sqrtCD - lnzqz0)
c
c--- Correct for stability dependence
	if(tair.eq.thsfc) then
c	 Neutral case
	   CD=CDN
	   CH2=CHN2
	   CE2=CEN
	   go to 9
	endif
c
c--- Obukhov length (Stull 1988, eq.9.7.5k, p.386)
        L2 = ustar*tair*uair / (kn*grav*(tair-thsfc))
        zeta = zref/L2
c
c--- Stability functions (Liu et al 1979, p.1723) become negative for
c--- low wind speeds, now using stability parameters f(RiB) from
c--- Louis (1979) and Louis et al. (1981)
c--- (Note A&M and Stull use different sign conventions for psim2)
c--- Monin-Obukhov similarity theory applied to the surface layer
c--- works only when the winds are NOT calm and u* is not zero
c
	if(Ri2.gt.0.) then        ! Stable surface layer case
           psim2 = -7.*zeta
           psih2 = psim2
           psiq = psih2
	else if(Ri2.lt.0.) then ! Unstable surface layer case
           x2=(1.-16.*zeta)**0.25
           psim2 = 2.*alog((1.+x2)/2.) + alog((1.+x2*x2)/2.) - 
     &             2.*atan(x2)+ pi/2.
           psih2 = 2.*alog((1.+x2*x2)/2.)
           psiq = psih2
	endif

c
c--- Get transfer coefficients from Andreas & Murphy (2.14-2.15)
c
        CD = kn**2 / (alog(zref/z0)-psim2)**2
        sqrtCD = CD**0.5
        CH2 = alphah2*kn*sqrtCD / (kn/sqrtCD - lnztz0-(psih2-psim2))
        CE2 = alphae*kn*sqrtCD / (kn/sqrtCD - lnzqz0-(psiq-psim2))
c---test 10 /08/2001
         CH2 =0.00175
         CE2 =0.00175   
c
c---  For ECMWF or NCEP analyses, convert 10 m winds to 2 m winds
c
    9   u_data = uair    ! Save wind speed at 10 m
c
     	fsens=rho2*cp*CH2*uair*(thsfc-tair)   !  This line used to be numbered
	flat= rho2*lv*CE2*uair*(qsfc-qair)    !  for goto 9
 100  format(8f15.9)
c
c---03/08/2001  calcul de zrchu----
c 
        zchu1=rho2*cp*CH2*uair
        zchu2=rho2*lv*CE2*uair
        uair = u_data    !  Return wind speed to 10 m value
	ctf=CH2
	cte=CE2
c
cTEA linearly interpolate
        if (usave.gt.0.) then
           fsens = fsens*(usave/uair)
           flat  = flat*(usave/uair)
           ctf = ctf*(usave/uair)
           cte = cte*(usave/uair)
        endif
c--- restrict flat to non-negative values
c test du 10/08/2001                  
C          flat  = max(flat,0.)
  	return
	end
