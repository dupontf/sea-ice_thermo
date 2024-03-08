	subroutine flx(hi2,thsfc,tair,qair,fsens,flat,qsfc,
     &			zchu1,zchu2,uair2,zref)
! inputs :hi2(m),thsfc(K),qsfc(kg/kg),tair(K),qair(kg/kg),uair(m/s),zref(m)
! outputs:fsens (W/m^2), flat (W/m^2)
!--- Computes turbulent fluxes of sensible and latent heat (z=10m).
!--- Based on Andreas, E.L., 1987: A Theory for the Scalar Roughness 
!--- and the Scalar Transfer Coefficients over Snow and Sea Ice,
!--- Boundary Layer Meteorology, v.38, 159-184.
! 
        implicit real(8) (a-h,o-z)
!       implicit real (a-h,o-z)        
        integer ii
!--- Constants:
!real aw,cp,cw2,epsi,epss,epsw,gd,grav,hilead,pci,one
!       real pcs,pi,qb2,qi,rd,rv,qs,qs0,rho2,rhoice,rhosf,rhosm,rhow
!       real sigma,t0, ai2, bi, ts2, ta2
!
!real alphae,alphah2,alrs,CEN,CHN2,cte,ctf
!       real kn,L2,lnztz0,lnzqz0,nu,psih2,psim2,psiq,uair
!       real qair,qsfc,Ri2,Rstar,sqrtCD,tair,thsfc,uair,usave,ustar
!       real Ri2,Rstar,sqrtCD,usave,ustar
!       real(8) tair,thsfc,uair2,hi2,qair,fsens,flat,qsfc
!       real(8) zchu1,zchu2,CH2,CE2,CD,CDN
!       real u_data,x2,z0,zeta,zq2,zref,zt
        dimension b0t(3),b1t(3),b2t(3), b0q(3),b1q(3),b2q(3)
        real(8) lv,kn,L2,lnztz0,lnzqz0,nu
        
!!!	include 'const.cmn'
!
!--- Values of coefficients for polynomials, Table I, p.177
	data b0t/1.250, 0.149, 0.317/
	data b1t/0.000,-0.550,-0.565/
	data b2t/0.000, 0.000,-0.183/
	data b0q/1.610, 0.351, 0.396/
	data b1q/0.000,-0.628,-0.512/
	data b2q/0.000, 0.000,-0.180/
!--03/08/2001      
        uair=uair2
!--- Constants
	kn=0.4
	nu=1.461e-5
	alphah2=1.0
	alphae=1.0
	ai2   = 21.8746
	bi    =-265.5
        grav  = 9.81          !  gravitational acceleration (m/s2)
        lv    = 2.501e+6      !  latent heat of vaporization (J/kg)
        rho2  = 1.275         !  density of dry air (kg/m3)
        cp    = 1005.         !  heat capacity of dry air (J/kg.K)
        zref  = 10.           !  (m)
!--- cst rajoutees le 02/08/2001---
        one   = 1.0
        pi    = 4.0 * tan(one)
!--- get surface sphum
        ts2=thsfc-273.16
        qsfc  = 0.622*6.11/1013.*exp(min(ai2*ts2/(ts2-bi),10.))
!
         
! FD	if(uair.le.0.) then
! FD	  fsens=0.
! FD	  flat=0.
! FD          ctf=0.
! FD          cte=0.
! FD	  return

! FD        else if (uair.lt.0.5) then
        uair = max( uair, 0.1) ! FD capping

        if (uair.lt.0.5) then
!TEA save value (will linearly interpolate down from 0.5 m/s)
          usave = uair
          uair = 0.5
        else
          usave = -1.0
	endif
!
	Ri2=grav*zref*(tair-thsfc)/(tair*uair*uair)
!
!--- Thickness dependent roughness length (Guest & Davidson 1991 JGR)
	if(hi2.lt.0.01) then
	   z0=8.0e-4
	elseif(hi2.ge.0.01.and.hi2.le.0.10) then
	   z0=4.5e-4
	elseif(hi2.gt.0.10.and.hi2.le.0.30) then 
	   z0=2.4e-3
	elseif(hi2.gt.0.30.and.hi2.le.2.00) then
	   z0=1.3e-3
	elseif(hi2.gt.2.) then
	   z0=2.0e-3
	endif
!
!--- Using reference height of 10 m, get u* from uair and z0 using (1)
!
	ustar = uair * kn / log(zref/z0)
!
!--- Get roughness Reynolds number R* = u* z0 / v  (p.163)
!
	Rstar = ustar * z0 / nu
!
!--- Get ln(zT/z0) and ln(zQ/z0) from (53) and Table I  Andreas (p.177)
!
	if(Rstar.le.0.135) then
	   ii=1
	elseif(Rstar.gt.0.135.and.Rstar.lt.2.5) then
	   ii=2
	elseif(Rstar.ge.2.5) then
	   ii=3
	endif

	alrs=log(Rstar)
	lnztz0 = b0t(ii) + b1t(ii)*alrs + b2t(ii)*alrs*alrs
	lnzqz0 = b0q(ii) + b1q(ii)*alrs + b2q(ii)*alrs*alrs
	zt = z0*exp(lnztz0)     ! Roughness length for temperature
	zq2 = z0*exp(lnzqz0)     ! Roughness length for q

!
!--- Get neutral drag coefficient CD from (10)  Andreas p. 162
!
	CDN = kn**2 / (log(zref/z0))**2
!
!--- Get neutral CH2 and CE2 from (11), (12) using alphah2=alphae=1.0
!(p.162)
!
	sqrtCD = CDN**0.5
	CHN2 = alphah2*kn*sqrtCD / (kn/sqrtCD - lnztz0)
	CEN = alphae*kn*sqrtCD / (kn/sqrtCD - lnzqz0)
!
!--- Correct for stability dependence
	if(tair.eq.thsfc) then
!	 Neutral case
	   CD=CDN
	   CH2=CHN2
	   CE2=CEN
	   go to 9
	endif
!
!--- Obukhov length (Stull 1988, eq.9.7.5k, p.386)
        L2 = ustar*tair*uair / (kn*grav*(tair-thsfc))
        zeta = zref/L2
!
!--- Stability functions (Liu et al 1979, p.1723) become negative for
!--- low wind speeds, now using stability parameters f(RiB) from
!--- Louis (1979) and Louis et al. (1981)
!--- (Note A&M and Stull use different sign conventions for psim2)
!--- Monin-Obukhov similarity theory applied to the surface layer
!--- works only when the winds are NOT calm and u* is not zero
!
	if(Ri2.gt.0.) then        ! Stable surface layer case
           psim2 = -7.*zeta
           psih2 = psim2
           psiq = psih2
	else if(Ri2.lt.0.) then ! Unstable surface layer case
           x2=(1.-16.*zeta)**0.25
           psim2 = 2.*log((1.+x2)/2.) + log((1.+x2*x2)/2.) - 
     &             2.*tan(x2)+ pi/2.
           psih2 = 2.*log((1.+x2*x2)/2.)
           psiq = psih2
	endif

!
!--- Get transfer coefficients from Andreas & Murphy (2.14-2.15)
!
        CD = kn**2 / (log(zref/z0)-psim2)**2
        sqrtCD = CD**0.5
        CH2 = alphah2*kn*sqrtCD / (kn/sqrtCD - lnztz0-(psih2-psim2))
        CE2 = alphae*kn*sqrtCD / (kn/sqrtCD - lnzqz0-(psiq-psim2))
!---test 10 /08/2001
         CH2 =0.00175
         CE2 =0.00175 
!
!---  For ECMWF or NCEP analyses, convert 10 m winds to 2 m winds
!
    9   u_data = uair    ! Save wind speed at 10 m
!
     	fsens=rho2*cp*CH2*uair*(thsfc-tair)   !  This line used to be numbered
	flat= rho2*lv*CE2*uair*(qsfc-qair)    !  for goto 9
 100  format(8f15.9)
!
!---03/08/2001  calcul de zrchu----
! 
        zchu1=rho2*cp*CH2*uair
        zchu2=rho2*lv*CE2*uair
        uair = u_data    !  Return wind speed to 10 m value
	ctf=CH2
	cte=CE2
!
!TEA linearly interpolate
        if (usave.gt.0.) then
           fsens = fsens*(usave/uair)
           flat  = flat*(usave/uair)
           ctf = ctf*(usave/uair)
           cte = cte*(usave/uair)
        endif
!--- restrict flat to non-negative values
! test du 10/08/2001                  
!          flat  = max(flat,0.)
  	return
	end
