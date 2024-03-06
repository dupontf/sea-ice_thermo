!
!  bloc "ice.com" : 'include' in the routines linked with the ice
!  modif : 25/09/98

! tfsn      Melting point temperature of the snow
! tfsg      Melting point temperature of the ice
! xkn       Conductivity of the snow
! xkg       Conductivity of the ice
! rcpn      Density times specific heat for the snow
! rcpg      Density times specific heat for the ice
! rhog      Density of the ice
! rhon      Density of the snow
! emig      Emissivity of the ice
! sglace    Salinity of the ice
! hmelt     Maximum melting at the bottom
! acrit(2)  Minimum fraction for leads
! hgcrit(2) Ice thickness for lateral accretion
! hgmin     Ice thickness corr. to max. energy stored in brine pocket
! hndif     Computation of temp. in snow or not 
! hgdif     Computation of temp. in ice or not 
! hglim     Minimum ice thickness
! amax      Maximum lead fraction
! uscomi    =1.0/(1.0-amax)
! beta      Numerical caracteritic of the scheme for diffusion in ice
! ddtb      Time step for ice thermodynamics (s)
! swiqst    Energy stored in brine pocket or not 
! parlat    Percentage of energy used for lateral ablation
! hakspl    Slope of distr. for Hakkinen-Mellor's lateral melting
! hibspl    Slope of distribution for Hibler's lateral melting
! exld      Exponent for leads-closure rate
! hakdif    Coefficient for diffusions of ice and snow
! hth       Threshold thickness for comp. of eq. thermal conductivity
! hnzst     Thickness of the surf. layer in temp. computation
! parsub    Switch for snow sublimation or not
! cnscg     ratio  rcpn/rcpg
! nbits     Number of time steps in Newton -Raphson procedure
! stefan    Stefan-Boltzman constant
! xsn       Sublimation heat for the snow
! vkarmn    von Karman constant
! cevap     Latent heat of evaporation of water
! zemise    Emissivity of water
! rhoesn    1/rhon
! firg      IR flux over the ice (only used for outputs)
! fcsg      Sensible heat flux over the ice (only used for outputs)
! fleg      Latent heat flux over the ice (only used for outputs)
! ts        Surface temperature of the ice
! tfu       Melting point temperature of sea water
! hnbq      Snow thickness
! hgbq      Ice thickness
! hgbqp     Ice production/melting
! albq      Leads fraction
! qstobq    Energy stored in the brine pockets
! fbbq      Heat flux at the ice base
! tbq       Temperature inside the ice/snow layer
! dmnbq     Variation of snow mass
! dmgbq     Variation of ice mass
! qlbq      heat balance of the lead (or of the open ocean)
! qcmbq     Energy needed to bring the ocean surface layer until its freezing 
!	    point (at a factor 2)
! thcm      part of the solar energy used in the lead heat budget
! fstrbq    Solar flux transmitted trough the ice
! ffltbq    Array linked with the max heat contained in brine pockets (?)
! fscmbq    Linked with the solar flux below the ice (?)
! fsbbq     Also linked with the solar flux below the ice (?)
! qfvbq     Array used to store energy in case of toral lateral ablation (?)
! xzo       rugosity of the ice (no more used)
! dmgwi     Variation of the mass of snow ice
! psbq      Surface air pressure
! tabq      Surface air temperature
! qabq      Surface air humidity
! vabq      Surface wind velocity
! hnplbq    Snow precipitation
! fevabq    Evaporation flux
! fsolcn    Solar flux at the ocean surface
! fsolg     Solar flux at the ice surface
! flecn     Latent heat flux at the ocean surface
! fcscn     Sensible heat flux at the ocean surface
! tenagx    Wind stress at the ice surface (x)
! tenagy    Wind stress at the ice surface (y)
! albg      03/08/2001 albedo obtenu du forcage simip2
! albege    Albedo of the snow or ice (only for outputs)
! tairox    Wind stress at the ocean surface (x)
! tairoy    Wind stress at the ocean surface (y)
! ratbqg    Longwave downward radiation flux over the ice
! ratbqo    Longwave downward radiation flux over the ocean
! cloud     Cloud fraction
! tdew      Air relative humidity
! albecn    Albedo of the ocean (only for outputs)
! tauc      Cloud optical depth
! runoff    river runoff
! sdvt      u*^2/(Stress/density)
! fcm1      Solar flux at the ocean surface
! fcm2      Non-solar flux at the ocean surface
! fwat      Freshwater flux (change of definition between the routines) 
! reslum    Relative absorption of solar radiation in each ocean level
! thcon_i_swi conductivity formula switch
!
!--COMMON blocs :
!
!------------------------------------------------------------------------------
      real(8)
     &  deltaz_i_phy(maxnlay) ,       !: thicknesses of the physical ice layers
     &  z_i_phy(maxnlay)      ,       !: cotes of the physical ice layers
     &  deltaz_s_phy(maxnlay) ,       !: thicknesses of the physical snow layers
     &  z_s_phy(maxnlay)              !: cotes of the physical snow layers

      COMMON /ice_grid/
     &  deltaz_i_phy,z_i_phy,deltaz_s_phy,z_s_phy

      real(8)
     &  tfsn,tfsg,xkn,xkg,rcpn,rcpg,rhog,rhon,
     &  emig,sglace,hmelt,acrit(2),hgcrit(2),hgmin,hndif,
     &  hgdif,hglim,amax,uscomi,beta,ddtb,swiqst,parlat,
     &  hakspl,hibspl,exld,hakdif,hth,hnzst,parsub,cnscg,nbits

      COMMON / ice_constants /
     &  tfsn,tfsg,xkn,xkg,rcpn,rcpg,rhog,rhon,
     &  emig,sglace,hmelt,acrit,hgcrit,hgmin,hndif,
     &  hgdif,hglim,amax,uscomi,beta,ddtb,swiqst,parlat,
     &  hakspl,hibspl,exld,hakdif,hth,hnzst,parsub,cnscg,nbits

      real(8)
     &  stefan,xsn,vkarmn,cevap,zemise,rhoesn
     
      COMMON / fluxsf /
     &  stefan,xsn,vkarmn,cevap,zemise,rhoesn
     
      real(8)
     &  firg(imax,jmax),fcsg(imax,jmax),fleg(imax,jmax)

      COMMON / comdia /
     &  firg,fcsg,fleg

      real(8)
! FD     &  ts(imax,jmax),tfu(imax,jmax),hnbq(imax,jmax),
     &  tfu(imax,jmax),hnbq(imax,jmax),
     &  hgbq(imax,jmax),hgbqp(imax,jmax),albq(imax,jmax),
     &  qstobq(imax,jmax),fbbq(imax,jmax),
     &  dmnbq(imax,jmax),dmgbq(imax,jmax),
     &  qlbq(imax,jmax),qcmbq(imax,jmax),thcm(imax,jmax),
     &  fstrbq(imax,jmax),ffltbq(imax,jmax),fscmbq(imax,jmax),
     &  fsbbq(imax,jmax),qfvbq(imax,jmax),xzo(imax,jmax),
     &  dmgwi(imax,jmax),total(imax,jmax)
c
      COMMON / comban /
! FD     &  ts,tfu,hnbq,
     &  tfu,hnbq,
     &  hgbq,hgbqp,albq,
     &  qstobq,fbbq,
     &  dmnbq,dmgbq,
     &  qlbq,qcmbq,thcm,
     &  fstrbq,ffltbq,fscmbq,
     &  fsbbq,qfvbq,xzo,
     &  dmgwi,total
c
      real(8)
     &  psbq(imax,jmax),tabq(imax,jmax),
     &  qabq(imax,jmax),vabq(imax,jmax),
     &  hnplbq(imax,jmax),fevabq(imax,jmax),fsolcn(imax,jmax),
     &  hnpbq(imax,jmax), 
     &  fsolg(imax,jmax),flecn(imax,jmax),fcscn(imax,jmax),
     &  tenagx(imax,jmax),tenagy(imax,jmax),albg(imax,jmax),
     &  albege(imax,jmax),tairox(imax,jmax),tairoy(imax,jmax),
     &  ratbqg(imax,jmax),ratbqo(imax,jmax),cloud(imax,jmax),
     &  tdew(imax,jmax),
     &  albecn(imax,jmax),tauc(imax,jmax),runoff(imax,jmax),
     &  sdvt(imax,jmax),fsolg2(imax,jmax)

      COMMON / comfor /
     &  psbq,tabq,
     &  qabq,vabq,
     &  hnplbq,fevabq,fsolcn,
     &  hnpbq, 
     &  fsolg,flecn,fcscn,
     &  tenagx,tenagy,albg,
     &  albege,tairox,tairoy,
     &  ratbqg,ratbqo,cloud,
     &  tdew,
     &  albecn,tauc,runoff,
     &  sdvt,fsolg2

      real(8)
     &  fcm1(imax,jmax),fcm2(imax,jmax),
     &  fwat(imax,jmax),
     &  reslum(imax,jmax,0:kmax+1)

      COMMON / comca /
     &  fcm1,fcm2,fwat,reslum

c global characteristics of the ice pack
      real(8)
     &  t_i(imax,jmax,maxnlay),t_s(imax,jmax,maxnlay),t_su(imax,jmax),
     &  t_bo(imax,jmax),ht_s(imax,jmax),ht_i(imax,jmax),
     &  s_i(imax,jmax,maxnlay)

      COMMON / ice_global /
     &  t_i,t_s,t_su,t_bo,ht_s,ht_i,s_i

      real(8)
     &  fc_int

      COMMON / heat_fluxes /
     &  fc_int

      integer :: layer,layer_a,nconv,numofday
      integer :: modul0, modul1a, modul2, modul3, modul4
      integer :: modul1b,modul5,nbot0,nbot1,ntop0,ntop1

c     Number of layers in the ice and snow
      integer :: nlayi0, nlayi1, nlays0, nlays1
      integer :: nlay_i, nlay_s
      integer :: thcon_i_swi

      COMMON/layers/ n_i, n_s, thcon_i_swi

      real(8) 
     &    zsim(0:maxnlay),tempsim(maxnlay),salsim(maxnlay),
     &             tempint(maxnlay),salint(maxnlay),
     &             thick0(maxnlay),
     &             thick1(maxnlay),
     &             hsold,hgold,hsnew,hgnew
 
      COMMON/simip/zsim,tempsim,salsim,
     &             tempint,salint,
     &             thick0,
     &             thick1,
     &             hsold,hgold,hsnew,hgnew
 
      real(8) zm0(0:maxnlay),zm1(0:maxnlay),qm0(maxnlay),
     &              qm1(0:maxnlay+2), sal_new_layer

      COMMON/vertres/zm0,zm1,qm0,qm1, sal_new_layer

      real(8) xlgm,xlgn,cpg,cpw,gammac,betak,tmut,tpw,
     &              cpoc,deltah, betak1, betak2, visc_br, beta_ocs

      COMMON/heateqcoe/xlgm,xlgn,cpg,cpw,gammac,betak,tmut,tpw,
     &              cpoc,deltah, betak1, betak2, visc_br, beta_ocs

      real(8) sal_read(11), hi_read(11), hgins, hnins,
     &              tsuins, oce_sal, oce_flx, num_sal, nday1,
     &              ipremjour, 
     &              i_sal
     
      COMMON/barrowconf/sal_read, hi_read, hgins, hnins,
     &              tsuins, oce_sal, oce_flx, num_sal, nday1,
     &              ipremjour, 
     &              i_sal
     
      real(8) numd_sn1, numd_sn2, numd_sn3, 
     &                  sn_prec_1, sn_prec_2

      COMMON/snowprecip/numd_sn1, numd_sn2, numd_sn3, 
     &                  sn_prec_1, sn_prec_2

      real(8) sf_mult, tabq_ano

      COMMON/tuneforcing/sf_mult, tabq_ano

      real(8) flu_beta, flu_bvtr, rad_io,
     &                frtr_si_phy, qsummer, d_br_mol, d_br_tur,
     &                ra_c, ra_smooth, e_tres, delta_cw, ini_swi, s_ini

      COMMON/fluidtpt/flu_beta, flu_bvtr, rad_io,
     &                frtr_si_phy, qsummer, d_br_mol, d_br_tur,
     &                ra_c, ra_smooth, e_tres, delta_cw, ini_swi, s_ini

      CHARACTER*2 gravdr
      CHARACTER*4 alb_char, pre_char, sal_char
      CHARACTER*12 name_file_bar_par
      CHARACTER*17 name_file_bar_for
       
      COMMON/chars/gravdr, alb_char, pre_char, sal_char, 
     &             name_file_bar_par, namefile_bar_for

      ! radiation transfer in sea ice
      REAL(8) ::
     &   rad_inot_s, rad_inot_i, rad_kappa_s, rad_kappa_i
 
      COMMON/icerad/
     &   rad_inot_s, rad_inot_i, rad_kappa_s, rad_kappa_i
     
!
!--fin du fichier "ice.com"
!-------------------------------------------------------------------------------
