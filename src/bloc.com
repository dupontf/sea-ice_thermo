c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  fichier "bloc.com" : incorpore par instruction 'include' dans les routines :
c      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN,
c      start, flucor, scale, slopez, slopes, scadew, scadns, scali,
c      uve, uvi, barot, uvb0et, uvbfet, uvm,
c      informe, defcst, defgrid, redforc, correct,
c      conti3d, etat, vdiffu, alph2dc, alphdkc, alphdec, raccord,
c      savrunb, redrunb, savrunc, redrunc,
c      ncdfout, moyen, streamv, meridflu, scadhv, checkwst, streamh, stream1h,
c      vague, local, binout, defgl, unigl, foroutp, lisstab, flowsurf.
c  inclus apres "type.com", "para.com".
c  modif : 02/02/98
 
c--blocs common :
 
c--Tableaux de variables en evolution (+ transfert entre routines) :
      common / iteration /
     &   pmix, numit, numspl, ninstb, nclin, nclmoy
 
      common / dynam1 / tpstot,
     &  fss(imax,jmax,nsmax), daeta(imax,jmax),
     &  eta(imax,jmax), ub(imax,jmax), vb(imax,jmax),
     &  u(imax,jmax,kmax), v(imax,jmax,kmax),
     &  scal(imax,jmax,kmax,nsmax), b(imax,jmax,kmax),
     &  bvf(imax,jmax,kmax), avsdz(imax,jmax,kmax),
     &  avudz(imax,jmax,kmax), w(imax,jmax,kmax+1),
     &  q2turb(imax,jmax,kmax+1), fqajc(imax,jmax,kmax)
     & ,vlturb(imax,jmax,kmax), avqdz(imax,jmax,kmax),
     &  tm2tur(imax,jmax,kmax)
 
c--Tableaux de flux et facteurs d'evolution et variables en evolution :
      common / dynam2 / deriv(nsmax),
     &  q(imax,jmax,kmax), fub(imax,jmax,kmax), fvb(imax,jmax,kmax),
     &  etaspl(imax,jmax), ubspl(imax,jmax), vbspl(imax,jmax),
     &  phizzz(imax,jmax,kmax+1,2), phihhh(imax,jmax,6),
     &  umoy(imax,jmax), vmoy(imax,jmax)
 
Csai  common / saison /
Csai &  tmens(imax,jmax,nmois),smens(imax,jmax,nseas),
Csai &  txmens(imax,jmax,nmois),tymens(imax,jmax,nmois),
Csai &  d2tmns(imax,jmax,nmois),d2smns(imax,jmax,nseas),
Csai &  d2txms(imax,jmax,nmois),d2tyms(imax,jmax,nmois),
Csai &  flxss(imax,jmax,nsmax,nmois),d2fxss(imax,jmax,nsmax,nmois),
Csai &  flxsur(imax,jmax,nsmax)
 
      common / bulk_forcing /
     &  scal0(kmax,nsmax), spvr, scalr(imax,jmax,kmax,nsmax),
     &  rappel(imax,jmax,kmax),rappes(imax,jmax,kmax,nsmax),ahrap,
     &  rappes1(1,1,kmax,26),rappes2(1,1,kmax,28),rappes3(1,1,kmax,29),
     &  rappes4(1,1,kmax,31), phimnx(imax,jmax,0:1,nsmax),
     &  phifu(imax,jmax), phifv(imax,jmax), phifs(imax,jmax,nsmax),
     &  phisu(imax,jmax), phisv(imax,jmax), phiss(imax,jmax,nsmax)
Cadh & ,flxus(imax,jmax,nsmax), flxvs(imax,jmax,nsmax)
     & ,phivs(imax,jmax,kmax,nsmax)
     & ,ust2s(imax,jmax),ust2b(imax,jmax)
 
      common / icouplage /
     &  master, icoupl, icoutp, itau_slow
 
Ccpl  common / tbforcing /
Ccpl &  tfmocn(imax,jmax,ntocn), ttoocn(imax,jmax,ntatm)
 
      common / cfmetric /
     &  cmx(imax,jmax,0:3), cmy(imax,jmax,0:3),
     &  smx(imax,jmax,0:3), smy(imax,jmax,0:3),
     &  cmxy(imax,jmax,0:3), smxy(imax,jmax,0:3),
     &  cmxdy(imax,jmax), cmydx(imax,jmax),
     &  fs2cor(imax,jmax), aire(imax,jmax),
Cfcc &  fcucor(imax,jmax), fcvcor(imax,jmax),
     &  covrai(imax,jmax)
 
      common / domain /
     &  dx, dy, unsdx, unsdy, uns2dx, uns2dy,
     &  z(kmax+1), dz(kmax), unsdz(kmax),
     &  zw(kmax+1), dzw(kmax+1), unsdzw(kmax+1),
     &  huy(imax,jmax), hux(imax,jmax), hu(imax,jmax),
     &  unshu(imax,jmax),
     &  tms(imax,jmax,kmax), tmu(imax,jmax,kmax)
 
      common / surfvol /
     &  ctmi(imax,jmax,kmax,0:1),
     &  zsurfs(kmax), zsurfo(kmax), zsurfv(kmax),
     &  zvols, zvolo, zvolv, zvolw, zsurf,
Ccpl &  zsurfsla(kmax,0:nltmax), zsurfola(kmax,0:nltmax),
Ccpl &  zsurfsba(kmax,0:nbsmax), zsurfoba(kmax,0:nbsmax),
Ccpl &  zvolsla(0:nltmax), zvolola(0:nltmax),
Ccpl &  zvolsba(0:nbsmax), zvoloba(0:nbsmax),
     &  unsvol
 
      common / limites /
     &  kniv(imax,jmax,-1:1),
     &  kfs(imax,jmax), kfu(imax,jmax), ks1, ks2, ku1, ku2,
     &  is1(jmax) , is2(jmax) , ims1, ims2, js1, js2,
     &  iu1(jmax) , iu2(jmax) , imu1, imu2, ju1, ju2,
     &  isf1(jmax), isf2(jmax), iuf1(jmax), iuf2(jmax),
     &  jcl1, jcl2, jeq, iberp, jberp, ibera, jbera,
     &  jdl1, jdl2, ijsdl, ijudl
 
      common / lcoins /
     &  i1coin(ncomax,kmax), i2coin(ncomax,kmax),
     &  i3coin(ncomax,kmax), i4coin(ncomax,kmax),
     &  n1coin(kmax), n2coin(kmax), n3coin(kmax), n4coin(kmax),
     &  ijslp(nlpmax), kslp(nlpmax), lslp(nlpmax),
     &  nxslp, nxyslp
 
      common / lerun /
     &  nstart,nend,nitrun,nlast,nsav,ninfo,nfr_out,
     &  nwjl,nwtal,nwm,nwa,nwtest,idyn,nit,premjour,TMIN,JIMP,indice1,
     &  indice2,indice3,indice4,lstab,kstart,kinput,koutpu,
     &  nitrap, ntmoy,lecture,
     &  nyear1 ! FD
 
      common / runpara /
     &  dts(kmax), dtu, dtb, cdbot, zlatz,
     &  ahs(kmax), ahu, ahe, alphxu, alphxv, alphyu, alphyv,
Ciso &  ai(kmax)  , slopemax(kmax),
Ciso &  aitd(kmax), slopmgm(kmax), afilt, ahh, avv,
     &  alphgr(nsmax), algrmn(nsmax), alphah(2), alphmi(kmax),
     &  alphaz(kmax), rifsmx, rifumx,
     &  avnu0(kmax), avnub(kmax), avk0(kmax), avkb(kmax),
     &  txiadu, txiads, txidfu, txidfs, txeflx(nsmax),
     &  bering, ajcmix, xslop, y, d, t, dtsd2,cdiv,
     &  year,day,zbath(kmax),dzbath(kmax),tenhrx(imax,jmax),
     &  tenhry(imax,jmax),vits(nsmax),vabqec(imax,jmax),depth,pCO2,
     &  depthrho,ACTT(kmax,30)
 
      common / coetur/
     &  vkappa,q2tmin,ghmax,ghmin,zlotur,vlmin,varfor,sqrghm,
     &  kajul
 
      character*6 refexp
      common / cerun /
     &  refexp
 
c--fin du fichier "bloc.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
