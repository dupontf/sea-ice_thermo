c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  fichier "reper.com" : incorpore par instruction 'include' dans les routines :
c     CLASS, GRADP, FEEDOM, UNIBIN, barot,
c     informe, defcst, defgrid, geogra, redforc,
c     ncdfout, moyen, streamh, streamv, meridflu,
c     checkwst, local, binout, sepgl, defgl, unigl.
c (commun a toutes les routines de traitement des resultats et de
c  definition du domaine) ; inclus apres "type.com", "para.com".
c-----
c  modif : 29/01/98
 
c--blocs common :
 
c-----
c- variables that must be transfered from defcst to defgrid and redforc
 
      common / kdefini /
     &  kbath(imax,jmax), kbath1(imax), kbath2(imax),
     &  kfond, kforc, mdforc
 
      common / deftra /
     &  rapp0(kmax), rapp1(kmax), unstyr, yforc(nsmax), unitfx(nsmax)
 
      character*40 filcor
      common / cfilcor /
     &  filcor
 
c-----
c- variables utilisees dans routine "informe" pour sortie sur fich. "evolu" :
      common / jinform /
     &  nvinfo, nvhsf, nferme, icheck, jcheck, kcheck
 
      common / sinform /
     &  zninfo, zmdeta,
     &  vinfor(ninfmx), scalwr(nsmax), vinfom(ninfmx)
 
      character*30 fmtw
      character*(nchsep) titvar
      common / cinform /
     &  fmtw, titvar(ninfmx)
 
c-----
c- variables utilisees dans la definition du domaine :
      common / distances /
     &  dlong, dlat, dniv,
     &  xlon1, ylat1, zniv1, xalon1, yalat1,
     &  dxaj, dyai, xaj1, yai1, xwpoln
 
      equivalence (dxwi , dlong), (dywj , dlat)
      equivalence ( xwi1, xlon1), ( ywj1, ylat1)
 
c- indices delimitant les grilles, zones, bassins et detroits :
      common / morceaux /
     &  ndhsf, jsep(imax), jnorth(imax), jgeogr(imax,jmax,0:1),
     &  jezon(0:nbsmax), iszon(jmax,0:nbsmax), iezon(jmax,0:nbsmax),
     &  jsbas(0:nbsmax), jebas(0:nbsmax), icl(-imax:imax+imax),
     &  isbas(jmax,0:nbsmax), iebas(jmax,0:kmax,0:nbsmax),
     &  ishsf(nhsfmx), iehsf(nhsfmx), jshsf(nhsfmx), jehsf(nhsfmx)
 
c- variables et coeff. utilisee pour Moyenne_Zonale(Vraie_Latitude)
      common / vraielat / yvrlat(imax,jmax),
     &  rgeogr(imax,jmax,0:1), ageogr(imax,jmax), bgeogr(imax,jmax)
 
c- titres associes :
      character*10 titzon
      character*20 tithsf
      common / cclieu /
     &  titzon(0:nbsmax), tithsf(nhsfmx)
 
c-----
c--variables utilisees pour regrouper les resultats sur a grille Globale :
C     common / globgrid /
C    & tmgl(imax,jmax,kmax,2), csgl(imax,jmax,4), surfgl(imax,jmax,2)
 
C     common / limglob /
C    & kfgl(imax,jmax), imdl1(jmax), imdl2
 
c--fin du fichier "reper.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
