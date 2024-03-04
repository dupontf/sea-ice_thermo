!-------------------------------------------------------------------------------
!
!  fichier "para.com" : incorpore par instruction 'include' dans les programmes
!   (et les routines des programmes) :
!      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN, TRSBATH.
!
 
!--parametres lies a la taille du domaine :
! [indispensable pour inclusion des fichiers bloc.com, reper.com, var??.com]
      parameter ( imax = 1 , jmax = 1 , kmax = 230 )
      parameter ( nsmax = 33 )
      parameter ( nbpt=imax*jmax )
      parameter ( ixjmax = imax*jmax , ijkmax = ixjmax*kmax )
 
 
!--Nombre maximum de coins (par niveaux et par type) :
      parameter ( ncomax = 100 )
!--Nombre maximum d'arretes (suivant X, Y, et pour les kmax Niv.)
!      (ordre de grandeur : imax*jmax)
      parameter ( nlpmax = 6000 )
 
!--parametres donnant le nombre de niveaux dans la glace
!     parameter ( nkb0 = 3 )
      parameter ( maxnlay  = 230 )
      parameter ( nlay_bio = 10 )
      parameter ( ntra_bio = 7  )
      parameter ( n_forc = 11 ) ! number of forcing fields
      parameter ( numout   = 84 )
      parameter ( jpl = 1 ) ! number of ice categories
 
!--determine le type de bassin : -1=rectangulaire, 0=ferme(rectangulaire + courb
!--1=cyclique , 2= 2 grilles separees qui sont raccordees, 3= grille en coordonn
      parameter ( ltest = -1 , jsepar = 50 )
 
!--parametres fixant le nombre de zones (= nb bassins), nb de detroits :
! [indispensable pour inclusion du fichier  reper.com]
      parameter ( nbsmax = 3 , nhsfmx = 10 )
!--parametres pour les sorties sur fichier "evolu" par la routine "informe" :
! [indispensable pour inclusion du fichier  reper.com]
      parameter ( ninfmx = 30 + nhsfmx + 3*kmax + 15)
      parameter ( nchinf = 5 , nchsep = nchinf+2 )
 
!--parametres lies a la frequence des donnees pour T , S , tau.x et tau.y
      parameter ( nmois = 12 , nseas = 4 )
 
!--parametres lies a la definition du tableau general utilise pour les sorties:
! [indispensable pour inclusion du fichier  var??.com]
      parameter( ltymax = 11 )
      parameter( krlmin = -4 - nsmax)
      parameter( nvmax = 99 , nv3dmx = 9+nsmax )
!- nv2dmx = Nb. Var. 2D (fixe+Ns+Ice) ; kv2dmx = Nb. Niv. reserves pour var. 2D
      parameter( nv2dmx = 8 + nsmax + 15 , kv2dmx = 10 + 3*nsmax + 6 )
      parameter( krlmax = krlmin + 1 + nv3dmx*kmax + kv2dmx )
!--parametre indiquant le rang "k" (ds le tableau general),
!   du 1er niveau occupe par la variable:
! [indispensable pour inclusion du fichier  vareq.com]
      parameter( krlu  =  0 )
      parameter( krlfc = krlu - nsmax - 4 )
      parameter( krlfs = krlfc + 1 )
      parameter( krlfs3= krlfc + 2 )
      parameter( krlfs4= krlfc + 3 )
      parameter( krlps = krlu - 4 )
      parameter( krlet = krlu - 3 )
      parameter( krlub = krlu - 2 )
      parameter( krlvb = krlu - 1 )
      parameter( krlv  = krlu + kmax   )
      parameter( krlt  = krlu + kmax*2 )
      parameter( krls  = krlu + kmax*3 )
      parameter( krls3 = krlu + kmax*4 )
      parameter( krls4 = krlu + kmax*5 )
      parameter( krlb  = krlu + kmax*(2+nsmax) )
      parameter( krln2 = krlu + kmax*(3+nsmax) )
      parameter( krlas = krlu + kmax*(4+nsmax) )
      parameter( krlau = krlu + kmax*(5+nsmax) )
      parameter( krlw  = krlu + kmax*(6+nsmax) )
      parameter( krltke= krlu + kmax*(7+nsmax) + 1 )
      parameter( krlajc= krlu + kmax*(8+nsmax) + 2 )
 
      parameter( krlusl= krlajc + kmax )
      parameter( krlvsl= krlusl + nsmax   + 2 )
      parameter( krlhac= krlusl + nsmax*2 + 4 )
      parameter( krleac= krlhac + 1 )
 
      parameter( krlhg = krleac + 1 )
      parameter( krlfq = krlhg  + 1 )
      parameter( krlqs = krlhg  + 2 )
      parameter( krlal = krlhg  + 3 )
      parameter( krlhn = krlhg  + 4 )
      parameter( krlts = krlhg  + 5 )
 
      parameter( krlvaf= krlt )
      parameter( krlvdf= krlvaf + kmax )
      parameter( krlaxt= krlvaf + kmax*2 )
      parameter( krlayt= krlvaf + kmax*3 )
      parameter( krlhat= krlvaf + kmax*4 )
      parameter( krlhdt= krlvaf + kmax*5 )
      parameter( krlvat= krlvaf + kmax*6 )
      parameter( krlvdt= krlvaf + kmax*7 + 1)
 
!--fin du fichier "para.com"
!-------------------------------------------------------------------------------
