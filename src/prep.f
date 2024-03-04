      PROGRAM PREP
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  PREProcess fichier text
c  traite un fichier text (Lignes de 80c au plus !) en substituant
c les n premiers caracteres des lignes qui commences par le mot selectione.
c  (dans la liste) par le mot correspondant (de la liste).
c--En entree : Fichier "prep.list" contient :
c    le mot selectione + les noms des fichiers a traiter
c    Nom du fichier de sortie = Nom du fichier d'entree//'.prep'
c  modif : 30/11/95
 
      parameter ( jmax = 100 , lgmot = 10)
      character*1 cc1
      character*10 cmot1(jmax), cmot2(jmax)
      character*10 ccsel, cblmot
      character*30 fmt
      character*50 filinp, filout
      character*80 line, linin, linout
 
c--initialisation :
      do 10 n=1,len(cblmot)
        cblmot(n:n) = ' '
 10   continue
 
c--ouverture et lecture de l'entete du fichier "prep.list"
      open(32, file='prep.list', status='old')
      read(32,*)
      read(32,*)
      read(32,*) nccmot, nbmots
      do 50 j=1,nbmots
        read(32,'(2(A10,5X))') cmot1(j), cmot2(j)
 50   continue
      read(32,*)
 
      do 60 jj=1,nbmots,8
        write(6,'(A)') 'Remplace les mots :'
        jje = min(nbmots,jj+8-1)
        write(6,'(10A)') cmot1(jj)(:nccmot),
     &               (', '//cmot1(j)(:nccmot),j=jj+1,jje)
        write(6,*) ' par :'
        write(6,'(10A)') cmot2(jj)(:nccmot),
     &               (', '//cmot2(j)(:nccmot),j=jj+1,jje)
 60   continue
      write(6,'(A)')'   en debut de ligne, dans les fichiers suivant :'
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Traitement du fichier d'entree "filinp" .                       |
c-----------------------------------------------------------------------
 
 200  continue
      read(32,'(A)',end=990 ) filinp
      ncfili = 0
      do 210 nn=1,len(filinp)
        if (filinp(nn:nn).ne.' ') ncfili = nn
 210  continue
      if (ncfili.eq.0) goto 990
 
      ncfili = min0(ncfili,45)
      filout = filinp(:ncfili)//'.prep'
 
      open(20, file=filinp(:ncfili), status='old')
      open(30, file=filout(:ncfili+5), status='unknown')
 
      write(6,'(A)') filinp(:ncfili)//' - debut du traitement'
 
      nmodif = 0
 250  continue
      read(20,'(A)',END=900) line
      numlin = numlin + 1
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) traitement d'une ligne du fichier d'entree .
c-----------------------------------------------------------------------
 
c- traitement d'une ligne :
      do 310 j=1,nbmots
        if ( line(:nccmot).eq.cmot1(j)(:nccmot) ) then
          nmodif = nmodif + 1
C         line(:nccmot) = cblmot(:nccmot)
          line(:nccmot) = cmot2(j)(:nccmot)
        endif
 310  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Ecriture de la ligne modifiee .
c-----------------------------------------------------------------------
 
      ncline = 1
      do 610 n=1,len(line)
        if(line(n:n).ne.' ') ncline = n
 610  continue
 
      write(30,'(A)') line(:ncline)
 
      goto 250
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
 900  continue
      if (nmodif.ge.1) then
        write(6,'(A,I8,A)') '  resultat dans '//filout(:ncfili+5)
     &                     //' :',  nmodif, ' modifs'
      else
        write(6,'(A,I8,A)') '  resultat dans '//filout(:ncfili+5)
      endif
 
      close(20)
      close(30)
 
      goto 200
 
 990  continue
      write(6,*) 'fin de la liste de fichier'
      close(32)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end
