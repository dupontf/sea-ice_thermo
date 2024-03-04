c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  fichier "type.com" : incorpore par instruction 'include' dans les programmes
c   (et les routines des programmes)	:
c      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN, TRSBATH.
c contient les definitions de type (sauf les chaines de caracteres).
c  modif : 08/01/95
 
c--declaration implicite de type (standard fortran) :
      implicit double precision (a-h,o-z)
Cray  implicit real (a-h,o-z)
 
c--declaration explicite de type :
      complex*16 acplex, bcplex, vcplex
Cray  complex*8  acplex, bcplex, vcplex
 
c--declaration explicite PVM de type :
Ccpl  integer PVMDEFAULT/0/
Ccpl  integer INTEGER4/3/
Ccpl  integer REAL8/6/
 
c--fin du fichier "type.com"
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
