C--VERSION:2005.04.08
C modified by Martin Vancoppenolle for use with LIM1D
C read of forcing

C  -----------------------------------------------------------------------
C             libUN : User level NetCDF READ / WRITE routines
C
C                     by Philippe Marbaix and Xavier Fettweis
C
C              Compatible with NetCDF version 3.x (or above).
C  -----------------------------------------------------------------------

C   User-frendly interface :
C   ------------------------

c   CF_INI_FILE   : Initialization of the netcf file  
c   CF_CREATE_DIM : Create axis/dimensions
c   CF_CREATE_VAR : Create variables
c   CF_CREATE_FILE: Write the netcdf file
c   CF_WRITE      : Write variables
c   CF_READ3D/2D  : Read variables
c   CF_OPEN       : Open  netcdf file
c   CF_CLOSE      : Close netcdf file
 
C   Main routines :
C   ---------------

c     UNscreate   : General file creation routine,
c                    defining multiple dimensions + attributes

c     UNwrite     : General variables writting routine
c                    (also updates 'range' attribute and variable if present)
c                   Note: Use UNlwrite to write 2D planes in 3D variables

c     UN(s)read   : Reading routine (grid coordinates + variable)

C   Complementary routines :
C   ------------------------

c     UNparam     : set optional parameters of libUN functions
c     UNwopen     : re-open file for writting
c     UNropen     : open file for reading
c     UNgtime     : Find time index for a given time value
c     UNgindx     : Generalization of UNgtime: find value in any 1D data.   
c     UNfindx     : modified version of UNgindx safe for non-monotonic data
c     UNclose     : close the NetCDF file
c     UNwratt     : Real attributes writting 
c     UNwcatt     : Characters attributes creation & writing

C   Double Precision :
C   ------------------

c     To be in double precision, type this 
c     > sed "s/REAL\*4/REAL\*8/g"      libUN.f  > libUN1.f
c     > sed "s/\_REAL/\_DOUBLE/g"      libUN1.f > libUN2.f
c     > sed "s/NF\_FLOAT/NF\_DOUBLE/g" libUN2.f > libUNd.f
c     > rm -f libUN1.f libUN2.f

C  -----------------------------------------------------------------------


C    +---------------------------+---------------------------------------+
C    +  Subroutine CD_INI_FILE : + Initialize the netcdf file            +
C    +---------------------------+---------------------------------------+

      SUBROUTINE CF_INI_FILE (filename, filetitle)

c     Input :
c     =======

c     filename  = name  of the netcdf file
c     filetitle = title in the netcdf file

      IMPLICIT NONE
 
      INCLUDE 'libUN.inc'

      CHARACTER *(*) filename,filetitle  

      CF_attnam(1) = 'actual_range'
      CF_attnum(1) = 2

      CF_varnbrtot =  0 ! Initialization
      CF_dimnbrtot = -1 ! Initialization

      CF_filenam   = filename
      CF_filetit   = filetitle

      END SUBROUTINE CF_INI_FILE


C    +-----------------------------+-------------------------------------+
C    +  Subroutine CF_CREATE_DIM : + Create dimensions/axis              +
C    +-----------------------------+-------------------------------------+

      SUBROUTINE CF_CREATE_DIM (dimname,dimunits,dimdim,vallues)

c     Input :
c     =======

c     dimname  = name of the axis/dimension
c     dimunits = units of the axis/dimension
c     dimdim   = dimensions of the axis/dimension
c     vallues  = vallues of the axis/dimension

      IMPLICIT NONE
 
      INCLUDE 'libUN.inc'

      CHARACTER *(*) dimname,dimunits
      
      INTEGER        dimdim,i
      REAL*4         vallues(dimdim)

      CF_dimnbrtot                 = CF_dimnbrtot + 1

      CF_dimnbrtot                 = max(0,CF_dimnbrtot)

      CF_dimnam(CF_dimnbrtot)      = dimname
      CF_dimnamuni(CF_dimnbrtot)   = dimunits    
      CF_dim(CF_dimnbrtot)         = dimdim
 
      do i = 1,dimdim
      CF_dimval(i,CF_dimnbrtot)    = vallues(i)
      enddo   
    
      END SUBROUTINE CF_CREATE_DIM

C    +-----------------------------+-------------------------------------+
C    +  Subroutine CF_CREATE_VAR : + Create variables			 +
C    +-----------------------------+-------------------------------------+

      SUBROUTINE CF_CREATE_VAR (varname,vartitle,varunits,varaxe4,
     .                          varaxe1,varaxe2,varaxe3)

c     Input :
c     =======

c     varname  = name of the variable
c     vartitle = title of the variable
c     varunits = units of the variable
c     varaxeX  = axes used by the variable (T,X,Y,Z)

      IMPLICIT NONE
 
      INCLUDE 'libUN.inc'

      CHARACTER *(*) varname,vartitle,varunits
      CHARACTER *(*) varaxe1,varaxe2,varaxe3,varaxe4
      
      CF_varnbrtot                 = max (0,CF_varnbrtot + 1)

      CF_varnam(CF_varnbrtot)      = varname
      CF_varnamdim(1,CF_varnbrtot) = varaxe1
      CF_varnamdim(2,CF_varnbrtot) = varaxe2
      CF_varnamdim(3,CF_varnbrtot) = varaxe3
      CF_varnamdim(4,CF_varnbrtot) = varaxe4
      CF_varnamuni(CF_varnbrtot)   = varunits
      CF_vardes(CF_varnbrtot)      = vartitle

      END SUBROUTINE CF_CREATE_VAR

C    +--------------------------------------+----------------------------+
C    +  Subroutine CF_CREATE_VAR_VIA_FILE : + Create variables           +
C    +--------------------------------------+----------------------------+

      SUBROUTINE CF_CREATE_VAR_VIA_FILE (filename)

c     Input :
c     =======

c     filename  = name of the file containing informations
c                 about the variables 

      IMPLICIT NONE
 
      INCLUDE 'libUN.inc'

      CHARACTER*200 filename
      
      CHARACTER*120 tmpvar

      OPEN(unit=999,status='old',file=filename)

980   CONTINUE
      READ (999,'(A120)',end=990) tmpvar

      IF (tmpvar(1:4).eq.'    ') THEN
       CF_varnbrtot = max (0,CF_varnbrtot + 1)
         READ (tmpvar,'(4x,5A9,A12,A50)')
     .         CF_varnam(CF_varnbrtot),
     .         CF_varnamdim(1,CF_varnbrtot),
     .         CF_varnamdim(2,CF_varnbrtot),
     .         CF_varnamdim(3,CF_varnbrtot),
     .         CF_varnamdim(4,CF_varnbrtot),
     .         CF_varnamuni(CF_varnbrtot),
     .         CF_vardes(CF_varnbrtot)
      ENDIF

      GOTO 980
990   CONTINUE

      END SUBROUTINE CF_CREATE_VAR_VIA_FILE

C    +------------------------------+------------------------------------+
C    +  Subroutine CF_CREATE_FILE : + Create the netcdf file             +
C    +------------------------------+------------------------------------+

      SUBROUTINE CF_CREATE_FILE (filename)

c     Input :
c     =======

c     filename  = name  of the netcdf file

      IMPLICIT NONE
 
      INCLUDE 'libUN.inc'

      CHARACTER *(*) filename  

      INTEGER        i,j,id

      INTEGER        UN1_dim(0:CF_dimnbrtot)

      REAL(4)        UN1_dimval(CF_dimmaxlen,0:CF_dimnbrtot)

      CHARACTER*31   UN1_dimnam(0:CF_dimnbrtot),
     .               UN1_dimnamuni(0:CF_dimnbrtot) 

      if(filename.ne.CF_filenam)then
       write(6,*) "ERROR: not "//CF_filenam
       stop
      endif
   
      DO i=0,CF_dimnbrtot
       UN1_dim(i)       = CF_dim(i)
       UN1_dimnam(i)    = CF_dimnam(i)
       UN1_dimnamuni(i) = CF_dimnamuni(i)
! FD debug
      write(*,*) i,TRIM(UN1_dimnam(i)),CF_dim(i)
       DO j=1,CF_dim(i)
        UN1_dimval(j,i) = CF_dimval(j,i)
       END DO
      END DO

! FD debug
!      write(*,*) CF_filenam
!      write(*,*) CF_filetit,CF_dimnbrtot,UN1_dim
!      write(*,*) CF_dimmaxlen, UN1_dimnam ,UN1_dimnamuni
!      write(*,*) UN1_dimval
!      write(*,*) CF_varmaxnbr,CF_varnbrtot,CF_varnam
!      write(*,*) CF_varnamdim,CF_varnamuni,CF_vardes
!      write(*,*) CF_attnbr,CF_attnam,CF_attnum,id
!      stop
      call UNscreate (CF_filenam,CF_filetit,CF_dimnbrtot,UN1_dim,
     .                CF_dimmaxlen, UN1_dimnam ,UN1_dimnamuni,
     .                UN1_dimval,
     .                CF_varmaxnbr,CF_varnbrtot,CF_varnam,
     .                CF_varnamdim,CF_varnamuni,CF_vardes,
     .                CF_attnbr,CF_attnam,CF_attnum,id)
      
      call UNclose   (id)   


      END SUBROUTINE CF_CREATE_FILE

C    +------------------------+------------------------------------------+
C    +  Subroutine CF_WRITE : + Writes variables                         +
C    +------------------------+------------------------------------------+


      SUBROUTINE CF_WRITE (FILEname, VARname , itime,
     &                    Ni,  Nj, Nlev, var)

c     Input :
c     =======

c     FILEname    = name of the netcdf file
c     VARname     = name of variables
c     itime       = index on time axis
c     Ni,Nj,Nlev  = X,Y,Z dimension
c     var         = array of vallues of the variable 

      IMPLICIT NONE

      INCLUDE 'libUN.inc'

      CHARACTER *(*) FILEname,VARname  
      INTEGER        itime
      INTEGER        Ni,  Nj, Nlev,fileid 
      REAL*4         var(Ni, Nj, Nlev)

      if(CF_filenamopened.ne.FILEname) then
      CALL UNwopen (FILEname,fileid)
      else
      fileid = CF_fileidopened
      endif

      CALL UNwrite (fileid,VARname ,itime,Ni,  Nj, Nlev, var)  

      if(CF_filenamopened.ne.FILEname) then
      call UNclose (fileid)
      endif

      END SUBROUTINE CF_WRITE

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine CF_READ1D : + Read variables                          +
C**  +-------------------------+-----------------------------------------+

      SUBROUTINE CF_READDIM (FILEname, VARdim , N)
      ! routine added by martin vancoppenolle to read 1D arrays

c     Input :
c     =======

c     FILEname    = name of the netcdf file
c     VARname     = name of variables
c     itime       = index on time axis
c     N           = X dimension

c     Output :
c     ========

c     var         = array of values of the variable 

      IMPLICIT NONE

      INCLUDE 'libUN.inc'

      CHARACTER *(*) FILEname,VARdim
      INTEGER        N
      INTEGER        dimID, Ierro, fileid
      CHARACTER*31   filetitle  

      if(CF_filenamopened.ne.FILEname) then
      CALL UNropen (FILEname,fileid,filetitle)
      else
      fileid = CF_fileidopened
      endif

      Ierro=NF_INQ_DIMID(fileid, VARdim, dimID)
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('READDIM', Ierro)

      Ierro=NF_INQ_DIM(fileid , dimID, VARdim, N)
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('READDIM', Ierro)

      END SUBROUTINE CF_READDIM


C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine CF_READ1D : + Read variables                          +
C**  +-------------------------+-----------------------------------------+

      SUBROUTINE CF_READ1D (FILEname, VARname , itime,
     .                      N, var)
      ! routine added by martin vancoppenolle to read 1D arrays

c     Input :
c     =======

c     FILEname    = name of the netcdf file
c     VARname     = name of variables
c     itime       = index on time axis
c     N           = X dimension

c     Output :
c     ========

c     var         = array of values of the variable 

      IMPLICIT NONE

      INCLUDE 'libUN.inc'

      CHARACTER *(*) FILEname,VARname
      CHARACTER*31   var_units,filetitle  
      INTEGER        N, itime,level
      REAL*4         var(N)

      INTEGER        i,j,fileid

      if(CF_filenamopened.ne.FILEname) then
      CALL UNropen (FILEname,fileid,filetitle)
      else
      fileid = CF_fileidopened
      endif

      CALL UNsread (fileid, VARname, itime, 1, 1, 1,
     &              N , 1 , 1,var_units, var)          

      if(CF_filenamopened.ne.FILEname) then
      call UNclose (fileid)
      endif

      END SUBROUTINE CF_READ1D

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine CF_READ2D : + Read variables                          +
C**  +-------------------------+-----------------------------------------+

      SUBROUTINE CF_READ2D (FILEname, VARname , itime,
     .                      Ni,  Nj, Nlev, var)


c     Input :
c     =======

c     FILEname    = name of the netcdf file
c     VARname     = name of variables
c     itime       = index on time axis
c     Ni,Nj,Nlev  = X,Y,Z dimension

c     Output :
c     ========

c     var         = array of vallues of the variable 

      IMPLICIT NONE

      INCLUDE 'libUN.inc'

      CHARACTER *(*) FILEname,VARname
      CHARACTER*31   var_units,filetitle  
      INTEGER        Ni,  Nj, Nlev,itime,level
      REAL*4         var(Ni, Nj)

      INTEGER        i,j,fileid

      if(CF_filenamopened.ne.FILEname) then
      CALL UNropen (FILEname,fileid,filetitle)
      else
      fileid = CF_fileidopened
      endif

      CALL UNsread (fileid, VARname, itime, Nlev, 1, 1,
     &              Ni , Nj , 1,var_units, var)          

      if(CF_filenamopened.ne.FILEname) then
      call UNclose (fileid)
      endif

      END SUBROUTINE CF_READ2D

C    +-------------------------+-----------------------------------------+
C    +  Subroutine CF_READ3D : + Read variables                          +
C    +-------------------------+-----------------------------------------+


      SUBROUTINE CF_READ3D (FILEname, VARname , itime,
     .                      Ni,  Nj, Nlev, var)

c     Input :
c     =======

c     FILEname    = name of the netcdf file
c     VARname     = name of variables
c     itime       = index on time axis
c     Ni,Nj,Nlev  = X,Y,Z dimension

c     Output :
c     ========

c     var         = array of vallues of the variable 

      IMPLICIT NONE

      INCLUDE 'libUN.inc'

      CHARACTER *(*) FILEname,VARname
      CHARACTER*31   var_units,filetitle  
      INTEGER        Ni,  Nj, Nlev,itime,level 
      REAL*4         var(Ni, Nj,Nlev)

      INTEGER        i,j,fileid

      if(CF_filenamopened.ne.FILEname) then
      CALL UNropen (FILEname,fileid,filetitle)
      else
      fileid = CF_fileidopened
      endif

      CALL UNsread (fileid, VARname, itime, 0, 1, 1,
     &              Ni , Nj , Nlev,var_units, var)          

      if(CF_filenamopened.ne.FILEname) then
      call UNclose (fileid)
      endif

      END SUBROUTINE CF_READ3D

C**  +------------------------+------------------------------------------+
C**  +  Subroutine CF_CLOSE : + Close the file                           +
C**  +------------------------+------------------------------------------+

      SUBROUTINE CF_CLOSE (FILEname)

      IMPLICIT NONE

      INCLUDE 'libUN.inc'
   
      CHARACTER*(*) FILEname

      if(FILEname.eq.CF_filenamopened)then
      call UNclose (CF_fileidopened)
      else
      print *,FILEname//" not opened"
      endif
  
      CF_filenamopened = ""
      CF_fileidopened  = 0

      END SUBROUTINE CF_CLOSE

C**  +-----------------------+-------------------------------------------+
C**  +  Subroutine CF_OPEN : + open the file                             +
C**  +-----------------------+-------------------------------------------+

      SUBROUTINE CF_OPEN (FILEname,FILEid)

      IMPLICIT NONE

      INCLUDE 'libUN.inc'
  
      INTEGER       FILEid 

      CHARACTER*(*) FILEname

      call UNwopen (FILEname,FILEid)
     
      CF_filenamopened = FILEname
      
      CF_fileidopened  = FILEid

      END SUBROUTINE CF_OPEN

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNscreate : +                                         +
C**  +-------------------------+                                         +
C**  +  * Purpose :                                                      +
C**  +     Create a NetCDF file, general version.                        +
C**  +     (Staggered grids + other extensions to UNcreate)              +
C**  +                                                                   +
C**  +  * How it works : calling routine must provide                    +
C**  +    -a list of dimensions                                          +
C**  +     (size of each dimens., names, units and values of coordinates)+
C**  +    -a list of variables                                           +
C**  +     (units, number of dimensions, names of selected dimensions)   +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +  -------                                                          +
C**  +                                                                   +
C**  +  General :                                                        +
C**  +   FILEnam          [char]: Name of the file to be created.        +
C**  +   title            [char]: Title attribute                        +
C**  +                                                                   +
C**  +  Dimensions:                                                      +
C**  +   TND                    : Total Number of SPATIAL dimensions     +
C**  +                            Notice : Set "time" to dimension No 0  +
C**  +   DFdim(0:TND)           : # discrete values for each dimension   +
C**  +                            Notice : DFdim(0).eq.0                 +
C**  +                            -> UNLIMITED TIME (coord. not defined) +
C**  +                               WARNING: In this case, the NetCDF   +
C**  +                               use a temporary space to duplicate  +
C**  +                               the file -> NOT RECOMMENDED         +
C**  +   MXdim                  : Maximum value of DFdim, = arrays size  +
C**  +   NAMdim(0:TND)    [char]: Name of dimensions, except time        +
C**  +   UNIdim(0:TND)    [char]: Units of dimensions (attribute)        +
C**  +   VALdim(MXdim,0:TND)[R4]: Values of coordinate for each dimension+
C**  +                                                                   +
C**  +  Variables:                                                       +
C**  +   Dvs                    : Variable's definitions array sizes,    +
C**  +   Nvs                    : Number of defined variables(Nvs.le.Dvs)+
C**  +   name_vs (Dvs)    [char]: name of variable.                      +
C**  +   unit_vs (Dvs)    [char]: physical units of variable (attribute) +
C**  +   Sdim_vs (4,Dvs)  [char]: name of Selected dims (in above list)  +
C**  +                            Blanked or '-' elements = not used     +
C**  +   lnam_vs (Dvs)    [char]: Long_name attribute (descript. of var.)+
C**  +                                                                   +
C**  +  List of real attributes to all variables:                        +
C**  +   Nra                    : Number of Real Attributes (.ge.1 !)    +  
C**  +   NAMrat(Nra)      [char]: NAMes of Real ATtributes  (''=none)    +
C**  +                            (initial value= 0; set it with UNwratt)+
C**  +   Nvals(Nra)             : Number of values of these attributes.  +
C**  +   ! Currently limited to 1 value (scalar) or 2 (2 elements vector)+
C**  +   ! EXCEPTION: Setting the last attribute name to '[var]_range'   +
C**  +                does create a variable (!) for level-by-level range+
C**  +                (very usefull for 3D + time fields)                +
C**  +                                                                   +
C**  +  NB : [char] variables may have any length.                       +
C**  +       blanks characters are NOT ALLOWED in any variable,          +
C**  +          except the "title".                                      +
C**  +          and the NetCDF variables defined here are always real*4  +
C**  +                                                                   +
C**  +  OUTPUT :                                                         +
C**  +  --------                                                         +
C**  +   FILEid                 : Index of the NetCDF file (remains open)+
C**  +-------------------------------------------------------------------+

      SUBROUTINE UNscreate (FILEnam, title,
     &      TND, DFdim, MXdim, NAMdim, UNIdim, VALdim, 
     &      Dvs, Nvs, name_vs, Sdim_vs, unit_vs, lnam_vs,
     &      Nra, NAMrat, Nvals,
     &      FILEid )

C +
      IMPLICIT NONE
 
      INCLUDE 'libUN.inc'

C +
      INTEGER icheck, MXND
C     ** Maximum number of dimensions 
      parameter (MXND = 100) 

C +   INPUT:      
C +   - - -
      CHARACTER *(*) FILEnam
      CHARACTER *(*) title  

      INTEGER        TND, DFdim(0:TND), MXdim
      CHARACTER *(*) NAMdim(0:TND)   
      CHARACTER *(*) UNIdim(0:TND)
      REAL*4         VALdim(MXdim,0:TND)

      INTEGER        Nvs, Dvs
      CHARACTER *(*) name_vs(Dvs)
      CHARACTER *(*) Sdim_vs(4,Dvs)
      CHARACTER *(*) unit_vs(Dvs)
      CHARACTER *(*) lnam_vs(Dvs)

      INTEGER        Nra
      CHARACTER *(*) NAMrat(Nra)
      CHARACTER*24   Host,Fdate
      CHARACTER*200  tmpchar
      INTEGER        Nvals(Nra)

C +   OUTPUT:    
C +   - - - -
      INTEGER FILEid 

C +   LOCAL:
C +   - - -
      INTEGER        VARSIZE
      EXTERNAL       VARSIZE
      CHARACTER*(30) tmpchr
      INTEGER        dimDID(0:MXND)
      INTEGER        dimVID(0:MXND), vsVID, vrVID
      INTEGER        dID(4), start(4), count(4), rdID(2)
      INTEGER        mimaID
      INTEGER        stride(4),imap(4)
      INTEGER        Ndim_vs
      INTEGER        ivs, igd, idi, ira, itmp
      INTEGER        Nlen
      INTEGER        dNlen(0:MXND)  
      INTEGER        Ierro, TTerr, ii,jj
      REAL*4         zero1(1), zero2(2)
      
      icheck= 0 !Debugging level


C*    0. Initialisations
C     ------------------
      IF (icheck.ge.1) WRITE(*,*) 'UNscreate : Begin'

C +   Routines which opens a file must reset libUN internals:
      CALL UNparam('RESET_PARAMS_',0.0_4)

      DO ii = 1,4
        stride(ii) = 1
      ENDDO
      zero1(1) = 0.
      zero2(1) = 0.
      zero2(2) = 0.
      TTerr = 0 !Total of error flags

      IF (TND .gt. MXND) THEN
        write(*,*)'UNscreate - Error: so much dimensions ?',TND
      END IF

C     Create a NetCDF file and enter define mode :
C     --------------------------------------------
      IF (icheck.ge.2) WRITE(*,*) 'FILEnam :', FILEnam

C     ** getting FILEnam [char] size :
      Nlen = VARSIZE(FILEnam)

      Ierro=NF_CREATE(FILEnam(1:Nlen), NF_CLOBBER , FILEid)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
C     ** identif.                       =>overwrite =error

C*    Time coordinate definition.
C     ---------------------------

C     ** Define dimension :    
      IF (icheck.ge.3) WRITE(*,*) '# time iters.:', DFdim(0)
      IF (DFdim(0).eq.0.) THEN
        Ierro=NF_DEF_DIM(FILEid , 'time', NF_UNLIMITED, dimDID(0))
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
        TTerr = TTerr + ABS(Ierro)
      ELSE
        Ierro=NF_DEF_DIM(FILEid , 'time', DFdim(0), dimDID(0))
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
        TTerr = TTerr + ABS(Ierro)
      END IF
      dNlen(0)= 4  ! 4 characters in the name 'time'...
      IF (NAMdim(0)(1:4).ne.'time') THEN
        WRITE(*,*) 'Sorry, NAMdim(0) must be ''time'' .'
        STOP
      END IF
       
C     ** Define variable for the time coordinate values :
      dID(1)    = dimDID(0)
      Ierro=NF_DEF_VAR(FILEid , 'time', NF_FLOAT,1 , dID,  dimVID(0))
C     **      ^^^^^^^^^^ FILEid  var name  type  dims  DIMid VARid  
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
      TTerr = TTerr + ABS(Ierro)
! FD add calendar definition
      Ierro=  NF_PUT_ATT_TEXT(FILEid , dimVID(0) ,'calendar',
     &                          7   ,'360_DAY')
      TTerr = TTerr + ABS(Ierro)


C     Spatial coordinates definitions : DIMS and VARs (locations).
C     ------------------------------------------------------------
C
      DO igd = 1,TND            !** BEGIN LOOP over all spatial dims
        IF (icheck.ge.3) WRITE(*,*) '  spatial dim:', NAMdim(igd)
        
C       ** getting NAMdim [char] size :
        Nlen = VARSIZE(NAMdim(igd)) 
        dNlen(igd) = Nlen  !For further use of NAMdim

        Ierro=NF_DEF_DIM(FILEid    , NAMdim(igd)(1:Nlen),
     &                     DFdim(igd),dimDID(igd))
C       **line1 ^^^^^^^^^^ FILEid    | dim name            
C       **line2            # values  | VARid
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
        TTerr = TTerr + ABS(Ierro)

        dID(1)     = dimDID(igd)  
        Ierro=NF_DEF_VAR(FILEid    , NAMdim(igd)(1:Nlen),
     &                     NF_FLOAT  ,    1  , dID     ,dimVID(igd))
C       **line1 ^^^^^^^^^^ FILEid    | dim name            
C       **line2            type      | #dims | dimsIDs | VARid 
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
        TTerr = TTerr + ABS(Ierro)

      END DO                    !** END   LOOP over all spatial dims

C     Special coordinate definition: MinMax (for [var]_range)
C     -------------------------------------------------------
      IF (NAMrat(Nra)(1:11).eq.'[var]_range') THEN

        Ierro=NF_DEF_DIM(FILEid, 'MinMax', 2, mimaID)
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
      ENDIF

C     Define the fields. 
C     ------------------

      DO ivs = 1,Nvs             !**BEGIN LOOP on var. num.
        IF (icheck.ge.3)
     &    WRITE (*,*) 'Defining variable ',name_vs(ivs)


C       Set space and time dimensions
C       - - - - - - - - - - - - - - -
C        ** Initialise number of dimensions :
         Ndim_vs= 0 

         DO idi = 1, 4           !** BEGIN LOOP on var dims.
         IF  (Sdim_vs(idi,ivs)(1:1).ne.' '
     &   .and.Sdim_vs(idi,ivs)(1:1).ne.'-') THEN !**skip undefined. 

C         ** getting Sdim_vs [char] size :
          Nlen =  VARSIZE(Sdim_vs(idi,ivs))

C         ** Searching for the dimension index from its name (Sdim_vs)    
          igd = 0
          DO WHILE (Sdim_vs(idi,ivs)(1:Nlen)
     &        .ne. NAMdim(igd)(1:dNlen(igd)) )
            IF (igd.eq.TND) THEN 
              write(*,*)'UNscreate-ERROR: Dimension not found:',
     &              Sdim_vs(idi,ivs)(1:Nlen)
              STOP
            END IF
            igd = igd + 1
          END DO               
C         ** Construct the dimensions id's for that variable (ivs):
          IF (icheck.ge.3)
     &       WRITE (*,*) 'using dimension ',NAMdim(igd), dimDID(igd)
          Ndim_vs      = Ndim_vs + 1
          dID(Ndim_vs) = dimDID(igd) 

        END IF
        END DO                   !** END   LOOP on var dims.

C       Define our special [var]_range field for 4D variables
C       - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF  (Ndim_vs.eq.4 
     &  .and.NAMrat(Nra)(1:11).eq.'[var]_range') THEN 

          Nlen = VARSIZE(name_vs(ivs))
          rdID(1)  = dID (3)  !(4D variable, 3th dim = level)
          rdID(2)  = mimaID   !(for min, max)
          tmpchr = name_vs(ivs)(1:Nlen)//'_range'
          itmp   = Nlen + 6
          Ierro =  NF_DEF_VAR(FILEid,tmpchr(1:itmp),
     &                        NF_FLOAT, 2, rdID, vrVID)
          IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
          TTerr = TTerr + ABS(Ierro)

        ENDIF

C       Define fields :
C       - - - - - - - -
        Nlen = VARSIZE(name_vs(ivs))
        Ierro=NF_DEF_VAR(FILEid , name_vs(ivs)(1:Nlen),
     &                     NF_FLOAT, Ndim_vs, dID     , vsVID)
C       **line1 ^^^^^^^^^^ FILEid | variable name
C       **line2            type   | #dims   | dimsIDs | VARid
        IF (Ierro.NE.NF_NOERR) 
     &      CALL HANDLE_ERR('UNscreate (field)', Ierro)
        TTerr = TTerr + ABS(Ierro)


C     Set the variable's attributes : 
C     -------------------------------

C       ** Units: 
C       - - - - - 
C       ** getting unit_vs [char] size :
        Nlen = VARSIZE(unit_vs(ivs))

        Ierro=  NF_PUT_ATT_TEXT(FILEid , vsVID ,'units',
     &                          Nlen   ,unit_vs(ivs)(1:Nlen))
c       **line1 ^^^^^^^^^^^^^^^ FILEid |var.id | attr.name
C       **line2                 length | attr.value
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
        TTerr = TTerr + ABS(Ierro)
        
C       ** Special case : units = sigma 
C       - - - - - - - - - - - - - - - -
C       In this case, CV convention advises to write the following
C        attribute :  positive = down
C
        Nlen = VARSIZE(lnam_vs(ivs))

        IF ( unit_vs(ivs)(1:Nlen) .EQ. '[sigma]'
     &  .OR. unit_vs(ivs)(1:Nlen) .EQ. 'sigma_level' ) THEN
          IF (icheck.ge.3) THEN
            WRITE(*,*) 'Unit = sigma -> setting positive attr'
          ENDIF   

          Ierro=  NF_PUT_ATT_TEXT(FILEid , vsVID ,'positive',
     &                            4      ,'down')
c         **line1 ^^^^^^^^^^^^^^^ FILEid |var.id | attr.name
C         **line2                 length | attr.value
          IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
        ENDIF

C       ** "long_name":
C       - - - - - - - -
        Nlen = VARSIZE(lnam_vs(ivs))

        IF (icheck.ge.3)
     &    WRITE (*,*) 'Write long_name ',lnam_vs(ivs)(1:Nlen)

        Ierro=NF_PUT_ATT_TEXT(FILEid , vsVID ,'long_name',
     &                          Nlen  ,lnam_vs(ivs)(1:Nlen)     )

        do jj=1,Nlen
          if(lnam_vs(ivs)(jj:jj).eq." ") lnam_vs(ivs)(jj:jj)="_"
          if(lnam_vs(ivs)(jj:jj).eq.".") lnam_vs(ivs)(jj:jj)="_"
          if(lnam_vs(ivs)(jj:jj).eq."(") lnam_vs(ivs)(jj:jj)="_"
          if(lnam_vs(ivs)(jj:jj).eq.")") lnam_vs(ivs)(jj:jj)="_"
          if(lnam_vs(ivs)(jj:jj).eq."/") lnam_vs(ivs)(jj:jj)="_"
        enddo

        Ierro=NF_PUT_ATT_TEXT(FILEid , vsVID ,'standard_name',
     &                          Nlen  ,lnam_vs(ivs)(1:Nlen)     )

        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
        TTerr = TTerr + ABS(Ierro)


C       ** From the list of real attributes (input argument) : 
C       - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
        DO ira = 1, Nra
        IF (NAMrat(ira)(1:1).ne.' ') THEN
        IF (NAMrat(ira)(1:11).eq.'valid_range') THEN 

C         ** The "valid_range" special attribute  :
          Ierro=NF_PUT_ATT_REAL(FILEid  ,vsVID ,'valid_range' ,
     &                          NF_FLOAT,2     , ValRange)
          TTerr = TTerr + ABS(Ierro)

        ELSE IF (NAMrat(ira)(1:11).ne.'[var]_range') THEN 
        
C         ** All "regular" attributes :
          Nlen = VARSIZE(NAMrat(ira))
          IF (Nvals(ira).eq.1) THEN
            Ierro=NF_PUT_ATT_REAL(FILEid,vsVID,NAMrat(ira)(1:Nlen),
     &                           NF_FLOAT,   Nvals  , zero1      )
            TTerr = TTerr + ABS(Ierro)
          ELSE IF (Nvals(ira).eq.2) THEN
            Ierro=NF_PUT_ATT_REAL(FILEid,vsVID,NAMrat(ira)(1:Nlen),
     &                           NF_FLOAT, Nvals  , zero2        )
            TTerr = TTerr + ABS(Ierro)
c
           END IF
        END IF
        END IF
        END DO

      END DO                     ! **END   LOOP on var. num.

C     Set 'unit' attribute for the dimensions:
C     ----------------------------------------

      DO igd = 0,TND         !** BEGIN LOOP over all spatial dims

C       ** getting NAMdim [char] size :
        Nlen = VARSIZE(UNIdim(igd))

        Ierro=NF_PUT_ATT_TEXT(FILEid , dimVID(igd),'units',
     &                          Nlen   , UNIdim(igd)        )

        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)

        Nlen = VARSIZE(NAMdim(igd))

        Ierro=NF_PUT_ATT_TEXT(FILEid , dimVID(igd),'long_name',
     &                          Nlen   , NAMdim(igd)        )

        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)

        Ierro=NF_PUT_ATT_TEXT(FILEid , dimVID(igd),'standard_name',
     &                          Nlen   , NAMdim(igd)        )

        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
        TTerr = TTerr + ABS(Ierro)

      ENDDO

C     Global attribute(s).
C     --------------------

C     ** Title (some general file descriptor) :
C     ** getting unit_vs [char] size :

      Nlen = VARSIZE(title)

      Ierro=NF_PUT_ATT_TEXT(FILEid ,NF_GLOBAL,'title',
     &                      Nlen    ,title(1:Nlen)       )

      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)


      Nlen = VARSIZE(CF_institution)

      Ierro=NF_PUT_ATT_TEXT(FILEid ,NF_GLOBAL,'institution',
     &                      Nlen    ,CF_institution)

      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)

c     CALL HostNm(Host, Ierro)

      tmpchar="libUN ("//CF_libUN_version//") - "//FDate()
c    &        " - "//Host

      Nlen = VARSIZE(tmpchar)

      Ierro=NF_PUT_ATT_TEXT(FILEid ,NF_GLOBAL,'history',
     &                      Nlen    ,tmpchar)     

      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)


      Nlen = VARSIZE(NF_INQ_LIBVERS())

      Ierro=NF_PUT_ATT_TEXT(FILEid ,NF_GLOBAL,'netcdf',
     &                      Nlen    ,NF_INQ_LIBVERS())


      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
      TTerr = TTerr + ABS(Ierro)


C     Leave define mode (!file remains open )
C     ---------------------------------------
      Ierro=NF_ENDDEF(FILEid)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)
      TTerr = TTerr + ABS(Ierro)


C     Writing of dimensions coordinates.
C     ----------------------------------

C     ** Time :
C     - - - - -

      start(1)= 1          !Vector of starting indexes values
      count(1)= DFdim(0)   !Vector of total # indexes values
        IF (icheck.ge.3)
     &    WRITE (*,*) 'Write coords for ',NAMdim(0),count(1)

C     ** Set 'imap' to write with NCVPTG; NCVPT could be enough ?
C     ** (imap tells NetCDF about the memory locations of var,
C     **  we choose NCVPTG because 
C     **  only a portion of VALdim is written.)
      imap(1) = 1
      imap(2) = 0                 ! Not used : write only 1 coord.

      Ierro=NF_PUT_VARM_REAL(FILEid ,dimVID(0), start        , count,
     &                         stride , imap    , VALdim(1,0)         )
C     **line 1 ^^^^^^^^^^^^^^^ ID file| id var. |read from...  |#data
C     **line 2                 step   |re-arrang|variable(beg.)
C     **                      (^^^^stride is not used)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)


C     ** Space coordinates :
C     - - - - - - - - - - - -

      DO igd = 1,TND          !** BEGIN LOOP over all spatial dims

        start(1)= 1
        count(1)= DFdim(igd)
        IF (icheck.ge.3)
     &    WRITE (*,*) 'Write coords for ',NAMdim(igd),count(1)


        Ierro=NF_PUT_VARM_REAL(FILEid ,dimVID(igd),start , count,
     &                           stride , imap      ,VALdim(1,igd))
C       **      ^^^^^^^^^^^^^^^^ see above
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNscreate', Ierro)

        TTerr = TTerr + ABS(Ierro)
 
      END DO                  !** END   LOOP over all spatial dims

C     Stop if an error occured.
C     -------------------------

      IF (TTerr.ne.0) THEN 
        STOP 'UNscreate : Sorry, an error occured.'
      ENDIF

C +
      RETURN
      END SUBROUTINE UNscreate

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNwrite :   +                                         +
C**  +-------------------------+                                         +
C**  +  * Writes a variable into a NetCDF file,                          +
C**  +    (the NetCDF file must have been created (or re-opened) and     +
C**  +     closed after all writing operations).                         +
C**  +  * Automatically updates attribute 'actual_range' if available    +
C**  +          "          "    special var. '[var]_range'    "          +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEid  : input file identifier (from UNcreate OR NetCDF open) +
C**  +    VARname : name given to the variable to write (must be in file)+
C**  +    itime   : No of time step to write to                          +
C**  +    Ni,Nj,Nlev: dimensions of 'var'                                +
C**  +              ! Nlev= 1 for 2D and 1D input variables.             + 
C**  +                Nj  = 1 for 1D input variables.                    +
C**  +              NB: can not write 1 level of 3D var only (->UNlwrite)+
C**  +                                                                   +
C**  +    var     : The variable to be writen                            +
C**  +                                                                   +
C**  +  REMARK :                                                         +
C**  +    Truncation of input data is permited:                          +
C**  +    If the dim of "var" > dim in the NetCDF file,                  +
C**  +    "var" is automatically truncted. However, this => WARNING      +
C**  +    message, UNLESS a specific truncation was "announced"          +
C**  +    in var:                                                        +
C**  +       To truncate the first dim to Li, let var(Ni,1,1) = Li       +
C**  +       To truncate the 2nd   dim to Lj, let var(1,Nj,1) = Lj       +
C**  +       ... (this has no effect exept cancel the "WARNING" message) +
C**  +-------------------------------------------------------------------+

      SUBROUTINE UNwrite (FILEid , VARname , itime,
     &                    Ni,  Nj, Nlev, var)

      IMPLICIT NONE

      INCLUDE 'libUN.inc'

      INTEGER icheck

      INTEGER Lvnam 
      PARAMETER (Lvnam=20)

C     ** input 
      INTEGER        FILEid 
      INTEGER        itime
      INTEGER        Ni,  Nj, Nlev 
      CHARACTER *(*) VARname 
      REAL*4         var(Ni, Nj, Nlev)

C     ** local :
      INTEGER    MXlv
      PARAMETER (MXlv=500) 
C                ^^^^Maximal # levels for a special output
      INTEGER  VARSIZE
      EXTERNAL VARSIZE
      INTEGER NVRi,  NVRj, NVRlev
      INTEGER Ierro, TTerr, Nvatts, vtype
      INTEGER dimID(4), dimSIZ(4), count(4)    
      INTEGER start(4),stride(4),imap(4)
      CHARACTER*(Lvnam) dimNAM(4) 
      CHARACTER*(Lvnam) recname
      CHARACTER*(30) tmpchr
      INTEGER varVID
      INTEGER VNlen, NDIMvar, NSDIvar, tiDI, itmp
      INTEGER iz, ii, jj, ll
      INTEGER iUNLIMDIM
      REAL*4 chkdim
      REAL*4 Arange(2),sValRange(2)
      REAL*4 Srange(MXlv,2)
      LOGICAL OkRange
      
      icheck= 0     !** 'debugging' level
      TTerr = 0     !** 'total number of errors
 
      IF (icheck.ge.1) WRITE(*,*) 'UNwrite : Begin'

C*    1. Get the variable field  and dims IDs
C     ----------------------------------------

      IF (icheck.ge.2) WRITE(*,*) 'FILEid  :', FILEid 

C     ** getting VARname  size :
      VNlen = VARSIZE (VARname)
      IF (icheck.ge.3) WRITE(*,*) 'VNlen  :', VNlen
      IF (icheck.ge.2) WRITE(*,*) 'VARname   :', VARname (1:VNlen)

C     ** variable field ID :
      Ierro=NF_INQ_VARID (FILEid, VARname (1:VNlen), varVID)

C     ** Cancel writing if an error occured : variable undefined ?
      IF (Ierro.ne.0.and.icheck.ge.1) THEN
         WRITE(*,*) 'UNwrite  Info  : Variable ',VARname(1:VNlen)
     &            ,' not found -> not written.' 
      END IF
      IF (Ierro.ne.0) GOTO 9999 !** UNwrite_end


C     ** Inquire about the number of dimensions in var :
C     **
! FD debug
      recname=''
      Ierro=NF_INQ_VAR(FILEid , varVID, recname, vtype,
     &                   NDIMvar,  dimID,  Nvatts)
C     **  line1          id/file  id/var  var name var type
C     **  line2          # dims   id/dims #attributes
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwrite', Ierro)

      IF (icheck.ge.2) WRITE(*,*) 'Ierro1. ', Ierro


C*    2. Dimensions : inquire about file + compare with input data.
C     -------------------------------------------------------------

C     2.1 Inquire dimensions names and sizes :
C +   - - - - - - - - - - - - - - - - - - - - -
      DO iz = 1,4
        dimSIZ(iz)=0
        dimNAM(iz)='       '
C       ** Set any unused dimension to "0" size / no name
      END DO 
      DO iz = 1,NDIMvar
        Ierro=NF_INQ_DIM(FILEid , dimID(iz), dimNAM(iz), dimSIZ(iz))
C       **                 id/file  id/dim     dimname      dimsize    
C       **                                     !output      output
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwrite', Ierro)
      END DO
      IF (icheck.ge.3) WRITE(*,*) 'NDIMvar  ',NDIMvar
      IF (icheck.ge.3) WRITE(*,*) 'Ierro 2.0',Ierro  

C     2.2 Set writing region according to field dimension : 2D or 3D
C +   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     ** Set horizontal dimensions (default, for most data) :
      count(1) = Ni
      count(2) = Nj
C +   ** Other default values:
      count(3) = 0
      count(4) = 0
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1

C +- ------3D+time variable in file-----------
      IF (NDIMvar.eq.4) THEN
C       ** 3D space + time: 
        NSDIvar = 3     ! # space dims
        tiDI    = 4     ! No. of the time dim
C       ** write 3D space:
        start(3) = 1    ! Start of index 3 in var (here = vert. levs)
        count(3) = Nlev ! # values of index 3 in var 
C       ** write one time step:
        start(4) = itime
        count(4) = 1
C +- ------3D *OR* 2D+time var in file--------
      ELSE IF (NDIMvar.eq.3) THEN
        IF (Nlev.EQ.1) THEN
C         ** 2D space + time (standard use of UNlib):
          NSDIvar = 2
          tiDI    = 3
C         ** ...write one time step:
          start(3) = itime
          count(3) = 1     
        ELSE
C         ** 3D (no time slice):
          NSDIvar = 3
          tiDI    = 0
C         ** ...write 3rd dimension:
          start(3) = 1    
          count(3) = Nlev
        ENDIF
C +- ------2D *OR* 1D+time var in file--------
      ELSE IF (NDIMvar.eq.2) THEN
        IF (Nj.EQ.1 .AND. dimNAM(2)(1:4).EQ.'time') THEN
C         ** Write a 1D vector at time= itime:
          NSDIvar = 1
          tiDI    = 2
          start(2) = itime
          count(2) = 1
        ELSE
C         ** Usual MAR 2D space (no time):
          NSDIvar = 2
          tiDI    = 0
        END IF
C +- ------1D *OR* 0D+time var in file--------
      ELSE IF (NDIMvar.eq.1) THEN
C       ** 1D space or time
        IF (Ni.eq.1) THEN
C         ** Write a single element (at itime)
          start(1) = itime
          count(1) = 1
          count(2) = 0
          NSDIvar = 0
          tiDI    = 1
        ELSE
C         ** Write a vector (use only "space" dim 1)
          NSDIvar = 1
          tiDI    = 0
          count(2)= 0
        END IF
      ELSE
         WRITE(*,*) 'UNwrite ERROR : data field dimension ?'
         STOP
      END IF

C     2.3 Compare file dimensions to input data.
C +   - - - - - - - - - - - - - - - - - - - - - -
C     ** Save variable size for use as "valid" size (-> range):
      NVRi   = Ni
      NVRj   = Nj
      NVRlev = Nlev
C     ** Space dimensions :
      IF (NSDIvar.GT.0) THEN
      DO iz = 1,NSDIvar
        IF      (dimSIZ(iz).gt.count(iz)) THEN
          write(*,*) 'UNwrite - WARNING: '
          write(*,*) ' Your field ',VARname,' has an empty part.'
          write(*,*) ' (for the dimension:',dimNAM(iz),')'
        ELSE IF (dimSIZ(iz).lt.count(iz)) THEN
C         ** Do display "warning" only if truncation
C            was not "correctly announced" (see header)
C            (NVR... => stop here when updating the range attribute)
          IF (iz.EQ.1) THEN 
            chkdim = var(Ni,1,1) 
            NVRi   = dimSIZ(1) 
          ELSE IF (iz.EQ.2) THEN 
            chkdim = var(1,Nj,1)
            NVRj   = dimSIZ(2)
          ELSE IF (iz.EQ.3) THEN 
            chkdim = var(1,1,Nlev)
            NVRlev = dimSIZ(3)
          ELSE  
            chkdim = 0.0
          ENDIF
          Ierro= NF_INQ_UNLIMDIM (FILEid, iUNLIMDIM) 
          IF (dimID(iz).NE.iUNLIMDIM) THEN
           IF (ABS(chkdim-dimSIZ(iz)).GT. 0.1 ) THEN
            write(*,*) 'UNwrite - WARNING: '
            write(*,*) ' Your field ',VARname,' will be truncated.'
            write(*,*) ' (for the dimension:',dimNAM(iz),')'
           ENDIF
           count(iz) = dimSIZ(iz)
          ENDIF
        END IF
      END DO
      END IF

C     ** Time dimension (when defined):
      IF (tiDI.ne.0) THEN
       IF (itime.gt.dimSIZ(tiDI)) THEN
         IF (icheck.ge.1) WRITE(*,*) 'Time limit, ID', dimID(tiDI) 
         Ierro= NF_INQ_UNLIMDIM (FILEid, iUNLIMDIM) 
         IF (dimID(tiDI).NE.iUNLIMDIM) THEN
            WRITE(*,*) 'UNwrite - ERROR:   '
            WRITE(*,*) ' Time index out of range '                        
            STOP
         ENDIF
        END IF
      END IF

      IF (icheck.ge.2) WRITE(*,*) 'Ierro2. ', Ierro
      IF (icheck.ge.2) WRITE(*,*) 'Dimension names :',dimNAM
      IF (icheck.ge.2) WRITE(*,*) 'dimSIZ :',dimSIZ
      IF (icheck.ge.2) WRITE(*,*) 'count  :',count
      IF (icheck.ge.2) WRITE(*,*) 'start  :',start
      IF (icheck.ge.2) WRITE(*,*) 'dimID  :',dimID 

C*    3. Write variable.
C     ------------------

C     ** Set 'imap' and WRITE with NCVPTG:
C     ** NOTE : since the arrays (grid_*) may be over-dimensionned,
C     **        we use the 'generalised' writing routine NCVPTG
C     ** (imap tells NetCDF about the memory locations of var)
      imap(1) = 1
      imap(2) = imap(1) * Ni      ! 1st dim of var = Ni 
      imap(3) = imap(2) * Nj      ! 2nd dim of var = Nj 
      imap(4) = 0                 ! (not used: 0 or 1 time step)   
      DO iz=1,4
        stride(iz)=1
      END DO
C     ** NOTE: stride is not used.

      Ierro=NF_PUT_VARM_REAL(FILEid , varVID  , start      , count,
     &                         stride , imap    , var(1,1,1) )
C     **  line1:              id/file | id/var  |read from...|#data
C     **  line2:              step    |re-arrang|variable(beg.)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwrite', Ierro)

      IF (icheck.ge.2) WRITE(*,*) 'Ierro3.2', Ierro

C*    4a. Update 'actual_range' attribute.               
C     ------------------------------------

C     If 'actual_range' available, get its current value:
C     - - - - - - - - - - - - - - - - - - - - - - - - - -

C     ** Get the old min and max values:
      Ierro=NF_GET_ATT_REAL(FILEid ,varVID ,'actual_range' ,
     &                        Arange )
c     **line1 ^^^^^^^^^^^^^^  FILEid |var.id | attr.name
C     **line2                 value

C     ** Cancel if an error occured : attribute undefined ?
      IF (Ierro.ne.0.and.icheck.ge.1) THEN
         WRITE(*,*) 'UNwrite  Info : attribute actual_range ' 
     &             ,' not found -> not written.'
      END IF
      IF (Ierro.ne.0) GOTO 9990 !** Next section

C     If 'valid_range' available, get its current value:
C     - - - - - - - - - - - - - - - - - - - - - - - - - -

C     ** Get the min/max valid range (outside = missing val):
      Ierro=NF_GET_ATT_REAL(FILEid ,varVID ,'valid_range' ,
     &                        sValRange)
      IF (Ierro.ne.0) THEN
         sValRange(1)=ValRange(1)
         sValRange(2)=ValRange(2)
      END IF

C     Update the min an max
C     - - - - - - - - - - - 

C     **If this is the first pass, initialise min and max:
      IF (      Arange(1).EQ. NF_FILL_REAL 
     .    .OR. (Arange(1).EQ. 0.0 .AND. Arange(2).EQ. 0.0) ) THEN
        OkRange = .false. 
      ELSE
        OkRange = .true.
      ENDIF

      DO ll=1, NVRlev
      DO jj=1, NVRj
      DO ii=1, NVRi 
        IF (  var(ii,jj,ll).GE.sValRange(1)
     &  .AND. var(ii,jj,ll).LE.sValRange(2)) THEN
           IF (OkRange) THEN
              Arange(1) = MIN(Arange(1), var(ii,jj,ll))
              Arange(2) = MAX(Arange(2), var(ii,jj,ll))
           ELSE        
              Arange(1) = var(ii,jj,ll)
              Arange(2) = var(ii,jj,ll)
              OkRange = .true.
           ENDIF
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      IF (icheck.ge.2) WRITE(*,*) 'Arange',Arange

C     Set attribute.
C     - - - - - - - -

      Ierro=NF_PUT_ATT_REAL(FILEid  ,varVID ,'actual_range' ,
     &                        NF_FLOAT,2      ,Arange)
c     **line1 ^^^^^^^^^^^^^^^ FILEid  |var.id | attr.name
C     **line2                 type    |len    | attr.value
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwrite', Ierro)
      TTerr = TTerr + ABS(Ierro)

C     ** Next section:
 9990 CONTINUE

C*    5. Update the optional '[var]_range' special variable.
C     ------------------------------------------------------
      IF (NDIMvar.eq.4.and.Nlev.lt.MXlv) THEN

C     If '[var]_range' available, get its current value:
C     - - - - - - - - - - - - - - - - - - - - - - - - - -

C     ** Get ID of variable [var]_range :
      tmpchr = VARname(1:VNlen)//'_range'
      itmp   = VNlen + 6
      Ierro=NF_INQ_VARID(FILEid, tmpchr(1:itmp), varVID)

C     ** Cancel if an error occured : undefined ?
      IF (Ierro.ne.0.and.icheck.ge.1) THEN
         WRITE(*,*) 'UNwrite  Info : [var]_range '
     &            ,' not found -> not written.'
      END IF
      IF (Ierro.ne.0) GOTO 9999 !** UNwrite_end

C     ** Get the old min and max values:
C     ** NOTE :
C     **        we use the 'generalised' reading routine NCVGTG
C     ** (imap tells NetCDF about the memory locations of var)
      imap(1) = 1
      imap(2) = imap(1) * MXlv   
      start(1)= 1
      start(2)= 1
      count(1)= Nlev
      count(2)= 2

C     ** (See UNread for explanations about NCVGTG)
      Ierro=NF_GET_VARM_REAL(FILEid, varVID, start, count,   
     &                         stride,  imap , Srange(1,1) )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwrite', Ierro)

C     Update the min an max
C     - - - - - - - - - - -
C     **If this is the first pass, initialise min and max:
C     **(Constant fields shall not be accounted for)
      DO ll=1, Nlev
        IF (Srange(ll,1).eq.Srange(ll,2)) THEN
          Srange(ll,1) = var(1,1,ll)
          Srange(ll,2) = var(1,1,ll) 
        ENDIF
      ENDDO

      DO jj=1, NVRj
      DO ii=1, NVRi
       DO ll=1, NVRlev
        Srange(ll,1) = MIN(Srange(ll,1), var(ii,jj,ll))
        Srange(ll,2) = MAX(Srange(ll,2), var(ii,jj,ll))
       ENDDO
      ENDDO
      ENDDO
      IF (icheck.ge.4) WRITE(*,*) 'Srange',Srange


C     Set special variable [var]_range
C     - - - - - - - - - - - - - - - - -
C     **(See UNread for explanations abtout NCVPTG)

      Ierro=NF_PUT_VARM_REAL(FILEid , varVID , start, count,
     &                         stride , imap   , Srange(1,1) )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwrite', Ierro)

      ENDIF  ! End Section 5.

C     UNwrite_end
C     -----------
      IF (icheck.ge.2) WRITE(*,*) 'Errors count:',TTerr
      IF (icheck.ge.2) WRITE(*,*) 'UNwrite : End'
 9999 CONTINUE
      RETURN
      END
C**
C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNlwrite :  +                                         +
C**  +-------------------------+                                         +
C**  +  * Writes a 2D horizontal LEVEL into a 3D+time NetCDF variable    +  
C**  +       OR  a 1D vector           into a 2D+time                    +
C**  +             --            ----         --                         +
C**  +    (SEE ALSO : UNwrite, for all dimensions - this a pecular case  +
C**  +     Note: 1D vectors are writen in the 1st dim of 2D+time)        +
C**  +                                                                   +
C**  +  * Automatically updates attribute 'actual_range' if available    +
C**  +          "          "    special var. '[var]_range'    "          +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEid  : input file identifier (from UNcreate OR NetCDF open) +
C**  +    VARname : name given to the variable to write (must be in file)+
C**  +    itime   : No of time step to write to                          +
C**  +    level   : No of level     to write to                          +
C**  +    Ni,  Nj : dimensions of 'var'...                               +
C**  +    var     : A 2D variable to be writen                           +
C**  +-------------------------------------------------------------------+

      SUBROUTINE UNlwrite (FILEid , VARname , itime,
     &                     ilev, Ni,  Nj, var)

      IMPLICIT NONE

      INCLUDE 'libUN.inc'

      INTEGER icheck

      INTEGER Lvnam 
      PARAMETER (Lvnam=20)

C     ** input 
      INTEGER FILEid 
      INTEGER itime, ilev
      INTEGER Ni,  Nj
      CHARACTER *(*) VARname 
      REAL*4 var(Ni, Nj)

C     ** local :
      INTEGER  VARSIZE
      EXTERNAL VARSIZE
      INTEGER Ierro, TTerr, Nvatts, vtype
      INTEGER dimID(4), dimSIZ(4), count(4)    
      INTEGER start(4),stride(4),imap(4)
      INTEGER iUNLIMDIM
      CHARACTER*(Lvnam) dimNAM(4) 
      CHARACTER*(Lvnam) recname
      CHARACTER*(30) tmpchr
      INTEGER varVID
      INTEGER VNlen, NDIMvar, NSDIvar, tiDI, ilDI, itmp
      INTEGER iz, ii, jj
      LOGICAL OkRange
      REAL*4 Arange(2), sValRange(2)
      REAL*4 Srange(2)
            
      icheck= 0     !** 'debugging' level
      TTerr = 0     !** 'total numbe of errors
 
      IF (icheck.ge.1) WRITE(*,*) 'UNlwrite : Begin'

C*    1. Get the variable field  and dims IDs
C     ----------------------------------------

      IF (icheck.ge.2) WRITE(*,*) 'FILEid  :', FILEid 

C     ** getting VARname  size :
      VNlen = VARSIZE (VARname)
      IF (icheck.ge.3) WRITE(*,*) 'VNlen  :',VNlen
      IF (icheck.ge.2) WRITE(*,*) 'VARname   :', VARname (1:VNlen)

C     ** variable field ID :
      Ierro=NF_INQ_VARID (FILEid, VARname (1:VNlen), varVID)

C     ** Cancel writing if an error occured : variable undefined ?
      IF (Ierro.ne.0.and.icheck.ge.1) THEN
         WRITE(*,*) 'UNlwrite  Info  : Variable ',VARname(1:VNlen)
     &            ,' not found -> not written.' 
      END IF
      IF (Ierro.ne.0) GOTO 9999 !** UNlwrite_end


C     ** Inquire about the number of dimensions in var :
C     **
      Ierro=NF_INQ_VAR(FILEid , varVID, recname, vtype,
     &                   NDIMvar,  dimID,  Nvatts)
C     **  line1          id/file  id/var  var name var type
C     **  line2          # dims   id/dims #attributes
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNlwrite', Ierro)

      IF (icheck.ge.2) WRITE(*,*) 'Ierro1. ', Ierro


C*    2. Dimensions : inquire about file + compare with input data.
C     -------------------------------------------------------------

C     2.1 Inquire dimensions names and sizes :
C +   - - - - - - - - - - - - - - - - - - - - -
      DO iz = 1,4
        dimSIZ(iz)=0
        dimNAM(iz)='       '
C       ** Set any unused dimension to "0" size / no name
      END DO

      DO iz = 1,NDIMvar
        Ierro=NF_INQ_DIM(FILEid , dimID(iz), dimNAM(iz), dimSIZ(iz))
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNlwrite', Ierro)
C       **           id/file   id/dim    dimname     dimsize    error
C       **                               !output     output
      END DO
      IF (icheck.ge.3) WRITE(*,*) 'NDIMvar  ',NDIMvar
      IF (icheck.ge.3) WRITE(*,*) 'Ierro 2.0',Ierro  

C     2.2 Set writing region according to field dimension :  3D
C +   - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     ** Set horizontal dimensions (all field dims):
      count(1) = Ni
      count(2) = Nj
      start(1) = 1
      start(2) = 1
C +- ------ 3D+time var in file--------
      IF (NDIMvar.eq.4) THEN
        NSDIvar = 2     ! # input space dims (for a 2D level)
        tiDI    = 4     ! No. of the time dim
C       ** write one level (set the level No) :
        start(3) = ilev ! Start of index 3 in var 
        count(3) = 1    ! # values of index 3 in var 
        ilDI     = 3
C       ** write one time step:
        start(4) = itime
        count(4) = 1
C +- ------ 2D+time var in file--------
      ELSE IF (NDIMvar.eq.3) THEN
        NSDIvar = 1     ! # input space dims (for a 1D vector)
        tiDI    = 3     ! No. of the time dim
C       ** write one "level" - here a 1D vector in the 1st dim.
        start(2) = ilev ! Start of index 2 in var
        count(2) = 1    ! # values of index 3 in var
        ilDI     = 2
C       ** write one time step:
        start(3) = itime
        count(3) = 1
      ELSE
         WRITE(*,*) 'UNlwrite ERROR : data field dimension ?'
         WRITE(*,*) '  NB: UNlwrite = only for (2 or) 3D +time.'
         STOP
      END IF

C     2.3 Compare file dimensions to input data.
C +   - - - - - - - - - - - - - - - - - - - - - -
C     ** Space dimensions :
      DO iz = 1,NSDIvar
        IF      (dimSIZ(iz).gt.count(iz)) THEN
          write(*,*) 'UNlwrite - WARNING: '
          write(*,*) ' Your field ',VARname,' has an empty part.'
          write(*,*) ' (for the dimension:',dimNAM(iz),')'
        ELSE IF (dimSIZ(iz).lt.count(iz)) THEN
          write(*,*) 'UNlwrite - WARNING: '
          write(*,*) ' Your field ',VARname,' will be truncated.'
          write(*,*) ' (for the dimension:',dimNAM(iz),')'
          count(iz) = dimSIZ(iz)
        END IF
      END DO

C     ** Space dimensions - check if requested level exists:
      IF (dimSIZ(ilDI).lt.ilev) THEN
        write(*,*) 'UNlwrite - ERROR: '
        write(*,*) ' The requested level =',ilev
        write(*,*) ' does not exist in the field ',VARname
        write(*,*) ' (for the dimension:',dimNAM(ilDI),')'
        STOP
      END IF

C     ** Time dimension (when defined):
      IF (tiDI.ne.0) THEN
       IF (itime.gt.dimSIZ(tiDI)) THEN
         IF (icheck.ge.1) WRITE(*,*) 'Time limit, ID', dimID(tiDI) 
         Ierro= NF_INQ_UNLIMDIM (FILEid, iUNLIMDIM) 
         IF (dimID(tiDI).NE.iUNLIMDIM) THEN
            WRITE(*,*) 'UNlwrite - ERROR:  '
            WRITE(*,*) ' Time index out of range '                        
            STOP
         ENDIF
        END IF
      END IF

      IF (icheck.ge.2) WRITE(*,*) 'Ierro2. ', Ierro
      IF (icheck.ge.2) WRITE(*,*) 'Dimension names :',dimNAM
      IF (icheck.ge.3) WRITE(*,*) 'dimSIZ :',dimSIZ
      IF (icheck.ge.3) WRITE(*,*) 'count  :',count
      IF (icheck.ge.3) WRITE(*,*) 'start  :',start
      IF (icheck.ge.3) WRITE(*,*) 'dimID  :',dimID 

C*    3. Write variable.
C     ------------------

C     ** Set 'imap' and WRITE with NCVPTG:
C     ** NOTE : since the arrays (grid_*) may be over-dimensionned,
C     **        we use the 'generalised' writing routine NCVPTG
C     ** (imap tells NetCDF about the memory locations of var)
      imap(1) = 1
      imap(2) = imap(1) * Ni      ! 1st dim of var = Ni 
      imap(3) = imap(2) * Nj      ! (not used: 1 level...)
      imap(4) = 0                 ! (not used: 0 or 1 time step)   
      DO iz=1,4
        stride(iz)=1
      END DO
C     ** NOTE: stride is not used.

      Ierro=NF_PUT_VARM_REAL (FILEid  , varVID  , start      , count,
     &                          stride  , imap    , var(1,1)          )
C     **  line1:                id/file | id/var  |read from...|#data
C     **  line2:                step    |re-arrang|variable(beg.)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNlwrite', Ierro)

      IF (icheck.ge.2) WRITE(*,*) 'Ierro3.2', Ierro

C*    4a. Update 'actual_range' attribute.               
C     ------------------------------------

C     If 'actual_range' available, get its current value:
C     - - - - - - - - - - - - - - - - - - - - - - - - - -

C     ** Get the old min and max values:
      Ierro=NF_GET_ATT_REAL(FILEid ,varVID ,'actual_range' ,
     &                        Arange )
c     **line1 ^^^^^^^^^^^^^^^ FILEid |var.id | attr.name
C     **line2                 value

C     ** Cancel if an error occured : attribute undefined ?
      IF (Ierro.ne.0.and.icheck.ge.1) THEN
         WRITE(*,*) 'UNlwrite  Info : attribute actual_range ' 
     &             ,' not found -> not written.'
      END IF
      IF (Ierro.ne.0) GOTO 9990 !** Next section

C     If 'valid_range' available, get its current value:
C     - - - - - - - - - - - - - - - - - - - - - - - - - -

C     ** Get the min/max valid range (outside = missing val):
      Ierro=NF_GET_ATT_REAL(FILEid ,varVID ,'valid_range' ,
     &                        sValRange)
      IF (Ierro.ne.0) THEN
         sValRange(1)=ValRange(1)
         sValRange(1)=ValRange(2)
      END IF

C     Update the min an max
C     - - - - - - - - - - - 

C     **If this is the first pass, initialise min and max:
      IF (      Arange(1).EQ. NF_FILL_REAL 
     .    .OR. (Arange(1).EQ. 0.0 .AND. Arange(2).EQ. 0.0) ) THEN
        OkRange = .false. 
      ELSE
        OkRange = .true.
      ENDIF

      DO jj=1, Nj
      DO ii=1, Ni
        IF (  var(ii,jj).GE.sValRange(1)
     &  .AND. var(ii,jj).LE.sValRange(2)) THEN
           IF (OkRange) THEN
              Arange(1) = MIN(Arange(1), var(ii,jj))
              Arange(2) = MAX(Arange(2), var(ii,jj))
           ELSE        
              Arange(1) = var(ii,jj)
              Arange(2) = var(ii,jj)
              OkRange = .true.
           ENDIF
        ENDIF
      ENDDO
      ENDDO
      IF (icheck.ge.2) WRITE(*,*) 'Arange',Arange

C     Set attribute.
C     - - - - - - - -

      Ierro=NF_PUT_ATT_REAL(FILEid  ,varVID ,'actual_range' ,
     &                        NF_FLOAT,2      ,Arange )
c     **line1 ^^^^^^^^^^^^^^^ FILEid  |var.id | attr.name
C     **line2                 type    |len    | attr.value
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNlwrite', Ierro)
      TTerr = TTerr + ABS(Ierro)

C     ** Next section:
 9990 CONTINUE
 

C*    5. Update the optional '[var]_range' special variable.
C     ------------------------------------------------------
      IF (NDIMvar.eq.4) THEN

C     If '[var]_range' available, get its current value:
C     - - - - - - - - - - - - - - - - - - - - - - - - - -

C     ** Get ID of variable [var]_range :
      tmpchr = VARname(1:VNlen)//'_range'
      itmp   = VNlen + 6
      Ierro=NF_INQ_VARID (FILEid, tmpchr(1:itmp), varVID)

C     ** Cancel if an error occured : undefined ?
      IF (Ierro.ne.0.and.icheck.ge.1) THEN
         WRITE(*,*) 'UNlwrite  Info : [var]_range '
     &             ,' not found -> not written.'
      END IF
      IF (Ierro.ne.0) GOTO 9999 !** UNlwrite_end

C     ** Get the old min and max values:
C     ** NOTE :
C     **        we use the 'generalised' reading routine NCVGTG
C     ** (imap tells NetCDF about the memory locations of var)
      imap(1) = 1
      imap(2) = 0                ! Not used (write only 1 lev)
      start(1)= ilev
      count(1)= 1   
      start(2)= 1
      count(2)= 2

C     ** (See UNread for explanations abtout NCVGTG)
      Ierro=NF_GET_VARM_REAL(FILEid, varVID, start      ,count,   
     &                         stride,  imap , Srange(1) )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNlwrite', Ierro)

C     Update the min an max
C     - - - - - - - - - - -
C     **If this is the first pass, initialise min and max:
C     **(Constant fields shall not be accounted for)
      IF (Srange(1).eq.Srange(2)) THEN
          Srange(1) = var(1,1)
          Srange(2) = var(1,1) 
      ENDIF

      DO jj=1, Nj
      DO ii=1, Ni
        Srange(1) = MIN(Srange(1), var(ii,jj))
        Srange(2) = MAX(Srange(2), var(ii,jj))
      ENDDO
      ENDDO
      IF (icheck.ge.4) WRITE(*,*) 'Srange',Srange


C     Set special variable [var]_range
C     - - - - - - - - - - - - - - - - -
C     **(See UNread for explanations abtout NCVPTG)

      Ierro=NF_PUT_VARM_REAL(FILEid , varVID , start        , count,
     &                         stride , imap   , Srange(1)  )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNlwrite', Ierro)

      ENDIF  ! End Section 5.

C     UNlwrite_end
C     -----------
      IF (icheck.ge.2) WRITE(*,*) 'Errors count:',TTerr
      IF (icheck.ge.2) WRITE(*,*) 'UNlwrite : End'
 9999 CONTINUE
      RETURN
      END
C**
C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNread :    +                                         +
C**  +-------------------------+                                         +
C**  +  * Reads a model variable from a NetCDF file,                     +
C**  +    and reads the coordinates of the grid upon wich it is defined. +
C**  +    (the NetCDF file must have been opened and must be closed      +
C**  +     after all reading operations). May read an x-y subregion.     +
C**  +                                                                   +
C**  +  INPUT :                                                          + 
C**  +    FILEid  : input file identifier (from NetCDF open)             +
C**  +    VARname  : name of the requested variable.                     +
C**  +    time : [integer*4] is the time index of the data field to read +
C**  +    level: [integer*4] (usefull for 3D-space fields only) :        +
C**  +                       if not=0 --> = no of the level              +
C**  +                                      -> output is 2D (l_dim = 1)  +
C**  +                       if  =0   --> read ALL levels                +
C**  +                                      -> output is 3D              +
C**  +    i_dbeg, j_dbeg      : horizontal indexes of requested region   +
C**  +                          in input data file                       + 
C**  +    i_dim, j_dim, l_dim : ...the dimensions of 'var',              + 
C**  +                       = the dimensions of the sub-region to read  + 
C**  +                       ! l_dim = 1 if level not=0                  + 
C**  +                       ! j_dim = 1 if var is 1D                    + 
C**  +  OUTPUT :                                                         + 
C**  +    varax1[i_dim] (real  )                                         + 
C**  +    varax2[j_dim]: Horizontal coordinates in the file (lat/lon,...)+ 
C**  +    varlev[l_dim]: vertical coordinate of the levels               + 
C**  +                   (! when level not=0, only varlev(1) is defined) + 
C**  +    var_units                 : physical units of var.             + 
C**  +    var[i_dim,j_dim,l_dim]    :                                    + 
C**  +                            data field values                      + 
C**  +                            (var must be defined, and is REAL  )   + 
C**  +                                                                   + 
C**  +-------------------------------------------------------------------+

      SUBROUTINE UNread
     &      (FILEid , VARname , time, level, i_dbeg, j_dbeg,
     &       i_dim   , j_dim   , l_dim    ,
     &       varax1  , varax2  , varlev,      
     &       var_units, var)

      IMPLICIT NONE
      INCLUDE 'libUN.inc'

      INTEGER icheck

      INTEGER Lvnam 
      PARAMETER (Lvnam=40)

C     ** input 
      INTEGER FILEid 
      INTEGER time, level, i_dbeg, j_dbeg
      INTEGER i_dim, j_dim, l_dim
      CHARACTER *(*) VARname 

C     ** output
      REAL*4    varax1(i_dim), varax2(j_dim), varlev(l_dim)
      CHARACTER *(*) var_units
      REAL*4    var (i_dim, j_dim, l_dim)

C     ** local :
      INTEGER  VARSIZE
      EXTERNAL VARSIZE
      REAL*4  varmin,varmax
      INTEGER Ierro, Nvatts, vtype
      INTEGER dimID(4), dimSIZ(4), dimREG(4)    
      INTEGER start(4),begREG(4),count(4),stride(4),imap(4)
      CHARACTER*(Lvnam) dimNAM(4) 
      CHARACTER*(Lvnam) dNAMver, dNAMtim
      CHARACTER*(Lvnam) recname
      CHARACTER*(10) Routine
      INTEGER ax1VID, ax2VID, verVID, timVID, varVID
      INTEGER VNlen, varNUMDIM
      INTEGER ii,jj,ll,z
      
      icheck= 0
C*    0. Initialisations
C     ------------------
      Routine= 'UNread'
      IF (icheck.ge.1) WRITE(*,*) 'UNread : Begin'

      DO ii = 1,4 
        stride(ii) = 1
        begREG(ii) = 1
        start (ii) = 1
      ENDDO
 
C*    1. Get the variable field  and dims IDs
C     ----------------------------------------

      IF (icheck.ge.3) WRITE(*,*) 'FILEid  :', FILEid 

C     ** getting VARname  size :
      VNlen = VARSIZE(VARname)
      IF (icheck.ge.3) WRITE(*,*) 'VNlen  :',VNlen
      IF (icheck.ge.2) WRITE(*,*) 'VARname   :', VARname (1:VNlen)

C     ** variable field ID :
      Ierro=NF_INQ_VARID (FILEid, VARname (1:VNlen), varVID)

C*    1b. Handle non-existing variables
C     ---------------------------------
      IF (Ierro.NE.NF_NOERR) THEN 
         IF (Ierro.EQ.NF_ENOTVAR .AND. iVarWarn.LE.1) THEN
            IF (iVarWarn.EQ.1) THEN
              write(*,*) 'WARNING (UNsread): variable not found:'
              write(*,*) '     ',varName
            ENDIF
            DO ll=1,l_dim
            DO jj=1,j_dim
            DO ii=1,i_dim
              var (ii,jj,ll)=VarRepl
            ENDDO
            ENDDO
            ENDDO
            RETURN  ! EXIT SUBROUTINE, read nothing
         ENDIF
         WRITE(*,*) 'Error reading variable: ', VARname(1:VNlen)
         CALL HANDLE_ERR('UNsread',Ierro)
      ENDIF

C     1c. Inquire about the number of dimensions in var
C     -------------------------------------------------

! FD debug
      recname=''
      Ierro=NF_INQ_VAR(FILEid   , varVID, recname, vtype,
     &                   varNUMDIM, dimID, Nvatts )
C     **  line1          id/file    id/var  var name  var type
C     **  line2          # dims    id/dims #attributes
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNsread', Ierro)

      IF (icheck.ge.3) WRITE(*,*) 'Ierro1. ', Ierro

C*    2. Dimensions : in the reading region and in the file.
C     ------------------------------------------------------

C     ** inquire dimensions names and sizes :
      DO z = 1,varNUMDIM
! FD debug
        dimNAM(z)=''
        Ierro=NF_INQ_DIM(FILEid , dimID(z), dimNAM(z), dimSIZ(z))
C       **                 id/file  id/dim    dimname    dimsize
C       **                                    !output    output
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR(Routine, Ierro)
      END DO

C     ** In this version, we read only a xy subregion of the file :
      dimREG(1) = i_dim
      dimREG(2) = j_dim
      begREG(1) = i_dbeg
      begREG(2) = j_dbeg
      IF (begREG(1).lt.1)  begREG(1) = 1
      IF (begREG(2).lt.1)  begREG(2) = 1
      
C     ** Set reading region according to field dimension : 2D or 3D
      IF (varNUMDIM.eq.4) THEN
C       ** for 3D fields :
        IF (level.gt.0) THEN
C       ** one level is read :
          dimREG(3) = 1 
          begREG(3) = level
          dNAMver   = dimNAM(3)
        ELSE
C       ** all levels are read :
          dimREG(3) = l_dim
          begREG(3) = 1    
          dNAMver   = dimNAM(3)
        END IF
C       ** one time step is read:
        dimREG(4) = 1 
        begREG(4)  = time 
        dNAMtim   = dimNAM(4)
      ELSE IF (varNUMDIM.eq.3) THEN
C       ** for 2D space fields + time:
C       ** one time step is read:
        dimREG(3) = 1      
        begREG(3) = time
        dNAMtim   = dimNAM(3)
        dimREG(4) = 0
        begREG(4) = 0
        dimNAM(4) = 'none'
      ELSE IF (varNUMDIM.eq.2) THEN
C       ** for 2D fields :
C       ** no time step is read:
        dimREG(3) = 0     
        begREG(3) = 0   
        dNAMtim   = 'none'   
        dimNAM(3) = 'none'
        dimREG(4) = 0
        begREG(4) = 0
        dimNAM(4) = 'none'
      ELSE IF (varNUMDIM.eq.1) THEN
C       ** for 1D variable :
C       ** not assumed to be on a XYZ grid,       
C       ** just read a vector  
        dimREG(1) = 1    ! this was added by Martin Vancop for 1d vectors
        begREG(1) = time ! this was added by Martin Vancop for 1d vectors
        dimREG(2) = 0
        begREG(2) = 0
        dimNAM(2) = 'none'
        dimREG(3) = 0
        begREG(3) = 0
        dimNAM(3) = 'none'
        dNAMtim   = 'none'
        dimREG(4) = 0
        begREG(4) = 0
        dimNAM(4) = 'none'
      ELSE
        WRITE(*,*) 'UNread ERROR : data field dimension ?'
        STOP
      END IF

      DO z = 1,varNUMDIM
        IF (begREG(z).gt.dimSIZ(z)) THEN
          write(*,*) 'UNread - ERROR   : requested area out      '
          write(*,*) '                   of file area.          '
          write(*,*) '  (for the dimension:' , dimNAM(z) , ')'
          STOP
        END IF
        IF (dimSIZ(z).lt.(dimREG(z)+begREG(z)- 1) ) THEN
          write(*,*) 'UNread - WARNING : empty portion in field, '
          write(*,*) '  requested region > file contents       '
          write(*,*) '  (for the dimension:' , dimNAM(z) , ')'
          dimREG(z) = dimSIZ(z) - begREG(z) + 1
        END IF
      END DO

      IF (icheck.ge.3) WRITE(*,*) 'Ierro2. ', Ierro
      IF (icheck.ge.2) WRITE(*,*) 'Dimension names :',dimNAM
      IF (icheck.ge.2) WRITE(*,*) 'dimSIZ :',dimSIZ
      IF (icheck.ge.2) WRITE(*,*) 'dimREG :',dimREG
      IF (icheck.ge.2) WRITE(*,*) 'begREG :',begREG
      IF (icheck.ge.3) WRITE(*,*) 'dimID  :',dimID 

C*    3. Get the variables IDs for the grid points locations. 
C     -------------------------------------------------------

      IF (varNUMDIM.ge.2) THEN
        Ierro=NF_INQ_VARID (FILEid, dimNAM(1), ax1VID)
        IF (Ierro.NE.NF_NOERR) THEN
          IF (Ierro.EQ.NF_ENOTVAR) THEN
            WRITE(*,*) 'Coordinate values not found:',dimNAM(1)
          ENDIF
          CALL HANDLE_ERR(Routine, Ierro)
        ENDIF
        Ierro=NF_INQ_VARID (FILEid, dimNAM(2), ax2VID)
        IF (Ierro.NE.NF_NOERR) THEN
          IF (Ierro.EQ.NF_ENOTVAR) THEN
            WRITE(*,*) 'Coordinate values not found:',dimNAM(2)
          ENDIF
          CALL HANDLE_ERR(Routine, Ierro)
        ENDIF
      ENDIF
      IF (varNUMDIM.ge.3) THEN
        Ierro=NF_INQ_VARID (FILEid, dNAMtim, timVID)
        IF (Ierro.NE.NF_NOERR) THEN
          IF (Ierro.EQ.NF_ENOTVAR) THEN
            WRITE(*,*) 'Coordinate values not found:',dNAMtim
          ENDIF
          CALL HANDLE_ERR(Routine, Ierro)
        ENDIF
      END IF
      IF (varNUMDIM.eq.4) THEN
        Ierro=NF_INQ_VARID (FILEid, dNAMver, verVID)
        IF (Ierro.NE.NF_NOERR) THEN
          IF (Ierro.EQ.NF_ENOTVAR) THEN
            WRITE(*,*) 'Coordinate values not found:',dNAMver
          ENDIF
          CALL HANDLE_ERR(Routine, Ierro)
        ENDIF
      END IF
C     **                      id/file  name    id/var

      IF (icheck.ge.3) WRITE(*,*) 'Ierro3. ', Ierro

C*    4. Get attributes.         
C     ------------------

      IF (varNUMDIM.ge.2) THEN   !Not for 1D vectors (special case)
C       ** units attribute 
        Ierro=NF_GET_ATT_TEXT (FILEid , varVID, 'units', 
     &                           var_units) 
        IF (Ierro.NE.NF_NOERR) THEN 
          IF (Ierro.EQ.NF_ENOTATT) THEN
            write(*,*) 'Note (UNread): units not found for'
            write(*,*) '     ',varName
            var_units=' '
          ELSE
            CALL HANDLE_ERR('UNread',Ierro)
          ENDIF
        ENDIF

        IF (icheck.ge.2) WRITE(*,*) 'var_units :', var_units
      ENDIF

C*    5. Get values.
C     --------------
C*    5.1 ...for the grid points locations.
C     -------------------------------------
      
C     ** Horizontal : always read, except for 1D vectors
      IF (varNUMDIM.ge.2) THEN  
        count(1)=dimREG(1)
        start(1)=begREG(1)
        Ierro=NF_GET_VARA_REAL(FILEid ,ax1VID,start,count,varax1)
C       **                       id/file id/var from  #data data
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR(Routine, Ierro)
        count(1)=dimREG(2)
        start(1)=begREG(2)
        Ierro=NF_GET_VARA_REAL(FILEid ,ax2VID,start,count,varax2)
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR(Routine, Ierro)
      ENDIF

C     ** vertical :  only for 3D fields.
      IF (varNUMDIM.eq.4) THEN
        start(1) =begREG(3)
        count(1) =dimREG(3)
        Ierro =  NF_GET_VARA_REAL(FILEid ,verVID,start,count,varlev)
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR(Routine, Ierro)
      END IF

      IF (icheck.ge.3) WRITE(*,*) 'Ierro5.1', Ierro

C*    5.2 ...for the the variable.
C     ----------------------------

C     ** Set 'imap' and READ with NCVGTG:
C     ** NOTE :                                                  
C     **        we use the 'generalised' reading routine NCVGTG 
C     ** (imap tells NetCDF about the memory locations of var) 
      imap(1) = 1
      imap(2) = imap(1) * i_dim  ! 1st dim of var = i_dim
      imap(3) = imap(2) * j_dim  ! 2nd dim of var = j_dim 
      imap(4) = 0                !  Should NEVER be used        
      Ierro=NF_GET_VARM_REAL(FILEid   ,  varVID ,begREG      , dimREG,
     &                         stride   ,   imap  ,var(1,1,1)          )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR(Routine, Ierro)

      IF (icheck.ge.3) WRITE(*,*) 'Ierro5.2', Ierro
      
C*    6. Check data
C     -------------
      IF (ireadchk.GE.1) THEN
       varmax = var (1,1,1)
       varmin = var (1,1,1)
       DO ll=1,l_dim
       DO jj=1,j_dim
       DO ii=1,i_dim
          var(ii,jj,ll)=var(ii,jj,ll)+0.0E0
C           This fixes underflow values but must compile with -fpe1
          varmax = MAX(var (ii,jj,ll),varmax)
          varmin = MIN(var (ii,jj,ll),varmin)
       ENDDO
       ENDDO
       ENDDO
       IF (varmin.LT.vReadMin .OR. varmax.GT.vReadMax) THEN
          write(*,*) 'WARNING (UNread): variable ', VARname
          write(*,*) '  is out of specified bounds;'
          write(*,*) '  min is:', varmin
          write(*,*) '  max is:', varmax
       ENDIF
      ENDIF
      
      IF (icheck.ge.2) WRITE(*,*) 'UNread : End' 

      END SUBROUTINE UNread

C**
C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNsread :   +                                         +
C**  +-------------------------+                                         +
C**  +  * Reads a model variable from a NetCDF file,                     +
C**  +    SIMPLIFIED VERSION of  UNread  : does NOT read coordinates.    +
C**  +                                                                   +
C**  +                                                                   +
C**  +  INPUT :                                                          + 
C**  +    FILEid  : input file identifier (from NetCDF open)             +
C**  +    VARname  : name of the requested variable.                     +
C**  +    time : [integer*4] is the time index of the data field to read +
C**  +    level: [integer*4] (usefull for 3D-space fields only) :        +
C**  +                       if not=0 --> = no of the level              +
C**  +                                      -> output is 2D (l_dim = 1)  +
C**  +                       if  =0   --> read ALL levels                +
C**  +                                      -> output is 3D              +
C**  +    i_dbeg, j_dbeg      : horizontal indexes of requested region   +
C**  +                          in input data file                       + 
C**  +    i_dim, j_dim, l_dim : ...the dimensions of 'var',              + 
C**  +                       = the dimensions of the sub-region to read  + 
C**  +                       ! l_dim = 1 if level not=0                  + 
C**  +                       ! j_dim = 1 if var is 1D                    + 
C**  +  OUTPUT :                                                         + 
C**  +    var_units                 : physical units of var.             + 
C**  +    var[i_dim,j_dim,l_dim]    :                                    + 
C**  +                            data field values                      + 
C**  +                            (var must be defined, and is REAL  )   + 
C**  +                                                                   + 
C**  +-------------------------------------------------------------------+

      SUBROUTINE UNsread
     &      (FILEid, VARname, time, level, i_dbeg, j_dbeg,
     &                                     i_dim , j_dim , l_dim,
     &       var_units, var)



      IMPLICIT NONE

C     ** input 
      INTEGER        FILEid 
      INTEGER        time, level, i_dbeg, j_dbeg
      INTEGER        i_dim, j_dim, l_dim
      CHARACTER *(*) VARname 

C     ** output
      CHARACTER *(*) var_units
      REAL*4         var (i_dim, j_dim, l_dim)
      REAL*4         varax1(i_dim), varax2(j_dim), varlev(l_dim)

      call UNread (FILEid , VARname , time, level, i_dbeg, j_dbeg,
     &       i_dim   , j_dim   , l_dim    ,
     &       varax1  , varax2  , varlev,      
     &       var_units, var)


      END SUBROUTINE UNsread

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNwcatt :   +                                         +
C**  +-------------------------+                                         +
C**  +  *Character Attributes creation and (over)writing                 +
C**  +    (the NetCDF file must be open, in data mode)                   +
C**  +  *WARNING: this routine (may?) use a temporary disk space         +
C**  +            equal to the file length (duplicate the file)          +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEid  : input file identifier (from UNcreate OR NetCDF open) +
C**  +    varnam  : name of variable to which attribute shall be attached+
C**  +              or 'GLOBAL_ATT'                                      +
C**  +    attnam  : name of writen attribute.                            +
C**  +    attval  : string to be assigned to attribute.                  +
C**  +              (never inclulde more than 3 consecutive blanks !)    +
c**  +                                                                   +
C**  +  Note : all arguments except FILEid  are strings of any length    +
C**  +-------------------------------------------------------------------+

      SUBROUTINE UNwcatt (FILEid , varnam, attnam, attval)

      INCLUDE 'libUN.inc'

C     **Input:

      INTEGER FILEid 
      CHARACTER*(*) varnam
      CHARACTER*(*) attnam
      CHARACTER*(*) attval

C     **Local:
      INTEGER  VARSIZE
      EXTERNAL VARSIZE
      INTEGER Nlen, Ierro, varVID, Vlen, TTerr
      INTEGER icheck
      icheck= 0     !** 'debugging' level

      IF (icheck.ge.1) WRITE(*,*) 'UNwcatt : Begin'

C*    Get the variable ID
C     -------------------

      IF (icheck.ge.2) WRITE(*,*) 'FILEid  :', FILEid 

C     ** getting varnam size :
      Nlen = VARSIZE(varnam)

C     ** Case of global attributes:
      IF (varnam(1:Nlen).EQ.'GLOBAL_ATT') THEN 
        varVID=NF_GLOBAL 

      ELSE

C     ** Get variable ID to which att is attached to:
        Ierro=NF_INQ_VARID (FILEid , varnam(1:Nlen), varVID)
        TTerr = ABS(Ierro)

C       ** Cancel writing if an error occured : variable undefined ?
        IF (Ierro.ne.0) THEN
           WRITE(*,*) 'UNwcatt -ERROR : Variable ',varnam(1:Nlen)
     &               ,' not found -> not written.'
        END IF
        IF (Ierro.ne.0) RETURN !** UNwcatt_end

      ENDIF

C     Switch to Define Mode, 
C       because attribute may be created or change size.
C     --------------------------------------------------
      Ierro=NF_REDEF (FILEid)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwcatt', Ierro)

C     Set attribute.
C     --------------

C     ** getting attnam [char] size :
      Nlen = VARSIZE(attnam)
C     ** getting attval [char] size :
      Vlen = VARSIZE(attval)

      Ierro=NF_PUT_ATT_TEXT(FILEid ,varVID ,attnam(1:Nlen),
     &                       Vlen  ,attval(1:Vlen)      )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwcatt', Ierro)
c     **line1^^^^ FILEid |var.id | attr.name
C     **line2     type   | len   | attr.value | flag
      TTerr = TTerr + ABS(Ierro)


C     Leave define mode (!file remains open )
C     ---------------------------------------
      Ierro=NF_ENDDEF(FILEid )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwcatt', Ierro)

      RETURN
      END 

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNwratt :   +                                         +
C**  +-------------------------+                                         +
C**  +  *Real   attributes writing  - ! Can not create new attrib !      +
C**  +    (the NetCDF file must be open)                                 +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEid  : input file identifier (from UNcreate OR NetCDF open) +
C**  +    varnam  : name given to the variable to write (must be in file)+
C**  +    attnam  : name of treated attribute.                           +
c**  +    Nvals   : Number of values of that attribute                   +
C**  +    atvalsi(Nvals) : Real   vector of values for attribute.        +
c**  +                                                                   +
C**  +-------------------------------------------------------------------+

C                    WARNING: this routine uses a temporary disk space
C                             equal to the file length (duplicate the file)
C                             (its use is NOT recommended)

      SUBROUTINE UNwratt (FILEid , varnam, attnam, Nvals, atvals)

      INCLUDE 'libUN.inc'

C     **Input:

      INTEGER FILEid , Nvals
      CHARACTER*(*) varnam
      CHARACTER*(*) attnam
      REAL*4        atvals(Nvals)

C     **Local:
      INTEGER  VARSIZE
      EXTERNAL VARSIZE
      INTEGER Nlen, Ierro, varVID
      INTEGER icheck, TTerr
      icheck= 0     !** 'debugging' level
      TTerr = 0

      IF (icheck.ge.1) WRITE(*,*) 'UNwratt : Begin'

C*    Get the variable ID
C     -------------------
      IF (icheck.ge.2) WRITE(*,*) 'FILEid  :', FILEid 

C     ** getting varnam size :
      Nlen = VARSIZE(varnam)

C     ** variable ID :
      Ierro=NF_INQ_VARID(FILEid , varnam(1:Nlen), varVID)
      TTerr = TTerr + ABS(Ierro)

C     ** Cancel writing if an error occured : variable undefined ?
      IF (Ierro.ne.0) THEN
         WRITE(*,*) 'UNwratt -ERROR : Variable ',varnam(1:Nlen)
     &            ,' not found -> not written.'
      END IF
      IF (Ierro.ne.0) GOTO 9999 !** UNwratt_end


C     Set attribute.
C     --------------

C     ** getting attnam [char] size :
      Nlen = VARSIZE(attnam)

      Ierro=NF_PUT_ATT_REAL(FILEid ,varVID ,attnam(1:Nlen),
     &           NF_FLOAT,nvals  ,atvals  )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNwratt', Ierro)
c     **line1^^^^FILEid |var.id | attr.name
C     **line2    type   | attr.value | flag
      TTerr = TTerr + ABS(Ierro)


 9999 continue
      RETURN
      END

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNwopen :   +                            libUN (0896) +
C**  +-------------------------+-----------------------------------------+
C**  +  * Open a NetCDF file for writing.                                +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEnam : file name                                            +
C**  +                                                                   +
C**  +  OUTPUT :                                                         +
C**  +    FILEid  : NetCDF file identifier ('logical unit')              +
C**  +---------------------------------------------------------------7++++
 
      SUBROUTINE UNwopen (FILEnam, FILEid )

      IMPLICIT NONE
      INCLUDE 'libUN.inc'

C     ** input
      CHARACTER*(*) FILEnam

C     ** output
      INTEGER FILEid

C     ** local :
      INTEGER Ierro
      INTEGER icheck

      icheck=0
      
C +   Routines which opens a file must reset libUN internals:
      CALL UNparam('RESET_PARAMS_',0.0_4)
      
C     ** Open NetCDF file, for read-only:
C     -----------------------------------
      Ierro=NF_OPEN(FILEnam,NF_WRITE,FILEid)
      IF (Ierro.NE.NF_NOERR) THEN
         WRITE(*,*) 'Error opening file: ', FILEnam            
         CALL HANDLE_ERR('UNwopen', Ierro)
      ENDIF


9999  continue
      RETURN
      END



C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNropen :   +                            libUN (0896) +
C**  +-------------------------+-----------------------------------------+
C**  +  * Open a NetCDF file for reading,                                +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEnam : file name                                            +
C**  +                                                                   +
C**  +  OUTPUT :                                                         +
C**  +    FILEid  : NetCDF file identifier ('logical unit')              +
C**  +    FILEtit : title of the NetCDF file                             +
C**  +              ! [CHAR], must be defined (length > length(title) !) +
C**  +---------------------------------------------------------------7++++

      SUBROUTINE UNropen (FILEnam, FILEid , FILEtit)

      IMPLICIT NONE
      INCLUDE 'libUN.inc'

C     ** input
      CHARACTER*(*) FILEnam

C     ** output
      INTEGER FILEid      
      CHARACTER*(*) FILEtit

C     ** local :
      INTEGER Ierro
      INTEGER icheck

      icheck=0
      
      IF (icheck.ge.2) WRITE(*,*) 'UNropen: Begin'
      IF (icheck.ge.2) WRITE(*,*) 'FILEnam: ', FILEnam

C +   Routines which opens a file must reset libUN internals:
      CALL UNparam('RESET_PARAMS_',0.0_4)

C     ** Open NetCDF file, for read-only:
C     -----------------------------------
      Ierro=NF_OPEN(FILEnam,NF_NOWRITE,FILEid)
      IF (Ierro.NE.NF_NOERR) THEN
         WRITE(*,*) 'Error opening file: ', FILEnam
         CALL HANDLE_ERR('UNropen', Ierro)
      ENDIF


C     ** Read title attribute, 
C     ------------------------

C     ** Read attribute:
! FD      Ierro=NF_GET_ATT_TEXT(FILEid, NF_GLOBAL, 'title', 
! FD     &             FILEtit)
! FD
! FDC     ** Display message if an error occured : 
! FDC     **  no title or title too long ? 
! FD      IF (Ierro.ne.0) THEN
! FD         WRITE(*,*) 'UNropen WARNING: no title or title too long' 
! FD      END IF
      IF (icheck.ge.2) WRITE(*,*) 'UNropen: End'

9999  continue
      RETURN
      END

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNgtime :   +                            libUN (0896) +
C**  +-------------------------+-----------------------------------------+
C**  +  * From a given value of desired 'time' coordinate,               +
C**  +    gets the coordinate index ('iteration no') + found time value  +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEid  : NetCDF file identifier (from UNropen)                +
C**  +    RQtime  : ReQuested time                                       +
C**  +                                                                   +
C**  +  OUTPUT :                                                         +
C**  +    RDtime  : The last time for wich RDtime .le. RQtime            +
C**  +    Ftime   : The next time value Following RDtime                 +
C**  +              (-1 if it would be after end-of-file)                +
C**  +    it      : The time index : RDtime = time(it)                   +
C**  +---------------------------------------------------------------7++++

      SUBROUTINE UNgtime (FILEid, RQtime, RDtime, Ftime, it) 

      IMPLICIT NONE
      INCLUDE 'libUN.inc'

      INTEGER Lvnam
      PARAMETER (Lvnam=20)

C     ** input
      INTEGER FILEid 
      REAL*4  RQtime

C     ** output
      REAL*4  RDtime, Ftime
      INTEGER it

C     ** local :
      INTEGER Ierro, timVID
      INTEGER timDID
      REAL*4  gtim
      INTEGER K, KHI, KLO, Kmax
      INTEGER Mindex(1)
      INTEGER icheck
      CHARACTER*(Lvnam) dimNAM(1)

      icheck= 0

C     ** Kmax= nb pas de temps dans le fichier, = dim(time):
C     ** - - - - - - - - - - - - - - - - - - - - - - - - - -
C     
      Ierro=NF_INQ_DIMID(FILEid, 'time', timDID)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgtime', Ierro)
C     **^^ Dimension'time' NetCDF index

      Ierro=NF_INQ_DIM(FILEid, timDID , dimNAM, Kmax  )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgtime', Ierro)
C     **         id/file  id/dim   dimname dimsize  error
C     **                           !output output

C     ** Read/Search the requested time step.
C     ** - - - - - - - - - - - - - - - - - - -

      Ierro=NF_INQ_VARID(FILEid, 'time',timVID)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgtime', Ierro)
C                                         **^^ Variable 'time' NetCDF index

      KLO=1
      KHI=Kmax

 1    IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2

C       ** Set the position of the needed time step:
        Mindex(1)= K
C       ** Get 1 time value (gtim = time(K)):
        Ierro=NF_GET_VAR1_REAL(FILEid, timVID, Mindex, gtim)
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgtime', Ierro)

        IF(gtim.GT.RQtime)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      it= KLO
C     ** read RDtime= time(KLO)
      Mindex(1)= KLO
      Ierro=NF_GET_VAR1_REAL(FILEid, timVID, Mindex, RDtime)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgtime', Ierro)
C     ** read Ftime= time(KHI)
      Mindex(1)= KHI
      Ierro=NF_GET_VAR1_REAL(FILEid, timVID, Mindex, Ftime)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgtime', Ierro)
 
C     ** IF the last available time step is before
C     **     the requested time, then KHI and KLO are the
C     **     two last available time step. Correct this :
      IF (RQtime.ge.Ftime) THEN
        RDtime= Ftime                   
        it = KHI
        Ftime= -1.0
      ENDIF

      RETURN
      END

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNgindx :   +                            libUN (0199) +
C**  +-------------------------+-----------------------------------------+
C**  +  * From a given value of a desired coordinate,                    +
C**  +    gets the coordinate index + found the coresp. coordinate value +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEid  : NetCDF file identifier (from UNropen)                +
C**  +    Cname   : The name of the coordinate                           +
C**  +    RQval   : The requested value for that coordinate              +
C**  +                                                                   +
C**  +  OUTPUT :                                                         +
C**  +    RDval   : The last value for wich RDval .le. RQval             +
C**  +    Fval    : The next val value Following RDval                   +
C**  +              (-1 if it would be after end-of-file)                +
C**  +    indx    : The val index : RDval = value_of_Cname(it)           +
C**  +---------------------------------------------------------------7++++

      SUBROUTINE UNgindx (FILEid, Cname, RQval, RDval, Fval, indx)

      IMPLICIT NONE
      INCLUDE 'libUN.inc'

      INTEGER Lvnam
      PARAMETER (Lvnam=20)

C     ** input
      INTEGER FILEid
      CHARACTER *(*) Cname
      REAL*4  RQval

C     ** output
      REAL*4  RDval, Fval
      INTEGER indx 

C     ** local :
      INTEGER  VARSIZE
      EXTERNAL VARSIZE
      REAL*4  gval
      INTEGER Ierro
      INTEGER varDID, VNlen, varVID, varNUMDIM
      INTEGER Nvatts, vtype
      INTEGER K, KHI, KLO, Kmax
      INTEGER Mindex(1), dimID(4)
      INTEGER icheck
      CHARACTER*(Lvnam) dimNAM(4)
      CHARACTER*13 recname

      icheck= 0

C     ** Kmax= nb pas de temps dans le fichier, = dim(val):
C     ** - - - - - - - - - - - - - - - - - - - - - - - - - -
C     ** get Cname string size :
      VNlen = VARSIZE (Cname)
C
C     ** get variable ID :
      Ierro=NF_INQ_VARID(FILEid , Cname (1:VNlen), varVID)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgindex', Ierro)
C
C     ** Inquire about the id of the dimension:
C     **
      Ierro=NF_INQ_VAR(FILEid , varVID, recname, vtype,
     &          varNUMDIM, dimID , Nvatts)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgindex', Ierro)
C     **  line1  id/file   id/var  var name  var type
C     **  line2   # dims   id/dims #attributes
      varDID = dimID(1) 
C     ^^^At last, the id of the relevant dimension.

      Ierro=NF_INQ_DIM(FILEid, varDID , dimNAM, Kmax  )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgindex', Ierro)
C     **         id/file  id/dim   dimname dimsize  error
C     **                           !output output
C     ** (Kmax is what we needed: size of the dimension)

C     ** Read/Search the requested val step.
C     ** - - - - - - - - - - - - - - - - - - -

      KLO=1
      KHI=Kmax

 1    IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2

C       ** Set the position of the needed val step:
        Mindex(1)= K
C       ** Get 1 val value (gval = val(K)):
        Ierro=NF_GET_VAR1_REAL(FILEid, varVID, Mindex, gval)
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgindex', Ierro)

        IF(gval.GT.RQval)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      indx= KLO
C     ** read RDval= val(KLO)
      Mindex(1)= KLO
      Ierro=NF_GET_VAR1_REAL(FILEid, varVID, Mindex, RDval)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgindex', Ierro)
C     ** read Fval= val(KHI)
      Mindex(1)= KHI
      Ierro=NF_GET_VAR1_REAL(FILEid, varVID, Mindex, Fval)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNgindex', Ierro)

C     ** IF the last available val step is before
C     **     the requested val, then KHI and KLO are the
C     **     two last available val step. Correct this :
      IF (RQval.ge.Fval) THEN
        RDval= Fval
        indx = KHI
        Fval= -1.0
      ENDIF

      RETURN
      END

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNfindx :   +                            (libUN  2003)+
C**  +-------------------------+-----------------------------------------+
C**  +  * Intended to replace UNgindx or UNgtime                         +
C**  +    From a given value of a desired coordinate,                    +
C**  +    gets the coordinate index + the coresp. coordinate value       +
C**  +    This version solves the issue of Dates at year change          +
C**  +    occuring because 1 jan is < 31 dec.  Not optimised.            +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEid  : NetCDF file identifier (from UNropen)                +
C**  +    Cname   : The name of the coordinate                           +
C**  +    RQval   : The requested value for that coordinate              +
C**  +                                                                   +
C**  +  OUTPUT :                                                         +
C**  +    RDval   : The file value closest to RQval                      +
C**  +    Fval    : The next value in the file                           +
C**  +              (-1 if after file end)                               +
C**  +              (This is mainly for compatibility with older version)+
C**  +    indx    : The val index : RDval = value_of_Cname(it)           +
C**  +              (-1 may be returned if the value can't be found)     +
C**  +---------------------------------------------------------------7++++

      SUBROUTINE UNfindx (FILEid, Cname, RQval, RDval, Fval, indx)

      IMPLICIT NONE
      INCLUDE 'libUN.inc'

      INTEGER Lvnam
      PARAMETER (Lvnam=20)

C     ** input
      INTEGER FILEid
      CHARACTER *(*) Cname
      REAL*4  RQval

C     ** output
      REAL*4  RDval, Fval
      INTEGER indx

C     ** local :
      INTEGER  VARSIZE
      EXTERNAL VARSIZE
      REAL*4  gval, bmatch, gdist
      INTEGER Ierro
      INTEGER varDID, VNlen, varVID, varNUMDIM
      INTEGER Nvatts, vtype
      INTEGER K, KHI, KLO, Kmax
      INTEGER Mindex(1), dimID(4)
      INTEGER icheck
      CHARACTER*(Lvnam) dimNAM(4)
      CHARACTER*13 recname

      icheck= 0

C     ** Kmax= nb pas de temps dans le fichier, = dim(val):
C     ** - - - - - - - - - - - - - - - - - - - - - - - - - -
C     ** get Cname string size :
      VNlen = VARSIZE (Cname)
C
C     ** get variable ID :
      Ierro=NF_INQ_VARID(FILEid , Cname (1:VNlen), varVID)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNfindex', Ierro)
C
C     ** Inquire about the id of the dimension:
C     **
      Ierro=NF_INQ_VAR(FILEid , varVID, recname, vtype,
     &          varNUMDIM, dimID , Nvatts)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNfindex', Ierro)
C     **  line1  id/file   id/var  var name  var type
C     **  line2   # dims   id/dims #attributes
      varDID = dimID(1)
C     ^^^At last, the id of the relevant dimension.

      Ierro=NF_INQ_DIM(FILEid, varDID , dimNAM, Kmax  )
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNfindex', Ierro)
C     **         id/file  id/dim   dimname dimsize  error
C     **                           !output output
C     ** (Kmax is what we needed: size of the dimension)

C     ** Read/Search the requested val step.
C     ** - - - - - - - - - - - - - - - - - - -

C     This is a workaround, not optimised as stated above.
C     We simply look at all values sequencially.
C
      bmatch=1.E10
      KLO=-1

      DO K=1,KMAX

C       ** Get 1 val value (gval = val(K)):
        Mindex(1)= K
        Ierro=NF_GET_VAR1_REAL(FILEid, varVID, Mindex, gval)
        IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNfindex', Ierro)

        gdist=ABS(gval-RQval)
        IF (gdist.LT.bmatch) THEN

         bmatch=gdist
         KLO=K

        ENDIF

      ENDDO

      indx= KLO

      KHI = min((KLO+1),KMAX)

C     ** read values...

      Mindex(1)= KLO
      Ierro=NF_GET_VAR1_REAL(FILEid, varVID, Mindex, RDval)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNfindex', Ierro)
C     ** read Fval= val(KHI)
      Mindex(1)= KHI
      Ierro=NF_GET_VAR1_REAL(FILEid, varVID, Mindex, Fval)
      IF (Ierro.NE.NF_NOERR) CALL HANDLE_ERR('UNfindex', Ierro)

      IF (KHI.EQ.KLO) THEN
        Fval= -1.0
      ENDIF

      IF (bmatch.GT.1.E9) THEN
        Fval= -1.0
        indx= -1
      ENDIF

      RETURN
      END

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNclose :   +                            libUN (0300) +
C**  +-------------------------+-----------------------------------------+
C**  +  * Close the desired file                                         +
C**  +    Created to suppress the need the directly call a netcdf        +
C**  +    routine from a program                                         +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +    FILEid  : NetCDF file identifier (from UNropen)                +
C**  +---------------------------------------------------------------7++++

      SUBROUTINE UNCLOSE(FILEid)

      IMPLICIT NONE
      INCLUDE 'libUN.inc'

      integer Ierro, FILEid

      Ierro=NF_CLOSE(FILEid)
      IF (Ierro.NE.NF_NOERR) THEN
         CALL HANDLE_ERR('UNclose', Ierro)
      ENDIF

      END
      
C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine UNparam :   +                            libUN (0202) +
C**  +-------------------------+-----------------------------------------+
C**  +  Changes some global libUN parameters                             +
C**  +  NB: default values are set at first libUN call                   +
C**  +                                                                   +
C**  +                                                                   +
C**  +  INPUT : pname   name of the parameters to set                    +
C**  +          pvalue  the requested new value                          +
C**  +                                                                   +
C**  +---------------------------------------------------------------7++++

      SUBROUTINE UNparam(pname,pvalue)

      IMPLICIT NONE

      INCLUDE 'libUN.inc'
 
      CHARACTER*(*) pname
      REAL*4  pvalue
      
      LOGICAL Lstart
      SAVE    Lstart
      DATA    Lstart /.true./

      IF      (pname.EQ.'RESET_PARAMS_') THEN
         IF (Lstart.OR.pvalue.GT.0.5) THEN
           vMissVal= 1.0E21   ! for missing values
           VarRepl = vMissVal ! for missing VARIABLES
           ValRange(1)= -vMissVal/10.
           ValRange(2)=  vMissVal/10.
           iVarWarn= 2
           vReadMin  = 0.0
           vReadMax  = 0.0
           ireadchk = 0
           Lstart   = .false.
          ENDIF
          
      ELSE IF (pname.EQ.'NOVAR_REPLACE') THEN
         VarRepl = pvalue  
         
      ELSE IF (pname.EQ.'NOVAR_WARNING') THEN
         iVarWarn= NINT(pvalue)
         
      ELSE IF (pname.EQ.'VALID_RANGE_MIN') THEN
         ValRange(1) = pvalue  

      ELSE IF (pname.EQ.'VALID_RANGE_MAX') THEN
         ValRange(2) = pvalue  
      
      ELSE IF (pname.EQ.'READOVER_WARN') THEN
         vReadMin  = - pvalue
         vReadMax  =   pvalue
         ireadchk = 1
      
      ELSE IF (pname.EQ.'READ_MIN_WARN') THEN
         vReadMin  =   pvalue
         ireadchk = 1
         
      ELSE IF (pname.EQ.'READ_MAX_WARN') THEN
         vReadMax  =   pvalue
         ireadchk = 1   
         
      ELSE 
         write(*,*) 'UNparam (libUN) Error: '       
         write(*,*) '  parameter undefined:', pname       

      ENDIF

      END
      
C**  +-------------------------+-----------------------------------------+
      SUBROUTINE UNversion(UNver,NCDFver)
C**  +-------------------------+-----------------------------------------+
      
      IMPLICIT NONE
      INCLUDE 'libUN.inc'

      CHARACTER*80 UNver,NCDFver

      UNver  = '2005.03.31'
      NCDFver= NF_INQ_LIBVERS()

      END

C**  +-------------------------------------------------------------------+
      FUNCTION VARSIZE(CHAvar)
C**  +-------------------------------------------------------------------+
      IMPLICIT NONE
      integer maxcha,iz,VARSIZE
      parameter (maxcha=512)
      character*(*)      CHAvar
      character*(maxcha) CHAtmp

      WRITE(CHAtmp,'(A)') CHAvar
      iz = 0
      do while ((CHAtmp(iz+1:iz+3).ne.'   ').and.(iz+3.le.maxcha))
        iz = iz + 1
      end do
      VARSIZE =  iz

      RETURN
      END


C**  +-------------------------------------------------------------------+
      SUBROUTINE HANDLE_ERR(LOCATION, STATUS)
C**  +-------------------------------------------------------------------+
      IMPLICIT NONE

      INCLUDE 'libUN.inc'

      character*(*) LOCATION
      integer STATUS
      IF (STATUS.NE.NF_NOERR) THEN
        WRITE(*,*) 'IN ROUTINE ', LOCATION
        WRITE(*,*) NF_STRERROR(STATUS)
        STOP 'Stopped'
      ENDIF
      END

C UN library: history of fixed bugs and updates.
C ----------------------------------------------
C
C                        961206 - UNgtime, trouble at end-of-file
C                        961218 - - all -, display 'artificial' errors
C                        970318 -   again, display 'artificial' errors
C                        971028 - (3 sub),'syntax'error on Cray computer
C                        971105 - Allowed variable "imap(1)", =8 for Cray
C                        980705 - "single element" extension to UNwrite.
C                        980709 - bug fixes (start) in UNwrite & UNlwrite
C                                 ("DATA" statement incorrectly used).
C                        980825 - Changed default "stride" to 1 for v3.x
C                        981222 - bug fix: allow UNwrite for unlim dims.
C                                 note that this should be tested.
C                        990110 - Added "UNgindx" = general. of UNgtime
C                               - Removed all "DATA" and all "//" in write
C                                 (the later should improve compatibility)
C                        990128 - UNwrite: added a "no warning" option.
C                        990323 - UNwrite: added 1D+time capability.
C                        990807 - UNwrite: added 3D-notime capability.
C  -----------------------------------------------------------------------
C                        000404 - Major upgrade: compatibility with
C                                 NetCDF v3.4
C                               - NOTE: Types other than REAL may be
C                                 accepted in UNread, but not tested
C  -----------------------------------------------------------------------
C                        000614 - Bug fixes: uninitialised error count
C                                 in UNwcatt, bug in UNclose. 
C                        000620 - Bug fix: UNropen (args. of get title fn)
C                        000713 - Bug fix: UNgtime (missing arg in a call)
C                                 (last tree caused by 000404 upgrade)
C  -----------------------------------------------------------------------
C                        000928 - UNlwrite: added 2D+time capability.
C                        001008 - All: CHARACTER*(*) declaration for units
C                                 and longer strings for intern. variables
C                        010417 - UNread: added var not found info
C                                 UNropen: added file not found info
C                        010715 - UNwrite + UNlwrite: 
C                                   fixed bug / unlimited time dim
C                        0107xx - UNwrite:  
C                                   missing values -> not in "range"
C                        020130 - All:    
C                                  .removed obsolete warnings about 
C                                   double precision in files.      
C                                  .added a version (libUN_dbl) with
C                                   REAL*8 as arguments - but still
C                                   creates REAL*4 in files.
C                        020526 - Added UNparam function,
C                                 which provide optional features such
C                                 as missing variable behavior control
C                        020808 - Very simple fix for underflows while 
C                                 reading some files; must use -fpe1
C                                 Fixed a bug -> out of range msg
C                        030121 - Enabled some non-standard NetCDF files
C                                 (missing units...) -> new warnings
C                                 rather then program stop.
C                        030215 - Added UNfindx for non-monotonic data
C                        030215 - Removed warning related to UNLIM dims
C                        030311 - Added VALID_RANGE attribute (option)
C                                 (if set, the range is accounted for
C                                 in the min/max set while writing vars)
C                        040902 - Improvements to "valid_range" attribute
C                               - Added attribute "positive=down"
C                                 if units are sigma or sigma_level
C                        050331 - Added "user friendly" interfaces
