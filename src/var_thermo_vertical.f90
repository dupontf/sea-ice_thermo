module var_thermo_vertical

implicit none

      ! arguments
      integer, parameter :: maxlay = 230
      integer nlice, nlsno
      double precision ti(0:maxlay)
      double precision ts(0:maxlay)
      double precision qi(1:maxlay)
      double precision qs(1:maxlay)
      double precision si(maxlay)
      double precision sinew(maxlay)
      double precision dzi(maxlay)
      double precision dzs(maxlay)
      double precision zi(0:maxlay)
      double precision zs(0:maxlay)

      double precision hi, hs
      double precision tsu, tbo
      double precision dhs
      double precision dhi_surf
      double precision dhi_bot
      double precision dh_sni
! heat flux
      double precision oceflx
      double precision netlw, dwnlw, pres
      double precision tair, qair, uair
      double precision tocn, uio   ! sst (degC) and relative ice-ocean velocity		[m/s]
      double precision fsens, flat
      double precision q0
      double precision fac_transmi, swrad
      double precision heatflx_bot, heatflx_top
      double precision tmelt
      double precision snowfall, fprecip
      double precision seasal
      double precision eskel
      double precision snosub
      double precision massmelt
      double precision fsalt
      double precision si_acc_new
      double precision frac_sni, fcons_sni
! conductivity flux
      double precision fcsu, fcbo
      double precision hfc_s(0:maxlay)
      double precision hfc_i(0:maxlay)
! interior absorption
      double precision swradab_i(maxlay)
      double precision swradab_s(maxlay)

! constants
      double precision :: &
        cp_ice    = 2.062d+03, &          ! sea ice specific heat
        cp_wat    = 3.99d+03, &           ! seawater specific heat
        emi   = 0.99d0, &             ! surface emissivity
        cond_sno    = 0.31d0, &             ! 0.31d0 ISPOL Olivier Lecomte communication ! ref value 0.31d0 ! snow thermal conductivity
        cond_ice    = 2.034d0, &            ! pure ice thermal conductivity
        mlfus   = 3.335d+05, &          ! massive latent heat, ice
        mlsub    = 2.834d+06, &          ! sublimation latent heat
        tp0     = 273.16d0, &          ! water triple point
        fracsal   = 0.054d0, &            ! rate between seawater freezing point and
        tmelt_sno   = 273.16d0, &           ! snow melting point
        tmelt_ice   = 273.16d0, &           ! sea ice melting point
        stefa = 5.6697d-08, &         ! stefa-boltzmann constant
        rhowat   = 1025.d0 , &           ! ocean mean density
        rhosno   = 330.d0, &             ! 355.d0 ISPOL Olivier Lecomte communication ! ref value 330.0  ! snow density
        rhoice   = 917.d0              ! sea ice density

! reference salinity
  double precision,save :: sref=34.80d0

        logical :: first=.true. ! required for the Huwald and FE models

! varying sigma coordinate
  integer ni,ns
  double precision, dimension  (maxlay) :: dziold
  double precision, dimension(0:maxlay) :: tiold, tinew, timid, ziold
  double precision ::   hminice=.001d0

  integer ith_cond ! conductivity formula switch



contains


!---------------------------------------------------------------------
! Freezing temperature of sea water at one standard atmosphere
! Ref: Winton (2000): A reformulated three-layer sea ice model 

double precision function Tfreeze1(s)

  implicit none

  double precision,intent(in) :: s ! Salinity
  double precision,parameter :: a=-0.054d0

  Tfreeze1=a*s

end function Tfreeze1

!---------------------------------------------------------------------
! Freezing temperature of sea water at one standard atmosphere
! Ref: (Millero, 1978) and Gill (1982): Atmosphere-Ocean Dynamics

double precision function Tfreeze2(s)

  implicit none

  double precision,intent(in) :: s ! Salinity

  double precision,parameter :: a=-0.0575d0
  double precision,parameter :: b=1.710523d-3
  double precision,parameter :: c=-2.154996d-4

  Tfreeze2=s*(a+b*sqrt(s)+c*s)

end function Tfreeze2


!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
FUNCTION func_ki(S,T)

implicit none
      DOUBLE PRECISION  func_ki, S, T

  select case (ith_cond)
   case(0)
     func_ki = func_ki0(S,T)
   case(1)
     func_ki = func_ki1(S,T)
  end select

END FUNCTION func_ki
!!-----------------------------------------------------------------------

!!-----------------------------------------------------------------------
!c     Ice thermal conductivity
!!-----------------------------------------------------------------------
!c         ki + beta S / T > 0 : ki should remain positive as T --> Tf = -mu*S

!c         ki > beta / mu  -->  beta < ki mu --> beta = (ki-0.02)*mu
!!-----------------------------------------------------------------------
!cc     Ice thermal conductivity
!cc         Ki should remain positive as T --> Tm = -muS
!cc                ki + beta S / Tm > 0
!cc                ki > beta / mu
!cc                beta < Ki mu --> beta = ( Ki - Kimin ) mu
!cc
!cc         Kimin can be that of air of ocean water, depending on assumption.
!cc         Kiair = 0.02 W/m/C and Kiwater = 0.33 W/m/C
!c!----------------------------------------------------------------------

FUNCTION func_ki0(S,T)

implicit none
      DOUBLE PRECISION  func_ki0, S, T, betanew
      DOUBLE PRECISION :: ki0 = 2.034d+0        ! thermal conductivity of fresh ice	[W/m/C]
      double precision :: mu = 0.054d0          ! empirical constant relating S and Tf	[C/psu]

      betanew = (ki0-0.02d0)*mu                 ! see HEADER
      func_ki0 = ki0 + betanew * S / MIN( T, -1d-20 )

END FUNCTION func_ki0

!------------------------------------------------------------------------------ 
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 
!  1) Thermal conductivity at the ice interfaces
!------------------------------------------------------------------------------ 
!
      ! Pringle et al., JGR 2007 formula
      ! 2.11 + 0.09 S/T - 0.011.T

!------------------------------------------------------------------------------ 


FUNCTION func_ki1(S,T)

implicit none
      DOUBLE PRECISION  func_ki1, S, T
      DOUBLE PRECISION :: &
        coeff1 = 0.09d0,  &            ! first th.cond. constant
        coeff2 = 0.011d0, &            ! second th.cond. constant
        zkimin   =  0.1d0
      func_ki1     = cond_ice + coeff1 * S / MIN( T, -1d-20 ) - coeff2 * T
      func_ki1     = MAX( func_ki1 , zkimin )

END FUNCTION func_ki1




!-----------------------------------------------------------------------
!     Heat capacity
!-----------------------------------------------------------------------

      FUNCTION func_cp(Tf,T1,T2)

      implicit none
      DOUBLE PRECISION  func_cp, T1, T2, Tf
      double precision :: TT

      TT = MIN ( T1, -1d-10 ) * MIN ( T2, -1d-10 )
      func_cp = MIN( cp_ice - mlfus * Tf / TT, 1d9 )

      END FUNCTION func_cp

!-----------------------------------------------------------------------
!     Energy of melt
!-----------------------------------------------------------------------

      FUNCTION func_qmelt(Tf,T)

      implicit none
      DOUBLE PRECISION  func_qmelt, T, Tf

      func_qmelt = cp_ice * ( Tf - T ) + mlfus * ( 1.d0 - func_bf (Tf ,T ) )

      END FUNCTION func_qmelt

!-----------------------------------------------------------------------
!     Energy of melt
!-----------------------------------------------------------------------

      FUNCTION func_qm(Tf,T)

      implicit none
      DOUBLE PRECISION  func_qm, T, Tf

      func_qm = cp_ice * ( Tf - T ) + mlfus * ( 1.d0 - func_bf (Tf ,T ) ) - cp_wat * Tf

      END FUNCTION func_qm

!-----------------------------------------------------------------------
!     Internal Energy
!-----------------------------------------------------------------------
!     E = cp_ice*(1-     bf  ) Ti + mlfus*(1 -    bf  ) + cp_wat*    bf   *Ti 
!       = cp_ice*(1- -mu*S/Ti) Ti + mlfus*(1- -mu*S/Ti) + cp_wat*(-muS/Ti)*Ti
!         cp_ice*(Ti -   Tf  )    + mlfus*(1 - Tf/Ti)   + cp_wat*Tf
!-----------------------------------------------------------------------

      FUNCTION func_El(Tf,T)

      implicit none
      DOUBLE PRECISION  func_El, Tf, T

      func_El = cp_ice * ( T - Tf ) - mlfus * ( 1.d0 - func_bf (Tf ,T ) ) + cp_wat * Tf

      END FUNCTION func_El

!-----------------------------------------------------------------------
!     Specific heat term in E
!-----------------------------------------------------------------------

      FUNCTION func_sh(Tf,T)

      implicit none
      DOUBLE PRECISION  func_sh, T, Tf

      func_sh = cp_ice*(1d0-func_bf(Tf,T)) + cp_wat*func_bf(Tf,T)

      END FUNCTION func_sh

!-----------------------------------------------------------------------
!     Latent heat term in E
!-----------------------------------------------------------------------

      FUNCTION func_lh(Tf,T)

      implicit none
      DOUBLE PRECISION  func_lh, Tf, T

      func_lh = -mlfus*(1d0-func_bf(Tf,T))

      END FUNCTION func_lh

!-----------------------------------------------------------------------
!     Brine fraction
!-----------------------------------------------------------------------

      FUNCTION func_bf(Tf,T)

      implicit none
      DOUBLE PRECISION  func_bf, Tf, T

      func_bf = Tf / MIN ( T, -1d-20 )

      END FUNCTION func_bf

!-----------------------------------------------------------------------
!     Latent heat
!-----------------------------------------------------------------------

      FUNCTION func_le(T)

      implicit none
      DOUBLE PRECISION  func_le, T
      DOUBLE PRECISION :: Lsub = 2.834d+6     ! spec. latent heat of sublim, ice/snow [J/kg]

      IF (T .LT. 0d0) THEN
         func_le = Lsub
      ELSE
         func_le = Lsub
      ENDIF

      END FUNCTION func_le


   !==============================================================================
   ! retrieve mean layer frome energy temperature
   !==============================================================================
   subroutine invert_energy_one(em,si,tt)
   implicit none
   double precision em,si
   double precision ztmelts,zaaa,zbbb,zccc,zdiscrim,tt
      !-------------------
      ! Ice temperatures
      !-------------------
             !Energy of melting q(S,T) [J.m-3]
             !Ice layer melt temperature
             ztmelts    =  -0.054d0*si
             !Conversion q(S,T) -> T (second order equation)
             zaaa       =  cp_ice
             zbbb       =  ( cp_wat - cp_ice ) * ztmelts - em - mlfus
             zccc       =  mlfus * ztmelts
             zdiscrim   =  SQRT( MAX(zbbb*zbbb - 4.d0*zaaa*zccc,0.d0) )
             tt = ( - zbbb - zdiscrim ) / ( 2.d0 *zaaa )
             tt = MIN( -0.01d0, MAX(-80d0, tt ) )

   end subroutine invert_energy_one

   subroutine invert_qm_one(em,si,tt)
   implicit none
   double precision em,si
   double precision ztmelts,zaaa,zbbb,zccc,zdiscrim,tt
      !-------------------
      ! Ice temperatures
      !-------------------
             !Energy of melting q(S,T) [J.m-3]
             !Ice layer melt temperature
             ztmelts    =  -0.054d0*si
             !Conversion q(S,T) -> T (second order equation)
             zaaa       =  cp_ice
             zbbb       =  ( - cp_ice ) * ztmelts - em - mlfus
             zccc       =  mlfus * ztmelts
             zdiscrim   =  SQRT( MAX(zbbb*zbbb - 4.d0*zaaa*zccc,0.d0) )
             tt = ( - zbbb - zdiscrim ) / ( 2.d0 *zaaa )
             tt = MIN( -0.01d0, MAX(-80d0, tt ) )

   end subroutine invert_qm_one

end module var_thermo_vertical
