module var_thermo_vertical

implicit none

      ! arguments
      integer, parameter :: maxlay = 230
      integer nlice, nlsno
      real(8), dimension(:), allocatable :: &
              ti, & !  ice temperature (0:nlice)
              ts, & ! snow temperature (0:nlsno)
              qi, & !  ice enthalpy (nlice)
              qs, & ! snow enthalpy (nlsno)
              si, & !  ice salinity (nlice)
              sinew, & !  ice salinity (nlice)
              dzi,& !  ice layer tickness (nlice)
              dzs,& ! snow layer tickness (nlsno)
              zi, & !  ice sigma levels (0:nlice)
              zs    ! snow sigma levels (0:nlice)
      real(8), dimension(:), allocatable :: &
              dziold, &
              tiold, tinew, timid, ziold
      real(8), dimension(:,:), allocatable :: &
              cs, &! bottom and top sigma level of each layer
              os   ! bottom and top 1-sigma level of each layer
! interior absorption
      real(8), dimension(:), allocatable :: &
              swradab_i, &! light penetration in ice
              swradab_s   ! light penetration in snow
! Brine information
      real(8), dimension(:), allocatable :: &
              brine_v, &! brine volume
              brine_u, &! brine flushing velocity
              brine_r   ! brine Rayleigh number

      real(8) hi, hs
      real(8) tsu, tbo
      real(8) dhs
      real(8) dhi_surf
      real(8) dhi_bot
      real(8) dh_sni
! heat flux
      real(8) oceflx
      real(8) netlw, dwnlw, pres
      real(8) tair, qair, uair
      real(8) tocn, uio   ! sst (degC) and relative ice-ocean velocity		[m/s]
      real(8) fsens, flat
      real(8) q0
      real(8) fac_transmi, swrad
      real(8) heatflx_bot, heatflx_top
      real(8) tmelt
      real(8) snowfall, fprecip
      real(8) seasal
      real(8) eskel
      real(8) snosub
      real(8) massmelt
      real(8) fsalt
      real(8) si_acc_new
      real(8) frac_sni, fcons_sni
! conductivity flux
      real(8) fcsu, fcbo
      real(8) hfc_s(0:maxlay) ! needed by LIM
      real(8) hfc_i(0:maxlay) ! needed by LIM

! constants
      real(8), parameter :: &
        cp_ice    = 2.062e+03_8, &          ! sea ice specific heat
        cp_wat    = 3.99e+03_8, &           ! seawater specific heat
        emi       = 0.99_8, &               ! surface emissivity
        cond_sno  = 0.31_8, &               ! 0.31d0 ISPOL Olivier Lecomte communication ! ref value 0.31d0 ! snow thermal conductivity
        cond_ice  = 2.034_8, &              ! pure ice thermal conductivity
        mlfus     = 3.335e+05_8, &          ! massive latent heat, ice
        mlsub     = 2.834e+06_8, &          ! sublimation latent heat
        tp0       = 273.16_8, &             ! water triple point
        fracsal   = 0.054_8, &              ! rate between seawater freezing point and
        tmelt_sno = 273.16_8, &             ! snow melting point
        tmelt_ice = 273.16_8, &             ! sea ice melting point
        stefa     = 5.6697e-08_8, &         ! stefa-boltzmann constant
        rhowat    = 1025.0_8 , &            ! ocean mean density
        rhosno    = 330.0_8, &              ! 355.d0 ISPOL Olivier Lecomte communication ! ref value 330.0  ! snow density
        rhoice    = 917.0_8                 ! sea ice density

! reference salinity
  real(8),save :: sref=34.80d0

        logical :: first=.true. ! required for the Huwald and FE models

! varying sigma coordinate
  integer ni,ns
  real(8) ::   hminice=.001d0

  real(8) ::   fsbr ! brine mass flux
  real(8) ::   fhbr ! brine heat flux (Wm-2)
  integer ith_cond ! conductivity formulation switch
  integer i_perm_eff, i_perm_for ! permeability formula switch
  integer i_Ra ! Rayleigh number formulation switch
  integer iliquid_cond ! liquidus interface switch

!--------------------------------------------------------------------------------------------------
! GN model constants (i_scheme = 2)
  real(8) :: Rc_GN     = 1.01_8     ! Critical Rayleigh number
  real(8) :: alpha_GN  = 1.56e-3_8  ! Brine flow (kg/m3/s)
  real(8) :: rho_br_GN = 1020._8    ! Brine density (kg/m3)
!--------------------------------------------------------------------------------------------------
! RJW model constants

  real(8) :: cl     = 4000000._8   ! Liquid volumetric heat capacity (J/m3)
  real(8) :: grav   = 9.81_8       ! Gravity (m/s^2)
  real(8) :: beta_s = 7.5e-4_8     ! Beta_S (g/kg)-1
  real(8) :: kl     = 0.523_8      ! W/m/K
  real(8) :: visc   = 1.8e-6_8     ! Kinematic viscosity
!  real(8) :: Rc_RJW = 40.0_8        ! Critical Rayleigh number (40 in the paper)
  real(8) :: Rc_RJW = 1.01_8        ! Critical Rayleigh number (40 in the paper)
  real(8) :: al_RJW = 0.03_8       ! tuned value?
!  real(8) :: al_RJW = 30.0_8       ! tuned value?
  real(8) :: w_pr   = -2.0e-8_8    ! Prescribed brine velocity (guessed value)

!--------------------------------------------------------------------------------------------------
! Vetal TCD 2013 Rayleigh constants
  real(8) :: kappa = 1.2e-7_8   ! Thermal diffusivity of brine
  real(8) :: mudyn = 1.9e-3_8   ! dynamical viscosity of brine
  real(8) :: B_S   = 0.81_8     ! Sensitivity of density to salinity (kg/m3/(g/kg))              

contains

subroutine allocate_thermo_1d

     allocate(ti(0:nlice), & !  ice temperature (0:nlice)
              ts(0:nlsno+1), & ! snow temperature (0:nlsno)
              qi(nlice),   & !  ice enthalpy (nlice)
              qs(nlsno),   & ! snow enthalpy (nlsno)
              si(nlice),   & !  ice salinity (nlice) ! sinew needed in LIM
              sinew(nlice),   & !  ice salinity (nlice) ! sinew needed in LIM
              dzi(nlice+nlsno),  & !  ice layer tickness (nlice)
              dzs(nlice+nlsno),  & ! snow layer tickness (nlsno)
              zi(0:nlice+nlsno), & !  ice sigma levels (0:nlice)
              zs(0:nlice+nlsno) )  ! snow sigma levels (0:nlice)
     allocate(dziold(nlice+nlsno),  &
              ziold(0:nlice+nlsno), &
              tiold(0:nlice+nlsno), &
              tinew(0:nlice+nlsno), &
              timid(0:nlice+nlsno))
     allocate( &
              cs(2,0:nlice+nlsno), & ! bottom and top sigma level of each layer
              os(2,0:nlice+nlsno))   ! bottom and top 1-sigma level of each layer
     allocate( &
              swradab_i(nlice), & ! light penetration in ice
              swradab_s(nlsno))   ! light penetration in snow
     allocate( &
              brine_v(nlice), &   ! brine volume
              brine_u(nlice), &   ! brine flushing velocity
              brine_r(nlice))     ! brine Rayleigh number

end subroutine allocate_thermo_1d

!---------------------------------------------------------------------
! Freezing temperature of sea water at one standard atmosphere
! Ref: Winton (2000): A reformulated three-layer sea ice model 

real(8) function Tfreeze1(s)

  implicit none

  real(8),intent(in) :: s ! Salinity
  real(8),parameter :: a=-0.054d0

  Tfreeze1=a*s

end function Tfreeze1

!---------------------------------------------------------------------
! Freezing temperature of sea water at one standard atmosphere
! Ref: (Millero, 1978) and Gill (1982): Atmosphere-Ocean Dynamics

real(8) function Tfreeze2(s)

  implicit none

  real(8),intent(in) :: s ! Salinity

  real(8),parameter :: a=-0.0575d0
  real(8),parameter :: b=1.710523d-3
  real(8),parameter :: c=-2.154996d-4

  Tfreeze2=s*(a+b*sqrt(s)+c*s)

end function Tfreeze2


!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
FUNCTION func_ki(S,T)

implicit none
      real(8)  func_ki, S, T

  select case (ith_cond)
   case(0)
     func_ki = func_ki0(S,T)
   case(1)
     func_ki = func_ki1(S,T)
   case(2)
     func_ki = 2.0_8
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
      real(8)  func_ki0, S, T, betanew
      real(8) :: ki0 = 2.034d+0        ! thermal conductivity of fresh ice	[W/m/C]
      real(8) :: mu = 0.054d0          ! empirical constant relating S and Tf	[C/psu]

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
      real(8)  func_ki1, S, T
      real(8) :: &
        coeff1 = 0.09d0,  &            ! first th.cond. constant
        coeff2 = 0.011d0, &            ! second th.cond. constant
        zkimin   =  0.1d0
      func_ki1     = cond_ice + coeff1 * S / MIN( T, -1d-20 ) - coeff2 * T
      func_ki1     = MAX( func_ki1 , zkimin )

END FUNCTION func_ki1




!-----------------------------------------------------------------------
!     Heat capacity (HB99)
!-----------------------------------------------------------------------

      FUNCTION func_cp(Tf,T1,T2)

      implicit none
      real(8)  func_cp, T1, T2, Tf
      real(8) :: TT

      TT = MIN ( T1, -1d-10 ) * MIN ( T2, -1d-10 )
      func_cp = cp_ice - mlfus * Tf / TT

      END FUNCTION func_cp

!-----------------------------------------------------------------------
!     Heat capacity (true derivative)
!-----------------------------------------------------------------------

      FUNCTION func_cph(Tf,T)

      implicit none
      real(8)  func_cph, T, Tf
      real(8) :: TT

      TT = MIN ( T, -1d-10 ) ** 2
      func_cph = cp_ice - mlfus * Tf / TT

      END FUNCTION func_cph

!-----------------------------------------------------------------------
!     derivative of Heat capacity with respect to temperature T1
!-----------------------------------------------------------------------

      FUNCTION func_cpdt(Tf,T1,T2)

      implicit none
      real(8)  func_cpdt, T1, T2, Tf
      real(8) :: TT

      TT = MIN ( T1**2, -1d-20 ) * MIN ( T2, -1d-10 )
      func_cpdt = mlfus * Tf / TT

      END FUNCTION func_cpdt

!-----------------------------------------------------------------------
!     Energy of melt
!-----------------------------------------------------------------------

      FUNCTION func_qmelt(Tf,T)

      implicit none
      real(8)  func_qmelt, T, Tf

      func_qmelt = cp_ice * ( Tf - T ) + mlfus * ( 1.d0 - func_bf (Tf ,T ) )

      END FUNCTION func_qmelt

!-----------------------------------------------------------------------
!     Energy of melt
!-----------------------------------------------------------------------

      FUNCTION func_qm(Tf,T)

      implicit none
      real(8)  func_qm, T, Tf

      func_qm = cp_ice * ( Tf - T ) + mlfus * ( 1.d0 - func_bf (Tf ,T ) ) - cp_wat * Tf

      END FUNCTION func_qm

!-----------------------------------------------------------------------
!     Internal Energy BL99
!-----------------------------------------------------------------------
!     E = cp_ice*(1-     bf  ) Ti + mlfus*(1 -    bf  ) + cp_wat*    bf   *Ti 
!       = cp_ice*(1- -mu*S/Ti) Ti + mlfus*(1- -mu*S/Ti) + cp_wat*(-muS/Ti)*Ti
!         cp_ice*(Ti -   Tf  )    + mlfus*(1 - Tf/Ti)   + cp_wat*Tf
!-----------------------------------------------------------------------

      FUNCTION func_El(Tf,T)

      implicit none
      real(8)  func_El, Tf, T

      func_El = cp_ice * ( T - Tf ) - mlfus * ( 1.d0 - func_bf (Tf ,T ) ) + cp_wat * Tf

      END FUNCTION func_El

!-----------------------------------------------------------------------
!     Internal Energy mushy Layer
!-----------------------------------------------------------------------
!     E = phi cp_wat Ti + (1-phi) (cp_ice Ti - mlfus )
!-----------------------------------------------------------------------

      FUNCTION func_El_mush(T,p)

      implicit none
      ! arguments
      real(8)  func_El_mush, T, p ! T stands for Ti and p for phi (the liquid fraction)
      ! locals
      real(8) op
      op = 1.0_8 - p
      func_El_mush = cp_wat * T * p + op * ( cp_ice * T - mlfus )

      END FUNCTION func_El_mush

!-----------------------------------------------------------------------
!     Specific heat capacity mushy Layer
!-----------------------------------------------------------------------
!     C = dE/dT = phi cp_wat + (1-phi) cp_ice
!-----------------------------------------------------------------------

      FUNCTION func_cp_mush(T,p)

      implicit none
      real(8)  func_cp_mush, T, p ! T stands for Ti and p for phi

      func_cp_mush = cp_wat * p + (1.0_8-p) * cp_ice

      END FUNCTION func_cp_mush

!-----------------------------------------------------------------------
!     derivate of enthalpy relative to liquid fraction
!-----------------------------------------------------------------------
!     dE/dp = cp_wat Ti - (cp_ice Ti -l mlfus )
!-----------------------------------------------------------------------

      FUNCTION func_dedp(T,p)

      implicit none
      real(8)  func_dedp, T, p ! T stands for Ti and p for phi

      func_dedp = cp_wat * T - cp_ice * T + mlfus

      END FUNCTION func_dedp

!-----------------------------------------------------------------------
!     interface for liquidus temperatur/salinity relationship
!-----------------------------------------------------------------------
      FUNCTION func_liqu_sa(T)

      implicit none
      real(8)  func_liqu_sa, T ! T stands for Ti and p for phi

      select case (iliquid_cond)
      case(0)
       func_liqu_sa = func_liqu_sa0(T)
      case(1)
       func_liqu_sa = func_liqu_sa1(T)
      case(2)
       func_liqu_sa = func_liqu_sa2(T)
      end select

      END FUNCTION func_liqu_sa
!!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!     Liquidus relation giving salinity in brine of mushy Layer
!-----------------------------------------------------------------------
!     C_L = -17.6 * T - 0.389 * T^2 - 0.00362 * T^3
!-----------------------------------------------------------------------

      FUNCTION func_liqu_sa0(T)

      implicit none
      real(8)  func_liqu_sa0, T ! T stands for Ti and p for phi

      func_liqu_sa0 = 100.0_8

      END FUNCTION func_liqu_sa0

      FUNCTION func_liqu_sa1(T)

      implicit none
      real(8)  func_liqu_sa1, T ! T stands for Ti and p for phi

      func_liqu_sa1 = - 1.0_8 / 0.054_8 * T

      END FUNCTION func_liqu_sa1

      FUNCTION func_liqu_sa2(T)

      implicit none
      real(8)  func_liqu_sa2, T ! T stands for Ti and p for phi

      func_liqu_sa2 = - ( 17.6_8 + ( 0.389_8 + 0.00362_8 * T ) * T ) * T

      END FUNCTION func_liqu_sa2


!-----------------------------------------------------------------------
!     interface for liquidus temperatur/salinity first derivate
!-----------------------------------------------------------------------
      FUNCTION func_liqu_sadt(T)

      implicit none
      real(8)  func_liqu_sadt, T ! T stands for Ti and p for phi

      select case (iliquid_cond)
      case(0)
       func_liqu_sadt = func_liqu_sadt0(T)
      case(1)
       func_liqu_sadt = func_liqu_sadt1(T)
      case(2)
       func_liqu_sadt = func_liqu_sadt2(T)
      end select

      END FUNCTION func_liqu_sadt
!-----------------------------------------------------------------------
!     Liquidus relation giving salinity in brine of mushy Layer
!-----------------------------------------------------------------------
!     dC_L/dT = -17.6 - 2 * 0.389 * T - 3 * 0.00362 * T^2
!-----------------------------------------------------------------------

      FUNCTION func_liqu_sadt0(T)

      implicit none
      real(8)  func_liqu_sadt0, T ! T stands for Ti and p for phi

      func_liqu_sadt0 = 0.0_8

      END FUNCTION func_liqu_sadt0

      FUNCTION func_liqu_sadt1(T)

      implicit none
      real(8)  func_liqu_sadt1, T ! T stands for Ti and p for phi

      func_liqu_sadt1 = - 1.0_8 / 0.054

      END FUNCTION func_liqu_sadt1

      FUNCTION func_liqu_sadt2(T)

      implicit none
      real(8)  func_liqu_sadt2, T ! T stands for Ti and p for phi

      func_liqu_sadt2 = - ( 17.6_8 + ( 2.0_8 * 0.389_8 + 3.0_8 * 0.00362_8 * T ) * T )

      END FUNCTION func_liqu_sadt2

!-----------------------------------------------------------------------
!     Specific heat term in E
!-----------------------------------------------------------------------

      FUNCTION func_sh(Tf,T)

      implicit none
      real(8)  func_sh, T, Tf

      func_sh = cp_ice*(1d0-func_bf(Tf,T)) + cp_wat*func_bf(Tf,T)

      END FUNCTION func_sh

!-----------------------------------------------------------------------
!     Latent heat term in E
!-----------------------------------------------------------------------

      FUNCTION func_lh(Tf,T)

      implicit none
      real(8)  func_lh, Tf, T

      func_lh = -mlfus*(1d0-func_bf(Tf,T))

      END FUNCTION func_lh

!-----------------------------------------------------------------------
!     Brine fraction
!-----------------------------------------------------------------------

      FUNCTION func_bf(Tf,T)

      implicit none
      real(8)  func_bf, Tf, T

      func_bf = Tf / MIN ( T, -1d-20 )

      END FUNCTION func_bf

!-----------------------------------------------------------------------
!     Latent heat
!-----------------------------------------------------------------------

      FUNCTION func_le(T)

      implicit none
      real(8)  func_le, T
      real(8) :: Lsub = 2.834d+6     ! spec. latent heat of sublim, ice/snow [J/kg]

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
   real(8) em,si
   real(8) ztmelts,zaaa,zbbb,zccc,zdiscrim,tt
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
   real(8) em,si
   real(8) ztmelts,zaaa,zbbb,zccc,zdiscrim,tt
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

   !-------------------------------------------
   ! Effective permeability
   !-------------------------------------------
   subroutine ice_permeability(ni,phib,perm,dz)
   implicit none
   ! arguments
   integer ni ! total number of ice layer
   real(8), dimension(ni) :: phib,perm,dz ! value of mid-point (center of ice layer)
   ! locals
   integer k1,k2
   real(8) zbvf_min, zbvf_sum, zh_sum, zperm
        
        if ( i_perm_eff == 1 ) then ! Minimum
        
            do k1 = 1,ni
               zbvf_min = minval( phib(1:k1) )
               if ( i_perm_for == 1 ) then     ! Freitag
                  perm(k1) = 1.995e-8_8 * zbvf_min**3.1_8
               elseif ( i_perm_for == 2 ) then ! Rees Jones and Worster -> this case leads to bizzare results
                  perm(k1) = 1.0e-8_8 * zbvf_min**3
               endif
            enddo
            
        elseif ( i_perm_eff == 2 ) then ! Harmonic Mean
            
            do k1 = 1,ni
               zbvf_sum = 0.0_8
               zh_sum   = 0.0_8
               do k2=1,k1
                  if ( i_perm_for == 1 ) then     ! Freitag
                       zperm = 1.995e-8_8 * phib(k2)**3.1_8
                  elseif ( i_perm_for == 2 ) then ! Rees Jones and Worster -> this case leads to bizzare results
                       zperm = 1.0e-8_8 * phib(k2)**3
                  endif
                  zbvf_sum = zbvf_sum + dz(k2) / zperm
                  zh_sum = zh_sum + dz(k2)
               enddo
               zbvf_sum = zbvf_sum / zh_sum
               perm(k1) = 1.0_8 / zbvf_sum
            enddo
            
        endif ! condition on i_perm_eff
   end subroutine ice_permeability


   !-------------------------------------------
   ! Local Rayleigh number
   !-------------------------------------------
   subroutine ice_rayleigh_local(Sbr,Sw,perm,Ra,z)
   implicit none
   ! arguments
   real(8) Sbr,  & ! brine salinity in considered ice layer
           Sw,   & ! brine at ice base
           perm, & ! effective permeability
           Ra,   & ! Rayleigh number
           z       ! heigh from ice base of the ice layer
   ! locals
   real(8) z1, z2

        if ( i_Ra == 1 ) then ! formulation of Rees-Jones and Worster (JGR2014)
        
           z1 = cl * grav * beta_s / ( kl * visc )
           z2 = Sbr - Sw
        
        elseif( i_Ra == 2 ) then ! formulation of Vancoppenolle et al (TCD2013)
        
           z1 = grav / ( kappa * mudyn )
           z2 = B_S * ( Sbr - Sw )
           
        endif
        
        Ra = z1 * z2 * z * perm

   end subroutine ice_rayleigh_local

   !---------------------------------------------------------------
   ! Critical depth, effective Rayleigh number, vertical velocity
   !---------------------------------------------------------------
   subroutine vertical_brine_velocity_RJW(ni,Ra,z,w_br,u_br)
   ! arguments
   integer ni ! total number of ice layer
   real(8), dimension(ni) :: Ra, z, u_br ! value of mid-point (center of ice layer)
   real(8), dimension(0:ni) :: w_br ! Darcy velocity (positive when convecting as new brine moves upward and replaces the part being flushed)
   ! locals
   integer k,kc
   real(8) Rae,zc,wamp

        ! RJW
        ! 1) if Ra is everywhere < Rc => no convection Rae =0
        ! 2) else: convection until z_c

        zc  = 0.0_8
        kc  = 0
        Rae = 0.0_8
        do k=1,ni
           if ( Ra(k) >= Rc_RJW ) then 
             kc = k
             zc = z(k) ! heigh at top of layer k
           endif
        enddo
        if (kc > 0) Rae = maxval( Ra(1:kc) ) - Rc_RJW
        wamp = Rae * al_RJW * kl / cl

        w_br(0:ni) = 0.0_8
        do k=1,kc
           w_br(k) = - wamp * ( z(k) - zc ) / zc**2
        enddo
        w_br(0) = 2.0_8 * w_br(1) - w_br(2) ! linear extrapolation for ice base
        do k=1,ni
           u_br(k) = w_br(k-1) - w_br(k)
        enddo

   end subroutine vertical_brine_velocity_RJW

   !---------------------------------------------------------------
   ! Critical depth, effective Rayleigh number, vertical velocity
   !---------------------------------------------------------------
   subroutine vertical_brine_velocity_GN(ni,Ra,dz,ub,wb)
   ! arguments
   integer ni ! total number of ice layer
   real(8), dimension(ni) :: Ra, dz, ub ! value of mid-point (center of ice layer)
   real(8), dimension(0:ni) :: wb ! Darcy velocity (positive when convecting 
                                  ! as new brine moves upward and replaces the part being flushed)
   ! locals
   integer k

       do k = 1,ni
          ub(k) = alpha_GN / rho_br_GN * max( 0.0_8, Ra(k) - Rc_GN ) * dz(k)
       enddo
       wb(ni)=0
       do k=ni,1,-1
          wb(k-1) = wb(k) + ub(k)
       enddo       

   end subroutine vertical_brine_velocity_GN

   !---------------------------------------------------------------
   ! Reject liquid fraction when fraction is too high
   !---------------------------------------------------------------
   subroutine reject_brine_velocity(ni,phib,dz,um,dt)
   ! arguments
   integer ni ! total number of ice layer
   real(8), dimension(ni) :: phib, dz, um ! value of mid-point (center of ice layer)
   real(8) dt
   ! locals
   integer k
   real(8) :: phibc = 0.8_8, &      ! critical liquid fraction above which drainage occurs
              taubc = 3600.0_8    ! restoring time for drainage in seconds (=1h)

       taubc = max(2.0_8*dt, taubc) ! make sure that the rate of flushing is not unsustainable

       um(1:ni) = 0.0_8
       do k = 1,ni
          if ( phib(k) > phibc ) um(k) = dz(k) / taubc
       enddo

   end subroutine reject_brine_velocity

end module var_thermo_vertical
