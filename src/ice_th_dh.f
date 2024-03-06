      SUBROUTINE ice_th_dh(nlay_s,nlay_i,kideb,kiut)
       !!------------------------------------------------------------------
       !!                ***         ROUTINE ice_th_dh         ***
       !! ** Purpose :
       !!           This routine determines variations of ice and snow thicknesses.
       !! ** Method  :
       !!           Ice/Snow surface melting arises from imbalance in surface fluxes
       !!           Bottom accretion/ablation arises from flux budget
       !!           Snow thickness can increase by precipitation and decrease by 
       !!              sublimation
       !!           If snow load excesses Archmiede limit, snow-ice is formed by
       !!              the flooding of sea-water in the snow
       !! ** Steps  
       !!           1) Compute available flux of heat for surface ablation
       !!           2) Compute snow and sea ice enthalpies
       !!           3) Surface ablation and sublimation
       !!           4) Bottom accretion/ablation
       !!           5) Case of Total ablation
       !!           6) Snow ice formation
       !!
       !! ** Inputs / Outputs
       !!
       !! ** External
       !!
       !! ** References : Bitz and Lipscomb, JGR 99
       !!                 Fichefet T. and M. Maqueda 1997, 
       !!                 J. Geophys. Res., 102(C6), 12609-12646   
       !!                 Vancoppenolle, Fichefet and Bitz, GRL 2005
       !!                 Vancoppenolle et al., OM08
       !!
       !! ** History  : 
       !!   original code    01-04 (LIM)
       !!   original routine
       !!               (05-2003) M. Vancoppenolle, Louvain-La-Neuve, Belgium
       !!               (05-2008) BIO-LIM
       !!
       !!------------------------------------------------------------------
       !! * Arguments
       !!------------------------------------------------------------------

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'

      ! Local Variables
      REAL(8) :: zrchu1(nbpt), zrchu2(nbpt), zqsat(nbpt), z_f_surf(nbpt) 
      REAL(8) :: zdeltah(maxnlay)
      LOGICAL l_write
      REAL(8) :: zgrr

      zqt_s_ini = 0.0
      zqt_s_fin = 0.0
      zdqt_s    = 0.0
      zqt_i_ini = 0.0
      zqt_i_fin = 0.0
      zdqt_i    = 0.0
      s_i_max   = 15.0

      ! Local Constants
      zeps = 1.0e-20
      l_write = .TRUE.

      IF ( l_write ) THEN
         WRITE(numout,*) ' ** ice_th_dh : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~ '
         WRITE(numout,*) 
      ENDIF
!
!------------------------------------------------------------------------------|
!  1) Calculate available heat for surface ablation         
!------------------------------------------------------------------------------|
!
      DO 20 ji = kideb, kiut
      
      z_f_surf(ji)  = fratsb(ji) + fleb(ji) + fcsb(ji) - fc_su(ji)
     &              + ab(ji)*fsolgb(ji)
      z_f_surf(ji)  = MAX(zero,z_f_surf(ji))
      z_f_surf(ji)  = z_f_surf(ji)*MAX(zero,
     &                SIGN(one,t_su_b(ji)-tfs(ji)))

      IF ( l_write ) THEN
         WRITE(numout,*) ' Available heat for surface ablation ... '
         WRITE(numout,*) 
         WRITE(numout,*) ' z_f_surf : ', z_f_surf(ji)
         WRITE(numout,*) ' fratsb   : ', fratsb(ji)
         WRITE(numout,*) ' fleb     : ', fleb(ji)
         WRITE(numout,*) ' fcsb     : ', fcsb(ji)
         WRITE(numout,*) ' fc_su    : ', fc_su(ji)
         WRITE(numout,*) ' ab       : ', ab(ji) 
         WRITE(numout,*) ' fsolgb   : ', fsolgb(ji)
         WRITE(numout,*)
         WRITE(numout,*) ' ht_i_b    : ', ht_i_b(ji)
         WRITE(numout,*) ' ht_s_b    : ', ht_s_b(ji)
         WRITE(numout,*) ' t_su_b   : ', t_su_b(ji)
         WRITE(numout,*)
      ENDIF

 20   CONTINUE

!
!------------------------------------------------------------------------------|
!  2) Snowfall and surface melt                                                |
!------------------------------------------------------------------------------|
!
      DO 40 ji = kideb, kiut

      ! total snow heat content for conservation
      zqt_s_ini = q_s_b(ji,1) * ht_s_b(ji)

      IF ( l_write ) THEN
         WRITE(numout,*) ' Surface ablation and sublimation ... '
         WRITE(numout,*) 
         WRITE(numout,*) ' zqt_s_ini : ', zqt_s_ini / ddtb
         WRITE(numout,*) ' ht_s_b    : ', ht_s_b(ji)
         WRITE(numout,*) ' q_s_b(1)  : ', q_s_b(ji,1)
         WRITE(numout,*)
      ENDIF

      !----------
      ! Snowfall
      !----------
      dh_s_prec(ji)    =  hnpbqb(ji)
      dh_s_melt(ji)    =  0.0
      zqprec           =  rhon * ( cpg * ( tpw - tabqb(ji) ) + xlgn )

      ! Conservation update
      zqt_s_ini        =  zqt_s_ini + zqprec*hnpbqb(ji)
      fprec            =  - zqprec * hnpbqb(ji) / ddtb

      IF ( l_write ) THEN
         WRITE(numout,*) ' snow falls! '
         WRITE(numout,*) ' dh_s_prec : ', dh_s_prec(ji)
         WRITE(numout,*) ' flux of h : ',
     &                   zqprec*hnpbqb(ji) / ddtb
         WRITE(numout,*) ' zqt_s_ini : ', zqt_s_ini / ddtb
         WRITE(numout,*)
      ENDIF

      !-----------
      ! Snow melt
      !-----------
      ! Melt of fallen snow
      zqfont_su        =  z_f_surf(ji) * ddtb 
      IF ( l_write ) WRITE(numout,*) ' snow melts! '
      IF ( l_write ) WRITE(numout,*) ' zqfont_su : ', zqfont_su / ddtb 

      zdeltah(1)       =  MIN( 0.0 , - zqfont_su / 
     &                    MAX( zqprec , zeps ) )
      zqfont_su        =  MAX( 0.0 , - dh_s_prec(ji) - zdeltah(1) ) * 
     &                    zqprec
      zdeltah(1)       =  MAX( - dh_s_prec(ji), zdeltah(1) )
      dh_s_melt(ji)    =  dh_s_melt(ji) + zdeltah(1)

      ! Melt of snow
      DO layer = 1, nlay_s !in case of melting of more than 1 layer
         zdeltah(layer) =  - zqfont_su / q_s_b(ji,layer)
         zqfont_su      = MAX( zero, - deltaz_s_phy(layer) - 
     &                    zdeltah(layer) ) * 
     &                    q_s_b(ji,layer)
         zdeltah(layer) = MAX( zdeltah(layer), - deltaz_s_phy(layer) )
         dh_s_melt(ji) =  dh_s_melt(ji) + zdeltah(layer) !resulting melt of snow    
      END DO
      dh_s_tot(ji)     =  dh_s_melt(ji) + dh_s_prec(ji)

      ! old and new snow thicknesses
      hsold            =  ht_s_b(ji)
      hsnew            =  ht_s_b(ji) + dh_s_tot(ji)

      ! if snow is still present zhn = 1, else zhn = 0
      zhn              =  1.0 - MAX( zero , SIGN( one , - hsnew ) )
      ht_s_b(ji)       =  MAX( zero , hsnew )

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Conservation test for snow
      zqt_s_fin = q_s_b(ji,1) * ht_s_b(ji)
      zdqt_s = zqt_s_fin - zqt_s_ini

      WRITE(numout,*)
      WRITE(numout,*) ' Conservation in snow... '
      WRITE(numout,*) ' dh_s_melt : ', dh_s_melt(ji)
      WRITE(numout,*) ' dh_s_prec : ', dh_s_prec(ji)
      WRITE(numout,*) ' ht_s_b    : ', ht_s_b(ji)
      WRITE(numout,*) ' zqt_s_ini : ', zqt_s_ini / ddtb
      WRITE(numout,*) ' zqt_s_fin : ', zqt_s_fin / ddtb
      WRITE(numout,*) ' zdqt_s    : ', zdqt_s / ddtb
      WRITE(numout,*) ' z_f_surf  : ', - z_f_surf(ji)
      IF ( zqt_s_fin.GT.0.0 ) THEN
         cons_err = ABS(zdqt_s / ddtb  + z_f_surf(ji) )
      ELSE
         cons_err = ABS(zqt_s_ini / ddtb + zdqt_s / ddtb )
      ENDIF
      WRITE(numout,*) ' Cons error, snow : ', cons_err
      WRITE(numout,*)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !------------------
      ! Ice surface melt
      !------------------
      IF ( l_write ) WRITE(numout,*) ' ice melts!  '

      zzf_surf = zqfont_su / ddtb
      zdqt_i            = 0.0
      dh_i_surf(ji)     =  0.0
      DO layer = 1, nlay_i
         zdeltah(layer) =  - zqfont_su / q_i_b(ji,layer)
         zqfont_su      =  MAX( zero , - deltaz_i_phy(layer) 
     &                  - zdeltah(layer) ) *  q_i_b(ji,layer)
         zdeltah(layer) =  MAX( zdeltah(layer) , - deltaz_i_phy(layer) )
         dh_i_surf(ji)  =  dh_i_surf(ji) + zdeltah(layer) !resulting melt of ice
         zdqt_i         = zdqt_i  + zdeltah(layer) * q_i_b(ji,layer) 
     &                  / ddtb
      END DO

      cons_err = ABS( zzf_surf + zdqt_i )

      IF ( l_write ) THEN
         WRITE(numout,*) ' Conservation in sea ice, surface '
         WRITE(numout,*) ' dh_i_surf: ', dh_i_surf(ji)
         WRITE(numout,*) ' ht_i_b   : ', ht_i_b(ji)
         WRITE(numout,*) ' zzf_surf : ', zzf_surf
         WRITE(numout,*) ' zdqt_i   : ', zdqt_i
         WRITE(numout,*) ' Cons error, ice : ', cons_err
         WRITE(numout,*)
      ENDIF
      
!
!------------------------------------------------------------------------------|
!  3) Sublimation at the surface                                               |
!------------------------------------------------------------------------------|
!
      !------------------
      ! Snow sublimation
      !------------------

      !------------------
      ! Ice sublimation
      !------------------

      !-------------------
      ! Snow condensation
      !-------------------

      ! 4.3) Snow/ice sublimation
      !
      ! If fleb is negative, snow condensates at the surface.
      !
      dh_s_subl(ji)    =  - parsub*fleb(ji)/(rhon*xsn)*ddtb
      dh_s_tot(ji)     =  dh_s_tot(ji) + dh_s_subl(ji)
      zdhcf            =  ht_s_b(ji) + dh_s_subl(ji) 
      ht_s_b(ji)       =  MAX(zero,zdhcf)
      dh_s_tot(ji)     =  ht_s_b(ji) - hsold
      dh_i_subl(ji)    =  - MAX(zero,-zdhcf)*rhon/rhog

      dh_i_surf(ji)    =  dh_i_surf(ji) + dh_i_subl(ji)

      hsnew            =  ht_s_b(ji)

      IF (ht_s_b(ji).le.0.0) THEN
         dh_s_tot(ji) =  MAX ( 0.0 , dh_s_tot(ji) )
      ENDIF

      IF ( l_write ) THEN
         WRITE(numout,*) ' Snow sublimation ... '
         WRITE(numout,*) ' '
         WRITE(numout,*) ' parsub : ', parsub
         WRITE(numout,*) ' dh_s_subl : ', dh_s_subl(ji)
         WRITE(numout,*) ' dh_i_subl : ', dh_i_subl(ji)
      ENDIF

 40   CONTINUE
!
!------------------------------------------------------------------------------|
!  4) Basal growth and melt                                                    |
!------------------------------------------------------------------------------|
!
      DO 50 ji = kideb, kiut

      IF ( l_write ) THEN
         WRITE(numout,*) ' Basal growth and melt ... '
         WRITE(numout,*) 
         WRITE(numout,*) ' fbbqb     : ', fbbqb(ji)
         WRITE(numout,*) ' fc_bo_i   : ', fc_bo_i(ji)
         WRITE(numout,*)
      ENDIF

      ! formation / melting of ice at the base is determined by the balance of
      ! the conductive heat flux in the ice (fc_bo_i), and the heat fluxes 
      ! from the ocean (fbbqb). Conductive heat flux is positive 
      ! downwards and fbbq is positive to the ice, i.e., upwards.
      ! Melt/formation rate is modulated by the enthalpy of the bottom ice layer. 
      ! accretion of ice at the base

      !--------------
      ! Basal growth
      !--------------
      ! Initial value (tested 1D, can be anything between 1 and 20)
      s_i_new      = 10.0
      num_iter_max = 10

      ! the growth rate (dh_i_bott) is function of the new ice
      ! heat content (q_i_b(nlay_i+1)). q_i_b depends on the new ice
      ! salinity (snewice). snewice depends on dh_i_bott
      ! it converges quickly, so, no problem

      ! Iterative procedure
      IF ( ( fc_bo_i(ji) + fbbqb(ji) ) .LT. 0.0 ) THEN

         IF ( l_write ) WRITE(numout,*) 
     &            ' Energy available for basal growth : ',
     &                  fc_bo_i(ji) + fbbqb(ji)

         DO iter = 1, num_iter_max
               tmelts             =   - tmut * s_i_new + tpw ! Melting point in K
               q_i_b(ji,nlay_i+1) = rhog*( cpg*(tmelts-t_bo_b(ji))
     &                 + xlgm*( 1.0-(tmelts-tpw)
     &                 / (t_bo_b(ji) - tpw) )
     &                 - cpw*(tmelts-tpw) ) 
               ! Bottom growth rate = - F*dt / q
               dh_i_bott(ji)    =  - ddtb*(fc_bo_i(ji) + fbbqb(ji) )
     &                             / q_i_b(ji,nlay_i+1)
         IF ( l_write ) WRITE(numout,*) 'qi',q_i_b(ji,nlay_i+1)
     &,dh_i_bott(ji) 
               !
               ! New ice salinity ( Cox and Weeks, JGR, 1988 )
               !
               ! zswi2  (1) if dh_i_bott/rdt .GT. 3.6e-7
               ! zswi12 (1) if dh_i_bott/rdt .LT. 3.6e-7 and .GT. 2.0e-8
               ! zswi1  (1) if dh_i_bott/rdt .LT. 2.0e-8
               !
               zgrr   = MIN( 1.0e-3, 
     &                  MAX ( dh_i_bott(ji) / ddtb , zeps ) )
               zswi2  = MAX( 0.0 , SIGN( 1.0 , zgrr - 3.6e-7 ) )
               zswi12 = MAX( 0.0 , SIGN( 1.0 , zgrr - 2.0e-8 ) ) * 
     &                  ( 1.0 - zswi2 )
               zswi1  = 1. - zswi2 * zswi12
               zfracs = zswi1  * 0.12 +   
     &                  zswi12 * ( 0.8925 + 0.0568 * 
     &                  LOG( 100.0 * zgrr ) ) +  
     &                  zswi2  * 0.26 /   
     &                  ( 0.26 + 0.74 * EXP ( - 724300.0 * zgrr ) )
               zfracs = 1.
               zds     = zfracs * oce_sal - s_i_new
               s_i_new = zfracs * oce_sal
         IF ( l_write ) WRITE(numout,*) 'si', s_i_new, zds, oce_sal

               ! the brine volume in the skeletal layer is equal to f
               e_skel  = zfracs

               ! salt flux due to initial salt entrapment
               fsbp = oce_sal * ( 1. - zfracs ) * dh_i_bott(ji) / ddtb *
     &                rhog / 1000.

!              WRITE(numout,*)
!              WRITE(numout,*) ' ddtb      : ', ddtb
!              WRITE(numout,*) ' zgrr      : ', zgrr
!              WRITE(numout,*) ' zswi12    : ', zswi12
!              WRITE(numout,*) ' zswi1     : ', zswi1
!              WRITE(numout,*) ' zswi2     : ', zswi2
!              WRITE(numout,*) ' dh_i_bott : ', dh_i_bott(ji)
!              WRITE(numout,*) ' oce_sal   : ', oce_sal
!              WRITE(numout,*) ' zfracs    : ', zfracs 
!              WRITE(numout,*) ' s_i_new   : ', s_i_new
!              WRITE(numout,*)

         END DO ! iter
      ENDIF ! fc_bo_i

      IF ( l_write ) THEN
         WRITE(numout,*) ' e_skel : ', e_skel
         WRITE(numout,*)
      ENDIF

      ! Final values
      IF ( (fc_bo_i(ji)+fbbqb(ji)) .LT. 0.0 ) THEN
      ! New ice salinity must not exceed 15 psu
         zoldsinew   = s_i_new
         s_i_new     = MIN( s_i_new , s_i_max )
         ! Metling point in K
         tmelts             =   - tmut * s_i_new + tpw
         ! New ice heat content (Bitz and Lipscomb, 1999)
         q_i_b(ji,nlay_i+1) = rhog*( cpg*(tmelts-t_bo_b(ji))
     &                 + xlgm*( 1.0-(tmelts-tpw)
     &                 / (t_bo_b(ji) - tpw) )
     &                 - cpw*(tmelts-tpw) ) 
         dh_i_bott(ji)    =  - ddtb*(fc_bo_i(ji) + fbbqb(ji) )
     &                    / q_i_b(ji,nlay_i+1)
         IF ( l_write ) WRITE(numout,*) ' dh_i_bott : ', dh_i_bott(ji)

      ENDIF 

      !-----------------
      ! Basal melt
      !-----------------
      IF ( ( fc_bo_i(ji) + fbbqb(ji) ) .GE. 0.0 ) THEN

         IF ( l_write ) WRITE(numout,*) ' Energy available for 
     &                                    basal melt   : ',
     &                  fc_bo_i(ji) + fbbqb(ji)

         zqfont_bo   = ddtb * ( fc_bo_i(ji) + fbbqb(ji) )
         zzf_base    = zqfont_bo / ddtb
         zdqt_i      = 0.0

         IF ( l_write ) WRITE(numout,*) ' zqfont_bo : ', zqfont_bo

         dh_i_bott(ji)     =  0.0
         DO layer = nlay_i, 1, -1
            zdeltah(layer) =  - zqfont_bo / q_i_b(ji,layer)
            zqfont_bo      =  MAX ( 0.0 , - deltaz_i_phy(layer) -
     &                        zdeltah(layer) )
     &                     *  q_i_b(ji,layer)
            dh_i_bott(ji)  =  dh_i_bott(ji) + zdeltah(layer)
            zdqt_i         =  zdqt_i + zdeltah(layer) * 
     &                        q_i_b(ji,layer) / ddtb
         END DO

         IF ( l_write ) WRITE(numout,*) ' dh_i_bott : ', dh_i_bott(ji)

         cons_err = ABS( zzf_base + zdqt_i )

         IF ( l_write ) THEN
            WRITE(numout,*) ' Conservation in sea ice, base '
            WRITE(numout,*) ' dh_i_bott: ', dh_i_bott(ji)
            WRITE(numout,*) ' ht_i_b   : ', ht_i_b(ji)
            WRITE(numout,*) ' zzf_base : ', zzf_base
            WRITE(numout,*) ' zdqt_i   : ', zdqt_i
!           WRITE(numout,*) ' Conservation error  ice surface : ', 
!    &                      cons_err
!           WRITE(numout,*)
         ENDIF

      ENDIF

      ! It can be than an internal temperature is greater than melt point
      ! then, see lim3 for correction

      ! new ice thickness
      zhgnew         = ht_i_b(ji) + dh_i_surf(ji) + dh_i_bott(ji)
      old_ht_i_b(ji) = ht_i_b(ji)

      ht_i_b(ji) = zhgnew

      !--------------------------------------
      ! Meltwater flow due to surface melt
      !--------------------------------------
      qsummer = ( - rhog * MIN ( dh_i_surf(ji) , 0.0 ) 
     &            - rhon * MIN ( dh_s_melt(ji) , 0.0 ) )  
      IF ( l_write ) THEN
         WRITE(numout,*) ' qsummer : ', qsummer
         WRITE(numout,*)
      ENDIF

 50   CONTINUE

!
!------------------------------------------------------------------------------|
!  5) Formation of snow-ice                                                    |
!------------------------------------------------------------------------------|
!
      ! When snow load excesses Archimede's limit, snow-ice interface goes down
      ! under sea-level, flooding of seawater transforms snow into ice
      ! dh_snowice is positive for the ice

      DO 70 ji = kideb, kiut

      dh_snowice(ji) = MAX( zero , ( rhon * ht_s_b(ji) + (rhog - rho0 )
     &                 * ht_i_b(ji)) / ( rhon + rho0 - rhog ) )
        
      zhgnew         = MAX( zhgnew , zhgnew + dh_snowice(ji) )
      zhnnew         = MIN( ht_s_b(ji) , ht_s_b(ji) - dh_snowice(ji) )

      ht_s_b(ji)  = zhnnew
      ht_i_b(ji)  = zhgnew

      IF ( l_write ) THEN
         WRITE(numout,*) ' dh_snowice : ', dh_snowice(ji)
         WRITE(numout,*)
         WRITE(numout,*) ' At the end of the routine ... '
         WRITE(numout,*) ' ht_s_b : ', ht_s_b(ji)
         WRITE(numout,*) ' ht_i_b : ', ht_i_b(ji)
      ENDIF

 70   CONTINUE

      RETURN

!------------------------------------------------------------------------------|
! Fin de la subroutine ice_th_dh
      END SUBROUTINE
