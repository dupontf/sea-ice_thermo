      SUBROUTINE ice_th_glohec(eti,ets,etilayer,kideb,kiut,jl, 
     &                         nlay_s, nlay_i)

      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_th_glohec *** 
      !!                 
      !! ** Purpose :  Compute total heat content for each category
      !!               Works with 1d vectors only
      !!
      !! history :
      !!  9.9  ! 07-04 (M.Vancoppenolle) original code
      !!-----------------------------------------------------------------------
      !! * Local variables

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'

      INTEGER, INTENT(in) ::  
     &   kideb, kiut,           ! bounds for the spatial loop
     &   jl                     ! category number

      REAL(8), DIMENSION (nbpt,jpl), INTENT(out) ::   
     &   eti, ets            ! vertically-summed heat content for ice /snow

      REAL(8), DIMENSION (nbpt,maxnlay), INTENT(out) ::   
     &   etilayer            ! heat content for ice layers

      REAL(8) ::  
     &   zdes,               ! snow heat content increment (dummy)
     &   zeps                ! very small value (1.e-10)

      INTEGER  ::  
     &   ji,jj,jk            ! loop indices

      !!-----------------------------------------------------------------------
      eti(:,:) = 0.0
      ets(:,:) = 0.0
      zeps     = 1.0e-10

      ! total q over all layers, ice [J.m-2]
      DO jk = 1, nlay_i
         DO ji = kideb, kiut
            etilayer(ji,jk) = q_i_b(ji,jk)  
     &                      * ht_i_b(ji) / nlay_i
            eti(ji,jl) = eti(ji,jl) + etilayer(ji,jk) 
         END DO
      END DO

      ! total q over all layers, snow [J.m-2]
      DO ji = kideb, kiut
         zdes = q_s_b(ji,1) * ht_s_b(ji) / nlay_s 
         ets(ji,jl) = ets(ji,jl) + zdes        
      END DO

      RETURN
!------------------------------------------------------------------------------      
      END SUBROUTINE ice_th_glohec

!==============================================================================

      SUBROUTINE ice_th_enmelt(kideb, kiut, nlay_s, nlay_i)
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_th_enmelt *** 
      !!                 
      !! ** Purpose :   Computes sea ice energy of melting q_i (J.m-3)
      !!
      !! ** Method  :   Formula (Bitz and Lipscomb, 1999)
      !!
      !! history : Martin Vancoppenolle, May 2007
      !!-------------------------------------------------------------------

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'

      REAL(8)                 ::       !: goes to trash
     &   ztmelts               ,       !: sea ice freezing point in K
     &   zeps 

      INTEGER                 ::     
     &   ji,                           !: spatial loop index
     &   jk                            !: vertical index

      LOGICAL                 ::
     &   ln_write = .FALSE.

      !!-------------------------------------------------------------------
      zeps = 1.0e-10
      IF ( ln_write ) THEN
         WRITE(numout,*) ' ** ice_th_enmelt : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*)
         WRITE(numout,*) ' kideb : ', kideb
         WRITE(numout,*) ' kiut  : ', kiut
         WRITE(numout,*) ' nlay_s: ', nlay_s
         WRITE(numout,*) ' nlay_i: ', nlay_i
      ENDIF

      ! Sea ice energy of melting
      DO jk = 1, nlay_i
         DO ji = kideb, kiut
            ztmelts      =   - tmut * s_i_b(ji,jk) + tpw 
            q_i_b(ji,jk) = rhog * ( cpg * ( ztmelts - t_i_b(ji,jk) )
     &                   + xlgm*( 1.0 - (ztmelts-tpw) / 
     &                     MIN((t_i_b(ji,jk)-tpw),-zeps) )  
     &                   - cpw      * ( ztmelts-tpw  ) ) 
         END DO !ji
      END DO !jk

      ! Snow energy of melting
      DO jk = 1, nlay_s
         DO ji = kideb,kiut
            q_s_b(ji,jk) = rhon * ( cpg  * ( tpw - t_s_b(ji,jk) ) + 
     &                     xlgm )
         END DO !ji
      END DO !jk

      WRITE(numout,*)

      RETURN
!------------------------------------------------------------------------------      
      END SUBROUTINE ice_th_enmelt

!==============================================================================

      SUBROUTINE ice_th_con_dif(kideb,kiut,nlay_s,nlay_i,jl)
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_th_con_dif *** 
      !!                 
      !! ** Purpose :   Test energy conservation after heat diffusion
      !!
      !! history :
      !!  9.9  ! 07-04 (M.Vancoppenolle) original code
      !!-------------------------------------------------------------------
      !! * Local variables

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'

      REAL(8)                  ::      !: ! goes to trash
     &   meance,                       !: mean conservation error
     &   max_cons_err,                 !: maximum tolerated conservation error
     &   max_surf_err                  !: maximum tolerated surface error

      INTEGER ::                     
     &   numce                         !: number of points for which conservation
                                       !  is violated
      INTEGER  ::  
     &   ji,jj,jk,                     !: loop indices
     &   zji, zjj

      LOGICAL  ::
     &   l_write
      !!---------------------------------------------------------------------

      max_cons_err =  0.1
      max_surf_err =  0.001
      l_write  = .TRUE.

      IF ( l_write ) THEN
         WRITE(numout,*) ' ** ice_th_con_dif : '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~ '
      ENDIF

      !--------------------------
      ! 1) Increment of energy
      !--------------------------
      ! vertically integrated
      DO ji = kideb, kiut
          dq_i(ji,jl) = qt_i_fin(ji,jl) - qt_i_in(ji,jl)   
     &                + qt_s_fin(ji,jl) - qt_s_in(ji,jl)
      END DO

      ! layer by layer
      DO ji = kideb, kiut
         DO jk = 1, nlay_i
            dq_i_layer(ji,jk) =  q_i_layer_fin(ji,jk) - 
     &                           q_i_layer_in(ji,jk)
         END DO
      END DO

      !-------------------------------------------
      ! 2) Atmospheric heat flux, ice heat budget
      !-------------------------------------------

      DO ji = kideb, kiut
          fatm(ji,jl) = ab(ji)*fsolgb(ji) + fratsb(ji)
     &                + fcsb(ji) + fleb(ji)
          sum_fluxq(ji,jl) = fc_su(ji) - fc_bo_i(ji)
     &                + fsolgb(ji)*(1-ab(ji)) - ftroce
      END DO

      IF ( l_write ) THEN
         DO ji = kideb, kiut
             WRITE(numout,*) ' fc_su   : ', fc_su(ji)
             WRITE(numout,*) ' fc_bo_i : ', fc_bo_i(ji)
             WRITE(numout,*) ' fsol*io : ', fsolgb(ji)*(1.-ab(ji))
             WRITE(numout,*) ' ftroce  : ', ftroce
         END DO
      ENDIF

      !-----------------------
      ! 3) Conservation error
      !-----------------------
      DO ji = kideb, kiut
          cons_error(ji,jl) = ABS( dq_i(ji,jl) / ddtb + 
     &                        sum_fluxq(ji,jl) )
      END DO

      numce = 0
      meance = 0.0
      DO ji = kideb, kiut
          IF ( cons_error(ji,jl) .GT. max_cons_err ) THEN
              numce = numce + 1
              meance = meance + cons_error(ji,jl)
          ENDIF
      IF (numce .GT. 0 ) meance = meance / numce

      IF ( cons_error(ji,jl) .GT. max_cons_err ) THEN
         WRITE(numout,*) ' Diffusion of heat - large cons. error '
         WRITE(numout,*) ' Maximum tolerated conservation error : ', 
     &                     max_cons_err
         WRITE(numout,*) ' Thickness category : ', jl
         WRITE(numout,*) ' DIF Mean conservation error ', 
     &                     meance, numit
      ENDIF

      ENDDO

      !---------------------------------------------------------
      ! 4) Surface error due to imbalance between Fatm and Fcsu
      !---------------------------------------------------------
      numce  = 0.0
      meance = 0.0

      DO ji = kideb, kiut
         surf_error(ji,jl) = ABS ( fatm(ji,jl) - fc_su(ji) )
         IF ( ( t_su_b(ji) .LT. tpw ) .AND. ( surf_error(ji,jl) .GT.  
     &          max_surf_err ) ) THEN
            numce = numce + 1 
            meance = meance + surf_error(ji,jl)
         ENDIF
         IF (numce .GT. 0 ) meance = meance / numce

         IF ( ( t_su_b(ji) .LT. tpw ) .AND. ( surf_error(ji,jl) .GT.  
     &             max_surf_err ) ) THEN
            WRITE(numout,*) ' Diffusion of heat - large surface error '
            WRITE(numout,*) ' Surf_error : ', surf_error(1,jl), 
     &                      ' (W/m2) '
            WRITE(numout,*) ' Maximum tolerated surface error : ', 
     &                     max_surf_err
            WRITE(numout,*) ' Thickness category : ', jl
            WRITE(numout,*) ' t_su_b: ', t_su_b(ji) 
            WRITE(numout,*) ' t_s_b : ', ( t_s_b(ji,layer), 
     &                      layer = 1, nlay_s )
            WRITE(numout,*) ' t_i_b : ', ( t_i_b(ji,layer), 
     &                      layer = 1, nlay_i )
         ENDIF

      ENDDO

      !-------------------
      ! 5) Layer by layer
      !-------------------
      IF ( l_write ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' Layer by layer ... '

         DO ji = kideb, kiut
            WRITE(numout,*) ' T_su : ', t_su_b(ji)
            WRITE(numout,*) ' *** Snow ***    h_s = ',
     &                  ht_s_b(ji) 

            IF ( ht_s_b(ji) .GT. 0.0 ) THEN
            zdqs = ( qt_s_fin(ji,jl) - qt_s_in(ji,jl) ) / ddtb
            WRITE(numout,*) ' Conservation error (W/m2): ', 
     &                   ABS( zdqs + fc_s(ji,0) - 
     &                   fc_s(ji,1) + radab_s(1) )
            ENDIF
         END DO

         DO ji = kideb, kiut
            WRITE(numout,*) ' *** Ice ***     h_i = ', ht_i_b(ji) 
            DO jk = 1, nlay_i
               WRITE(numout,*) ' LAYER : ', jk
               WRITE(numout,*) ' Conservation error (W/m2): ', 
     &            ABS( dq_i_layer(ji,jk) / ddtb + fc_i(ji,jk-1) - 
     &            fc_i(ji,jk) + radab_phy_i(jk) + radab_alg_i(jk) )
            END DO
         END DO

      ENDIF
!
!------------------------------------------------------------------------------      
      END SUBROUTINE ice_th_con_dif
!==============================================================================

      SUBROUTINE ice_th_con_dh(kideb,kiut,nlay_s,nlay_i,jl)
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_th_con_dh ***  
      !!                 
      !! ** Purpose :   Test energy conservation after ice growth and melt
      !!
      !! history :
      !!  9.9  ! 07-04 (M.Vancoppenolle) original code
      !!-------------------------------------------------------------------
      !! * Local variables

      INCLUDE 'type.com'
      INCLUDE 'para.com'
      INCLUDE 'const.com'
      INCLUDE 'bloc.com'
      INCLUDE 'ice.com'
      INCLUDE 'thermo.com'

      REAL(8)                  ::      !: ! goes to trash
     &   meance,                       !: mean conservation error
     &   max_cons_err,                 !: maximum tolerated conservation error
     &   max_surf_err                  !: maximum tolerated surface error

      INTEGER ::                     
     &   numce                         !: number of points for which conservation
                                       !  is violated
      INTEGER  ::  
     &   ji,jj,jk,                     !: loop indices
     &   zji, zjj
      LOGICAL  ::
     &   l_write

      !!---------------------------------------------------------------------

      l_write  = .TRUE.
      max_cons_err =  0.1
      max_surf_err =  0.001

      WRITE(numout,*) ' ** ice_th_con_dh : '
      WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~ '
            
      !------------------------
      ! 1) Increment of energy
      !------------------------
      DO ji = kideb, kiut ! vertically integrated
          dq_i(ji,jl) = qt_i_fin(ji,jl) - qt_i_in(ji,jl)   
     &                + qt_s_fin(ji,jl) - qt_s_in(ji,jl)
      END DO

      DO ji = kideb, kiut ! layer by layer
         DO jk = 1, nlay_i
            dq_i_layer(ji,jk) =  q_i_layer_fin(ji,jk) - 
     &                           q_i_layer_in(ji,jk)
         END DO
      END DO

      !-------------------------------------------
      ! 2) Atmospheric heat flux, ice heat budget
      !-------------------------------------------
      DO ji = kideb, kiut
          fatm(ji,jl) = fsolgb(ji) + fratsb(ji)
     &                + fcsb(ji) + fleb(ji)
          sum_fluxq(ji,jl) = fatm(ji,jl) + fbbqb(ji)
     &                     - ftroce + fsnic + fprec
      END DO

      IF ( l_write ) THEN
         DO ji = kideb, kiut
            WRITE(numout,*) ' fsolgb : ' , fsolgb(ji)
            WRITE(numout,*) ' fratsb : ' , fratsb(ji)
            WRITE(numout,*) ' fcsb   : ' , fcsb(ji)
            WRITE(numout,*) ' fleb   : ' , fleb(ji)
            WRITE(numout,*)
            WRITE(numout,*) ' fatm   : ' , fatm(ji,jl)
            WRITE(numout,*) ' fbbqb  : ' , fbbqb(ji)
            WRITE(numout,*) ' ftroce : ' , ftroce
            WRITE(numout,*) ' fprec  : ' , fprec
            WRITE(numout,*) ' fsnic  : ' , fsnic
            WRITE(numout,*)
            WRITE(numout,*) ' sum_fluxq : ', sum_fluxq(ji,jl)
            WRITE(numout,*) ' dq_i      : ', dq_i(ji,jl) / ddtb
         END DO
      ENDIF

      !-----------------------
      ! 3) Conservation error
      !-----------------------
      IF ( l_write ) WRITE(numout,*) ' kideb, kiut, jl : ', 
     &   kideb, kiut, jl

      DO ji = kideb, kiut
         cons_error(ji,jl) = ABS( dq_i(ji,jl) / ddtb + 
     &                        sum_fluxq(ji,jl) )
      END DO

      numce = 0
      meance = 0.0
      DO ji = kideb, kiut
         IF ( cons_error(ji,jl) .GT. max_cons_err ) THEN
            numce = numce + 1
            meance = meance + cons_error(ji,jl)
         ENDIF
         IF (numce .GT. 0 ) meance = meance / numce

         IF ( cons_error(ji,jl) .GT. max_cons_err ) THEN
            WRITE(numout,*) ' Growth and melt - large cons. error '
            WRITE(numout,*) ' Maximum tolerated conservation error : ', 
     &                        max_cons_err
            WRITE(numout,*) ' DH Mean conservation error ', 
     &                        meance, numit
         ENDIF
      ENDDO

!------------------------------------------------------------------------------      
      END SUBROUTINE ice_th_con_dh
