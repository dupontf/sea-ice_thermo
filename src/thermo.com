!		
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  Thermo.com is incorparated by an instruction include in ice_th.f , 
!      fontbc.f and acrlbq.f. It comprises the commons associated to
!      thermodynamic ice computation
!
! Correspondance between the variables
! qlbqb   qlbq
! qcmbqb  qcmbq
! thcmb   thcm
! fstbqb  fstrbq
! fltbqb  ffltbq
! fscbqb  fscmbq
! fsolgb  fsolg
! ratbqb  ratbqg
! psbqb   psbq
! tabqb   tabq
! qabqb   qabq
! vabqb   vabq
! qfvbqb  qfvbq
! tsb     ts
! tfub    tfu
! hnpbqb  zhnpbq
! hnbqb   hnbq
! hgbqb   hgbq
! albqb   albq
! qstbqb  qstobq
! fbbqb   fbbq
! tbqb    tbq
! dmgbqb  dmgbq
! dmnbqb  dmnbq
! qlbbqb  qlbsbq
! cldqb   cloud
! dmgwib  dmgwi
! npb     number of points where computations has to be done
! npac    correspondance between the points
! fratsb  firg
! fcsb    fcsg
! fleb    fleg
! albgb   albg   06/08/2001
! sal     salic  08/03/2002

      COMMON/combq/qlbqb(nbpt),qcmbqb(nbpt),thcmb(nbpt),fstbqb(nbpt),
     &             fltbqb(nbpt),fscbqb(nbpt),fsolgb(nbpt),fsolgb2(nbpt),
     &             ratbqb(nbpt),psbqb(nbpt),tabqb(nbpt),
     &             qabqb(nbpt),vabqb(nbpt),qfvbqb(nbpt),tsb(nbpt),
     &             tfub(nbpt),hnpbqb(nbpt),hnbqb(nbpt),hgbqb(nbpt),
     &             albqb(nbpt),qstbqb(nbpt),fbbqb(nbpt),
     &             dmgbqb(nbpt),dmnbqb(nbpt),qlbbqb(nbpt),cldqb(nbpt),
     &             dmgwib(nbpt),npb(nbpt),npac(nbpt),tfs(nbpt),
     &             albgb(nbpt),focea(nbpt),fsup(nbpt), qsfcb(nbpt), 
     &             tdewb(nbpt)

! heat diffusion equation
      common/diff/ht_s_b(nbpt),ht_i_b(nbpt),t_su_b(nbpt),t_bo_b(nbpt),
     &            t_i_b(nbpt,maxnlay),t_s_b(nbpt,maxnlay),
     &            s_i_b(nbpt,maxnlay),h_i(maxnlay),h_s(maxnlay),
     &            sn_i_b(maxnlay),
     &            ab(nbpt),s_i_mean, q_i_b(nbpt,0:maxnlay+2), 
     &            q_sal_i(nbpt,0:maxnlay+2),
     &            q_s_b(nbpt,0:maxnlay+2),
     &            old_ht_i_b(nbpt), e_i_b(maxnlay)
! vertical grid
      common/vgrid/z_i(0:maxnlay),z_s(0:maxnlay),
     &             dh_s_subl(nbpt), dh_s_prec(nbpt), dh_s_font(nbpt),
     &             dh_s_snic(nbpt),
     &             dh_i_subl(nbpt), dh_i_font(nbpt), dh_i_bott(nbpt),
     &             dh_i_snic(nbpt),dh_snowice(nbpt),
     &             dh_s_melt(nbpt),dh_i_melt(nbpt),
     &             dh_i_surf(nbpt),dh_s_tot(nbpt)
!
      common/comdbq/fratsb(nbpt),fcsb(nbpt),
     &              fleb(nbpt),dvsbqb(nbpt),dvbbqb(nbpt),dvlbqb(nbpt),
     &              dvnbqb(nbpt), fc_su(nbpt), fc_bo_i(nbpt), 
     &              fatm(nbpt,jpl), fc_i(nbpt,0:maxnlay),
     &              fc_s(nbpt,0:maxnlay)

! for computation of temperatures after accr/abl 
      common/ comtem /enthal(0:maxnlay),
     &		      aaa,bbb,ccc,discrim,
     &                enthalsi
      common/ brintrsp /q_br(maxnlay), dq_i_brf(maxnlay)
      common/ bitzlip /ftrice,tmelts,ftroce
      common/ vgrid /himax,hsmax
      common/ conserv /e_i_o,e_s_o,de_i,de_s,isnow,imelt,
     &                 z_e_i(maxnlay),z_f_i(0:maxnlay),
     &                 e_t_o,ti_old(maxnlay),si_old(maxnlay),
     &                 s_i_mmean(maxnlay),h_i_mmean,
     &                 dq_i(nbpt,jpl), sum_fluxq(nbpt,jpl),
     &                 cons_error(nbpt,jpl), surf_error(nbpt,jpl),
     &                 dq_i_layer(nbpt, maxnlay),
     &                 qt_i_in(nbpt,jpl), qt_s_in(nbpt,jpl),
     &                 q_i_layer_in(nbpt,maxnlay),
     &                 qt_i_fin(nbpt,jpl), qt_s_fin(nbpt,jpl),
     &                 q_i_layer_fin(nbpt,maxnlay), fprec, fsnic

      common/ salt / beta_sal, s_i_new, s_i_snic, e_skel, q_summer,
     &               diff_br(maxnlay), rayleigh(maxnlay), fsb, fsbp

      ! units is a flux, including radab
      common/radiation/radab_phy_i(maxnlay),radab_s(maxnlay),
     &                 radab_phy_s(maxnlay),radtr_i(0:maxnlay),
     &                 radab_alg_i(maxnlay),radtr_s(0:maxnlay)

      ! par, energy available at the top of each layer
      common/biology/chla_i(maxnlay), par(maxnlay)

      ! Remapping
      INTEGER snswi,icsuswi,icboswi,snicswi
      INTEGER snind,icsuind,icboind,snicind,limsum


!-- End of file 'thermo.com'
!-----------------------------------------------------------------------------!
