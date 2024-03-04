!
!		COMMONS FOR SEA ICE BIOGEOCHEMISTRY
!		===================================
!
!
! z_i_bio       : vertical coordinate of the middle of bio-layers
! deltaz_i_bio  : thickness of the bio-layers
! c_i_bio       : concentration of tracers in the brines
! cbu_i_bio     : bulk concentration of tracers in the ice
! cbun_i_bio    : new bulk concentration of tracers in the ice (dummy rray)
! t_i_bio       : temperature in the ice, interpolated on the bio grid
! s_i_bio       : salinity in the ice, interpolated on the bio grid
! e_i_bio       : relative brine volume, interpolated on the bio grid
! c_w_bio       : concentration of tracer in seawater

      COMMON /bioarrays/ z_i_bio(nlay_bio), deltaz_i_bio(nlay_bio),
     &   c_i_bio(ntra_bio,nlay_bio), cbu_i_bio(ntra_bio,nlay_bio),
     &   cbun_i_bio(ntra_bio,nlay_bio),
     &   ct_i_bio(ntra_bio),
     &   t_i_bio(nlay_bio), s_i_bio(nlay_bio), 
     &   e_i_bio(nlay_bio), c_w_bio(ntra_bio),
     &   ch_bo_bio(ntra_bio), ch_si_bio(ntra_bio),
     &   c_s_bio(ntra_bio), c_skel_bio(ntra_bio),
     &   f_su_tra(ntra_bio), f_bo_tra(ntra_bio),
     &   f_bogr(ntra_bio), f_sigr(ntra_bio),
     &   fdiff(ntra_bio,0:nlay_bio),
     &   diff_br_bio(nlay_bio), fcb(ntra_bio), fcbp(ntra_bio), 
     &   fcsi(ntra_bio), fcbm(ntra_bio), fcsu(ntra_bio), 
     &   fcb_max(ntra_bio)

      COMMON /biophyparams/
     &   frtr_si_bio      ,       !: fraction of sea water conc trapd in snow ice (bio)
     &   stickfac         ,       !: sticking factor for biological tracers
     &   par_fsw          ,       !: PAR / FSW conversion factor
     &   astar_alg        ,       !: specific absorption coeff (m-1 / (mg chla m-3))
     &   fdet_alg         ,       !: fraction of detrital absorption compared to algal absorption
     &   pp_presc         ,       !: prescribed primary production
     &   diff_da          ,       !: diffusivity for diatoms
     &   diff_nut         ,       !: diffusivity for nutrients
     &   h_bio                    !: thickness of the biological layer ('SL' & 'BAL' only)

      REAL(8), DIMENSION(ntra_bio) ::
     &   mt_i_bio_init, mt_i_bio_final,
     &   mt_i_bio

      REAL(8), DIMENSION(ntra_bio,nlay_bio) ::
     &   m_i_bio_init, m_i_bio_final


      COMMON/bioconserv/
     &   mt_i_bio_init, mt_i_bio_final,
     &   mt_i_bio,
     &   m_i_bio_init, m_i_bio_final

      LOGICAL ::
     &   flag_diff(ntra_bio)     , !: flag which describes diffusability of a tracer
     &   flag_active(ntra_bio)   , !: describe whether a tracer is assumed active or not
     &   flag_adsorb(ntra_bio)   , !: adsorbed or not ?
     &   ln_trbo                 , !: activate tracer basal entrapment
     &   ln_trsi                 , !: activate tracer surface entrapment
     &   ln_trdiff               , !: activate tracer diffusion
     &   ln_trremp               , !: activate tracer remapping
     &   ln_lim_dsi              , !: activate DSi limitation
     &   ln_lim_no3              , !: activate DSi limitation
     &   ln_lim_lig              , !: activate light limitation
     &   ln_lim_tem              , !: activate temperature inhibition
     &   ln_lim_sal              , !: activate brine salinity inhibition
     &   ln_lys                  , !: activate lysis
     &   ln_rem                  , !: activate remineralization
     &   ln_syn                  , !: activate diatom synthesis
     &   ln_dormant                !: activate diatom dormancy

      INTEGER(4) ::
     &   lim_sal_swi               !: switch for the type of brine salinity inhibition
     &   i_flux                    !: type of flux in the model (SL and BA only)

      COMMON /bioswi/
     &   flag_diff               ,
     &   flag_active             ,
     &   flag_adsorb             ,
     &   ln_trbo                 ,
     &   ln_trsi                 ,
     &   ln_trdiff               ,
     &   ln_trremp               ,
     &   ln_lim_dsi              ,
     &   ln_lim_no3              ,
     &   ln_lim_lig              ,
     &   ln_lim_tem              ,
     &   ln_lim_sal              ,
     &   ln_lys                  ,
     &   ln_rem                  ,
     &   ln_syn                  ,
     &   ln_dormant              ,
     &   lim_sal_swi             ,
     &   i_flux

      CHARACTER(len=3) ::
     &   biotr_i_nam(ntra_bio)  , !:
     &   biotr_i_typ(ntra_bio)    !:

      CHARACTER(len=8) ::         
     &   biotr_i_uni(ntra_bio)

      CHARACTER(len=2) ::         
     &   c_mod
         
      COMMON /biochar/
     &   biotr_i_nam, 
     &   biotr_i_typ,
     &   biotr_i_uni,
     &   c_mod

      COMMON /biomass/
     &   chla_i_bio(nlay_bio)     !: chlorophyll a concentration

      REAL(8) ::
     &   syn_bio(nlay_bio)       , !: diatom synthesis
     &   lys_bio(nlay_bio)       , !: diatom lysis
     &   rem_bio(nlay_bio)       , !: diatom remineralization
     &   lsr_bio(nlay_bio)       , !: lysis
     &   lim_lig(nlay_bio)       , !: light limitation
     &   lim_dsi(nlay_bio)       , !: DSi limitation
     &   lim_no3(nlay_bio)       , !: DSi limitation
     &   lim_tem(nlay_bio)       , !: temperature limitation
     &   lim_sal(nlay_bio)         !: salt limitation

      COMMON /biosources/
     &   syn_bio                 , !: diatom synthesis
     &   lys_bio                 , !: diatom lysis
     &   rem_bio                 , !: remineralization
     &   phi_bio(nlay_bio)       , !: photosynthesis
     &   resp_bio(nlay_bio)      , !: respiration
     &   excr_bio(nlay_bio)      , !: excretion
     &   lsr_bio                 , !: lysis
     &   lim_lig                 , !: light limitation
     &   lim_dsi                 , !: dsi limitation
     &   lim_no3                 , !: dsi limitation
     &   lim_tem                 , !: temperature limitation
     &   lim_sal                   !: salt limitation

      COMMON /biorad/             
     &   radab_alg_i_bio(nlay_bio) , !: Absorbed radiation (flux, w.m-2)
     &   pur_bio(nlay_bio)         , !: Photosynthetically usable radiation
     &   par_bio(nlay_bio)           !: Photosynthetically usable radiation
                                     !: (flux micromol quanta m-2 s-1)
      !---------------------------------------
      ! Constants and parameters of the model
      !---------------------------------------
      REAL(8) ::
     &   mumax_bio                      , !: maximum specific growth (s-1)
     &   kmax_bio                       , !: maximum specific photosynthesis (s-1)
     &   maint_bio                      , !: maintenance rate (s-1)
     &   klys_bio                       , !: autolysis rate (s-1)
     &   ek_bio                         , !: light adaptation micrmol quanta m-2 s-1
     &   khs_sr_bio                     , !: half sat for SR sytnhesis (-)
     &   ecg_bio                        , !: energetic cost growth (-)
     &   ke_bio                         , !: EPS excre ct (s-1)
     &   khs_fe_bio                     , !: half sat for Fe uptake (micrmol m-3)
     &   khs_si_bio                     , !: half sat for si uptake (mmol m-3)
     &   khs_n_bio                      , !: half sat for N uptake (mmol m-3)
     &   kmin_bio                       , !: fraction of Si loss that is remineralized (-)
     &   kdis_bio                       , !: frustule dissolution rate (s-1)
     &   kads_bio                       , !: adsorption constant for iron (s-1)
     &   rg_bio                         , !: temperature coefficient for diatom synthesis (deg C-1)
     &   lim_sal_wid                    , !: width of the limitation function (case 4)
     &   lim_sal_smax                   , !: salinity at which limitation function is 1 (case 4)
     &   chla_c                         , !: Chlorophyll-to-carbon ratio
     &   si_c                           , !: Si cell quota
     &   no3_c                            !: N cell quota

      ! Parameters
      REAL(8) :: 
     &   fe_c        = 0.04             , !: Fe cell quota
     &   fe_c_eps    = 0.08               !: Fe adsorption on EPS

      COMMON /bioparams/
     &   mumax_bio, kmax_bio, maint_bio, klys_bio, ek_bio, khs_sr_bio, 
     &   ecg_bio, ke_bio, khs_fe_bio, khs_si_bio, khs_n_bio, kmin_bio, 
     &   kdis_bio, kads_bio, rg_bio, chla_c, si_c, no3_c, lim_sal_wid, 
     &   lim_sal_smax

