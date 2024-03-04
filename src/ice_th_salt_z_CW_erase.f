      subroutine ice_th_salt_z(nlay_i,i_deb,i_fin)
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!! ** Purpose :
!!        This routine computes new salinities in the ice
!!
!! ** Method  : Vertical salinity profile computation 
!!              from segregation of salts, expulsion of brine
!!              gravity drainage, flushing
!!           
!! ** Steps
!!
!! ** Arguments
!!
!! ** Inputs / Outputs
!!
!! ** External
!!
!! ** References : Cox and Weeks, JGR, 1988
!!
!! ** History  : 
!!    (06-2003) Martin Vancoppenolle, Louvain-la-Neuve, Belgique
!!    (06-2004) Martin Vancoppenolle, Louvain-La-Neuve, Belgium
!!
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|

      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'thermo.com'

      dimension zs_i_old(maxnlay),z_brine_vol(maxnlay)
      dimension ds_i_exp(maxnlay)
      dimension ds_i_grd(maxnlay),ds_i_flu(maxnlay)
      dimension zds_i_exp(maxnlay)
      dimension zgrad_t(maxnlay),zt_i_bc(maxnlay)
! salt segregation variables
      dimension zold_cote(0:maxnlay),znew_cote(0:maxnlay)
      dimension zold_msal(0:maxnlay),znew_sal(0:maxnlay+2)
! brine expulsion variables
      dimension zs_b_old(maxnlay), zs_b_new(maxnlay)
      dimension zr_b_old(maxnlay), zr_b_new(maxnlay)
      dimension z_alpha(3,4)
! artificial flooding variables
      dimension ds_i_flo(maxnlay)
! flushing variables
      dimension ztrid(maxnlay,3), zindterm(maxnlay)
      dimension zdiagbis(maxnlay), zindtbis(maxnlay)
      dimension zsnew(maxnlay)

      dimension zalpha__i(maxnlay)
      dimension tube(nlay_i,2)   
      integer   ntubes, t, switch

! means
! local parameters
      zmu = 0.054
      zeps = 1.0e-6

      i_switch_gra = 2       ! switch for gravity drainage

c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c 0) Preliminary quantities
c------------------------------------------------------------------------------|

      do 10 ji = i_deb, i_fin

      do layer = 1, nlay_i
         zs_i_old(layer)    = s_i_b(ji,layer)                  !old salinities
         si_old(layer)      = zs_i_old(layer)
         zt_i_bc(layer)     = t_i_b(ji,layer) - 273.15         !celsius temps
         z_brine_vol(layer) = -zmu*s_i_b(ji,layer) / 
     &                        min(-zeps,zt_i_bc(layer))         !brine volume
                                                               !in non-dim units
      end do
      WRITE(*,*) ' zs_i_old    : ', ( zs_i_old(layer), 
     &            layer = 1, nlay_i )
      WRITE(*,*) ' zt_i_bc     : ', ( zt_i_bc(layer), 
     &            layer = 1, nlay_i )
      WRITE(*,*) ' s_i_b       : ', ( s_i_b    (ji,layer), 
     &            layer = 1, nlay_i )
      WRITE(*,*) ' t_i_b       : ', ( t_i_b    (ji,layer), 
     &            layer = 1, nlay_i )
      WRITE(*,*) ' z_brine_vol : ', ( z_brine_vol(layer), 
     &            layer = 1, nlay_i )

      zhiold = ht_i_b(ji) - dh_i_bott(ji) !old ice thickness in case of bottom freezing
      zformrate  =  max(dh_i_bott(ji)/ddtb,0.0) ! bottom formation rate of sea-ice

 10   continue


c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c 1) Salt segregation during ice growth 
c------------------------------------------------------------------------------|
c Cox and Weeks, JGR, 1988
      do 20 ji = i_deb, i_fin
      zsalmix   = 34.0 ! salinity of mixed layer to be changed in 3-D version
     !!!! to be changed

! computation of segregation coefficient zfracs  (fraction of salt retained
! in seawater during ice growth)
!     zfracs     =  0.12
!     if (zformrate.gt.2.0e-08) then
!        zfracs  =  0.8925 + 0.0568*log(100.0*zformrate)
!     endif
!     if (zformrate.gt.3.6e-07) then
!        zfracs  =  0.26 / (0.26 + 0.74*exp(-724300.0*zformrate))
!     endif

      if (zformrate.gt.0.0) then
!        zs_new_lay =  zfracs*zsalmix !increase of bulk salinity
! new method (iterative for ice growth)
         zs_new_lay = s_i_new
         sal_new_layer = zs_new_lay
         write(3157,*) zfracs, zformrate, dh_i_bott(ji)*zfracs*zsalmix 
      endif 

 20   continue
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c 2) Brine expulsion if ice coldens 
c------------------------------------------------------------------------------|
c Cox and Weeks, 1988      
! coefficients for brine salinity in function of temperature
      do 30 ji = i_deb, i_fin

      z_alpha(1,1) = -3.9921
      z_alpha(1,2) = -22.7
      z_alpha(1,3) = -1.0015
      z_alpha(1,4) = -0.019956

      z_alpha(2,1) = 206.24
      z_alpha(2,2) = -1.8907
      z_alpha(2,3) = -0.060868
      z_alpha(2,4) = -0.0010247

      z_alpha(3,1) = -4442.1
      z_alpha(3,2) = -277.86
      z_alpha(3,3) = -5.501
      z_alpha(3,4) = -0.03669

      do layer = 1, nlay_i

!switches for old brine salinity      
         if (ti_old(layer).gt.250.25) then
            izswi_ex = 1
         else if (ti_old(layer).gt.229.15) then
            izswi_ex = 2
         else
            izswi_ex = 3
         endif
         
!old brine salinity and density
         ztioldc = ti_old(layer) - 273.15 !celsius temperatures
         zs_b_old(layer) = z_alpha(izswi_ex,1)
     &   + z_alpha(izswi_ex,2)*ztioldc
     &   + z_alpha(izswi_ex,3)*ztioldc*ztioldc
     &   + z_alpha(izswi_ex,4)*ztioldc*ztioldc*ztioldc
!MV 2005 change brine salinity formula
         zs_b_old(layer) = -ztioldc/0.054
         zs_b_old(layer) = max(zs_b_old(layer),zeps)
         zr_b_old(layer) = 1000.0 + 0.8*zs_b_old(layer)

!switches for new brine salinity      
         if (t_i_b(ji,layer).gt.250.25) then
            izswi_ex = 1
         else if (t_i_b(ji,layer).gt.229.15) then
            izswi_ex = 2
         else
            izswi_ex = 3
         endif

!new brine salinity and density
         zs_b_new(layer) = z_alpha(izswi_ex,1)
     &   + z_alpha(izswi_ex,2)*zt_i_bc(layer)
     &   + z_alpha(izswi_ex,3)*zt_i_bc(layer)*zt_i_bc(layer)
     &   + z_alpha(izswi_ex,4)*zt_i_bc(layer)*zt_i_bc(layer)
     &     *zt_i_bc(layer)

!MV 2005 change brine salinity formula
         zs_b_new(layer) = -zt_i_bc(layer)/0.054

         zs_b_new(layer) = max(zs_b_new(layer),zeps)
         zr_b_new(layer) = 1000 + 0.8*zs_b_new(layer)

         zdummy = ( ( (zs_b_new(layer)/max(zs_b_old(layer),zeps)) )** 
     &   (-0.09051254089)  )
         zdummy = zdummy*zr_b_new(layer)/zr_b_old(layer)
         zdummy = zdummy*exp(0.8/917.0*
     &            (zs_b_old(layer)-zs_b_new(layer)))
         zds_i_exp(layer) = zs_i_old(layer)*(zdummy-1.0)

         zds_i_exp(layer) = min(zds_i_exp(layer),0.0)
      end do

!drainage in underlying layer
! 1) only drainage downwards 
      do layer = 2, nlay_i
!        ds_i_exp(layer) = zds_i_exp(layer) - zds_i_exp(layer-1) ! salt drained
!                                !is reported in underlying layer
         if (t_su_b(ji).gt.273.14) then
            ds_i_exp(layer) = 0.0
         endif
      end do

! 2) drainage up and down only half of zds_i_exp(nlay_i) is put in the mixed layer
      ds_i_exp(1)      = ( zds_i_exp(1)-zds_i_exp(2) )/2.0
      ds_i_exp(nlay_i) = zds_i_exp(nlay_i) - zds_i_exp(nlay_i-1)/2.0
      do layer = 2, nlay_i - 1 
         ds_i_exp(layer) = zds_i_exp(layer)
     &                   - ( zds_i_exp(layer-1) + zds_i_exp(layer+1) ) 
     &                   / 2.0
      end do
!! normal expulsion
!     do layer = 1, nlay_i
!        ds_i_exp(layer) =  zds_i_exp(layer)
!     end do

! no expulsion in summer
      do layer = 1, nlay_i
         if (t_su_b(ji).gt.273.14) then
            ds_i_exp(layer) = 0.0
         endif
      end do

 30   continue

c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c 3) Gravity drainage 
c------------------------------------------------------------------------------|
c Cox and Weeks, 1988

      do 40 ji = i_deb, i_fin

! temperature gradient
      z_h_lay    = ht_i_b(ji)/real(nlay_i) 
      zgrad_t(1) = 2.0*( t_i_b(ji,1) - (ht_s_b(ji)*t_i_b(ji,1) +
     &             z_h_lay*t_s_b(ji,1))/(z_h_lay + ht_s_b(ji)))/z_h_lay

      do layer = 2, nlay_i - 1
         zgrad_t(layer) = ( t_i_b(ji,layer+1) - t_i_b(ji,layer-1) ) /
     &                    z_h_lay
      end do
      zgrad_t(nlay_i) = - 2.0*(t_i_b(ji,nlay_i)-t_bo_b(ji))/z_h_lay

      izswiper = 1 ! this 'permeability' switch equals one as long
                   ! as the brine volume of each layer is gt 0.05
                   ! if not, all the upper layers will have their
                   ! gravity drainage equal to 0
      do layer = nlay_i, 1, -1
         write(*,*) zgrad_t(layer)/100,layer, numit
         ! this switch equals 1 if gravity drainage and 0 if not
         ! gravity drainage occurs if brine volume of the layer is gt 0.05
         izswigrd = max(0.0,sign(1.0,z_brine_vol(layer)-0.05))
         ! permeability switch
         if (z_brine_vol(layer).le.0.05) then
            izswiper = 0
         endif
         ! temperature gradient must be upwards ....!!!
         IF (zgrad_t(layer) .LE. 0.0 ) THEN
            izswigrd = 0.0
         ENDIF
         

!        !gravity drainage parameters
         zdelta          = 1.68d-7
c        zdelta = zdelta - 3.0*zdelta/10.0
         zdelta = zdelta - 6.5*zdelta/10.0
!        zdelta = zdelta*sqrt(max((1.0-ht_i_b(ji)/1.20),0.0))
         zeta            = 20.0d0
         ds_i_grd(layer) = zdelta*(1.0-zeta*z_brine_vol(layer))*
     &                     zgrad_t(layer)*ddtb*real(izswigrd)
     &                     *real(izswiper)
         ds_i_grd(layer) = min(0.0, ds_i_grd(layer))

! Another idea
!        z_s_temp_gd = zsalmix*z_brine_vol(layer)
!        ds_i_grd(layer) = real(izswiper)*(z_s_temp_gd-s_i_b(ji,layer))

!        write(*,*) numit, layer, (s_i_b(ji,layer)-ds_i_grd(layer))
!    &                    /z_brine_vol(layer)
!        if ((s_i_b(ji,layer)-ds_i_grd(layer))
!    &       /z_brine_vol(layer).le.(34.0)) then
!           ds_i_grd(layer) = 0.0
!        endif

         if (t_su_b(ji).gt.273.14) then
            ds_i_grd(layer) = 0.0
         endif
         
      end do

! drainage in underlying layer
!     do layer = 2, nlay_i
!        ds_i_grd(layer) = ds_i_grd(layer) - ds_i_grd(layer-1) ! salt drained
!                                !is reported in underlying layer
!     end do

 40   continue

c------------------------------------------------------------------------------|
c 4) Flushing 
c------------------------------------------------------------------------------|
      do 50 ji = i_deb, i_fin

      WRITE(*,*) ' ice_th_salt_z : '
      WRITE(*,*) ' ~~~~~~~~~~~~~   '
      WRITE(*,*) ' t_su_b : ', t_su_b(1)
      WRITE(*,*) ' tpw    : ', tpw
      izswiflu = max(0.0,sign(1.0,t_su_b(ji)-tpw)) !0 si hiver 1 si ete
      WRITE(*,*) ' izsiwflu : ', izswiflu 
                                                   !occurs if surface is melting
c minimum brine volume
      zbvmin = 1.0 !minimum brine volume       
      zpor_tres = 0.050 !Porosity treshold
      zbeta     = 0.30 !tunes the amount of water from the melt put on the reservoir
      zpor_tres = flu_bvtr
      zbeta = flu_beta
      WRITE(*,*) ' zpor_tres, zbeta : ', zpor_tres,zbeta
      WRITE(*,*) ' z_brine_vol : ', ( z_brine_vol(layer), 
     &            layer = 1, nlay_i )

      izswiflu = max(0.0,sign(1.0,t_su_b(ji)-tpw)) !0 si hiver 1 si ete
      do layer = 1, nlay_i
         zbvmin = min(z_brine_vol(layer),zbvmin) !minimum brine volume
         if (z_brine_vol(layer).lt.zpor_tres) then    !if brine volume of one layer
            izswiflu = 0                         !is smaller than 0.05 then
         endif                                   !no flushing occurs
      end do

c     write(*,*) '****************************************'
c     write(*,*) 'Brine volume minimum', zbvmin
c     write(*,*) 'Flushing ', numit, 'yes ?',izswiflu

c mean brine density and viscosity
      zrhobmean = 0.0 !mean brine density
      zvisbmean = 1.7e-3 !mean brine viscosity
      do layer = 1, nlay_i
         zrhobmean = zrhobmean+zr_b_new(layer)/real(nlay_i)
         zvisbmean = zvisbmean
      end do
c flux of brine through all layers is determined by minimum brine volume
c units = kg.m-2.sec-1
c Eicken 2004 permeability
      zcoeff11=4.708e-14
      zcoeff12=0.07690*1000.0
      zcoeff21=3.738e-11
      zcoeff22=0.007265*1000.0
      if (zbvmin.le.0.096) then
         zalpha = zcoeff11*exp(zcoeff12*zbvmin)
      else
         zalpha = zcoeff21*exp(zcoeff22*zbvmin)
      endif
c     zalpha    = 1.0e-12
c     zalpha    = 1.0e-8             !permeability of sea-ice (m2)

      zbrinespd = zalpha*zrhobmean*9.81/zvisbmean !Darcy Speed
      zbrinespdr = zbrinespd/zbvmin !Real speed
c     write(*,*) 'Porosity',zalpha,'Darcy Speed',zbrinespd
c     write(*,*) 'Real speed', zbrinespdr,'Brinevol',zbvmin
c     write(*,*) 'brine speed', zbrinespd, 'viscosity',zvisbmean,
c    &           'density', zrhobmean
      zbrineflx = zrhobmean*zbrinespd
c     write(*,*) 'brineflx',zbrineflx
      !(speed of brine flow [m.s-1] times density [kg.m-3] gives a mass flux
      ! [kg.m-2.s-1], which multiplied by ddtb gives mass of brine which can
      ! be flowed out the ice by gravity per square meter)
      zbrineflx = zbrineflx*real(izswiflu) !what can sustain the ice by porosity

c meltwater reservoir content and salinity
c meltwater inflow
      zdh_i     = min(0.0,dh_i_surf(ji))
      zdh_s     = min(0.0,dh_s_tot(ji))
      z_melt_inflow= zbeta*( -rhog*min(dh_i_surf(ji),0.0)   !equals Q*/time step/
     &                       -rhon*min(dh_s_tot(ji) ,0.0) ) !units = [kg.m-2]
      if ((numit.ge.517).and.(numit.le.537)) then
         ztotalmelt = z_melt_inflow*real(izswiflu) + ztotalmelt
         write(84,*) numit,'Deltah_snow', dh_s_tot(ji)
      else
         ztotalmelt = 0.0
      endif 

!     write(*,*) 'ZTOTAL', numit, ztotalmelt
      write(8774,*) '*****************'
      write(8774,*) '****',numit,'****'
      write(8774,*) 'melt inflow',z_melt_inflow
      write(8774,*) 'swi flushin',izswiflu
      write(8774,*) ' Darcy flux ', zbrineflx,' avail flux ', 
     &                              z_melt_inflow/ddtb
      write(8774,*) ' Q*deltat/rho ', z_melt_inflow/zr_b_new(1)
      write(8774,*) ' e*deltah ', z_brine_vol(1)*ht_i_b(ji)/real(nlay_i)

c meltwater reservoir content
c     z_wat_res0 = z_wat_res
c     z_wat_res  = z_wat_res + z_melt_inflow
c     z_brine_outflow = zbrineflx*ddtb

c     if (z_brine_outflow.le.z_wat_res) then
c        z_wat_res = z_wat_res - z_brine_outflow
c     else
c        z_brine_outflow = z_wat_res
c        z_wat_res   = z_wat_res - z_brine_outflow
c     endif
c     z_wat_res  = max(z_wat_res,0.0)

c meltwater reservoir salinity       
c     z_sal_res0= z_sal_res
c     if (z_wat_res.gt.0.0) then
c     z_sal_res = zbeta*rhog*max(-dh_i_surf(ji),0.0)
c    &            *(s_i_b(ji,1)-z_sal_res0)/
c    &            z_wat_res+z_sal_res0
c     else
c     z_sal_res = 0.0
c     endif

c another fashion to calculate brine : amount of flushed water limited by hydrostatic equilibrium
      z_brine_outflow = z_melt_inflow
c everything is flowed out
c and amount of flushing is limited by hydrostatic head

c salinity variations due to flushing in a continuum tube where brine flow is 
c dictated by minimum permeability, and where source water comes from the surface
      z_sal_res = 0.0

c Euler explicit is conservative but unstable 
c 1st and 2nd order Euler implicit are stable but not conservative
c 1) terms of the tridiagonal system
      zmasssalt0 = 0.0
      zmasssalt1 = z_brine_outflow*zs_b_new(nlay_i)/1000.0
c    &             /z_brine_vol(nlay_i)
c    &             /zr_b_new(nlay_i)/ht_i_b(ji)*real(nlay_i)
      if (izswiflu.eq.1) WRITE(8774,*) '*** HERE ***'
      do layer = 1, nlay_i
         zmasssalt0 = zmasssalt0 + 
     &          917.0*ht_i_b(ji)/real(nlay_i)*s_i_b(ji,layer)/1000.0
         write(8774,*) 'ZSBNEW BEFOR', zs_b_new(layer), layer 

         zalphai = -z_brine_outflow/zr_b_new(layer)
     &             /z_brine_vol(layer)
     &             /ht_i_b(ji)*real(nlay_i)

         write(8774,*) 'ZALPHAI',zalphai
         ztrid(layer,2)  = 1 - zalphai
         ztrid(layer,1)  = zalphai
         ztrid(layer,3)  = 0.0
         zindterm(layer) = zs_b_new(layer)
      end do
      ztrid(1,1) = 0.0

c 2) first transformation
      zindtbis(1) =  zindterm(1)
      zdiagbis(1) =  ztrid(1,2)

      do numeq = 2, nlay_i
         zdiagbis(numeq)  =  ztrid(numeq,2) - ztrid(numeq,1)*
     &                       ztrid(numeq-1,3)/zdiagbis(numeq-1)
         zindtbis(numeq)  =  zindterm(numeq) - ztrid(numeq,1)*
     &                       zindtbis(numeq-1)/zdiagbis(numeq-1)
      end do

c 3) back to the brine salinities
      zs_b_new(nlay_i)    =  zindtbis(nlay_i)/zdiagbis(nlay_i)
      do numeq = nlay_i - 1, 1, -1
         zs_b_new(numeq)  =  (zindtbis(numeq) - ztrid(numeq,3)*
     &                        t_i_b(ji,numeq+1))/zdiagbis(numeq)
      end do
      
c 4) back to the ice salinities
      do layer = 1, nlay_i
         write(8774,*) 'ZSBNEW AFTER', zs_b_new(layer), layer 
         zsnew(layer) = zs_b_new(layer)*z_brine_vol(layer)
         write(8774,*) 'ZSNEW', zsnew(layer), layer 
         ds_i_flu(layer) = real(izswiflu)*(zsnew(layer)-s_i_b(ji,layer))
         zmasssalt1 = zmasssalt1 + 
     &           917.0*ht_i_b(ji)/real(nlay_i)*(zsnew(layer))/1000.0
      end do
      zdiffmass = zmasssalt1 - zmasssalt0
      write(8774,*) ' Mass of salt conservation ? '
      write(8774,*) ' Before ', zmasssalt0
      write(8774,*) ' After  ', zmasssalt1
      write(8774,*) ' Differ ', zdiffmass

c 5) Correction to ice salinities in order to get mass conservation
      zcorrsal = - zdiffmass/real(nlay_i)*real(izswiflu)
      write(8774,*) 'Correction in kg of salt', zcorrsal
      zcorrsal = - zcorrsal / 917. / ht_i_b(ji) * real(nlay_i) * 1000.0
      write(8774,*) 'Correction in ppt       ', zcorrsal
      do layer = 1, nlay_i
         ds_i_flu(layer) = ds_i_flu(layer) + zcorrsal 
      end do
      if (ds_i_flu(layer) .gt. 0.0) then
      write(*,*) 'ALERTE, flushing positif',numit,ds_i_flu(layer),layer
      endif
c     ds_i_flu(layer) = max(ds_i_flu(layer),0.0)
c     layer = 1
c     ds_i_flu(layer) = real(izswiflu)*z_brine_outflow*
c    &                  (z_sal_res-zs_b_new(layer))/
c    &                  (zr_b_new(layer)*ht_i_b(ji)/real(nlay_i)) 
c     write(8774,*) layer,' S_B ',zs_b_new(layer),' RH_B ',zr_b_new(layer)
c     write(8774,*) ' ds_flu ',ds_i_flu(layer)
c     write(8774,*) ' Q*Delta_t/rho ', z_brine_outflow/zr_b_new(layer)
c     write(8774,*) ' E*Delta_h ', z_brine_vol(layer)*ht_i_b(ji)/
c    &                           real(nlay_i)
c     do layer = 2, nlay_i
c        ds_i_flu(layer) = real(izswiflu)*z_brine_outflow*
c    &                     (zs_b_new(layer-1)-zs_b_new(layer))/
c    &                     (zr_b_new(layer)*ht_i_b(ji)/real(nlay_i)) 
c     write(8774,*) layer,'S_B ',zs_b_new(layer),'RH_B ',zr_b_new(layer)
c     write(8774,*) 'ds_flu ',ds_i_flu(layer)
c        write(84,*) layer,'FLU',ds_i_flu(layer)
c salt diffusivity for the last layer
c     end do
c     zdiffus = 1.0e-3
c     zsal_bot = 10.0 !salinity of the bottom ice
c     ds_i_flu(nlay_i) = real(izswiflu)*zdiffus*
c    & (zsal_bot-2*zs_b_new(nlay_i)+zs_b_new(nlay_i-1))
c    & /ht_i_b(ji)/ht_i_b(ji)*real(nlay_i)*real(nlay_i)
c     write(8774,*) nlay_i,'S_B ',zs_b_new(nlay_i),'RH_B ',
c    &              zr_b_new(nlay_i)
c     write(8774,*) 'ds_flu ',ds_i_flu(nlay_i)
c     write(*,*) real(izswiflu)*zdiffus*
c    & (zsal_bot-2.0*zs_b_new(nlay_i)+zs_b_new(nlay_i-1))
c    & /ht_i_b(ji)/ht_i_b(ji)*real(nlay_i)*real(nlay_i)

c     write(*,*) zs_b_new(nlay_i-1),zs_b_new(nlay_i)

c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c 4b) Artificial term due to flooding 
c------------------------------------------------------------------------------|
c     if ((dh_s_melt(ji).eq.0.0).and.(ht_s_b(ji).ne.0.0)) then
c        ds_i_flo(1) = real(nlay_i)*6.679e-3
c     else
         ds_i_flo(1) = 0.0
c     endif
      do layer = 2, nlay_i
         ds_i_flo(layer) = 0.0
      end do

 50   continue
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c 5) Updates ice salinity and computes the means
c------------------------------------------------------------------------------|
      do 60 ji = i_deb, i_fin

! vertical means
      zs_i_mean   = 0.0
      ds_seg_mean = 0.0
      ds_exp_mean = 0.0
      ds_grd_mean = 0.0
      ds_flu_mean = 0.0

      do layer = 1, nlay_i
         zs_i_mean   = zs_i_mean   + s_i_b(ji,layer)/real(nlay_i) ! vertical mean salinity
      end do
      ds_seg_mean = max(dh_i_bott(ji),0.0)/max(dh_i_bott(ji)
     &               ,0.0000000000001)*
     &               (sal_new_layer*max(dh_i_bott(ji),0.0)+
     &               zhiold*zs_i_mean)/ht_i_b(ji)-
     &               zs_i_mean

      zs_i_mean   = 0.0

      do layer = 1, nlay_i
         s_i_b(ji,layer) = s_i_b(ji,layer)
     &                   + ds_i_exp(layer) + ds_i_grd(layer)
     &                   + ds_i_flu(layer) 
! Flooding added end of november 2004
c    &                   + ds_i_flo(layer)
         s_i_b(ji,layer) = max(s_i_b(ji,layer),0.1)
         zs_i_mean   = zs_i_mean   + s_i_b(ji,layer)/real(nlay_i) ! vertical mean salinity
         ds_exp_mean = ds_exp_mean + ds_i_exp(layer)/real(nlay_i)
         ds_grd_mean = ds_grd_mean + ds_i_grd(layer)/real(nlay_i)
         ds_flu_mean = ds_flu_mean + ds_i_flu(layer)/real(nlay_i)
      end do
! monthly profiles
      do layer = 1, nlay_i
         s_i_mmean(layer) = s_i_mmean(layer) + s_i_b(ji,layer)/15.0
      end do
      h_i_mmean = h_i_mmean + ht_i_b(ji)/15.0

 60   continue
       
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c 6) Output 
c------------------------------------------------------------------------------|
      do 70 ji = i_deb, i_fin
 751  format(10f12.3)
 752  format(10f8.3)
 753  format(11f8.3) 
 754  format(10f8.3)
 755  format(4f8.3)

      write(7501,751) zs_i_mean, ds_seg_mean, ds_exp_mean, ds_grd_mean,
     &                ds_flu_mean

      write(7502,752) (s_i_b(ji,layer),layer=1,nlay_i) 
      write(7505,754) (ds_i_exp(layer),layer=1,nlay_i)
      write(7506,754) (ds_i_grd(layer),layer=1,nlay_i)
      write(7507,754) (ds_i_flu(layer),layer=1,nlay_i)
c reservoir terms
      write(7508,755) z_wat_res, z_sal_res, z_melt_inflow, 
     &                z_brine_outflow

      if (mod(numit-908,15).eq.0) then
         write(7503,753) (s_i_mmean(layer),layer=1,nlay_i),h_i_mmean
         do layer = 1, nlay_i
            s_i_mmean(layer) = 0.0
         end do
         h_i_mmean = 0.0
      endif

 70   continue     


c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c fin de la subroutine
      end subroutine 
