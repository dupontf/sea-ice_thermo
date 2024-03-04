c
      dimension vicmom(imax,jmax,35)
      equivalence ( vicmom(1,1,1) , sxg(1,1) ) 
      equivalence ( vicmom(1,1,2) , syg(1,1) ) 
      equivalence ( vicmom(1,1,3) , sxxg(1,1) ) 
      equivalence ( vicmom(1,1,4) , syyg(1,1) ) 
      equivalence ( vicmom(1,1,5) , sxyg(1,1) ) 
      equivalence ( vicmom(1,1,6) , sxn(1,1) ) 
      equivalence ( vicmom(1,1,7) , syn(1,1) ) 
      equivalence ( vicmom(1,1,8) , sxxn(1,1) ) 
      equivalence ( vicmom(1,1,9) , syyn(1,1) ) 
      equivalence ( vicmom(1,1,10), sxyn(1,1) ) 
      equivalence ( vicmom(1,1,11), sxa(1,1) ) 
      equivalence ( vicmom(1,1,12), sya(1,1) ) 
      equivalence ( vicmom(1,1,13), sxxa(1,1) ) 
      equivalence ( vicmom(1,1,14), syya(1,1) ) 
      equivalence ( vicmom(1,1,15), sxya(1,1) ) 
      equivalence ( vicmom(1,1,16), sxc0(1,1) ) 
      equivalence ( vicmom(1,1,17), syc0(1,1) ) 
      equivalence ( vicmom(1,1,18), sxxc0(1,1) ) 
      equivalence ( vicmom(1,1,19), syyc0(1,1) ) 
      equivalence ( vicmom(1,1,20), sxyc0(1,1) ) 
      equivalence ( vicmom(1,1,21), sxc1(1,1) ) 
      equivalence ( vicmom(1,1,22), syc1(1,1) ) 
      equivalence ( vicmom(1,1,23), sxxc1(1,1) ) 
      equivalence ( vicmom(1,1,24), syyc1(1,1) ) 
      equivalence ( vicmom(1,1,25), sxyc1(1,1) ) 
      equivalence ( vicmom(1,1,26), sxc2(1,1) ) 
      equivalence ( vicmom(1,1,27), syc2(1,1) ) 
      equivalence ( vicmom(1,1,28), sxxc2(1,1) ) 
      equivalence ( vicmom(1,1,29), syyc2(1,1) ) 
      equivalence ( vicmom(1,1,30), sxyc2(1,1) ) 
      equivalence ( vicmom(1,1,31), sxst(1,1) ) 
      equivalence ( vicmom(1,1,32), syst(1,1) ) 
      equivalence ( vicmom(1,1,33), sxxst(1,1) ) 
      equivalence ( vicmom(1,1,34), syyst(1,1) ) 
      equivalence ( vicmom(1,1,35), sxyst(1,1) ) 
c
      common/comma2/sxg(imax,jmax),syg(imax,jmax),
     &              sxxg(imax,jmax),syyg(imax,jmax),
     &              sxyg(imax,jmax),
     &              sxn(imax,jmax),syn(imax,jmax),
     &              sxxn(imax,jmax),syyn(imax,jmax),
     &              sxyn(imax,jmax),
     &              sxa(imax,jmax),sya(imax,jmax),
     &              sxxa(imax,jmax),syya(imax,jmax),
     &              sxya(imax,jmax),
     &              sxc0(imax,jmax),syc0(imax,jmax),
     &              sxxc0(imax,jmax),syyc0(imax,jmax),
     &              sxyc0(imax,jmax),
     &              sxc1(imax,jmax),syc1(imax,jmax),
     &              sxxc1(imax,jmax),syyc1(imax,jmax),
     &              sxyc1(imax,jmax),
     &              sxc2(imax,jmax),syc2(imax,jmax),
     &              sxxc2(imax,jmax),syyc2(imax,jmax),
     &              sxyc2(imax,jmax),
     &              sxst(imax,jmax),syst(imax,jmax),
     &              sxxst(imax,jmax),syyst(imax,jmax),
     &              sxyst(imax,jmax)
c
c
