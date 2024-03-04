c
c		COMMONS FOR ICE DYNAMICS.
c		=========================
c
c
c nlminn   Minimum index for computation of ice drift (NH)
c nlmaxn   Maxiumum index for computation of ice drift (NH)
c nlmins   Minimum index for computation of ice drift (SH)
c nlmaxs   Maxiumum index for computation of ice drift (SH)
c zepsd1   First tolerance parameter
c zepsd2   Second tolerance parameter
c usdt     Inverse of the time step
c alpha    Coefficient for semi-implicit coriolis
c bound    Type of boundary conditions
c dm       Diffusion constant for dynamics
c om       Relaxation constant
c resl     Maximum value for the residual of relaxation
c nbitdf   Number of iterations for free drift
c nbiter   Number of sub-time steps for relaxation
c nbitdr   Max. number of iterations for relaxation
c cw       Drag coefficient for oceanic stress
c rhoco    rho0*cw
c rhoco2   rhoco*rhoco
c angvg    Turning angle for oceanic stress
c sangvg   sin(angvg)
c cangvg   cos(angvg)
c pstarh   First bulk-rheology parameter/2
c c        Second bulk-rhelogy parameter
c zetamn   Minimun value for viscosity
c creepl   Creep limit
c usecc2   1.0/(ecc*ecc)
c sber     Test if transport of ice at Bering or not
c iameth   Method for advection
c zfn      coriolis * 2
c wght     weight of the 4 neighbours to compute averages
c akappa   first group of metric coefficients
c alambd   second group of metric coefficients
c bkappa   third group of metric coefficients
c npo1i0   number of points where there is an ocenic velocity but 
c          not an ice velocity
c ipo1i0   index i of these points
c jpo1i0   index j of these points
c ug       Ice velocity (x)
c vg       Ice velocity (y)
c uo       Ocean velosity used in ice dynamics (x)
c vo       Ocean velosity used in ice dynamics (y)
c uost     Fixed ocean velocity (x)
c vost     Fixed ocean velocity (y)
c hnm      Mean snow thickness
c hgm      Mean ice thickness
c uvdif    Diffusion velocity for scalars
c ren      Reynolds number for the grid
c gridsz   Grid size for diffusion constant
c dxs1     Lenght of the grid squares sides (x)
c dxs2     Lenght of the grid squares sides (y)
c dxc1     Lenght of the grid squares centres (x)
c dxc2     Lenght of the grid squares centres (y)
c zindfa   Mask for scalars
c area     Surface of a grid square
c dfhu     Modified diffusivity coefficient (x)
c dfhv     Modified diffusivity coefficient (y)
c
      common /latitd/ nlminn,nlmaxn,nlmins,nlmaxs
c
      common /ctbqd1/ zepsd1,zepsd2,usdt,alpha,bound,dm,om,resl,
     &                nbitdf,nbiter,nbitdr
c
      common /ctbqd2/ cw,rhoco,rhoco2,angvg,sangvg,cangvg,pstarh,c,
     &                zetamn,creepl,usecc2,sber,iameth
c
      common /comgeo/ zfn(imax,jmax),wght(imax,jmax,2,2),
     &                akappa(imax,jmax,2,2),alambd(imax,jmax,2,2,2,2),
     &                bkappa(imax,jmax,2,2)

      common /o1i0/   npo1i0,ipo1i0(5),jpo1i0(5)
c
      common /combqd/ ug(imax,jmax),vg(imax,jmax),
     &                uo(imax,jmax),vo(imax,jmax),
     &                uost(imax,jmax),vost(imax,jmax),
     &                hnm(imax,jmax),hgm(imax,jmax)
c
      common /transf/ uvdif,ren,gridsz
c
      common /comadv/ dxs1(imax,jmax),dxs2(imax,jmax),
     &                dxc1(imax,jmax),dxc2(imax,jmax),
     &                zindfa(imax,jmax),area(imax,jmax)
c
      common /comdff/ dfhu(imax,jmax),dfhv(imax,jmax)
c
