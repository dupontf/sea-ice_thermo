objcnvt=cliob.o defcst.o defgrid.o\
 flx.o forcing_nc.o gather.o \
 ocesla.o pdecli.o rappeloc.o scater.o shine.o\
 ice_th.o ice_th_diff.o ice_th_dh.o\
 PrGlobal.o biocesla.o CO2flux.o Hsb.o inibiol.o phyto.o\
 rappelbio.o resbio.o swamco.o vitsedim.o\
 ice_rad.o\
 ice_phy_remap.o ice_phy_relay.o ice_phy_grid.o\
 ice_sal_diff.o ice_sal_diff_CW.o ice_sal_column.o ice_sal_conserv.o\
 ice_bio_column.o ice_bio_conserv.o\
 ice_bio_ini.o ice_bio_diff.o ice_bio_remap.o\
 ice_bio_grid.o ice_bio_interp_phy2bio.o ice_bio_interp_bio2phy.o\
 ice_bio_sms.o\
 ice_output.o\

objlim=biolim.o defcst.o defgrid.o\
 var_thermo_vertical.o \
 ice_thermo_lim.o \
 ice_thermodynamic_FV.o\
 ice_thermodynamic_FE.o\
 flx.o forcing_nc.o gather.o init_nc.o\
 ocesla.o pdecli.o scater.o shine.o\
 ice_rad.o\
 ice_th.o ice_th_diff.o ice_th_dh.o\
 ice_phy_remap.o ice_phy_relay.o ice_phy_grid.o\
 ice_sal_diff.o ice_sal_diff_CW.o ice_sal_column.o ice_sal_conserv.o\
 ice_bio_column.o ice_bio_conserv.o\
 ice_bio_ini.o ice_bio_diff.o ice_bio_remap.o\
 ice_bio_grid.o ice_bio_interp_phy2bio.o ice_bio_interp_bio2phy.o\
 ice_bio_sms.o\
 ice_output.o\
 ice_cons.o\
 libUN.o\


# inibiol.o Ccphy.o Hsb.o Organ.o biology.o Graz.o Nutrnt.o PrGlobal.o
# profil.o resu.o mid.o flux.o glace.o resekman.o price.o pricet.o resflux.o\
# conser.o freq.o gfreq.o ekmant.o inter.o


#libncdf=/opt/intel/ifort-9.1/lib/libnetcdf_ifc.a
#libncdf=/astr/home/vancop/LIM3_P01/modipsl/lib/netcdf-64bit/libnetcdf.a
# libncdf=-lnetcdf -lnetcdff
libncdf=-lnetcdf -lnetcdff

## compil=/opt/fortran90/bin/f90
compil=gfortran -g -fno-align-commons -I/usr/include #-fbounds-check #-check all
link=$(compil)
# typical quick one
#fast=-xHOST -O3 -no-prec-div -align all -convert big_endian
fast= #-O3
extra=-fdefault-real-8
#fast= -O3 -no-prec-div -align all -convert big_endian
# check everything - !!! changes the results
#fast=-align all -convert big_endian -w -vec_report0 -zero -g -O0 -fpe0 -fp-model precise -fp-stack-check -debug all -ftrapuv -traceback -C
# checks stuff but not everything without changing the results
#fast=-align all -convert big_endian -C

CPPFLAGS=-DTHERMO_LIM
CPPFLAGS=-DTHERMO_FV
#CPPFLAGS=-DTHERMO_FE

#---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

$(objcnvt) : type.com para.com bloc.com const.com dynami.com moment.com ice.com\
 reper.com thermo.com

#---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#para.com: para.Com prep.rac
#	prepf para.Com para.com
#
#bloc.com: bloc.Com prep.rac
#	prepf bloc.Com bloc.com
#
#type.com: type.Com prep.rac
#	prepf type.Com type.com


#---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

siob: $(objcnvt)
	$(compil) $(fast) $(extra) -o runsiob $(objcnvt)

biolim: $(objlim)
	$(compil) $(fast) $(extra) -o biolim $(objlim) $(libncdf)

%.o: %.f90
	$(compil) $(fast) $(extra) -c $^ -o $@

.f.o:
	$(compil) $(fast) $(extra) -c $<

.F.o:
	$(compil) ${CPPFLAGS} $(fast) $(extra) -c $<

clean:
	rm -f *.o *.mod biolim
