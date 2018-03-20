RM = /bin/rm -f
UTILIB = ./public_utilities.a

# to use the g77 fortran compiler uncomment the following line
# FC =  g77 -ffixed-line-length-132

# to use the gfortran fortran compiler uncomment the following line
#FFLAGS   = -std=legacy -ffixed-line-length-none
#FFLAGS   = -std=legacy -ffixed-line-length-none -fcheck=all -Wall
FFLAGS    = -std=legacy -ffixed-line-length-none -fcheck=all -Wunused -Wuninitialized
FC        =  gfortran $(FFLAGS)

# to use the gfortarn fortran compiler in a 64 bit machine uncomment the following line
# FC =  gfortran -ffixed-line-length-none -frecord-marker=4

.f.o:
	$(FC) -c $<

add_filters: add_filters.f filter.dec $(UTILIB)
	$(FC) add_filters.f $(UTILIB) -o add_filters

seefilter: seefilter.f filter.dec $(UTILIB)
	$(FC) seefilter.f $(UTILIB) -o seefilter

sun_color: sun_color.f filter.dec $(UTILIB)
	$(FC) sun_color.f $(UTILIB) -o sun_color

galaxevpl: galaxevpl.f $(UTILIB)
	$(FC) galaxevpl.f $(UTILIB) -o galaxevpl

zeropoint: zeropoint.f $(UTILIB)
	$(FC) zeropoint.f $(UTILIB) -o zeropoint

csp_galaxev: csp_galaxev.f csp.dec $(UTILIB)
	$(FC) csp_galaxev.f $(UTILIB) -o csp_galaxev

big_csp_galaxev: ./big_csp/big_csp_galaxev.f $(UTILIB)
	cd ./big_csp ; make big_csp_galaxev ; cd ..

add_bursts: add_bursts.f $(UTILIB)
	$(FC) add_bursts.f $(UTILIB) -o add_bursts

vel_disp: vel_disp.f $(UTILIB)
	$(FC) vel_disp.f $(UTILIB) -o vel_disp

downgrade_resolution: downgrade_resolution.f $(UTILIB)
	$(FC) downgrade_resolution.f $(UTILIB) -o downgrade_resolution

ascii_ised: ascii_ised.f $(UTILIB)
	$(FC) ascii_ised.f $(UTILIB) -o ascii_ised

bin_ised: bin_ised.f $(UTILIB)
	$(FC) bin_ised.f $(UTILIB) -o bin_ised

cm_evolution: cm_evolution.f $(UTILIB)
	$(FC) cm_evolution.f $(UTILIB) -o cm_evolution

flux_evolution: flux_evolution.f $(UTILIB)
	$(FC) flux_evolution.f $(UTILIB) -o flux_evolution

mm_evolution: mm_evolution.f $(UTILIB)
	$(FC) mm_evolution.f $(UTILIB) -o mm_evolution

zmag: zmag.f $(UTILIB)
	$(FC) zmag.f $(UTILIB) -o zmag

interpolate_ised: interpolate_ised.f $(UTILIB)
	$(FC) interpolate_ised.f $(UTILIB) -o interpolate_ised

clean:
	$(RM) public_utilities.a $(OBJFILES) ./big_csp/*.o
	$(RM) bin_ised ascii_ised vel_disp csp_galaxev add_bursts downgrade_resolution
	$(RM) zeropoint galaxevpl seefilter add_filters cm_evolution mm_evolution zmag
	$(RM) sun_color interpolate_ised
	$(RM) big_csp_galaxev 

all:
	make public_utilities.a \
	bin_ised ascii_ised vel_disp csp_galaxev add_bursts downgrade_resolution \
	zeropoint galaxevpl seefilter add_filters cm_evolution mm_evolution zmag \
	sun_color interpolate_ised

#	make big_csp_galaxev

#---------------------------------------------------------------------------
#commands to build library public_utilities.a
OBJFILES = chaext.o filter_n.o largo.o sort2.o trapz1.o linear.o iplace.o \
	rminim.o rmaxim.o percent.o rf_color.o name_sed.o gw_indices.o \
	irec_check.o gw_ix_sed.o gw_ix_sed_lick_system.o clyman.o ymag2.o \
	b912.o b4000.o ew_balmer.o file_header.o dattim.o t.o \
	zx.o locate.o trapz2.o trapz3.o trapz_lr.o gaussian_v_disp.o \
	gauss.o chabimf.o gw_etal_dat.o gw_thesis_dat.o gc_d4000.o \
	sfr.o sfr_0_b.o expand_time_steps.o convolve_new.o copyright.o \
	usrsfr.o xreverse.o sort.o trapez.o delete_files.o \
	gal_mass.o sdss_color.o read_ised.o nread.o evol_ised.o \
	ke_2filt_correction.o dismod.o midpnt.o dl.o fmag.o cosmol_c.o \
	lick_system.o degrade_resolution.o trapzq.o usort.o smthpl.o \
	bh_ix_sed.o degrade_to_constant_resolution.o sigma_indous.o \
	fwhm_lick_ids.o chelium.o pgas.o nargu.o ke_nfilt_correction.o \
	mixfilenames.o file_w_ages.o no_bl.o interpolate_ssp.o s500.o \
	intrp_ised.o

public_utilities.a: $(OBJFILES)
	ar ruv public_utilities.a $?
	ranlib public_utilities.a
