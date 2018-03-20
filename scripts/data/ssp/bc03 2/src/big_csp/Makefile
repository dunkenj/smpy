PU	  = ../public_utilities.a

FFLAGS	  = -std=legacy -ffixed-line-length-none -fcheck=all -Wunused -Wuninitialized
FC	  = gfortran $(FFLAGS)
CFLAGS	  = -c

RM = /bin/rm -f

.f.o:
	$(FC) -c $<

#commands to build library big_csp.a
OBJFILES = convolve_tx.o expand_time_steps.o sfr_0_b.o usrsfr.o

# 	$(FC) big_csp_galaxev.f -mcmodel=medium $(OBJFILES) $(PU) -o ../big_csp_galaxev

big_csp_galaxev: big_csp_galaxev.f csp.dec $(OBJFILES) $(PU)
	$(FC) big_csp_galaxev.f $(OBJFILES) $(PU) -o ../big_csp_galaxev
