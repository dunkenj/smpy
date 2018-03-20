	SUBROUTINE ke_2filt_correction(kf,x,y,n,z,h0,q0,lun,imag)

c	Compute k and e correction in filter NF

c	Array declarations
	include 'filter.dec'
	include 'cosmo.dec'
	include 'SSP_0.dec'
	integer kf(2)
	real k1_corr_z,k1_ABcorr_z,k2_corr_z,k2_ABcorr_z,k12_corr_0,k12_ABcorr_0
	real k1_corr_0,k1_ABcorr_0,k2_corr_0,k2_ABcorr_0,k12_corr_z,k12_ABcorr_z
	real x(n),y(n),y0(imw)
	data icall/0/
	save f1,f2,y0

c	Report cosmology, selected filters, etc.
	if (icall.eq.0) then
		icall=1
		write (lun,100) h,omega,omega_lambda,q,clambda,tu,ttg,zf
		call meaning(1,lun,kf(1))
		write (lun,'(a,19x,2(3x,4i9),3x,3i9)') '# Filter:',(kf(1),i=1,11)
		write (lun,101)
		write (lun,102)
		write (lun,103)
		if (imag.eq.2) then
			write (lun+1,100) h,omega,omega_lambda,q,clambda,tu,ttg,zf
			call meaning(1,lun+1,kf(2))
       			write (lun+1,'(a,19x,2(3x,4i9),3x,3i9)') '# Filter:',(kf(2),i=1,11)
			write (lun+1,101)
			write (lun+1,102)
			write (lun+1,103)

			write (lun+2,100) h,omega,omega_lambda,q,clambda,tu,ttg,zf
			call meaning(2,lun+2,kf(2))
			write (lun+2,107) '# Filters:',(kf(1),kf(2),i=1,9)
			write (lun+2,104)
			write (lun+2,105)
			write (lun+2,106)
		endif
100     format ('# Cosmology: Ho = ',f4.0,2x,'Omega =',f5.2,2x,'Omega_lambda =',f5.2,2x,'qo = ',f7.4,2x,
     *          'Lambda =',1pE10.3,2x,'tu = ',0pf5.2,' Gyr',2x,'tg = ',f5.2,' Gyr', 2x,' zf = ',f6.2/'#')
101	format ('#  z      LTT     tz      dm         M_rf     M_ne     M_ev     m_ev',
     *          '      M_AB_rf  M_AB_ne  M_AB_ev  m_AB_ev     e+k_cor  k_cor_ev k_cor_ne')
102	format ('#         Gyr     Gyr     mag        mag      mag      mag      mag',
     *          '        AB_mag   AB_mag   AB_mag   AB_mag       mag      mag      mag')
103	format ('# (1)     (2)     (3)     (4)        (5)      (6)      (7)      (8)',
     *          '         (9)      (10)     (11)     (12)        (13)     (14)     (15)')
104	format ('#  z      LTT     tz      dm         C_rf     C_ne     C_ev',
     *          '      C_AB_rf  C_AB_ne  C_AB_ev     e+k_cor  k_cor_ev k_cor_ne')
105	format ('#         Gyr     Gyr     mag        mag      mag      mag',
     *          '        AB_mag   AB_mag   AB_mag       mag      mag      mag')
106	format ('# (1)     (2)     (3)     (4)        (5)      (6)      (7)',
     *          '         (8)      (9)      (10)        (11)     (12)     (13)')
107	format (a,23x,3(2x,i3,'-',i3),3x,3(2x,i3,'-',i3),3x,3(2x,i3,'-',i3))
	endif

c	Compute cosmological distance modulus.
	dm=dismod(h0,q0,z)

c	Compute light travel time from z=0 to z
	tl=tu-t(h0,q0,z,clambda)
	tz=ttg-tl

c	If z=0, store sed to compute e+k and k-corrections
	if (z.eq.0.) then
		do i=1,n
		y0(i)=y(i)
		enddo
c		Compute zero points in filters kf(1) and kf(2) in the Vega system
		f1=vega_0p_n(kf(1))
		if (imag.eq.2) f2=vega_0p_n(kf(2))
	endif

c	Compute absolute magnitude in filter kf(1) (Vega system)
	f1_ev   = f1 - 2.5*alog10(filter_n(kf(1),x,y ,n,z,kerr))
	f1_rf   = f1 - 2.5*alog10(filter_n(kf(1),x,y ,n,0.,kerr))
	f1_ne   = f1 - 2.5*alog10(filter_n(kf(1),x,y0,n,z,kerr))
	f1_ne_0 = f1 - 2.5*alog10(filter_n(kf(1),x,y0,n,0.,kerr))
	f1_ev_0 = f1_rf

c	Compute absolute AB magnituds in filter kf(1)
	ABf1_ev   = AB_mag(kf(1),x,y ,n,z)
	ABf1_rf   = AB_mag(kf(1),x,y ,n,0.)
	ABf1_ne   = AB_mag(kf(1),x,y0,n,z)
	ABf1_ne_0 = AB_mag(kf(1),x,y0,n,0.)
	ABf1_ev_0 = ABf1_rf

c	Compute apparent magnitude
	f1mag    = f1_ev   + dm
	f1magAB  = ABf1_ev + dm

c	Compute k-correction
c	   For z=0 sed
	      k1_corr_0   = f1_ne   - f1_ne_0
	      k1_ABcorr_0 = ABf1_ne - ABf1_ne_0
c	   For evolved sed
	      k1_corr_z   = f1_ev   - f1_rf
	      k1_ABcorr_z = ABf1_ev - ABf1_rf

c	Compute (e+k)-correction
	ek1_corr   = f1_ev   - f1_ne_0
	ek1_ABcorr = ABf1_ev - ABf1_ne_0

c	Verify if second filter has been entered
	if (imag.eq.2) then
c		Compute absolute magnitude in filter kf(2) (Vega system)
		f2_ev   = f2 - 2.5*alog10(filter_n(kf(2),x,y ,n,z,kerr))
		f2_rf   = f2 - 2.5*alog10(filter_n(kf(2),x,y ,n,0.,kerr))
		f2_ne   = f2 - 2.5*alog10(filter_n(kf(2),x,y0,n,z,kerr))
		f2_ne_0 = f2 - 2.5*alog10(filter_n(kf(2),x,y0,n,0.,kerr))
		f2_ev_0 = f2_rf

c		Compute absolute AB magnituds in filter kf(2)
		ABf2_ev   = AB_mag(kf(2),x,y ,n,z)
		ABf2_rf   = AB_mag(kf(2),x,y ,n,0.)
		ABf2_ne   = AB_mag(kf(2),x,y0,n,z)
		ABf2_ne_0 = AB_mag(kf(2),x,y0,n,0.)
		ABf2_ev_0 = ABf2_rf

c		Compute apparent magnitude
		f2mag    = f2_ev   + dm
		f2magAB  = ABf2_ev + dm
	
c		Compute k-correction
c		   For z=0 sed
		      k2_corr_0   = f2_ne   - f2_ne_0
		      k2_ABcorr_0 = ABf2_ne - ABf2_ne_0
c		   For evolved sed
		      k2_corr_z   = f2_ev   - f2_rf
		      k2_ABcorr_z = ABf2_ev - ABf2_rf

c		Compute (e+k)-correction
		ek2_corr   = f2_ev   - f2_ne_0
		ek2_ABcorr = ABf2_ev - ABf2_ne_0

c		Compute color kf(1) - kf(2) in Vega system
		c12_ev   = f1_ev   - f2_ev
		c12_rf   = f1_rf   - f2_rf
		c12_ne   = f1_ne   - f2_ne
		c12_ne_0 = f1_ne_0 - f2_ne_0
		c12_ev_0 = f1_ev_0 - f2_ev_0

c		Compute color kf(1) - kf(2) in AB system
		ABc12_ev   = ABf1_ev   - ABf2_ev
		ABc12_rf   = ABf1_rf   - ABf2_rf
		ABc12_ne   = ABf1_ne   - ABf2_ne
		ABc12_ne_0 = ABf1_ne_0 - ABf2_ne_0
		ABc12_ev_0 = ABf1_ev_0 - ABf2_ev_0

c		Compute color k-correction
		k12_corr_0   = c12_ne   - c12_ne_0
		k12_corr_z   = c12_ev   - c12_ev_0
		k12_ABcorr_0 = ABc12_ne - ABc12_ne_0
		k12_ABcorr_z = ABc12_ev - ABc12_ev_0

c		Compute color (e+k)-corrections
		ek12_corr   = c12_ev   - c12_ne_0
		ek12_ABcorr = ABc12_ev - ABc12_ne_0
	endif

c	Write results
	if (z.eq.0.) then
		write (lun,108) z,tl,tz,dm,f1_rf,f1_ne,f1_ev,f1mag,
     *                          ABf1_rf,ABf1_ne,ABf1_ev,f1magAB,ek1_ABcorr,k1_ABcorr_z,k1_ABcorr_0
		if (imag.eq.2) then
			write (lun+1,108) z,tl,tz,dm,f2_rf,f2_ne,f2_ev,f2mag,
     *                          ABf2_rf,ABf2_ne,ABf2_ev,f2magAB,ek2_ABcorr,k2_ABcorr_z,k2_ABcorr_0
			write (lun+2,109) z,tl,tz,dm,c12_rf,c12_ne,c12_ev,
     *                          ABc12_rf,ABc12_ne,ABc12_ev,ek12_ABcorr,k12_ABcorr_z,k12_ABcorr_0
		endif
	else
		write (lun,110) z,tl,tz,dm,f1_rf,f1_ne,f1_ev,f1mag,
     *                          ABf1_rf,ABf1_ne,ABf1_ev,f1magAB,ek1_ABcorr,k1_ABcorr_z,k1_ABcorr_0
		if (imag.eq.2) then
			write (lun+1,110) z,tl,tz,dm,f2_rf,f2_ne,f2_ev,f2mag,
     *                          ABf2_rf,ABf2_ne,ABf2_ev,f2magAB,ek2_ABcorr,k2_ABcorr_z,k2_ABcorr_0
			write (lun+2,111) z,tl,tz,dm,c12_rf,c12_ne,c12_ev,
     *                          ABc12_rf,ABc12_ne,ABc12_ev,ek12_ABcorr,k12_ABcorr_z,k12_ABcorr_0
		endif
	endif
108	format ('#',f5.3,3f8.3,3x,4f9.4,3x,4f9.4,3x,3f9.4)
109	format ('#',f5.3,3f8.3,3x,3f9.4,3x,3f9.4,3x,3f9.4)
110	format (f6.3,3f8.3,3x,4f9.4,3x,4f9.4,3x,3f9.4)
111	format (f6.3,3f8.3,3x,3f9.4,3x,3f9.4,3x,3f9.4)
	return
	end

	SUBROUTINE MEANING(IFILE,LUN,NF)

c	Prints meaning of different quantities in output files
	if (ifile.eq.1) then
		write (lun,'( a)')     '# Meaning of printed variables (see Tables 5 and 6 of bc03.ps for details):'
		write (lun,'( a)')     '#   z        = redshift (z=0 entry in table corresponds to a distance of 10 pc)'
		write (lun,'( a)')     '#   ltt      = light travel time from z=0 to z in Gyr'
		write (lun,'( a)')     '#   tz       = age of galaxy at redshit = z'
		write (lun,'( a)')     '#   dm       = cosmological distance modulus in magnitude units'
		write (lun,'(a,i3,a)') '#   M_rf     = M_rf(z) in filter ',nf,' (Vega system)'
		write (lun,'(a,i3,a)') '#   M_ne     = M_ne(z) in filter ',nf,' (Vega system)'
		write (lun,'(a,i3,a)') '#   M_ev     = M_ev(z) in filter ',nf,' (Vega system)'
		write (lun,'(a,i3,a)') '#   m_ev     = apparent magnitude M_ev(z) + dm(z) in filter ',nf,' (Vega system)'
		write (lun,'(a,i3,a)') '#   M_AB_rf  = AB M_rf(z) in filter ',nf,' (AB system)'
		write (lun,'(a,i3,a)') '#   M_AB_ne  = AB M_ne(z) in filter ',nf,' (AB system)'
		write (lun,'(a,i3,a)') '#   M_AB_ev  = AB M_ev(z) in filter ',nf,' (AB system)'
		write (lun,'(a,i3,a)') '#   m_AB_ev  = apparent AB magnitude M_ev(z) + dm(z) in filter ',nf,' (AB system)'
		write (lun,'(a,i3,a)') '#   e+k_cor  = (e+k)-correction in filter ',nf
		write (lun,'(a,i3,a)') '#   k_cor_ev = k-correction computed for sed at redshift z in filter ',nf
		write (lun,'(a,i3,a)') '#   k_cor_ne = k-correction computed for sed at redshift 0 in filter ',nf
		write (lun,'( a)')     '#'
	elseif (ifile.eq.2) then
		write (lun,'( a)')     '# Meaning of printed variables (see Tables 5 and 6 of bc03.ps for details):'
		write (lun,'( a)')     '#   z        = redshift (z=0 entry in table corresponds to a distance of 10 pc)'
		write (lun,'( a)')     '#   ltt      = light travel time from z=0 to z in Gyr'
		write (lun,'( a)')     '#   tz       = age of galaxy at redshit = z'
		write (lun,'( a)')     '#   dm       = cosmological distance modulus in magnitude units'
		write (lun,'( a)')     '#   C_rf     = color M1_rf(z) - M2_rf(z) (Vega system)'
		write (lun,'( a)')     '#   C_ne     = color M1_ne(z) - M2_ne(z) (Vega system)'
		write (lun,'( a)')     '#   C_ev     = color M1_ev(z) - M2_ev(z) (Vega system)'
		write (lun,'( a)')     '#   C_AB_rf  = AB color M1_rf(z) - M2_rf(z) (AB system)'
		write (lun,'( a)')     '#   C_AB_ne  = AB color M1_ne(z) - M2_ne(z) (AB system)'
		write (lun,'( a)')     '#   C_AB_ev  = AB color M1_ev(z) - M2_ev(z) (AB system)'
		write (lun,'( a)')     '#   e+k_cor  = color (e+k)-correction'
		write (lun,'(a,i3,a)') '#   k_cor_ev = color k-correction computed for sed at redshift z'
		write (lun,'(a,i3,a)') '#   k_cor_ne = color k-correction computed for sed at redshift 0'
		write (lun,'( a)')     '#'
	endif
	return
	end
