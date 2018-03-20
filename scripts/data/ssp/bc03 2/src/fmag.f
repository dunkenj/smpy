	Real*4 FUNCTION FMAG(z,kf,fmagAB,iwrite)

c	Returns magnitude in filter kf at redshift Z for pre-chosen cosmological model

c	Array declarations required also in calling program
	character file*128,filtid*64
	real k_corr_0,k_ABcorr_0,k_corr_z,k_ABcorr_z,xk_buff(2)
	common /mag/ file,icall,h,omega,omega_lambda,clambda,q,tu,tgal,tl,tz,zf,dm,filtid,
     *	f_ev,f_rf,f_ne,f_ne_0,ABf_ev,ABf_rf,ABf_ne,ABf_ne_0,k_corr_0,k_ABcorr_0,k_corr_z,k_ABcorr_z,
     *	ek_corr,ek_ABcorr,xk_buff

c	Meaning of variables in call to function fmag (see Tables 5 and 6 of bc03.ps for details):
c	  z		= redshift to compute apparent magnitude
c	  kf		= number of filter chosen to compute galaxy magnitude
c	  fmag		= apparent magnitude M_ev(z) + dm(z) in filter kf (Vega system)
c	  fmagAB	= apparent AB magnitude M_ev(z) + dm(z) in filter kf (AB system)
c	  iwrite	= print screen output (iwrite > 0)

c	Meaning of variables in common /mag/
c	  file		= BC GALAXEV *.ised model file name
c	  icall		= 0 when calling function for the first time
c	  h		= Hubble constant (Km/s/Mpc)
c	  omega		= density parameter
c	  omega_lambda	= density parameter associated to cosmological constant
c	  clambda     	= cosmological constant
c	  q		= deceleration parameter
c	  tu		= age of the universe in this cosmology
c	  tgal		= age of galaxies
c	  tl		= light travel time from z=0 to z
c	  tz		= age of galaxies at redshift z
c	  zf		= redshift of galaxy formation
c	  dm		= cosmological distance modulus in magnitude units
c	  filtid	= filter identification
c	  f_ev   	= M_ev(z) in filter kf (Vega system)
c	  f_rf   	= M_rf(z) in filter kf (Vega system)
c	  f_ne   	= M_ne(z) in filter kf (Vega system)
c	  f_ne_0 	= M_ne(0) in filter kf (Vega system)
c	  ABf_ev   	= AB M_ev(z) in filter kf (AB system)
c	  ABf_rf   	= AB M_rf(z) in filter kf (AB system)
c	  ABf_ne   	= AB M_ne(z) in filter kf (AB system)
c	  ABf_ne_0 	= AB M_ne(0) in filter kf (AB system)
c	  k_corr_z    	= k-correction in filter kf computed from sed at redshift z (Vega system)
c	  k_corr_0    	= k-correction in filter kf computed from sed at redshift 0 (Vega system)
c	  k_ABcorr_z  	= k-correction in filter kf computed from sed at redshift z (AB system)
c	  k_ABcorr_0  	= k-correction in filter kf computed from sed at redshift 0 (AB system)
c	  ek_corr   	= (e+k)-correction in filter kf (Vega system)
c	  ek_ABcorr 	= (e+k)-correction in filter kf (AB system)

c	Arrray declarations internal to fmag
	include 'read_ised.dec'
	include 'filter.dec'
	real yz(imw),y0(imw)
	save kl,last,y0
	data jcall/0/

c	Read *.ised file if first call to routine
	if (icall.eq.0) then
c		Read filter file
		write (6,*)
		write (6,*) 'Reading files and setting up zero points...'
		if (iread.eq.0) call filter0
c		Read galaxy sed
1		call read_ised(file,tgal,kl,kerr)
		if (kerr.ne.0) goto 2
		jcall=0
c		Compute evolution vs z
		call evol_ised(kl,0.,h,q,clambda,yz,jcall,last)
		write (6,*)
		icall=1
c		Store z=0 sed in array y0
		do i=1,iw
		y0(i)=yz(i)
		enddo
		goto 3
2		write (6,*)
c		Ask for *.ised file name
		write (6,'(x,3a,$)') 'BC_GALAXEV model sed in file [',file(1:largo(file)),'] = '
		read (5,'(a)',end=10) file
		call chaext(file,'ised',nn)
		goto 1
	endif

c	Compute cosmological distance modulus
3	dm=dismod(h,q,z)

c	Compute light travel time from z=0 to z
	tl=tu-t(h,q,z,clambda)
	tz=tgal*1.E-9-tl

c	Compute zero point in filter kf (Vega system)
	f0=vega_0p_n(kf)

c	Load filter id for filter kf
	filtid=fid(kf)

c	Get evolved sed at redshift z = yz
	call evol_ised(kl,z,h,q,clambda,yz,jcall,last)

c	Compute absolute magnitude in filter kf (Vega system)
c	  Magnitude of evolved sed at redshift z
	     f_ev   = f0 - 2.5*alog10(filter_n(kf,w,yz,iw,z,lerr))
c	  Magnitude of evolved sed at redshift 0
	     f_rf   = f0 - 2.5*alog10(filter_n(kf,w,yz,iw,0.,lerr))
c	  Magnitude of z = 0 sed at redshift z
	     f_ne   = f0 - 2.5*alog10(filter_n(kf,w,y0,iw,z,lerr))
c	  Magnitude of z = 0 sed at redshift 0
	     f_ne_0 = f0 - 2.5*alog10(filter_n(kf,w,y0,iw,0.,lerr))

c	Compute absolute AB magnitude in filter kf
c	  Magnitude of evolved sed at redshift z
	     ABf_ev   = AB_mag(kf,w,yz,iw,z)
c	  Magnitude of evolved sed at redshift 0
	     ABf_rf   = AB_mag(kf,w,yz,iw,0.)
c	  Magnitude of z = 0 sed at redshift z
	     ABf_ne   = AB_mag(kf,w,y0,iw,z)
c	  Magnitude of z = 0 sed at redshift 0
	     ABf_ne_0 = AB_mag(kf,w,y0,iw,0.)

c	Compute apparent magnitudes
	fmag    = f_ev   + dm
	fmagAB  = ABf_ev + dm

c	Compute k-corrections:
c	  For z=0 sed
	      k_corr_0   = f_ne   - f_ne_0
              k_ABcorr_0 = ABf_ne - ABf_ne_0
c	  For evolved sed
	      k_corr_z   = f_ev   - f_rf
              k_ABcorr_z = ABf_ev - ABf_rf

c	Compute (e+k)-corrections
	ek_corr   = f_ev   - f_ne_0
        ek_ABcorr = ABf_ev - ABf_ne_0

c	Print results if required
	if (iwrite.gt.0) then
		write (6,'(/x,a)') 'Meaning of printed variables (see Tables 5 and 6 of bc03.ps for details):'
		write (6,'(a)')    '  z       = redshift'
		write (6,'(a)')    '  ltt     = light travel time from z=0 to z in Gyr'
		write (6,'(a)')    '  tz      = age of galaxy at redshit = z'
		write (6,'(a)')    '  dm      = cosmological distance modulus in magnitude units'
		write (6,'(a)')    '  M_rf    = M_rf(z) in filter kf'
		write (6,'(a)')    '  M_ne    = M_ne(z) in filter kf'
		write (6,'(a)')    '  M_ev    = M_ev(z) in filter kf'
		write (6,'(a)')    '  m_ev    = apparent magnitude M_ev(z) + dm(z) in filter kf'
		write (6,'(a)')    '  e+kcor  = (e+k)-correction in filter kf'
		write (6,'(a)')    '  k_cor_z = k-correction in filter kf computed from sed at redshift z'
		write (6,'(a)')    '  k_cor_0 = k-correction in filter kf computed from sed at redshift 0'
		write (6,'(a)')    '  Z.P.    = Zeropoint definition, either Vega or AB system'
		write (6,'(/x,a,i4,2a)') 'In filter No.',kf,' = ',filtid(1:largo(filtid))
		write (6,'(3x,2a)') 'z   ltt(Gyr) tz(Gyr)   dm      M_rf     M_ne     M_ev     m_ev    e+kcor   k_cor_z  k_cor_0  Z.P.'
		write (6,101) z,tl,tz,dm,f_rf,f_ne,f_ev,fmag,ek_corr,k_corr_z,k_corr_0,ABf_rf,ABf_ne,ABf_ev,fmagAB,ek_ABcorr,k_ABcorr_z,k_ABcorr_0
101		format (f6.3,3f8.3,7f9.4,3x,'Vega'/30x,7f9.4,4x,'AB')
		write (6,'(x,140a)') ('-',i=1,61),'^',('-',i=1,37)
        endif
	return
10	stop
	end
