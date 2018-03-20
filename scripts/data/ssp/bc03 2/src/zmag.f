	Program zmag

c	Returns galaxy magnitude in a given filter at redshift z for chosen cosmological model

c	Array declarations required in program calling function FMAG
	character file*128,filtid*64
	real k_corr_0,k_ABcorr_0,k_corr_z,k_ABcorr_z,xk_buff(2)
        common /mag/ file,icall,h,omega,omega_lambda,clambda,q,tu,tgal,tl,tz,zf,dm,filtid,
     *  f_ev,f_rf,f_ne,f_ne_0,ABf_ev,ABf_rf,ABf_ne,ABf_ne_0,k_corr_0,k_ABcorr_0,k_corr_z,k_ABcorr_z,
     *  ek_corr,ek_ABcorr,xk_buff

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
c	  tgal		= age of galaxy today
c	  tl		= light travel time from z=0 to z
c	  tz		= age of galaxy at redshift z
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
c         k_corr_z      = k-correction in filter kf computed from sed at redshift z (Vega system)
c         k_corr_0      = k-correction in filter kf computed from sed at redshift 0 (Vega system)
c         k_ABcorr_z    = k-correction in filter kf computed from sed at redshift z (AB system)
c         k_ABcorr_0    = k-correction in filter kf computed from sed at redshift 0 (AB system)
c	  ek_corr   	= (e+k)-correction in filter kf (Vega system)
c	  ek_ABcorr 	= (e+k)-correction in filter kf (AB system)

c	Array declarations specific to this program
	character name*128,aux*24
	real xaux(3)
	data icall/0/

c	Default values
c	Default model
	data file/'bc2003_hr_m62_chab_ssp.ised'/
c	Cosmology:
	data h/70./,omega/0.30/,omega_lambda/0.70/,z/0./,zl/0./
c	Default filter
	kf=15

c	Compute parameters derived from predefined parameters
c	Obtain cosmological constant and paramter q
1	clambda=cosmol_c(h,omega,omega_lambda,q)

c       Age of the universe and maximum age for galaxy (rounded off)
        tu=t(h,q,0.,clambda)
        tgal=int(tu)*1.E9
        zf=zx(tgal*1.E-9,h,q,clambda)

c	Write predefined parameters
2	write (6,*)
	write (6,*) 'Predefined parameters:'
	write (6,*)
	write (6,'(2a)')       ' Set 1 -> BC03 model in file:       ',file(1:largo(file))
	write (6,'(a,f6.2,a)') ' Set 2 -> Hubble constant H0:       ',h,' Km/s/Mpc'
	write (6,'(a,f6.2)')   '          Omega:                    ',omega
	write (6,'(a,f6.2)')   '          Omega_lambda:             ',omega_lambda
	write (6,'(a,1pe10.2)')'            Lambda:                 ',clambda
	write (6,'(a,f6.2)')   '            Parameter q:            ',q
	write (6,'(a,f7.3,a)') '            Age of the universe:    ',tu,' Gyr'
	write (6,'(a,f7.3,a)') ' Set 3 -> Age of galaxy today:      ',tgal*1.E-9,' Gyr'
	write (6,'(a,f7.3,a)') '            z of galaxy formation:  ',zf

3	write (6,*)
	write (6,'(x,a)') 'Enter redshift (z) and filter to compute magnitude'
	write (6,'(x,a)') '       To modify parameters in Set N, enter z = -N'
	write (6,'(/x,a,f6.3,a,i3,a,$)') 'Enter z, filter number [',z,',',kf,'] = '
	read (5,'(a)',err=31,end=10) aux
c	Decode z and kf
	if (largo(aux).gt.0) then
		i=index(aux,',')
		if (i.eq.0) then
			read (aux,*,err=31) z
		elseif (i.eq.1) then
			read (aux(2:),*,err=31) xp
		elseif (i.gt.1) then
			read (aux,*,err=31) z,xp
		endif
	else
		xp=kf
	endif
c	read (5,'(2f10.0)',err=31,end=10) z,xp
	goto 32
31	write (6,'(x,a)') 'Error in selection. Please try again'
	goto 3
32	if (xp.gt.0) kf=nint(xp)

	if (z.lt.0) then
		if (z.eq.-1.) then
c			Ask for *.ised file name
			write (6,'(/x,3a,$)') 'BC_GALAXEV model sed in file [',file(1:largo(file)),'] = '
			read (5,'(a)',end=10) name
			if (largo(name).gt.0) then
        			call chaext(name,'ised',nn)
        			file=name
			endif
			z=zl
			icall=0
			goto 2
		elseif (z.eq.-2.) then
c			Ask for Ho, Omega, and Omega_lambda
			write (6,'(/x,a,f4.0,a,f5.3,a,$)') 'Enter cosmological parameters'
4			write (6,'(x,a,f4.0,a,f5.3,a,f5.3,a,$)')  'Enter Ho [',h,'], Omega [',omega,'], Omega_lambda [',omega_lambda,'] = '
			call nread(xaux,nc,*5,*10)
			if (nc.gt.0.) then
				h=xaux(1)
				omega=xaux(2)
				omega_lambda=xaux(3)
			endif
			z=zl
			goto 1
5			write (6,'(x,a)') 'Error in cosmological parameters. Please try again'
			goto 4
		elseif (z.eq.-3.) then
c			Ask for galaxy age
			write (6,'(/x,a,f6.2,a  )') 'Age of this universe  = tu =',tu,' Gyr'
6			write (6,'(x,a,f6.2,a,$)') 'Enter age of galaxy [',tgal*1.e-9,' Gyr] = '
			call qread(ttg,nc,*7,*10)
			write (6,'(x,a,$)')        '        at redshift [0.000] = '
			call qread(zg,mc,*7,*10)
			if (mc.gt.0) then
				tl=tu-t(h,q,zg,clambda)
			else
				tl=0.
			endif
			if (nc.gt.0.) then
				tgal=(tl+ttg)*1.E9
				zf=zx(tgal*1.E-9,h,q,clambda)
			endif
			z=zl
			icall=0
			goto 2
7			write (6,'(x,a)') 'Error in galaxy age. Please try again'
			goto 6
		endif
	else
c		Compute required magnitude
		iwrite=1
		xmag = FMAG(z,kf,fmagAB,iwrite)
		zl=z
		goto 2
	endif
10      end
