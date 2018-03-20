	SUBROUTINE ke_Nfilt_correction(kf,x,y,n,z,h0,q0,lun,nmag)

c	Compute apparent magnitude in nmag filters contained in array kf

c	Array declarations
	include 'filter.dec'
	include 'cosmo.dec'
	character*1024 aux
	integer kf(nmag)
	real x(n),y(n),f0(50),f1mag(50),f1magAB(50)
	data icall/0/
	save f0

c	Report cosmology, selected filters, etc.
	if (icall.eq.0) then
		icall=1
		write (lun,100)      h,omega,omega_lambda,q,clambda,tu,ttg,zf
		write (lun  ,'(a)') '# Meaning of printed variables (see Tables 5 and 6 of bc03.ps for details):'
		write (lun  ,'(a)') '#   z        = redshift (z=0 entry in table corresponds to a distance of 10 pc)'
		write (lun  ,'(a)') '#   ltt      = light travel time from z=0 to z in Gyr'
		write (lun  ,'(a)') '#   tz       = age of galaxy at redshit = z'
		write (lun  ,'(a)') '#   dm       = cosmological distance modulus in magnitude units'
		write (lun  ,'(a)') '#   m_ev     = apparent magnitude M_ev(z) + dm(z) in each filter (Vega system)'
		write (lun  ,'(a)') '#'
		write (lun+1,100)    h,omega,omega_lambda,q,clambda,tu,ttg,zf
		write (lun+1,'(a)') '# Meaning of printed variables (see Tables 5 and 6 of bc03.ps for details):'
		write (lun+1,'(a)') '#   z        = redshift (z=0 entry in table corresponds to a distance of 10 pc)'
		write (lun+1,'(a)') '#   ltt      = light travel time from z=0 to z in Gyr'
		write (lun+1,'(a)') '#   tz       = age of galaxy at redshit = z'
		write (lun+1,'(a)') '#   dm       = cosmological distance modulus in magnitude units'
		write (lun+1,'(a)') '#   m_ev     = apparent magnitude M_ev(z) + dm(z) in each filter (AB system)'
		write (lun+1,'(a)') '#'
		write (lun,'(a,19x,3x,50i9)') '# Filter:',(kf(i),i=1,nmag)
		write (lun+1,'(a,19x,3x,50i9)') '# Filter:',(kf(i),i=1,nmag)
		write (lun,101) ('   m_ev  ',i=1,nmag)
		write (lun+1,101) ('   m_ev  ',i=1,nmag)
		write (lun,102) ('   mag   ',i=1,nmag)
		write (lun+1,102) ('   ABmag ',i=1,nmag)
		write (aux,103) (i+4,i=1,nmag)
		l=largo(aux)
		aux(l:l)=' '
		l=largo(aux)
		write (lun,'(a)') aux(1:l)
		write (lun+1,'(a)') aux(1:l)
	endif
100     format ('# Cosmology: Ho = ',f4.0,2x,'Omega =',f5.2,2x,'Omega_lambda =',f5.2,2x,'qo = ',f7.4,2x,
     *          'Lambda =',1pE10.3,2x,'tu = ',0pf5.2,' Gyr',2x,'tg = ',f5.2,' Gyr', 2x,' zf = ',f6.2/'#')
101	format ('#  z      LTT     tz      dm      ',50a) 
102	format ('#         Gyr     Gyr     mag     ',50a)
c103	format ('# (1)     (2)     (3)     (4)   ',50(5x,'(',i2,')'))
103	format ('# (1)     (2)     (3)     (4)        (',50(i2,')     ('))

c	Compute cosmological distance modulus.
	dm=dismod(h0,q0,z)

c	Compute light travel time from z=0 to z
	tl=tu-t(h0,q0,z,clambda)
	tz=ttg-tl

c	If z=0, compute zero points
	if (z.eq.0.) then
c		Compute zero points in filters in array kf(i)
		do i=1,nmag
		f0(i)=vega_0p_n(kf(i))
		enddo
	endif

c	Compute magnitudes at this redshift
	do i=1,nmag
c	Compute absolute magnitude in filters in array kf(i) (Vega system)
	f1_ev   = f0(i) - 2.5*alog10(filter_n(kf(i),x,y,n,z,kerr))
c	Compute absolute AB magnitude in filters in array kf(i)
	ABf1_ev = AB_mag(kf(i),x,y,n,z)
c	Compute apparent magnitudes
	f1mag(i)    = f1_ev   + dm
	f1magAB(i)  = ABf1_ev + dm
	enddo

c	Write results
	if (z.eq.0.) then
		write (lun,  108) z,tl,tz,dm,(f1mag(i),i=1,nmag)
		write (lun+1,108) z,tl,tz,dm,(f1magAB(i),i=1,nmag)
	else
		write (lun,  109) z,tl,tz,dm,(f1mag(i),i=1,nmag)
		write (lun+1,109) z,tl,tz,dm,(f1magAB(i),i=1,nmag)
	endif
108	format ('#',f5.3,3f8.3,3x,100f9.4)
109	format (f6.3,3f8.3,3x,100f9.4)
	return
	end
