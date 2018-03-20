	FUNCTION FILTER_N(I,X,Y,N,Z,KERR)

!	Use FILTR_MASTER function to compute required quantity

!	Array declaration
	real x(n),y(n)

!	Returns number of photons through Ith filter.
	filter_n = filtr_master(i,x,y,n,z,0,kerr)
	return

	ENTRY FILTER(I,X,Y,N,Z,KERR)
!	Returns total flux through Ith filter.
	filter   = filtr_master(i,x,y,n,z,1,kerr)
	return

	ENTRY F_MEAN(I,X,Y,N,Z,KERR)
!	Returns mean flux >>> (Fnu) <<< in ith filter (/Hz)
!	Use FILTR_MASTER function to compute required quantity
	f_mean  = filtr_master(i,x,y,n,z,2,kerr)
	return
	end

	FUNCTION FILTR_MASTER(I,X,Y,N,Z,IK,KERR)

!	Optimized in memory usage by G. Bruzual. 03-NOV-2008, Merida
!	Adapted and revised by G. Bruzual.       11-AUG-2003, Merida
!	Last version written by G. Bruzual.      10-May-1999, Merida
!	Original version written by G. Bruzual.  21-May-1983, Durham

!	Returns number of photons through I filter.
!	The s.e.d. is Y(X) with N data points, assumed flambda vs lambda.
!	Z is the redshift to be applied to the s.e.d.

!	Improved version that samples correctly the s.e.d. and the
!	filter response function at all points available in these arrays.

!	The new version of the F_MEAN routine requires the filter
!	response function to be interpolated at each point in the
!	sed which is not a point in the filter response. All points in
!	the filter response are also used. The format of the binary
!	file was changed to reduce its size. Filters are stored
!	sequentially in r(i,k), k=1 => wavelength, k=2 => response
!	function. Information about starting point and number of
!	points per filter is kept in file, as well as filter id label.

!	For a flat = 1 filter response function, and Y in the standard Lo/A
!	units used by bc_galaxev, the actual number of photons is 1.926E41xfilter_n

	INCLUDE 'filter.dec'
	real x(n),y(n),LINEAR,rfa(naux),rfb(naux)

!	Extract filter I from binary filter file
        if (i <= 0) then
!		Introduced to allow computation of K-correction with routine ~/is/k_correct.f
		filtr_master=1.
                return
	endif

!	Extract Ith filter in position K of filter arrays
	call build_filter_array(i,x,y,n,z,k)

!	Interpolate sed at shifted wavelength of filter
!	Select quantity to compute
	filtr_master=0.
	kerr=0.
	z1=1.+z
	m=0
	l=0
	do j=pos_i(k),pos_f(k)
	m=m+1
	if (ik <= 0) then
!		Compute filter_n = number of photons through filter
!		Compute number of photons, instead of total flux (Fukugita et al. 1996)
		xf(m)=xlam(j)
		rf(m)=z1*xlam(j)*rlam(j)*LINEAR(xf(m),x,y,n,l)

	elseif (ik <= 1) then
!		Compute filter = flux through filter
		xf(m)=xlam(j)
		rf(m)=rlam(j)*LINEAR(xf(m),x,y,n,l)

	elseif (ik <= 2) then
!		Compute f_mean = mean flux Fnu through filter
!		Returns mean flux >>> (Fnu) <<< in ith filter (/Hz)
!		As: f_nu=int(dnu Fnu Rnu/h*nu)/int(dnu Rnu/h*nu)
!		ie: f_nu=int(dlm Flm Rlm lm / c)/int(dlm Rlm/lm)
		xf(m) =xlam(j)*z1	! wavelength in detector''s frame
		rfa(m)=rlam(j)*xf(m)*LINEAR(xlam(j),x,y,n,l)
		rfb(m)=rlam(j)/xf(m)

	else
		write (6,*) 'FILTR_MASTER: unknown option:',ik
		stop
	endif
	enddo

!	Compute integral below filter. Factor (1+z) already included in lambda scale.
	if (ik <= 1) then
		filtr_master=TRAPZ1(xf,rf,m)
	else
		f_mean=TRAPZ1(xf,rfa,m)/TRAPZ1(xf,rfb,m)
		f_mean=f_mean/2.997925e+18
!		F(lambda)*dlambda = F[lambda/(1+z)]*dlambda/(1+z)
		f_mean=f_mean/z1
		filtr_master=f_mean
	endif

!	Check result
	if (filtr_master <= 0.) then
		kerr=1
		if (ierr(i).eq.0) then
			ierr(i)=1
!			ierr(i)=0
!			write (6,*) '--- Negative or zero flux through filter',i,filtr_master,' ---'
!			write (6,*) '---   Error reported only once. It may occur more than once. ---'
!			write (6,*) '--- Negative or zero flux through filter',i,filtr_master,' Error reported only once. ---'
      		endif
	endif
	return
	end

	SUBROUTINE BUILD_FILTER_ARRAY(I,X,Y,N,Z,K)

!	Extract filter I from binary filter file and adds points
!	corresponding to wavelength points in sed.

	INCLUDE 'filter.dec'
	real x(n),y(n),linear
        data iread/0/,jread/0/,ireset/1/,zlast/-100./
        data ierr/mxft*0/,zzlast/-100./,nlast/0/

!	Read filter file
	if (iread.eq.0) Call FILTER0

!	Check filter number
	if (i.gt.mxft) then
		write (6,*) 'Filter No.',i,' not available.'
		stop
	endif

!	Check z value
	z1=1.+z
	if (z.ne.zzlast.or.ireset.gt.0.or.n.ne.nlast) then
!		write (6,*)
! 		write (6,*) '--- Reset filter arrays...z,zzlast,ireset,n,nlast =',z,zzlast,ireset,n,nlast
! 		write (6,*) '--- Wavelength array in use: n, w(1), w(n)        =',n,x(1),x(n)
		ireset=0
		zzlast=z
		nlast=n
		nfilt=0
		ntot=0
		pos_i(1)=1
		do k=1,mxft
		ifilt(k)=-10
		ierr(k) =0
		enddo
	endif

!	Extract filter
!	Check if filter has been stored already
	do k=1,nfilt
	if (ifilt(k).eq.i) return
	enddo

!	Extract ith-filter. Shift by (1+z)
	m=0
	do k=ni(i),nl(i)
	m=m+1
	xf(m)=r(k,1)/z1
	rf(m)=r(k,2)
	enddo

!	Add wavelength points in the sed (x,y) in the galaxy restframe
	l=0
	do k=1,n
!	xz=x(k)/z1	! This statement is wrong. x(i) is already in the galaxy rest frame
	xz=x(k)		! This is correct. Emma C. Lake noticed this error in March 2017. !!!!!!!!!!!
	if (xz >= xf(np(i))) then
		goto 1
	elseif (xz >= xf(1)) then
		m=m+1
		xf(m)=xz
!		Interpolate shifted filter at shifted wavelength of sed
		rf(m)=LINEAR(xz,xf,rf,np(i),l)
	endif
	enddo

!	Sort (xf,rf) arrays according to xf
1	call sort2(m,xf,rf)
!	write (6,*) '                      ',m,' total points filters + SED shifted to z =',z
!	write (6,*)
!	write (6,*)   '--- Rebuilding filter array for filter number:',i,nl(i)-ni(i)+1,ni(i),nl(i),m,xf(1),xf(m)

!	Store sorted arrays
	do k=1,m
	ntot=ntot+1
	if (ntot.gt.nbuff) then
		write (6,*) 'STOP. Size of filter buffer exceeded:',ntot,nbuff
		stop
	endif
	xlam(ntot)=xf(k)
	rlam(ntot)=rf(k)
	enddo

!	Store filter No. and ending point in sorted arrays
	nfilt=nfilt+1
	ifilt(nfilt)=i
	pos_f(nfilt)=ntot
	pos_i(nfilt+1)=ntot+1
	k=nfilt
	kt=pos_f(k)-pos_i(k)+1
	if (kt.gt.naux) then
		write (6,*) 'STOP. Size of filter buffer exceeded:',kt,naux
		stop
	endif
!	write (6,'(i2,i4,3i8,3x,a)') nfilt,i,kt,pos_i(k),pos_f(k),(fid(i)(1:57))
	return
	end

	SUBROUTINE FILTER0

!	Reads filter response functions from file FILTERBIN.RES

	INCLUDE 'filter.dec'
	character filtfile*256
	close (81)

!	SUN Unixf77: Get file name from environment variable FILTERS
	call getenv('FILTERS',filtfile)
	open (81,file=filtfile,form='unformatted',status='old')

!	VAX VMS Fortran: Read file assigned to FILTERS
!	open (81,name='FILTERS',form='unformatted',status='old')

	l=largo(filtfile)
	write (6,'(/x,2a)') 'Reading Filter File: ',filtfile(1:l)
!	read (81,err=1) nf,np,r
!	New format
	read (81,err=1) nf,(ni(i),i=1,nf),(nl(i),i=1,nf),(np(i),i=1,nf),
     &		(fid(i),i=0,nf),ltot,(r(i,1),i=1,ltot),(r(i,2),i=1,ltot)
	close (81)
	iread=1
	write (6,'(i4,a,i3,a,$)') nf,' filters defined, out of ',mxft,' maximum'
!'
	write (6,'(x,a)') '    ...done'
	return
1	stop 'Program exits because of error reading file FILTERBIN.RES'
	end

	FUNCTION COLOR_N(I,J,X,Y,N,Z)

!	Uses number of photons inside filter (Fukugita et al. 1996)

!	Returns color for filter pair (I,J) for s.e.d. in (X,Y) with N
!	data points observed at redshift Z. The zeropoint is not added
!	and must be added in the calling program.

	INCLUDE 'filter.dec'
	dimension a(2),iz(2),x(n),y(n)
	iz(1)=i
	iz(2)=j
	do k=1,2
	a(k)=filter_n(iz(k),x,y,n,z,kerr)
	if (a(k).le.0.) goto 1
	enddo
	color_n=-2.5*alog10(a(1)/a(2))
	return
1	color_n=990.
	return
	end

	FUNCTION ZEROP_N(I,J)

!	Uses number of photons inside filter (Fukugita et al. 1996)

!	Returns zeropoint for filter pair (I,J). The A0 V stellar s.e.d.
!	is read from formatted file whose name is A0VSED the first time
!	the function is used.

	INCLUDE 'filter.dec'
	if (jread.eq.0) call readA0V
	zerop_n=-COLOR_N(i,j,xa0v,ya0v,jread,0.)
!	ireset=1
	return
	end

	SUBROUTINE readA0V

!	Reads A0V sed in file A0VSED

	INCLUDE 'filter.dec'
	character sedfile*80
	data xa0v,ya0v,za0v/ma0v*1.,ma0v*1.,ma0v*1./

!	Get file name from environment variable A0VSED
	if (jread.gt.0) return
	call getenv('A0VSED',sedfile)
	open (81,file=sedfile,form='formatted',status='old')
!	write (6,'(/x,2a)') 'Reading A0V Reference SED: ',sedfile(1:largo(sedfile))
	read (81,'(a)',err=3) sedfile
!	write (6,'(x,a)') sedfile (1:78)
	do n=1,10000
	read (81,*,end=1,err=3) xa0v(n),ya0v(n)
	enddo
1	close (81)
	jread=n-1
!	write (6,'(i5,a,$)') jread,' data points'
!	write (6,'(x,a)') '    ...done'
	return
3	stop 'Program exits because of error reading A0 V s.e.d. file'
	end

	FUNCTION VEGA_0P_N(MF)

!	Uses number of photons inside filter (Fukugita et al. 1996)

!	Returns zero point with respect to VEGA sed
!	Vega sed is first expressed in inits of Lsun/A, as in the *.ised files
!	The A0 V stellar s.e.d. is read from file A0VSED the first time the function is used.

	INCLUDE 'filter.dec'
	real yaux5(1250)
	data kread5/0/,yaux5/1250*0./

!	Vega apparent V magnitude
	vmag=0.03
!	Vega bolometric correction (old value)
!	bc=-0.1164
!	Vega bolometric correction (from Lang, Astrophysical Data, Table 9.3, p. 118)
	bc=-0.25
!	Vega bolometric magnitude
	vegabol=vmag+bc
!	write (6,*) 'vmag,bc,vegabol =', vmag,bc,vegabol

	if (jread.eq.0) call readA0V

	if (kread5.eq.0) then
		kread5=1
!		Compute total flux below Vega sed
		tot=trapz1(xa0v,ya0v,jread)
!		Desired total flux corresponding to mbol = vegabol in units of Lsun
		stot=10.**(-0.4*(vegabol-4.75))
!		Scale sed to desired stot, final sed in units of Lsun/A as *.ised files
		scl=stot/tot
		do i=1,jread
		yaux5(i)=scl*ya0v(i)
!		write (95,*) xa0v(i),yaux5(i)
		enddo
	endif

!	Compute vega zero for filter mf
	vega_0p_n=2.5*alog10(filter_n(mf,xa0v,yaux5,jread,0.,kerr))
!	write (6,*) 'VEGA in filter',mf,vega_0p_n
!	ireset=1
	return
	end

	FUNCTION VEGA_ABMAG(MF)

!	Uses number of photons inside filter (Fukugita et al. 1996)

!	Returns ABmag of star Vega in filter MF
!	The A0 V stellar s.e.d. is read from file A0VSED the first time the function is used.
!	Vega sed is first expressed in inits of Lsun/A, as in the *.ised files

	INCLUDE 'filter.dec'
	real yaux5(1250)
	data kread5/0/,yaux5/1250*0./

!	Vega apparent V magnitude
	vmag=0.03
!	Vega bolometric correction (old value)
!	bc=-0.1164
!	Vega bolometric correction (from Lang, Astrophysical Data, Table 9.3, p. 118)
	bc=-0.25
!	Vega bolometric magnitude
	vegabol=vmag+bc
!	write (6,*) 'vmag,bc,vegabol =', vmag,bc,vegabol

	if (jread.eq.0) call readA0V

	if (kread5.eq.0) then
		kread5=1
!		Compute total flux below Vega sed
		tot=trapz1(xa0v,ya0v,jread)
!		Desired total flux corresponding to mbol = vegabol in units of Lsun
		stot=10.**(-0.4*(vegabol-4.75))
!		Scale sed to desired stot, final sed in units of Lsun/A as *.ised files
		scl=stot/tot
		do i=1,jread
		yaux5(i)=scl*ya0v(i)
!		write (95,*) xa0v(i),yaux5(i)
		enddo
	endif

!	Compute vega zero for filter mf
	vega_vega  = vega_mag(mf,xa0v,yaux5,jread,0.)
	vega_abmag =   ab_mag(mf,xa0v,yaux5,jread,0.)
	write (6,*) 'VEGA in filter',mf,vega_vega,vega_abmag
!	ireset=1
	return
	end

	FUNCTION AB_MAG(MF,X,Y,MP,Z)

!	Compute absolute AB magnitude in filter MF (via S. Charlot)

!	Array declarations
	real x(mp),y(mp)

!	dl = 10 pc in units of Mpc
        dl = 1.e-5

!	This factor is SQRT(4*pi*(3.0856E24)^2/Lsun)
        AB0 = 5. * alog10(1.7684e+08 * dl)

!	Compute AB absolute magnitude
        AB_mag = AB0 - 2.5*alog10(f_mean(mf,x,y,mp,z,kerr)) - 48.6

	return
	end

	FUNCTION VEGA_MAG(MF,X,Y,MP,Z)

!	Compute absolute VEGA magnitude in filter MF

!	Array declarations
	real x(mp),y(mp)
	data v0/0./,vmag/-99./,fv/0./

!	Compute V magnitude Vega zero point
	if (v0 <=0.) then
		v0=vega_0p_n(15)
	endif

!	Compute Vmag only if MF>0
	if (mf > 0) then
		fv = filter_n(15,x,y,mp,0.,kerr)	! flux through filters V (15)
		if (fv > 0.) then
			fv   = 2.5*alog10(fv)		! log V flux
			vmag = v0 - fv			! Vmag for this sed
		else
			vega_mag = -99.
			return
		endif
	endif

!	Compute flux through filters MF
	nf = abs(mf)
	fx = filter_n(nf,x,y,mp,0.,kerr)
!	Compute magnitude
	if (fx >0.) then
!		Compute V-MF zero point
		zp = zerop_n(15,nf)
!		Compute V-MF color
		color = zp - fv + 2.5*alog10(fx)
!		Compute Vega_mag in filter mf
		vega_mag = vmag - color
	else
		vega_mag = 99.99
	endif
	return
	end

	FUNCTION ST_MAG(I,X,Y,N,Z,FLUX)

!	Compute absolute ST magnitude in filter I

!	Source: http://www.stsci.edu/hst/acs/analysis/zeropoints

!	STmag and ABmag: Both systems define the absolute physical flux density for a point source.
!	The conversion is chosen so that the magnitude at V corresponds roughly to that in the Johnson system.
!	In the STmag system, the flux density is expressed per unit wavelength, while
!	in the ABmag system, the flux density is expressed per unit frequency.
!	The definitions are:
!	        STmag = -2.5 Log F_lam -21.10
!	        ABmag = -2.5 Log F_nu - 48.60 
!	where F_nu is expressed in erg cm-2 s-1 Hz-1, and F_lam in erg cm-2 s-1 Ang-1.
!	An object with a constant flux distribution F_nu = 3.63 x 10-20 erg cm-2 s-1 Hz-1 at all wavelengths
!	will have ABmag=0 at all wavelengths, and similarly an object with F_lam = 3.63 x 10-9 erg cm-2 s-1 Ang-1 will have STmag=0.

!	Conversion from BC03 flux units:
!	  The stellar sed''s in BC03 are in units of Lsun/A.
!	  To obtain fluxes in ergs/sec/cm**2/A, multiply by Lsun = 3.826E33 ergs/sec, and divide by
!	  the surface area of a sphere of radius 10 pc =  area(10pc) = 1.19E40 cm2. Then, use the
!	  conversion factor:
!		 f = (3.826/1.19)*10.**(33-40) = 3.215E-07, and
!		fm = -2.5 log10 (f) = 16.232

!	Array declarations
	real x(n),y(n)
	data fm/16.232/

!	Compute effective flux in filter MF
        Flux   = effective_flux(i,x,y,n,z,xeff)
	if (flux > 0.) then
        	ST_mag = fm - 2.5*alog10(Flux) - 21.10
	else
        	ST_mag = -99.
	endif
	return
	end

	FUNCTION AREAFILT(N,XEFF)

!	Returns area below filter number n and effective wavelength

	INCLUDE 'filter.dec'

	real x(naux),y(naux),z(naux)
	if (iread.eq.0) Call FILTER0
	m=0
	do i=ni(n),nl(n)
	m=m+1
	x(m)=r(i,1)
	y(m)=r(i,2)
	z(m)=y(m)*x(m)
	enddo
	areafilt=trapz1(x,y,np(n))
	xeff=trapz1(x,z,np(n))/areafilt
	return
	end

	SUBROUTINE FILTERID(USERID)

!	Returns array with filter identification
	include 'filter.dec'
	character*64 userid(0:mxft)

!	Read filter file
	if (iread.eq.0) Call FILTER0
	do i=0,nf
	userid(i)=fid(i)
	enddo
	return
	end

	FUNCTION IFILTER()
!	Checks which filter file is in use.
	character filtfile*80
	call getenv('FILTERS',filtfile)
	ifilter=index(filtfile,'FILTERBIN.RES')
!	write (6,*) filtfile
	return
	end

	SUBROUTINE YNORM(I,X,Y,N,Z,Y0)

!	Normalizes Y(X) such that the flux seen through the Ith filter at
!	redshift Z is = Y0

	INCLUDE 'filter.dec'
	real x(n),y(n),a(4000)
	data a/4000*1./
	aflux=filter(i,x,y,n,z,kerr)
	area =filter(i,x,a,n,z,kerr)
	ya=aflux/area
!	write (6,*) 'ya = ',ya
	ya=y0/ya
	do j=1,n
	y(j)=ya*y(j)
	enddo
	return
	end

	FUNCTION EFFECTIVE_FLUX(I,X,Y,N,Z,XEFF)

!	Computes effective flux seen through the Ith filter at redshift Z (monochormatic flux)

	real x(n),y(n)
       	aflux=filter(i,x,y,n,z,kerr)
!      	area =areafilt(i,xeff)
       	area =areafilt(i,xeff)/(1.+z)		! Factor (1.+z) added by GBA. May 2016
       	effective_flux=aflux/area
       	return
       	end

	FUNCTION EFFECTIVE_NPHOT(I,X,Y,N,Z,XEFF)

!	Computes effective flux seen through the Ith filter at redshift Z (monochormatic flux)

	real x(n),y(n)
	anphot=filter_n(i,x,y,n,z,kerr)
!	area  =areafilt(i,xeff)
	area  =areafilt(i,xeff)/(1.+z)		! Factor (1.+z) added by GBA. May 2016
	effective_nphot=anphot/area
       	return
	end

	FUNCTION SUN_COLOR_N(I,J)

!	Returns color of solar spectrum for filter pair (I,J)

!	Variables
	parameter (msun=1238)
	character sedfile*80
	real xsun(msun),ysun(msun)
	data xsun,ysun,nsun/msun*0.,msun*0.,0/

!	Read solar spectrum
	if (nsun.eq.0) then
!		Get file name from environment variable A0VSED
		call getenv('SUNSED',sedfile)
		open (81,file=sedfile,form='formatted',status='old')
		write (6,'(/x,2a)') 'Reading solar spectrum: ',sedfile(1:largo(sedfile))
		read (81,'(a)',err=3) sedfile
		write (6,'(x,a)') sedfile (1:78)
		do n=1,10000
		read (81,*,end=1,err=3) xsun(n),ysun(n)
		enddo
1		close (81)
		nsun=n-1
		write (6,'(i5,a,$)') nsun,' data points'
		write (6,'(x,a)') '    ...done'
	endif

!	Compute color in desired bands
	SUN_COLOR_N = ZEROP_N(I,J) + COLOR_N(I,J,XSUN,YSUN,NSUN,0.)

	return
3	stop 'Program exits because of error reading solar spectrum file'
	end
