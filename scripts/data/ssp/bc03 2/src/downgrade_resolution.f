	PROGRAM DOWNGRADE_RESOLUTION

c	Downgrades spectral resolution of _hr_ *.ised files in the optical range

c	Applies either a gaussian smoothing, rebinning or velocity dispersion

c	Array declaration
        parameter (jb=3)
	include 'csp.dec'
        character name*256,save*256,aux*256
	logical stelib
        real t(jts),x(imw),y(imw),w(imw),z(imw),shigh(imw)
        real snr(0:jts),pnr(0:jts),bh(0:jts),sn(0:jts),wd(0:jts),rm(0:jts)
        real str(0:jts),sf(0:jts),evf(0:jts)
!	real bol(0:jts)
	data save/' '/

c	Check if correct filter file is in use.
	j=ifilter()
	if (j.eq.0) then
		write (6,'(2a)') char(7),'Please assign correct filter file'
		write (6,'(2a)') char(7),'Use command: stdfilt'
		stop
	endif

c	Ask for file name
1	call copyright(6)
2	write (6,'(x,3a,$)') 'BC_GALAXEV *.ised in file [',save(1:largo(save)),'] = '
	read (5,'(a)',end=10) save
	if (largo(save).gt.0) name=save
	call chaext(name,'ised',nm)
	open (11,file=name,form='unformatted',status='old',err=9)
c	Read basic parameters from input file
	read (11) ks,(t(i),i=1,ks),ml,mu,iseg,
     &  (xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),
     &  totm,totn,avs,io,tauo,id,tcut,ttt,ttt,ttt,id,id,igw,stelib
	read (11) n,(x(i),i=1,n)

c	Check if _lr_, _hrs_, or _hr_ model
	if (n < 2100 ) then		! BaSeL lr model
!	if (n.lt.6900.and.n.ne.5112) then
!		Wrong kind of model (5112 allows for Miles models)
		write (6,*)
		write (6,*) 'The file you have entered does not correspond to a high resolution model. Please try again.'
		write (6,*)
		close (11)
		goto 2
	endif

c	Choose smoothing option
	write (6,*)
	write (6,'(x,a)') 'Options to downgrade spectral resolution of sed''s in range from 3300 to 9300 A:'
	write (6,'(x,a)')  'Choose: 1 = Rebin sed''s'
	write (6,'(x,a)')  '        2 = Apply gaussian broadening function'
	write (6,'(x,a)')  '        3 = Stellar velocity dispersion'
	write (6,'(x,a)')  '        4 = From IndoUS to STELIB resolution'
3	write (6,'(x,a,$)')'Choice: '
	read (5,'(i10)',err=3,end=10) iop
	write (6,*)
	save=name
	i=index(name,'.')-1
	if (iop.eq.1) then
c		Rebin sed
		write (6,'(x,a)')   'Sample bin widths =  5 A to rebin to the Pickles library wavelength scale'
		write (6,'(x,a)')   '                    20 A to rebin to the BaSeL wavelength scale'
4		write (6,'(x,a,$)') 'Enter bin width   = '
		read (5,*,err=4,end=10) dw
		name=name(1:i) // '_rb_xxxx'
		if (dw.lt.10.) then
			write (name(i+5:),'(i1,a)') int(dw),'A  ' 
		elseif (dw.lt.100.) then
			write (name(i+5:),'(i2,a)') int(dw),'A ' 
		elseif (dw.lt.1000.) then
			write (name(i+5:),'(i3,a)') int(dw),'A' 
		endif
	elseif (iop.eq.2) then
c		Apply gaussian broadening function
		write (6,'(x,a)')   'Sample FWHM = 10 A to approximate the Pickles library spectral resolution'
		write (6,'(x,a)')   '              20 A to approximate the BaSeL spectral resolution'
		write (6,'(x,a)')   '              -1 A to approximate the Lick/IDS spectral resolution (Worthey and Ottaviani 1997)'
5		write (6,'(x,a,$)') 'Enter FWHM  = '
		read (5,*,err=5,end=10) fw
		name=name(1:i) // '_gb_xxxx'
		if (fw.eq.-1.) then
			write (name(i+5:),'(a)') 'lick' 
		elseif (fw.lt.10.) then
			write (name(i+5:),'(i1,a)') int(fw),'A  ' 
		elseif (fw.lt.100.) then
			write (name(i+5:),'(i2,a)') int(fw),'A ' 
		elseif (fw.lt.1000.) then
			write (name(i+5:),'(i3,a)') int(fw),'A' 
		endif
	elseif (iop.eq.3) then
c		Compute effects of stellar velocity dispersion on spectral resolution
6		write (6,'(x,a,$)') 'Enter sigma for velocity dispersion (km/s) = '
		read (5,'(f10.0)',err=6,end=10) sigma
		name=name(1:i) // '_vd_sigma_xxxx'
		if (sigma.ge.1000.) then
			write (name(i+11:),'(i4)') nint(sigma)
		elseif (sigma.ge.100.) then
			write (name(i+11:),'(i3,a)') nint(sigma),' '
		elseif (sigma.ge.10.) then
			write (name(i+11:),'(i2,a)') nint(sigma),'  '
		elseif (sigma.ge.0.) then
			write (name(i+11:),'(i1,a)') nint(sigma),'   '
		else
			write (6,*) 'Invalid value of sigma: ',sigma
			goto 6
		endif
	elseif (iop.eq.4) then
c		Downgrades from IndoUS to Stelib spectral resolution
c		Compute sigma as a function of wavelength
		m=n
		do ir=1,n
		w(ir)=x(ir)
		shigh(ir)=sigma_indous(x(ir))
		enddo
		name=name(1:i) // '_stelib_res'
	else
		write (6,*) 'Wrong option. Please try again.'
		goto 3
	endif

c	Open output files
	write (6,*)
	write (6,'(x,3a,$)') 'Generic output file name [',name(1:largo(name)),'] = '
	read (5,'(a)',end=10) aux
	if (largo(aux).gt.0) name=aux
c	Open ascii files for output
	lun=20
	call open_color_files(lun,name)
c	Open binary file and write basic data
	call chaext(name,'ised',nn)
	open (lun,file=name,form='unformatted',status='unknown')
       	write (lun) ks,(t(i),i=1,ks),ml,mu,iseg,
     &		    (xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),
     &		    totm,totn,avs,io,tauo,id,tcut,ttt,ttt,ttt,id,id,igw,stelib
	igw=0
	stelib=.true.

c	Do required task
	do k=1,ks
	read (11) n,(y(i),i=1,n),inx,(y(i),i=n+1,n+inx)
	if (iop.eq.1) then
		call rebin_sed(x,y,n,dw,w,z,m)
	elseif (iop.eq.2) then
		call gbfun_sed(x,y,n,fw,w,z,m)
	elseif (iop.eq.3) then
		call vdsig_sed(x,y,n,sigma,w,z,m)
	elseif (iop.eq.4) then
c		Transform to Stelib resolution
c		FWHM = 3. A
		d=3.
		c=2.*sqrt(2.*alog(2.))
		slow=d/c
		call degrade_to_constant_resolution(x,y,n,slow,shigh,z)
	endif
	if (k.eq.1) write (lun) m,(w(i),i=1,m),inx,(y(i),i=n+1,n+inx)
	write (lun) m,(z(i),i=1,m),inx,(y(i),i=n+1,n+inx)
        call rf_color(1,t(k),w,z,m,lun,1.E6,str(k),sf(k),evf(k),snr(k),pnr(k),bh(k),sn(k),wd(k),rm(k),toff(k),bolms(k),gasms(k),galms(k))
        if (k.gt.1) call percent(6,k,ks,'DOWNGRADE_RESOLUTION ' // name(1:nargo(name)))
	enddo

c	Add remaining records
	read (11,end=8) m,(bflx(i),i=1,m)
	write (lun) m,(bflx(i),i=1,m)
	read (11,end=8) m,(strm(i),i=1,m)
	write (lun) m,(strm(i),i=1,m)
	read (11,end=8) m,(sf(i),i=1,m)
	write (lun) m,(sf(i),i=1,m)
	read (11,end=8) m,(evfl(i),i=1,m)
	write (lun) m,(evfl(i),i=1,m)
	read (11,end=8) m,(snbr(i),i=1,m)
	write (lun) m,(snbr(i),i=1,m)
	read (11,end=8) m,(pnbr(i),i=1,m)
	write (lun) m,(pnbr(i),i=1,m)
	read (11,end=8) m,(bhtn(i),i=1,m)
	write (lun) m,(bhtn(i),i=1,m)
	read (11,end=8) m,(sntn(i),i=1,m)
	write (lun) m,(sntn(i),i=1,m)
	read (11,end=8) m,(wdtn(i),i=1,m)
	write (lun) m,(wdtn(i),i=1,m)
	read (11,end=8) m,(rmtm(i),i=1,m)
	write (lun) m,(rmtm(i),i=1,m)
	read (11,end=8) m,(toff(i),i=1,m)
	write (lun) m,(toff(i),i=1,m)
	read (11,end=8) m,(bolms(i),i=1,m)
	write (lun) m,(bolms(i),i=1,m)
	read (11,end=8) m,(gasms(i),i=1,m)
	write (lun) m,(gasms(i),i=1,m)
	read (11,end=8) m,(galms(i),i=1,m)
	write (lun) m,(galms(i),i=1,m)
8	close (11)
	close (lun)
c	Delete unwanted files
	call delete_files(name,9000,1)
	i=index(name,'.')-1
	write (6,'(/x,2a)') 'Finished building model: ',name(1:i)
	name=save
	goto 1
9       write (6,'(x,3a)') 'File ',name(1:largo(name)),' not found. Please try again.'
        goto 2
10      end

	SUBROUTINE REBIN_SED(X,Y,N,DW,W,Z,M)

c	Rebin STELIB sed (X,Y,N) in optical range to desired width (DW A)
c	Rebinned sed is returned in (W,Z,M)

c	Array declaration
	real x(n),y(n),w(n),z(n)

c	Rebinning in wavelength range from w0 to near 9290 A
c	Setup parameters
	if (dw.le.5) then
		w0=3325.
		wf=9295.
		wl=3322.
	else
		w0=3330.
		wf=9290.
		wl=3322.
	endif
	if (n > 8000) then		! xmiless hr sed
		w0 = 1000.
		wl = 1000.-dw
	endif

c	Up to wl, keep original sed
	m=0
	do j=1,n
	if (x(j) <= wl) then
		m=m+1
		w(m)=x(j)
		z(m)=y(j)
	endif	
	enddo

c	Build rebinned sed from w0 to 9290 A
	do j=1,n
	if (w0 <= wf) then
		m=m+1
		w(m)=w0
		z(m)=trapzq(x,y,n,w0-dw/2.,w0+dw/2.,ierr)/dw
		w0=w0+dw
	endif
	enddo

c	Above last w0, keep original sed
	do j=1,n
	if (x(j) >= w0) then
		m=m+1
		w(m)=x(j)
		z(m)=y(j)
	endif
	enddo
	return
	end

	SUBROUTINE GBFUN_SED(X,Y,N,FWHMS,W,Z,M)

c	Degrades resolution of sed in (X,Y,N) by applying a gaussian broadening function
c	The FWHM of the gaussian is assumed to be constant (= FWHM) with wavelength.
c	Resulting sed is returned in (W,Z,M)

c	Array declaration
	real x(n),y(n),w(n),z(n)

c	Degrade the sed
	if (fwhms.eq.-1.) then
c		Transform to Lick/IDS spectral resolution
		do i=1,n
		w(i)=x(i)
		enddo
		m=n
		call lick_system(x,y,n,z)
	else
		call degrade_resolution(x,y,n,fwhms,w,z,m)
	endif

!	Check for type of sed
	w0=3322.
	if (n > 8000) w0 = 1000.		! xmiless hr sed

c	Up to w0, keep original sed
	m=0
	do j=1,n
	if (x(j) <= w0) then
		m=m+1
		w(m)=x(j)
		z(m)=y(j)
	endif	
	enddo

c	Above 9290 A keep original sed
	w0=9290.
c	Find position K corresponding to w0 in array X
	call locate(x,n,w0,k)
	m=k+1
	do j=1,n
	if (x(j).gt.w0) then
		m=m+1
		w(m)=x(j)
		z(m)=y(j)
	endif
	enddo
	return
	end

	SUBROUTINE VDSIG_SED(X,Y,N,SIGMA,W,Z,M)

c	Degrades resolution of sed in (X,Y,N) introducing the stellar velocity dispersion
c	Resulting sed is returned in (W,Z,M)

c	Array declaration
	parameter (kw=25000)
	real x(n),y(n),w(kw),z(kw)

c	Compute new sed
	do i=1,n
	w(i)=x(i)
	z(i)=y(i)
	enddo
	m=n
	call gaussian_v_disp(w,z,m,sigma)

!	Check for type of sed
	w0=3330.
	if (n > 8000) w0 = 1000.		! xmiless hr sed

c	Up to w0, keep original sed
	if (x(1) <= w0) then
		m=0
		do j=1,n
		if (x(j).le.w0) then
			m=m+1
			w(m)=x(j)
			z(m)=y(j)
		endif	
		enddo
	endif

c	Above 9290 A keep original sed
	w0=9290.
	if (x(n).gt.w0) then
c		Find position k corresponding to w0 in array X
		call locate(x,n,w0,k)
		m=k+1
		do j=1,n
		if (x(j).gt.w0) then
			m=m+1
			w(m)=x(j)
			z(m)=y(j)
		endif
		enddo
	endif
	return
	end
