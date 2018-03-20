	PROGRAM CM_EVOLUTION

c	Computes apparent magnitude vs z, color vs z, and e+k and k-correction
c	in user selected filters for BC GALAXEV evolving model sed''s
c	The user selects the filter(s) to use from the list in filters.log
c	The user enters values for the cosmological parameters

	include 'read_ised.dec'
	include 'filter.dec'
	include 'cosmo.dec'
	include 'SSP_13.dec'

	character file*96,name*96,ext*24,ext1*4,ext2*4
	integer kf(2)
	real yz(imw),xaux(24)
	data lun/26/,file/' '/,kf/0,0/
c	Cosmology:
	data h/70./,omega/0.30/,omega_lambda/0.70/

c	Check if correct filter file is in use.
	j=ifilter()
	if (j.eq.0) then
		write (6,'(2a)') char(7),'Please assign correct filter file'
		stop
	endif

c	Welcome user
	write (6,*)
	write (6,'(x,a)') 'Galaxy Spectral Evolution Library (GALAXEV)'
	write (6,'(x,a)') 'UNIX Version (C) 1995-2015 - G. Bruzual and S. Charlot - All Rights Reserved'

c	Ask for Ho, Omega, and Omega_lambda
	write (6,'(/x,a,f4.0,a,f5.3,a,$)') 'Enter cosmological parameters'
3	write (6,'(x,a,f4.0,a,f5.3,a,f5.3,a,$)')  'Enter Ho [',h,'], Omega [',omega,'], Omega_lambda [',omega_lambda,'] = '
	call nread(xaux,nc,*9,*10)
	if (nc.gt.0.) then
		h=xaux(1)
		omega=xaux(2)
		omega_lambda=xaux(3)
	endif
c	Omega is now entered bu the user
c	omega=1.-omega_lambda

c	Obtain cosmological constant and q
	clambda=cosmol_c(h,omega,omega_lambda,q)

c	Age of the universe and maximum age for galaxies (rounded off)
	tu=t(h,q,0.,clambda)
	tg=int(tu)*1.E9

c	Ask for galaxy age
	write (6,*)
	write (6,'(x,a,f6.2,a  )') 'Age of this universe  = tu =',tu,' Gyr'
2	write (6,'(x,a,f6.2,a,$)') 'Enter age of galaxy today = tg [',tg*1.e-9,' Gyr] = '
	call qread(ttg,nc,*8,*10)
	if (nc.gt.0.) tg=ttg*1.E9
	ttg=tg*1.E-9
	zf=zx(ttg,h,q,clambda)

c	Report cosmology
	write (6,'(/x,a)') 'Cosmological model in use:'
	write (6,100) h,omega,omega_lambda,q,clambda
100	format (' Ho = ',f4.0,3x,'Omega =',f5.2,3x,'Omega_lambda =',f5.2,3x,'qo = ',
     *            f7.4,3x,'Lambda =',1pE10.3)
	write (6,101) tu,ttg,zf
101	format (' tu = ',f5.2,' Gyr',3x,'tg = ',f5.2,' Gyr', 3x,' zf = ',f6.2)

c	Ask for *.ised file name
	write (6,*)
1	write (6,'(x,3a,$)') 'BC_GALAXEV model sed in file [',file(1:largo(file)),'] = '
	read (5,'(a)',end=10) name
	call chaext(name,'ised',nn)
	file=name
c	Read *.ised file
	call read_ised(file,tg,kl,kerr)
	if (kerr.ne.0) goto 1

c	Read filter file
	call filter0
c	Read A0Vsed
	call readA0V

c	Ask for filter number
	write (6,'(/x,a)')   'Enter up to 2 filters separated by commas.'
4	write (6,'( x,a,$)') 'Compute magnitude/color through filter numbers = '
	read (5,'(2i10)',err=7,end=10) kf(1),kf(2)
	if (kf(1).le.0.or.kf(1).gt.nf) then
		write (6,'(x,a,i5,a)') 'Error in filter number',kf(1),'. Please try again'
		goto 4
	else
		imag=1
		write (6,'(/x,a)') 'Selected filters: '
		write (6,'(x,i5,x,a)') kf(1),fid(kf(1))
	endif
	if (kf(2).lt.0.or.kf(2).gt.nf) then
		write (6,'(x,a,i5,a)') 'Error in filter number',kf(2),'. Please try again'
		goto 4
	elseif (kf(2).gt.0) then
		imag=2
		write (6,'(x,i5,x,a)') kf(2),fid(kf(2))
	endif

c	Build output file name and write header into it
	ext1='F000'
	mf=kf(1)
	if (mf.lt.10) then
		write (ext1(4:4),'(i1)') mf
	elseif (mf.lt.100) then
		write (ext1(3:4),'(i2)') mf
	elseif (mf.le.nf) then
		write (ext1(2:4),'(i3)') mf
	endif
	ext='magnitude_' // ext1
	call chaext(name,ext,nn)
	open (lun,file=name,form='formatted',status='unknown')
	call file_header(lun,file,1)
	write (6,'(/x,2a)') 'Output written to file(s): ',name(1:largo(name))
	if (imag.eq.2) then
		mf=kf(2)
		ext2='F000'
		if (mf.lt.10) then
			write (ext2(4:4),'(i1)') mf
		elseif (mf.lt.100) then
			write (ext2(3:4),'(i2)') mf
		elseif (mf.le.nf) then
			write (ext2(2:4),'(i3)') mf
		endif
		ext='magnitude_' // ext2
		call chaext(name,ext,nn)
		write (6,'(28x,a)') name(1:largo(name))
		open (lun+1,file=name,form='formatted',status='unknown')
		call file_header(lun+1,file,1)
		ext='color_' // ext1 // '_' // ext2
		call chaext(name,ext,nn)
		write (6,'(28x,a)') name(1:largo(name))
		open (lun+2,file=name,form='formatted',status='unknown')
		call file_header(lun+2,file,1)
	endif

c	Compute magnitudes
	write (6,'(/x,a)') 'Computing magnitudes/colors...'
	icall=0
	do kz=1,120000
	z=zk(kz,zf)
	call evol_ised(kl,z,h,q,clambda,yz,icall,last)
	if (last.ne.0) goto 10
       	call ke_2filt_correction(kf,w,yz,iw,z,h,q,lun,imag)
	if (z.eq.zf) goto 10
	enddo
	close (lun)
	goto 10
7	write (6,'(x,a)') 'Error in filter number. Please try again'
	goto 4
8	write (6,'(x,a)') 'Error in galaxy age. Please try again'
	goto 2
9	write (6,'(x,a)') 'Error in parameters. Please try again'
	goto 3
10	end
