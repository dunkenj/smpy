	PROGRAM MM_EVOLUTION

c	Computes apparent magnitude vs z through several filters
c	The user selects the filter(s) to use from the list in filters.log
c	The user enters values for the cosmological parameters

	include 'read_ised.dec'
	include 'filter.dec'
	include 'cosmo.dec'
	include 'SSP_13.dec'

c	Variables
	character file*128,name*128,ext*24,argu*1024
	integer kf(50)
	real yz(imw),xaux(24)
	data lun/26/,file/' '/,kf/50*0/
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
	write (6,'(x,a)') 'UNIX Version (C) 1995-2013 - G. Bruzual and S. Charlot - All Rights Reserved'

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

c	Ask for filter numbers
4	write (6,'(/x,a)')   'Enter up to 50 filters separated by commas.'
	write (6,'( x,a,i3,a,$)') 'Compute magnitude through filter numbers (1:',nf,') = '
	read (5,'(a)',end=10) argu
	ngu = nargu(argu)
	read (argu,*,err=4) (kf(j),j=1,ngu)
	do i=1,ngu
	if (kf(i).le.0.or.kf(i).gt.nf) then
		write (6,'(x,a,i5,a)') 'Error in filter number',kf(i),'. Please try again'
		goto 4
	endif
	enddo
	write (6,'(/x,a)') 'Selected filters: '
	do i=1,ngu
	write (6,'(x,i5,x,a)') kf(i),fid(kf(i))
	enddo

c	Build output file name and write header into it
	ext='multi_mag_vega'
	call chaext(name,ext,nn)
	open (lun,file=name,form='formatted',status='unknown')
	call file_header(lun,file,1)
	write (6,'(/x,2a)') 'Output written to file(s): ',name(1:largo(name))
	ext='multi_mag_AB'
	call chaext(name,ext,nn)
	open (lun+1,file=name,form='formatted',status='unknown')
	call file_header(lun+1,file,1)
	write (6,'( x,2a)') '                           ',name(1:largo(name))

c	Compute magnitudes
	write (6,'(/x,a)') 'Computing magnitudes/colors...'
	icall=0
	do kz=1,120000
	z=zk(kz,zf)
	call evol_ised(kl,z,h,q,clambda,yz,icall,last)
	if (last.ne.0) goto 10
       	call ke_nfilt_correction(kf,w,yz,iw,z,h,q,lun,ngu)
	if (z.eq.zf) goto 10
	enddo
	close (lun)
	goto 10
!7	write (6,'(x,a)') 'Error in filter number. Please try again'
!	goto 4
8	write (6,'(x,a)') 'Error in galaxy age. Please try again'
	goto 2
9	write (6,'(x,a)') 'Error in parameters. Please try again'
	goto 3
10	end
