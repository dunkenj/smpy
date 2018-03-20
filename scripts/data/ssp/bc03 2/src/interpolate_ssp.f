	SUBROUTINE INTERPOLATE_SSP(NAME)

!	Interpolates BC03 model SSP at metallicity Z/Zo indicated in name (after character ',')

!	Variables
	parameter (jb=10)
	parameter (mmod=16)
	include 'csp.dec'
	integer bcopt
	real w(imw),h(imw)
        logical stelib
	character*(*) name,files(mmod)*128
	real zmod(mmod),zx(mmod)
!	data zmod/0.0001, 0.0004,  0.004, 0.008,  0.02,  0.05,  0.10/   						! BC03
!	data zmod/0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.010, 0.014, 0.017, 0.02, 0.03, 0.04/	! CB13

!	Check if interpolation has been requested
	i = index(name,',')
	call s500('c','',0.)
	if (i <= 0) then
!		write (501,*) 0.
!		write (501,*) 1.
!		write (501,*) name(:largo(name))
!		write (501,*) 0.
!		write (501,*) 'empty'
		call s500('f','',0.)
		call s500('f','',1.)
		call s500('a',name(:largo(name)),0.)
		call s500('f','',0.)
		call s500('a','empty',0.)
		return
	endif
	name = name(i+1:)
	ngu = nargu(name)
	if (ngu == 1) then
                read (name,*) z
		idth=0
	elseif (ngu == 2) then
                read (name,*) z,idth
	endif
	call s500('f','',z)

!	Check models in use. Define file names
	bcopt = 0
	call system('/bin/ls *ssp.ised > fort.202')
	read (202,'(a)') files(1)
	read (files(1),'(2x,i4)') bcopt
	close (202)
	call system('/bin/rm fort.202')
	if (bcopt == 2003) then
!		Use BC 2003 models
		nmod      = 6
		zsun      = 0.02
		zmod( 1)  = 0.0001
		zmod( 2)  = 0.0004
		zmod( 3)  = 0.004
		zmod( 4)  = 0.008
		zmod( 5)  = 0.02
		zmod( 6)  = 0.05
!		zmod( 7)  = 0.10
		if (index(files(1),'_stelib_') > 1) then
			files( 1) = 'bc2003_z0001_chab_hr_stelib_ssp.ised'	! = bc2003_hr_stelib_m22_chab_ssp.ised = bc2003_hr_m22_chab_ssp.ised
			files( 2) = 'bc2003_z0004_chab_hr_stelib_ssp.ised'	! = bc2003_hr_stelib_m32_chab_ssp.ised = bc2003_hr_m32_chab_ssp.ised
			files( 3) = 'bc2003_z004_chab_hr_stelib_ssp.ised'	! = bc2003_hr_stelib_m42_chab_ssp.ised = bc2003_hr_m42_chab_ssp.ised
			files( 4) = 'bc2003_z008_chab_hr_stelib_ssp.ised'	! = bc2003_hr_stelib_m52_chab_ssp.ised = bc2003_hr_m52_chab_ssp.ised
			files( 5) = 'bc2003_z020_chab_hr_stelib_ssp.ised'	! = bc2003_hr_stelib_m62_chab_ssp.ised = bc2003_hr_m62_chab_ssp.ised
			files( 6) = 'bc2003_z050_chab_hr_stelib_ssp.ised'	! = bc2003_hr_stelib_m72_chab_ssp.ised = bc2003_hr_m72_chab_ssp.ised
!			files( 7) = 'bc2003_z100_chab_hr_stelib_ssp.ised'	! = bc2003_hr_stelib_m82_chab_ssp.ised = bc2003_hr_m82_chab_ssp.ised
			bcopt     = 1
		elseif (index(files(1),'_xmiless_') > 1) then
			files( 1) = 'bc2003_z0001_chab_hr_xmiless_ssp.ised'	! = bc2003_hr_xmiless_m22_chab_ssp.ised
			files( 2) = 'bc2003_z0004_chab_hr_xmiless_ssp.ised'	! = bc2003_hr_xmiless_m32_chab_ssp.ised
			files( 3) = 'bc2003_z004_chab_hr_xmiless_ssp.ised'	! = bc2003_hr_xmiless_m42_chab_ssp.ised
			files( 4) = 'bc2003_z008_chab_hr_xmiless_ssp.ised'	! = bc2003_hr_xmiless_m52_chab_ssp.ised
			files( 5) = 'bc2003_z020_chab_hr_xmiless_ssp.ised'	! = bc2003_hr_xmiless_m62_chab_ssp.ised
			files( 6) = 'bc2003_z050_chab_hr_xmiless_ssp.ised'	! = bc2003_hr_xmiless_m72_chab_ssp.ised
!			files( 7) = 'bc2003_z100_chab_hr_xmiless_ssp.ised'	! = bc2003_hr_xmiless_m82_chab_ssp.ised
			bcopt     = 2
		else
			write (6,*) 'Unknown series of models: ',files(1)(:largo(files(1)))
			stop
		endif
	elseif (bcopt == 2013 .or. bcopt == 2016) then
!		Use CB 2013 models
		nmod      = 16
		zsun      = 0.017
		zmod( 1)  = 0.0000
		zmod( 2)  = 0.0001
		zmod( 3)  = 0.0002
		zmod( 4)  = 0.0005
		zmod( 5)  = 0.001
		zmod( 6)  = 0.002
		zmod( 7)  = 0.004
		zmod( 8)  = 0.006
		zmod( 9)  = 0.008
		zmod(10)  = 0.010
		zmod(11)  = 0.014
		zmod(12)  = 0.017
		zmod(13)  = 0.020
		zmod(14)  = 0.030
		zmod(15)  = 0.040
		zmod(16)  = 0.060
		if (index(files(1),'_stelib_') > 1) then
			files( 1) = 'cb2016_z0000_chab_hr_stelib_ssp.ised'	! = cb2013_s2_z0000y230n_chab_hr_stelib_ssp.ised
			files( 2) = 'cb2016_z0001_chab_hr_stelib_ssp.ised'	! = cb2013_s2_z0001y249n_chab_hr_stelib_ssp.ised
			files( 3) = 'cb2016_z0002_chab_hr_stelib_ssp.ised'	! = cb2013_s2_z0002y249n_chab_hr_stelib_ssp.ised
			files( 4) = 'cb2016_z0005_chab_hr_stelib_ssp.ised'	! = cb2013_s2_z0005y249n_chab_hr_stelib_ssp.ised
			files( 5) = 'cb2016_z001_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z001y250n_chab_hr_stelib_ssp.ised
			files( 6) = 'cb2016_z002_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z002y252n_chab_hr_stelib_ssp.ised
			files( 7) = 'cb2016_z004_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z004y256n_chab_hr_stelib_ssp.ised
			files( 8) = 'cb2016_z006_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z006y259n_chab_hr_stelib_ssp.ised
			files( 9) = 'cb2016_z008_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z008y263n_chab_hr_stelib_ssp.ised
			files(10) = 'cb2016_z010_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z010y267n_chab_hr_stelib_ssp.ised
			files(11) = 'cb2016_z014_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z014y273n_chab_hr_stelib_ssp.ised
			files(12) = 'cb2016_z017_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z017y279n_chab_hr_stelib_ssp.ised
			files(13) = 'cb2016_z020_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z020y284n_chab_hr_stelib_ssp.ised
			files(14) = 'cb2016_z030_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z030y302n_chab_hr_stelib_ssp.ised
			files(15) = 'cb2016_z040_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z040y321n_chab_hr_stelib_ssp.ised
			files(16) = 'cb2016_z060_chab_hr_stelib_ssp.ised' 	! = cb2013_s2_z060y356n_chab_hr_stelib_ssp.ised
			bcopt     = 3
		elseif (index(files(1),'_xmiless_') > 1) then
			files( 1) = 'cb2016_z0000_chab_hr_xmiless_ssp.ised'	! = cb2013_s2_z0000y230n_chab_hr_xmiless_ssp.ised
			files( 2) = 'cb2016_z0001_chab_hr_xmiless_ssp.ised'	! = cb2013_s2_z0001y249n_chab_hr_xmiless_ssp.ised
			files( 3) = 'cb2016_z0002_chab_hr_xmiless_ssp.ised'	! = cb2013_s2_z0002y249n_chab_hr_xmiless_ssp.ised
			files( 4) = 'cb2016_z0005_chab_hr_xmiless_ssp.ised'	! = cb2013_s2_z0005y249n_chab_hr_xmiless_ssp.ised
			files( 5) = 'cb2016_z001_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z001y250n_chab_hr_xmiless_ssp.ised
			files( 6) = 'cb2016_z002_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z002y252n_chab_hr_xmiless_ssp.ised
			files( 7) = 'cb2016_z004_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z004y256n_chab_hr_xmiless_ssp.ised
			files( 8) = 'cb2016_z006_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z006y259n_chab_hr_xmiless_ssp.ised
			files( 9) = 'cb2016_z008_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z008y263n_chab_hr_xmiless_ssp.ised
			files(10) = 'cb2016_z010_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z010y267n_chab_hr_xmiless_ssp.ised
			files(11) = 'cb2016_z014_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z014y273n_chab_hr_xmiless_ssp.ised
			files(12) = 'cb2016_z017_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z017y279n_chab_hr_xmiless_ssp.ised
			files(13) = 'cb2016_z020_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z020y284n_chab_hr_xmiless_ssp.ised
			files(14) = 'cb2016_z030_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z030y302n_chab_hr_xmiless_ssp.ised
			files(15) = 'cb2016_z040_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z040y321n_chab_hr_xmiless_ssp.ised
			files(16) = 'cb2016_z060_chab_hr_xmiless_ssp.ised' 	! = cb2013_s2_z060y356n_chab_hr_xmiless_ssp.ised
			bcopt     = 4
		else
			write (6,*) 'Unknown series of models: ',files(1)(:largo(files(1)))
			stop
		endif
	else
		write (6,*) 'Unknown model set:',bcopt
		stop
	endif

!	Express Z in solar units
	do i=1,nmod
	zx(i)=zmod(i)/zsun
	if (z == zx(i)) then
		name = files(i)
!		write (501,*) 1.
!		write (501,*) name(:largo(name))
!		write (501,*) 0.
!		write (501,*) 'empty'
		call s500('f','',1.)
		call s500('a',name(:largo(name)),0.)
		call s500('f','',0.)
		call s500('a','empty',0.)
		return
	endif
	enddo

!	Check extreme values
	if ( z <= zx(1) ) then
		name = files(1)
!		write (501,*) 1.
!		write (501,*) name(:largo(name))
!		write (501,*) 0.
!		write (501,*) 'empty'
		call s500('f','',1.)
		call s500('a',name(:largo(name)),0.)
		call s500('f','',0.)
		call s500('a','empty',0.)
		return
	elseif ( z >= zx(nmod) ) then
		name = files(nmod)
!		write (501,*) 1.
!		write (501,*) name(:largo(name))
!		write (501,*) 0.
!		write (501,*) 'empty'
		call s500('f','',1.)
		call s500('a',name(:largo(name)),0.)
		call s500('f','',0.)
		call s500('a','empty',0.)
		return
	endif

!	Find bracketing Z's
	do i=1,nmod-1
	if ( z > zx(i) .and. z <= zx(i+1) ) then
		write (6,*) 'Using z:',i,z,zx(i),zx(i+1)
		z1 = alog10(zx(i))
		z2 = alog10(zx(i+1))
		zo = alog10(z)
		w1 = (z2-zo)/(z2-z1)
		w2 = 1. - w1
		i1 = i
		i2 = i+1
!		write (501,*) w1
!		write (501,*) files(i1)(:largo(files(i1)))
!		write (501,*) w2
!		write (501,*) files(i2)(:largo(files(i2)))
		call s500('f','',w1)
		call s500('a',files(i1)(:largo(files(i1))),0.)
		call s500('f','',w2)
		call s500('a',files(i2)(:largo(files(i2))),0.)
	endif
	enddo
!	Build filename of temporary .ised file
	if (bcopt == 1) then
		name = 'bc2003_hr_stelib_mxx_chab_ssp.ised'
	elseif (bcopt == 2) then
		name = 'bc2003_hr_xmiless_mxx_chab_ssp.ised'
	elseif (bcopt == 3) then
!		name = 'cb2013_s2_mxx_chab_hr_stelib_ssp.ised'
		name = 'cb2016_mxx_chab_hr_stelib_ssp.ised'
	elseif (bcopt == 4) then
!		name = 'cb2013_s2_mxx_chab_hr_xmiless_ssp.ised'
		name = 'cb2016_mxx_chab_hr_xmiless_ssp.ised'
	else
		write (6,*) 'Stopping on bcopt =',bcopt
	endif
	call gen_file_name(name,z,idth)
	open  (3,file=name,form='unformatted',status='unknown')

!	Perform interpolation
	write (6,*) 'Reading file: ',files(i1)
	open  (1,file=files(i1),form='unformatted',status='old')
!	Read basic parameters from SSP file
	read  (1) nsteps,(tb(i),i=0,nsteps-1),ml,mu,iseg,(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),totm,totn,avs,jo,tauo,id,tcut,ttt,ttt,ttt,id,id,igw,stelib
!	Modify id record
	i=index(id,'X=')
	write (id(i:),'(a,f6.3)') 'SSP interpolated to Z/Zo = ',z
	write (3) nsteps,(tb(i),i=0,nsteps-1),ml,mu,iseg,(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),totm,totn,avs,jo,tauo,id,tcut,ttt,ttt,ttt,id,id,igw,stelib
	read  (1) inl,(w(i),i=1,inl)
	write (3) inl,(w(i),i=1,inl)
	write (6,*) 'Reading file: ',files(i2)
	write (6,*) 'Weights =',w1,w2
	write (6,*) 'Writing file: ',name(:largo(name))
	open  (2,file=files(i2),form='unformatted',status='old')
!	Read basic parameters from SSP file
	read  (2)
	read  (2)
	do n=0,nsteps-1
	read  (1) inl,(w(i),i=1,inl)
	read  (2) inl,(h(i),i=1,inl)
	write (3) inl,(w1*w(i)+w2*h(i),i=1,inl)
	enddo
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
	read  (1,end=1) nsteps,(w(i),i=1,nsteps)
	read  (2,end=1) nsteps,(h(i),i=1,nsteps)
	write (3)       nsteps,(w1*w(i)+w2*h(i),i=1,nsteps)
1	close (1)
	close (2)
	close (3)
	return
	end

	subroutine gen_file_name(name,z,idth)

!	Generates unique and new output file name 

	character*(*) name,aux*6,bux*128

	i = index(name,'mxx')+3
	write (bux,'(a,i2,a)') '_t',idth,'_' // name(i+1:largo(name))
	if (bux(3:3) == ' ') bux(3:3) = '0'
	do j=0,10000
	write (aux,'(a,i5)') 'z',int(z*10**4)+j
	do k=1,6
	if (aux(k:k) == ' ') aux(k:k) ='0'
	enddo
	name = name(:i) // aux // bux
	open(99,file=name,status='new',iostat=iopenstatus)
	if (iopenstatus == 0) then
		close (99)
		return
	endif
	enddo
	return
	end
