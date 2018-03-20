	PROGRAM GALAXEVPL

!	Writes a formatted file with galaxy spectral energy distributions
!	stored in files *.ised

!	parameter (nw=25000,no=250,js=20000)
	parameter (nw=55000,no=250,js=221)
	character namech*256,ext*12,name2*256,save*256,argu*1024,id*80,xrgu*1024	! ,headcol*12
	integer n(no)
	real c(no),t(js),cnu(nw),aux(no),bux(js),h(nw),w(nw),f(nw,js),p(6)
	data save/' '/,aux/no*0./

!	Get file name from argument list or ask user for file name.
!	Extension "ised" added by code.
	ja = iargc()
	jn = ja
	jm = ja
	jw = ja
1	call copyright(6)
	write (6,*)
	write (6,'(1x,a)') 'GALAXEVPL: This program will write a formatted file with galaxy s.e.d.''s'
	write (6,'(1x,2a)') 'read from an unformatted *.ised file.'
11	if (jn.gt.0) then
		call getarg(1,namech)
	else
		write (6,'(/1x,3a,$)') 'Input file name [',save(1:margo(save)),'] (CTRL/D exits) = '
		read (5,'(a)',end=10) namech
		if (largo(namech).eq.0) namech=save
	endif
!	call chaext(namech,'ised',mmm)
	call chext(namech,'ised')
	save=namech

!	Open and read file
	write (6,'(/1x,2a)') 'Reading file: ',namech(1:largo(namech))
	close (1)
	open (1,file=namech,status='old',form='unformatted',err=25)
!	read (1) ks,(t(i),i=1,ks)
	read (1) ks,(t(i),i=1,ks),xml,xmu,iseg,(p,i=1,iseg),totm,totn,avs,io,tau,id
!	write (6,*) xml,xmu,iseg,totm,totn,avs,io,tau
	write (6,*) 'Model ID: ',id
!	write (lun) nsteps+1,0.,(tb(i),i=1,nsteps),ml,mu,iseg,(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),totm,totn,avs,io,tau,id,tau,tau,1.,1.,id,id,-2,stelib
!	write (6,*) ks,(t(i),i=1,ks)
	read (1) iw,(w(i),i=1,iw)
	write (6,100)
     &  'In this file there are',ks,' galaxy s.e.d.''s, from ',t(1),' to ',t(ks),' yr'
100	format (/1x,a,i7,a,1pe9.2,a,e9.2,a)
	write (6,102) w(1),w(iw),iw
102	format (/' Each sed covers from lambda = ',F6.1,' to ',1pe12.5,' A in ',I6,' steps.')

!	Get age of sed''s to be listed from argument list or ask user
2	if (jm.gt.1) then
		call getarg(2,argu)
	else
		write (6,'(/1x,a,$)') 'Enter age (in Gyr) of up to 250 sed''s (separated by commas) = '
		read (5,'(a)',end=10) argu
	endif
!	Found a problem with gfortran in reading more arguments that have been actually entered.
!	Use function nargu(argu) to count how many time steps have been entered
	if (argu == '-all' .or. argu == '-vall') then
		ngu = min(ks,no)
		do j=1,ngu
!		aux(j)=t(j)*1.E-09
		aux(j)=j
		enddo
		aux(1)=-1
	else
		ngu = nargu(argu)
		read (argu,*,err=27) (aux(j),j=1,ngu)
	endif

!	Find sed''s in model file with age closest to values requested
	ic=0
	if (aux(1).lt.0) then
		aux(1)=-aux(1)
		tfac=1.
	else
		tfac=1.E+09
	endif
	do i=1,ngu
	if (aux(i).eq.0.) then
		ic=ic+1
		aux(ic)=t(1)
		if (tfac.eq.1) aux(ic)=1
	elseif (aux(i).gt.0.) then
		ic=ic+1
		aux(ic)=aux(i)*tfac
	endif
	enddo
!	Sort array aux and suppress duplicate entries
        call usort(ic,aux)
!	Find corresponding record numbers
	do i=1,ic
	if (tfac.gt.1.) then
		do k=1,ks
		bux(k)=abs(aux(i)-t(k))
		enddo
		closest=rminim(bux,1,ks,kmin)
	else
		kmin=nint(aux(i))
	endif
	n(i)=kmin
	enddo

!	Get wavelength range to copy to ascii file from argument list or ask user
!	ip = Number of flux values to be printed
3	ip = iw
!	if (argu == '-all') then
!		w1=0.
!		w2=0.
!		w0=0.
!		f0=0.
!		z=0.
!	elseif (jw.gt.2) then
	if (jw > 2) then
		call getarg(3,argu)
!		read (argu,'(5f10.0)',err=26) w1,w2,w0,f0,z
		read (argu,*,err=26) w1,w2	!,w0,f0,z
	else
		write (6,*)
		write (6,*) 'Enter desired wavelength range = [W1,W2] (default: full range).'
		write (6,*) '     o If you want all s.e.d.''s scaled to flux = F0 at lambda = W0, enter the'
		write (6,*) '          desired values (default: no scaling)'
		write (6,*) '     o If you want the output as Fnu vs. lambda, enter W1 with a minus sign.'
		write (6,*) '     o If you want all s.e.d.''s scaled so that F0 is the flux measured through'
		write (6,*) '          filter N at redshift Z, enter W0 = -N'
		write (6,*)
		write (6,'(1x,a,$)') 'W1,W2,W0,F0,Z = '
		read (5,'(a)',end=10) xrgu
		w1=0.
		w2=0.
		w0=0.
		f0=0.
		z=0.
		if (largo(xrgu) > 0) then
			ngu = nargu(xrgu)
			if (ngu.eq.1) read (xrgu,*,err=26) w1
			if (ngu.eq.2) read (xrgu,*,err=26) w1,w2
			if (ngu.eq.3) read (xrgu,*,err=26) w1,w2,w0
			if (ngu.eq.4) read (xrgu,*,err=26) w1,w2,w0,f0
			if (ngu.eq.5) read (xrgu,*,err=26) w1,w2,w0,f0,z
		endif
!		read (5,'(5e20.0)',err=26,end=10) w1,w2,w0,f0,z
	endif

!	Perform required action in wavelength range and model flux
	i1=1
	i2=ip
	i0=0
	nu=0
	j0=0
	if (w0 < 0.) then
		j0=int(abs(w0))
		w0=0.
		g0=f0
		f0=0.
	endif
	if (w1.lt.0) nu=1
	w1=abs(w1)
	do i=1,ip
	cnu(i)=1.
	if (nu.ne.0) cnu(i)=w(i)*w(i)*1.E-8/3.E10
	if (w(i) <= w0) i0=i
	if (w(i) <= w1) i1=i
!	if (w(i) <= w2) i2=i+1
	if (w(i) <= w2) i2=i
	enddo
	if (i1 >= i2) goto 26
	do j=1,ks
	read (1) ip,(h(i),i=1,ip)
	if (j0 > 0) call ynorm(j0,w,h,ip,z,g0)
	do i=1,ip
	f(i,j)=h(i)
	enddo
	enddo
!	If monochromatic normalization requested
	do i=1,ic
	if (f0.gt.0.) then
		c(i)=f0/f(i0,n(i))/cnu(i0)
	else
		c(i)=1.
	endif
	enddo

!	Get output file name from argument list or ask user
	if (ja > 3) then
		call getarg(4,xrgu)
		if (index(xrgu,'default')>0) then
			ja=3
			argu='default'
		endif
	endif
	if (ja > 3) then
		call getarg(4,namech)
	else
		if (n(1).lt.10) then
			write (ext,'(i1)') n(1)
		elseif (n(1).lt.100) then
			write (ext,'(i2)') n(1)
		elseif (n(1).lt.1000) then
			write (ext,'(i3)') n(1)
		elseif (n(1).lt.10000) then
			write (ext,'(i4)') n(1)
		else
			ext='glxvpl'
		endif
		if (argu == 'default') then
			call chext(namech,ext)
		elseif (argu == '-all') then
			if (ja > 2) then
				call getarg(3,namech)
			else
				call chext(namech,'all')
			endif
		elseif (argu == '-vall') then
			if (ja > 2) then
				call getarg(3,namech)
			else
				call chext(namech,'vall')
			endif
		else
			call chext(namech,ext)
			write (6,'(/1x,3a,$)') 'Default output file name [',trim(namech),'], or enter your own = '
			read (5,'(a)',end=10) name2
			if (largo(name2) > 0) namech=name2
		endif
	endif

!	Open file and write requested columns
	close (2)
	open (2,file=namech,form='formatted',status='unknown')
	if (argu == '-vall') then
		write (6,'(x,a)') 'List next to the flux the value from column J of file I:'
		write (6,'(x,a)') '            1 = *.1color'
		write (6,'(x,a)') '            2 = *.2color'
		write (6,'(x,a)') '            3 = *.3color'
		write (6,'(x,a)') '            4 = *.4color'
		write (6,'(x,a)') '            5 = *.5color'
		write (6,'(x,a)') '            6 = *.9color'
		write (6,'(x,a)') '            7 = *.1ABmag'
4		write (6,'(x,a,$)') 'Enter I,J = '
		read (5,'(2i10)',err=4,end=10) i,j
		if     (i == 1) then
			ext = '1color'
		elseif (i == 2) then
			ext = '2color'
		elseif (i == 3) then
			ext = '3color'
		elseif (i == 4) then
			ext = '4color'
		elseif (i == 5) then
			ext = '5color'
		elseif (i == 6) then
			ext = '9color'
		elseif (i == 7) then
			ext = '1ABmag'
		else
			write (6,*) 'Option I =',i,' not implemented'
			stop
		endif
		name2=namech
		call chext(name2,ext)
		open(666,file=name2)
!		Read header
		do l=1,30
		read (666,*)
		enddo
		do l=1,220
		read (666,*) (aux(l),k=1,j),bux(l)
		enddo
		close (666)
		write (6,*)
		write (6,105) 'Output file name = ',trim(namech),'Input file name  = ',trim(save)
		write (6,1051) '3rd column lists = ',trim(name2),', Columns =',j,j+1
		write (2,105) 'Output file name = ',trim(namech),'Input file name  = ',trim(save)
		write (2,1051) '3rd column lists = ',trim(name2),', Columns =',j,j+1
		write (2,'(a)') '#'
!		write (2,'(a)') '# Lambda(A)  Flux        Prop(J)    Prop(J+1)  Age(yr)   Rec'
		write (2,'(a)') '#  Age(yr)   Lambda(A)   Flux       Prop(J)    Prop(J+1)   Rec'
!		Add column for t=0, the same as for the first record
		do j=1,ic
		do i=i1,i2
		if (j == 1) then
			write (2,1091) t(n(1)),w(i),c(1)*cnu(i)*f(i,n(1)),aux(1),bux(1),j-1
		else
			write (2,1091) t(n(j)),w(i),c(j)*cnu(i)*f(i,n(j)),aux(j-1),bux(j-1),j-1
		endif
		enddo
		enddo
	else
		write (6,*)
		write (2,105) 'Input file name  = ',trim(save),'Output file name = ',trim(namech)
!		write (2,105) 'Output file name = ',trim(namech),'Input file name  = ',trim(save)
		write (2,107) (  n(j) ,j=1,ic)
!		write (2,108) (t(n(j)),j=1,ic)
		write (2,106) (   j+1 ,j=1,ic)
		write (2,1081)('      at    ',j=1,ic)
		write (2,1082) (t(n(j)),j=1,ic)
!		write (2,1083) (headcol(t(n(j))),j=1,ic)
		write (6,*)
		write (6,105) 'Input file name  = ',trim(save),'Output file name = ',trim(namech)
!		write (6,105) 'Output file name = ',trim(namech),'Input file name  = ',trim(save)
!		write (6,107) (  n(j) ,j=1,ic)
!!		write (6,108) (t(n(j)),j=1,ic)
!		write (6,106) (   j+1 ,j=1,ic)
!		write (6,1081)('      at    ',j=1,ic)
!		write (6,1082) (t(n(j))*1.E-6,j=1,ic)
!		write (6,1082) (t(n(j)),j=1,ic)
		write (6,'(/1x,i4,2a)') ic,' columns written to file: ',trim(namech)
		do i=i1,i2
		write (2,109) w(i),(c(j)*cnu(j)*f(i,n(j)),j=1,ic)
		enddo
	endif
105	format ('# ',2a)
1051	format ('# ',3a,2i4)
106	format ('# Column',250I12)
107	format ('# Record',250I12)
!108	format ('# Age(yr)     ',1p250E12.3)
1081	format ('# SED(Lo/A) ',250a)
1082	format ('# Lambda(A) ',1p250E12.3)
!1083	format ('# Wavelength',250a)
109	format (1PE12.6,1P250E12.4)
1091	format (1P3E11.4,2e12.4,i4)
	write (6,'(/1x,a)') '-------------------------------------------------------------------'
	if (iargc().gt.0) stop
	j=0
	goto 1

!	Error messages
25	write (6,'(1x,3a/)') 'File ',namech(1:largo(namech)),' not found. Please try again.'
	jn=0
	goto 11
26	write (6,'(1x,a/)') 'Error entering wavelengths. Please try again.'
	jw=0
	goto 3
27	write (6,'(1x,a/)') 'Error entering age of model. Please try again.'
	jm=0
	goto 2
10	end

	CHARACTER*12 FUNCTION HEADCOL(T)
!	Writes age in suitabel format for fits file column header
	headcol=''
	if (t == 0) then
		write (headcol,100) 0.,0
	else
!		Get exponent
		l=int(alog10(t))
!		Process age
		tx = t*10.**(-l)
		if (l<10) then
			write (headcol,100) tx,l
		else
			write (headcol,101) 10.*tx,l-1
		endif
	endif
100	format(2x,'t',f5.3,'E',i1,'yr')
101	format(2x,'t',f5.2,'E',i1,'yr')
	i = index(headcol,'.')
	headcol(i:i) = 'p'
	return
	end
