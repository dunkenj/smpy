	SUBROUTINE SDSS_COLOR(t,x,y,inw,lun,bolflux)

!	computes magnitude and colors in AB system

!	Array declarations
	parameter (nc=20,lb=50,mxft=300)
	character genfile*96,envfile*96,rfcolorfile*96,fid(0:mxft)*64
	integer no(lb),n1(nc),n2(nc),kerr(nc)
	real x(inw),y(inw),col(nc),fx(lb),ab(0:nc)
	real gmag
	common /kpercent/ ipcall,iplast,iphead,jc
	common /r_filter/ iread,jread,ireset
	data icall/0/,last/0/,nb/0/
	data nb,no,n1,n2,ng,nv,mc,kb,v0/0,lb*0,nc*0,nc*0,4*0,0./
	data kerr/nc*0/
!	Add V and g filters to list of desired filters
	integer glog
	data glog/121/			! SDSS g filter
	common /f1001/ mc,ab

!	Check for zero flux sed
	if (t.eq.0..or.bolflux.le.0.) return

!	Check if number of points has changed
	if (inw.ne.last) then
!		icall=0
!		write (6,*) 'Resetting sdss:',inw,last
!		ireset=1
		last=inw
	endif

	if (icall.eq.0) then
		icall=1

!		Format of file modified 07/01/2003:
!		To avoid recompiling this routine the arrays n1,n2 of
!		up to 12 elements each are now read from a file.
!		Get file name from environment variable RF_COLORS_ARRAYS
		envfile='RF_COLORS_ARRAYS'
		call getenv(envfile,rfcolorfile)
		close (1)
		open (1,file=rfcolorfile,status='old',form='formatted',err=2)
!		write (6,*)
!		write (6,*) 'List of filters in file: ',rfcolorfile(1:largo(rfcolorfile))
!		Skip first lines
		do i=1,100
		read (1,'(a)') genfile
		if (index(genfile,'AB system').gt.0) goto 5
		enddo
!5		write (6,*) genfile(1:largo(genfile))
5		kb=0
!		Read filter pairs
		do i=1,nc
		read (1,'(2i4)',end=10) n1(i),n2(i)
		kb=kb+1
		fx(kb)=n1(i)
		kb=kb+1
		fx(kb)=n2(i)
		enddo
10		mc=i-1
		close (1)
!		Make sure that SDSS g filter = filter glog in filter file, are included
		kb=kb+1
		fx(kb)=glog
!		Sort filters in numerical order
		call sort(kb,fx)
!		Find independent filters in arrays n1 and n2 and store in array no
		nb=1
		no(1)=fx(1)
		do i=2,kb
		if (fx(i).gt.fx(i-1)) then
			nb=nb+1
			no(nb)=nint(fx(i))
		endif
		enddo
!		Write filter ID
		call filterid(fid)
!		write (6,'(x,a)') 'Selected filters:'
!		do i=1,nb
!		write (6,'(i4,'': '',i4,3x,a)') i,no(i),fid(no(i))
!		enddo
!		write (6,*) 'Selected colors:'
!		write (6,'(20i4)') (n1(i),i=1,mc)
!		write (6,'(20i4)') (n2(i),i=1,mc)
!		write (6,*)

!		Fill arrays with filter numbers
		do i=1,mc
		do j=1,nb
		if (n1(i).eq.no(j)) n1(i) = j
		if (n2(i).eq.no(j)) n2(i) = j
		if (no(j).eq.glog)  ng    = j
		enddo
		enddo

	endif

!	Compute flux through each of nb filters
	do i=1,nb
	if (kerr(i).eq.0) then
		fx(i)=f_mean(no(i),x,y,inw,0.,kerr(i))
	else
		fx(i)=0.
	endif
	enddo

!	Compute colors in file RF_COLORS.sdss on  AB  mag system
	do i=1,mc
	col(i)=-2.5*alog10(fx(n1(i))/fx(n2(i)))
	enddo

!	Compute absolute g AB magnitude for a 1 Mo galaxy
!	It is -27.5 magnitudes brighter for a 1E11 Mo galaxy
!	10 pc in units of Mpc
	dl=1.e-5
	gmag=-2.5*alog10(fx(ng))-48.6
	gmag=gmag+(5. * alog10(1.7684e+08 * dl))	! this factor is SQRT(4*pi*(3.0856E24)^2/Lsun)

!	Compute flux inside GALEX filters
	nfuv = 139	! FUV Galex filter
	ffuv = f_mean(nfuv,x,y,inw,0.,kerruv)
	nnuv = 140	! NUV Galex filter
	fnuv = f_mean(nnuv,x,y,inw,0.,kerruv)

!	Compute flux inside a 100 A square filtered centered at 1500 A
	n1500 = 235     ! filter number
	f1500 = f_mean(n1500,x,y,inw,0.,kerruv)

!	Write colors to file
	tl=alog10(t)
!	write (lun+13,101) tl,bolmag,gmag,(col(i),i=1,mc),ffuv,fnuv,f1500
	write (lun+13,101) tl,gmag-col(1),gmag,(gmag-col(i),i=2,mc),ffuv,fnuv,f1500
101	format (f10.6,5f10.4,3x,7f10.4,3x,2f10.4,1p3e13.4)
!	write (1001,102) gmag-col(1),gmag,(gmag-col(i),i=2,mc)
!102	format (20f8.4)

!	Store AB mags
	ab(0) = gmag-col(1)
	ab(1) = gmag
	do i=2,mc
	ab(i) = gmag-col(i)
	enddo

!	Prepare filters for next call to filter_n
!	ireset=1

	return
2	write (6,*) 'Error opening file: ',rfcolorfile(1:largo(rfcolorfile))
	stop
	end
