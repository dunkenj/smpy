        FUNCTION USRSFR(TP)
 
c	Returns SFR computed according to usr defined SFR(t) read from file.
 
c	Array declarations
	character namex*256,buf*128
	include 'jb.dec'
	include 'csp.dec'
	real xbuf(10)
	data isfrc/0/
 
c	Read file if first time that function is called
	if (tp < 0.) isfrc=0
	if (isfrc.eq.0) then
		if (tp == -1.) then
                	read (5,'(a)',end=10) namex
		else
!			Read sfr from file fort.500 in the Chen et al. option
!			namex = 'fort.500'
			namex = 'internal buffer'
			isfrc=0
			do i=8,nsfrp
			call r500(i,'a',buf,zap)
			if (i<0) goto 1
			if (buf(1:1) .ne. '#') then
				isfrc = isfrc + 1
                		read (buf,*,err=11) time(isfrc),usr_sfr(isfrc)
			endif
			enddo
		endif
                close (1)
                open (1,file=namex,status='old',err=10)
c		Check type of file
		if (index(namex,'4color').gt.0) then
c			Read columns 1 and 10 of *.4color file
			do i=1,30
			read (1,*)
			enddo
			do i=1,nsfrp
			read (1,*,end=1,err=11) xbuf
			time(i)=10.**xbuf(1)
			usr_sfr(i)=xbuf(10)
			isfrc=i
			enddo
		else
			isfrc=0
                	do i=1,nsfrp
			read (1,'(a)',end=1) buf
			if (buf(1:1) .ne. '#') then
!                		read (1,*,end=1,err=11) time(i),usr_sfr(i)
				isfrc = isfrc + 1
                		read (buf,*,err=11) time(isfrc),usr_sfr(isfrc)
			endif
	                enddo
		endif
1		close (1)
		write (6,*) isfrc,' data points read from file ',namex(1:largo(namex))
 
!		Add, if needed, a first point at age = 0, and a second point just before first usr point.
		if (time(1) > 0.) then
			do i=1,isfrc
			time(isfrc+3-i)    = time(isfrc+1-i)
			usr_sfr(isfrc+3-i) = usr_sfr(isfrc+1-i)
			enddo
			time(1) = 0.
			time(2) = 0.99995*time(3)
			usr_sfr(1)=0.
			usr_sfr(2)=0.
			isfrc=isfrc+2
		endif
		tcut=time(isfrc)
		write (6,*) isfrc,' data points in user defined SFR'
		write (6,*) 'Tcut = ',tcut

!		Check if tcut < 20 Gyr and add extra points
!		if (tp == -1. .and. tcut < tb(nsteps-1)) then
!			isfrc=isfrc+1
!			time(isfrc)=10.**(alog10(tcut)+0.00001)
!			usr_sfr(isfrc)=0.
!			write (6,*) isfrc,time(isfrc),usr_sfr(isfrc)
!		endif

!		Double check for repeated time steps
		do i=3,isfrc
		if (time(i) <= time(i-1)) then
			time(i) = 10.**(alog10(time(i-1))+0.00001)
		endif
		enddo

!		Add extra points to array tb
		if (isfrc <= jts-its) then
!			Add all points in array usr_sfr
			call add_points
		else
!			Find singular points in array usr_sfr
			call singular_points
		endif
!		write (6,*) isfrc,time(isfrc)

!		Total stellar mass
		tstelm=trapz1(time,usr_sfr,isfrc)
		write (6,*) 'Total mass formed in stars (Mo) =',tstelm
!		do i=1,isfrc
!		tstelm=trapz1(time,usr_sfr,i)
!		write (900,*) i,time(i),usr_sfr(i),tstelm
!		enddo
	endif
 
!	Interpolate SFR at time TP
	if (tp.lt.0.) then
		usrsfr=0.
	elseif (tp.eq.0.) then
		usrsfr=usr_sfr(1)
	elseif (tp.gt.time(isfrc)) then
		usrsfr=0.
	else
		call locate (time,isfrc,tp,i1)
		a1=(tp-time(i1))/(time(i1+1)-time(i1))
		usrsfr=a1*usr_sfr(i1+1)+(1.-a1)*usr_sfr(i1)
!		write (6,*) i1,tp,time(i1),time(i1+1),a1,usrsfr
	endif
        return
10	write (6,'(1x,2a)') 'File not found: ',namex(1:largo(namex))
	usrsfr=-100.
	return
11	write (6,'(1x,2a)') 'Error reading file: ',namex(1:largo(namex))
	usrsfr=-200.
	return
	end
 
	SUBROUTINE SINGULAR_POINTS
 
!	Finds singular points (peaks, maximum) in array usr_sfr(time)
 
!	Array declarations
	include 'jb.dec'
	include 'csp.dec'
	integer is(nsfrp)
 
!	Find maximum in SFR and points 1.5 times above smoothed sfr
	ic=0
	tmax=-1.e20
	imax=0
	do i=1,isfrc
	tadd(i)=usr_sfr(i)
	enddo
	call smthpl(tadd,isfrc,5)
	do i=1,isfrc
	r=usr_sfr(i)/tadd(i)
	if (r.gt.1.5) then
		do j=-2,2
		ic=ic+1
		is(ic)=i+j
		enddo
	endif
c	elseif (usr_sfr(i).gt.tmax) then
	if (usr_sfr(i).gt.tmax) then
		tmax=usr_sfr(i)
		imax=i
	endif
	enddo
	if (imax.eq.1.or.imax.eq.isfrc) then
c		Look for changes in slope
		call sfr_slope(imax)
	endif
	do i=-3,3
	ic=ic+1
	is(ic)=imax+i
	enddo
 
c	Store array tb in array tadd
	jadd=0
	do i=0,nsteps-1
	jadd=jadd+1
	tadd(jadd)=tb(i)
	enddo
 
c	Add time steps not included in tb(i)
	do i=1,ic
	if (is(i).ge.1) then
		tnew=time(is(i))
		dtmin=1.E10
		do n=0,nsteps-1
		dt=abs(tb(n)-tnew)
		if (dt.lt.dtmin) then
			dtmin=dt
			nmin=n
		endif
		enddo
c		write (6,*) i,tnew,tb(nmin),dtmin,dtmin/tb(nmin)
		if (dtmin/tb(nmin).gt.1.E-4) then
			jadd=jadd+1
			tadd(jadd)=tnew
		endif
	endif
	enddo
	eps=0.001
	call usort(jadd,tadd,eps)
 
c	Find time steps not included in array tb
	n=jadd
	do j=1,n
	do i=0,nsteps-1
	if (tadd(j).eq.tb(i)) then
		tadd(j)=-tadd(j)
	endif
	enddo
	enddo
	jadd=0
	do j=1,n
	if (tadd(j).gt.0.) then
		jadd=jadd+1
		tadd(jadd)=tadd(j)
	endif
	enddo
	if (jadd.gt.jts) then
		write (6,*) 'SINGULAR_POINTS: Too many time steps to add',jadd,'.  Maximum =',jts
		stop
	else
c		write (6,*) 'SINGULAR_POINTS: ',jadd,' points to add'
	endif
	return
	end
 
	SUBROUTINE SFR_SLOPE(I)
 
c	Finds points where derivative of array usr_sfr(time) changes abruptly
 
c	Array declarations
	include 'jb.dec'
	include 'csp.dec'
 
c	Find change of slope
	s1=(usr_sfr(3)-usr_sfr(2))/(alog10(time(3))-alog10(time(2)))
	do i=4,isfrc
	slope=(usr_sfr(i)-usr_sfr(i-1))/(alog10(time(i))-alog10(time(i-1)))
	d=abs((slope-s1)/s1)
	if (d.gt.1) return
	s1=slope
	enddo
	return
	end
 
	SUBROUTINE ADD_POINTS
 
c	Add points in array time(i) to array tb(i)
 
c	Array declarations
	include 'jb.dec'
	include 'csp.dec'
 
c	Store array tb in array tadd
	jadd=0
	do i=0,nsteps-1
	jadd=jadd+1
	tadd(jadd)=tb(i)
	enddo
 
c	Store array time in array tadd
	do i=1,isfrc
	jadd=jadd+1
	tadd(jadd)=time(i)
	enddo
 
c	Sort array tadd
	call usort(jadd,tadd)
 
c	Find time steps in array tadd not included in array tb
	n=jadd
	do j=1,n
	do i=0,nsteps-1
	if (tadd(j).eq.tb(i)) then
		tadd(j)=-tadd(j)
	endif
	enddo
	enddo
	jadd=0
	do j=1,n
	if (tadd(j).gt.0.) then
		jadd=jadd+1
		tadd(jadd)=tadd(j)
	endif
	enddo
	if (jadd.gt.jts) then
		write (6,*) 'ADD_POINTS: Too many time steps to add',jadd,'.  Maximum =',jts
		stop
	else
c		write (6,*) 'ADD_POINTS: ',jadd,' points to add'
	endif
	return
	end
