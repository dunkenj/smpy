	PROGRAM SEEFILTER

c	Writes an ASCII file with the response function of the
c	selected filters.

	include 'filter.dec'
	real x(5000),y(5000),z(5000)

	call FILTER0

	l1=0
	l2=0
1	write (6,100) 'Extract filter number = '
	read (5,101,err=1,end=10) n
100	format (1x,a,$)
101	format (i10)
	nt=np(n)
	l1=l1+3
	l2=l2+nt+2
	write (28,'(''# '',a,i4,2a)') 'Filter',n,': ',fid(n)(1:largo(fid(n)))
	write (28,'(''# '',i4,a,i4,a,i4)') nt,' data points. Lines',l1,' to',l2
	nt=0
	do i=ni(n),nl(n)
	nt=nt+1
	x(nt)=r(i,1)
	y(nt)=r(i,2)
	z(nt)=y(nt)*x(nt)
	write (28,102) x(nt),y(nt)
102	format (1p2e14.5)
	enddo
	area=trapz1(x,y,nt)
	weff=trapz1(x,z,nt)/area
	write (6,'(x,a,i4,2a)') 'Filter',n,': ',fid(n)
	write (6,'(x,a,1p2e12.4)') 'Range covered        = ',x(1),x(nt)
	write (6,'(x,a,1pe12.4)') 'Area below filter    = ',area
	write (6,'(x,a,1pe12.4/)') 'Effective wavelength = ',weff
	write (6,*) '------------------------------------------------------'
	write (9,'(x,a,i4,2a)') 'Filter',n,': ',fid(n)
	write (9,'(x,a,1p2e12.4)') 'Range covered        = ',x(1),x(nt)
	write (9,'(x,a,1pe12.4)') 'Area below filter    = ',area
	write (9,'(x,a,1pe12.4/)') 'Effective wavelength = ',weff
	write (9,*) '------------------------------------------------------'
	l1=l2
	goto 1
10	write (6,*) 'Output in files fort.28 and fort.9'
	end

