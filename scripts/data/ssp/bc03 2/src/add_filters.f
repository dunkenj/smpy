	PROGRAM ADD_FILTERS

c	Modified 18-JUN-2003: G. Bruzual, CIDA
c		Changed input format to free format
c		First record is filter ID, with a '#' as the first character
c		Then follow 2 columns, lambda(A), R(lambda)
c		No need to know number of points in each response function

c	Modified 11-MAY-1999: G. Bruzual, CIDA
c		Filter response functions are stored in a binary file
c		as published, with no addition of extra points by interpolation.
c		The new version of the FILTER routine requires the filter
c		response function to be interpolated at each point in the
c		sed which is not a point in the filter response. All points in
c		the filter response are also used. The format of the binary
c		file was changed to reduce its size. Filters are stored
c		sequentially in r(i,k), k=1 => wavelength, k=2 => response
c		function. Information about starting point and number of
c		points per filter is kept in file, as well as filter id label.

c	Modified 09-JUN-1983: G. Bruzual
c		Program was modified on June-09-83 to include interpolation
c		of the filter response functions every 25 A only for those
c		filters for which the resulting number of points is .le. 250
c		and .gt. previous number of data points. A 25 A wavelength
c		step was found appropriate to compute synthetic colors with
c		the available stellar and galaxy s.e.d.''s.
c		For the IR filters a 50 A step is used.
c		For the L filter a 75 A step is used.

c	Written: 21-MAY-1983: G. Bruzual, Durham.

c	Produces a binary file with filter response functions.
c	The response functions are assumed to be in a formatted file.
c	The response functions are read into array r(i,j,k),
c	where j identifies the filter, i the data point, and k=1 => wavelength
c	scale, whereas k=2 => response function
c	Assumes up to mxft filters. The maximum number of points per
c	filter is mxwf, declared as parameters in filter.dec.

c	Array declaration
	include 'filter.dec'
	character ans,name*96,aux*128
	integer mp(0:mxft)
	real w1(naux),c1(naux)
	fid(0)=' '

c	Ask for input file
	write (6,'(x,a,$)') 'Filter response functions in (formatted) file = '
	read (5,'(a)',end=2) name
	write (6,'(x,a,$)') 'Normalize response functions to area = 1 (y or [n]) ? '
	read (5,'(a)',end=2) ans
	iarea=0
	if (ans.eq.'y'.or.ans.eq.'Y') iarea=1
	close (1)
	open (1,file=name,form='formatted',status='old')

c	Open filter log file
	close (3)
	open (3,file='filters.log',status='unknown')

c	Process input file
	nf=0
	ltot=0
	ktot=0
c	Read file
	do i=1,mxwp
	read (1,'(a)',end=1) aux
	if (aux(1:1).eq.'#') then
		mp(nf)=ktot
		nf=nf+1
		fid(nf)=aux(3:)
		ktot=0
	else
		ktot=ktot+1
		ltot=ltot+1
		read (aux,*) r(ltot,1),r(ltot,2)
		if (r(ltot,2).lt.0.) r(ltot,2)=0.
	endif
	enddo
1	mp(nf)=ktot
c	Define number of points per filter
	do i=1,nf
	np(i)=mp(i)
	if (i.eq.1) then
		ni(i)=1
	else
		ni(i)=ni(i-1)+np(i-1)
	endif
	nl(i)=ni(i)+np(i)-1
c	write (8,'(8i8)') i,np(i),ni(i),nl(i),ltot

c	Normalize filter, if requested
	if (iarea.gt.0) then
!		write (6,*) 'Normalizing to area = 1 filter',i
		l=0
		do k=ni(i),nl(i)
		l=l+1
		w1(l)=r(k,1)
		c1(l)=r(k,2)
		enddo
		area=trapz1(w1,c1,np(i))
		do k=ni(i),nl(i)
		r(k,2)=r(k,2)/area
		enddo
	endif
	write (6,101) i,fid(i),np(i)
	write (3,101) i,fid(i),np(i)
101	format (i4,1x,a,i6)
	enddo

c	Write output file
	write (6,'(/x,a,$)') 'Write binary file (y or [n]) ? '
	read (5,'(a)',end=2) ans
	if (ans.ne.'Y'.and.ans.ne.'y') stop 'Filter index in file filters.log'
	close (1)
	write (6,'(x,a,$)') 'Enter output (binary) file name = '
	read (5,'(a)',end=2) name
	open (1,file=name,form='unformatted',status='unknown')
c	write (1) nf,np,r
c	New format (11-MAY-1999)
	write (1) nf,(ni(i),i=1,nf),(nl(i),i=1,nf),(np(i),i=1,nf),
     &		(fid(i),i=0,nf),ltot,(r(i,1),i=1,ltot),(r(i,2),i=1,ltot)
	write (6,'(x,a)') 'Filter index in file: filters.log'
	stop 'Binary file written.'
2	end
