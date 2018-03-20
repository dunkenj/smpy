	PROGRAM INTERPOLATE_ISED

!	Writes to standard output interpolated *.ised file at requested ages

!	Array declarations
	parameter (nw=25000,no=250)
	character*1024 name,argu
	real tx(no),x(nw),y(nw),z(nw,no)

!	Get argument list
        jn = iargc()
	if (jn > 1) then
!		Get file name from argument list
		call getarg(1,name)
		call chaext(name,'ised',m)
!		Read requested ages
                call getarg(2,argu)
                ngu = nargu(argu)
                read (argu,*,err=1) (tx(j),j=1,ngu)
!		Call interpolation routine and store sed's
		do j=1,ngu
		call intrp_ised(name,tx(j),x,y,n)
		do i=1,n
		z(i,j) = y(i)
		enddo
		enddo
		write (6,100) name(1:largo(name))
		write (6,101) (j+1,j=1,ngu)
		write (6,102) (tx(j),j=1,ngu)
		write (6,103) ('    Flux   ',j=1,ngu)
		do i=1,n
		write (6,104) x(i),(z(i,j),j=1,ngu)
		enddo
100		format ('# File: ',a)
101		format ('# Column  ',250I11)
102		format ('# Age(yr)     ',1p250E11.3)
103		format ('# Lambda(A)   ',250a)
104		format (1pe14.6,250e11.4)
		stop
	endif
1	write (6,*) 'Command line usage:'
	write (6,*) 'interpolate_ised   filename    t(1),t(2),t(3),...          (enter ages t(1),t(2),t(3),...  in yr)'
	end
