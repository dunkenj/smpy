	subroutine dattim(ifile,k,a)

c	Prints current time and date at center of page

	character a,date*24
	call fdate(date)
	if (a.eq.'s'.or.a.eq.'S') then
		write (ifile,100) date
100		format ('#',22x,a)
	elseif (a.eq.'l'.or.a.eq.'L') then
		write (ifile,101) date
101		format ('#',53x,a)
	endif
	return
	end
