	Program ZEROPOINT

c	Finds zeropoint for filter pair (I,J)

	write (6,*) 'Return zero point for any filter pair.'
	write (6,*) 'To compute color, use functions zerop_n and color_n'
	write (6,*)

1	write (6,'(x,a,$)') 'Enter filter pair = i,j = '
	read (5,*,end=2,err=1) i,j
	zp=zerop_n(i,j)
	write (6,100) zp
	goto 1
100	format (35x,'Zero Point = ',f8.5/)
2	end
