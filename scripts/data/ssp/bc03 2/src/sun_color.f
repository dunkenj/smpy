	Program sun_color

c	Computes color of Kurucz solar spectrum in filter pair (I,J)

	write (6,*) 'Returns color of Kurucz solar spectrum model in a given filter pair'
	write (6,*)

1	write (6,'(x,a,$)') 'Enter filter pair = i,j = '
	read (5,*,end=2,err=1) i,j
	sc=sun_color_n(i,j)
	write (6,100) sc
	goto 1
100	format (35x,'Sun color = ',f8.5/)
2	end
