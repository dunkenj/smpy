	SUBROUTINE NO_BL(A)

c	Elimina espacios en blanco (todos) en A

	character*(*) a

	if (largo(a).eq.0) return
1	i=index(a(1:largo(a)),' ')
	if (i.eq.0) return
	a(i:)=a(i+1:)
	goto 1
	end
