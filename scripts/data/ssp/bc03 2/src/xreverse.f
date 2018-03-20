	SUBROUTINE XREVERSE(X,N)

c	Reverses order of real array X

	real*4 x(n),y(1000000)
	nt=n+1
	do k=1,n
	nt=nt-1
	y(nt)=x(k)
	enddo
	do k=1,n
	x(k)=y(k)
	enddo
	return
	end

	SUBROUTINE XREVERSE8(X,N)

c	Reverses order of real*8 array X

	real*8 x(n),y(1000000)
	nt=n+1
	do k=1,n
	nt=nt-1
	y(nt)=x(k)
	enddo
	do k=1,n
	x(k)=y(k)
	enddo
	return
	end
