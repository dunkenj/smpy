	INTEGER FUNCTION IPLACE(X0,XX,N1,N2)

!	Finds L such that XX(L) < X0 < XX(L+1)
!       Re written by GBA June 8, 2012

!       Variables
	real xx(n2)

!	Procedure
	x(i)=xx(n1+i-1)
	n=n2-n1+1
	if ( x(1) > x(2) ) then
		i=2
	else
		i=1
	endif
	l=0
	do while ( l < n-1 )
	l = l + 1
	if ( i == 1 ) then
		if ( x(l) <= x0 .and. x0 < x(l+1) ) then
			iplace=l+n1-1
			return
		endif
	else
		if ( x(l) >= x0 .and. x0 > x(l+1) ) then
			iplace=l+n1-1
			return
		endif
	endif
	enddo
	if ( x0 == x(n) ) then
		iplace=l+n1-1
	else
		iplace=0
	endif
	return
	end

!      INTEGER FUNCTION IPLACE(X0,XX,N1,N2)
!c     FINDS L SUCH THAT   XX(L).LT.X0.LT.XX(L+1)
!      DIMENSION XX(N2)
!      X(I)=XX(N1+I-1)
!      N=N2-N1+1
!      I=1
!      IF (X(1).GT.X(2)) I=2
!      L=0
!    1 L=L+1
!      IF (L.EQ.N) GOTO 4
!      GOTO (2,3),I
!    2 IF (X(L).LE.X0.AND.X0.LT.X(L+1)) GOTO 5
!      GOTO 1
!    3 IF (X(L).GE.X0.AND.X0.GT.X(L+1)) GOTO 5
!      GOTO 1
!    4 IF (X0.NE.X(N)) GOTO 6
!    5 IPLACE=L+N1-1
!      RETURN
!    6 IPLACE=0
!c	write (6,*) 'IPLACE: X0,N1,N2,X(1),X(N)',x0,n1,n2,x(1),x(n)
!c	write (6,*) 'x(2) =',x(2)
!      RETURN
!      END
