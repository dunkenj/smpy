	REAL FUNCTION LINEAR (X0,X,Y,N,I)

!	Interpolates linearly the function y(x) at x=x0.
!	Re written by GBA June 8, 2012

!	Variables
	real x(n),y(n)

!	Several options
	if ( i >= 0 ) then
		i = iplace(x0,x,max0(1,i),n)
	else
		i = -i
	endif
	if ( i == 0 ) then
		linear = 0.
!		linear = -1.e10
		ierr = linerr(x0,x(1),x(n),n)
	elseif ( x0 == x(n) ) then
		linear = y(n)
	else
		linear = y(i) + (y(i+1)-y(i))*(x0-x(i))/(x(i+1)-x(i))
	endif
	return
	end

	FUNCTION LINERR(X0,X1,XN,N)
!	Reports error in interpolation
	data ierr/0/
	if (ierr.eq.0) then
		ierr=1
		linerr=1
		write (6,'(10x,a)') ' ---------------------------------------------------------------------------------'
		write (6,1) X0,X1,XN,N
		write (6,'(10x,a)') ' ----------- Error reported only once, but it may occur more than once. ----------'
1		format (10x,' --- LINEAR:  X0 = ',1PE10.3,' is outside X range ---',2E10.3,I10)
!		stop
	endif
	return
	end

!      REAL FUNCTION LINEAR (X0,X,Y,N,I0)
!c     INTERPOLATES LINEARLY THE FUNCTION Y(X) AT X=X0.
!      REAL X(N),Y(N)
!	data ierr/0/
!      IF (I0) 21,20,20
!   20 I=IPLACE(X0,X,MAX0(1,I0),N)
!      GOTO 22
!   21 I=-I0
!   22 I0=I
!      IF (I.EQ.0) GOTO 2
!      IF (X0.EQ.X(I).OR.X(I).EQ.X(I+1)) GOTO 1
!      LINEAR=Y(I) + (Y(I+1)-Y(I))*(X0-X(I))/(X(I+1)-X(I))
!      RETURN
!    1 LINEAR=Y(I)
!      RETURN
!    2 if (ierr.eq.0) then
!	ierr=1
!c	write (6,*) 'I,I0 =',i,i0
!	write (6,3) X0,X(1),X(N),N
!	write (6,*) '---   Error reported only once. It may occur more than once. ---'
!      endif
!    3 FORMAT (' --- LINEAR:  X0 = ',1PE10.3,' is outside X range ---',2E10.3,I10)
!      LINEAR=0.
!      RETURN
!      END
