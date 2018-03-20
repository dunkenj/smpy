	FUNCTION TRAPZQ(X,Y,N,X1,X2,IERR)

c	Finds area below y(x) from x=X1 to x =X2
c	X1 and X2 are not necessarily values of x(i)
c	but must fulfill X1>=x(1) and X2<=x(n)

c	Starts search for X1 at last point found

c	Array declaration
	real x(n),y(n)
	data i1/1/

	trapzq=0.
	ierr=0.
c	Check for values outside valid range
	if (x1.lt.x(1)) then
		write (6,*) 'TRAPZQ: X1 outside valid range:',x1,x(1)
		ierr=1
		return
c		stop
	endif
	if (x2.gt.x(n)) then
		write (6,*) 'TRAPZQ: X2 outside valid range:',x2,x(n)
		ierr=2
		return
c		stop
	endif

c	Search array from last value of i1
c	Reset i1
	if (x1.lt.x(i1)) i1=1

c	Find first element in array X above X1
	do i=i1,n
	if (x(i).ge.x1) then
		i1=i
c		Find by interpolation values of Y1=y(X1)
		Y1 = ( (x(i1)-x1)*y(i1-1) + (x1-x(i1-1))*y(i1) ) / (x(i1)-x(i1-1))
c		Add area from x(i1-1) to X1
		trapzq=trapzq+(x(i1)-x1)*(Y1+y(i1))/2.
		goto 1
	endif
	enddo

c	Find last element in array X below X2
1	do i=i1,n
	if (x(i).ge.x2) then
		i2=i-1
c		Find by interpolation values of Y2=y(X2)
		Y2 = ( (x(i2+1)-x2)*y(i2) + (x2-x(i2))*y(i2+1) ) / (x(i2+1)-x(i2))
c		Add area from x(i2) to X2
		trapzq=trapzq+(x2-x(i2))*(Y2+y(i2))/2.
		goto 2
	endif
	enddo

c	Compute area from x(i1) to x(i2)
2	do i=i1+1,i2
	trapzq=trapzq+(x(i)-x(i-1))*(y(i)+y(i-1))/2.
	enddo 

	return
	end
