	FUNCTION TRAPZ2(X,Y,N,X1,X2,IERR)

c	Finds area below y(x) from x=X1 to x =X2
c	X1 and X2 are not necessarily values of x(i)
c	but must fulfill X1>=x(1) and X2<=x(n)

c	Array declaration
	real x(n),y(n)
	data jerr1,jerr2/2*0/

	TRAPZ2=0.
	ierr=0.

!	Check for values outside valid range
	if (x1 < x(1)) then
		if (jerr1.eq.0) then
			write (6,*) 'TRAPZ2: X1 outside valid range:',x1,x(1),x(n),n
			jerr1=1
		endif
		ierr=1
		return
!		stop
	elseif (x1 == x(1)) then
		i1 = 1
	else
!		Find first element in array X above X1
		do i=1,n
		if (x(i) > x1) then
			i1=i
!			Find by interpolation values of Y1=y(X1)
			Y1 = ( (x(i1)-x1)*y(i1-1) + (x1-x(i1-1))*y(i1) ) / (x(i1)-x(i1-1))
!			Add area from x(i1-1) to X1
			trapz2=trapz2+(x(i1)-x1)*(Y1+y(i1))/2.
			goto 1
		endif
		enddo
	endif
1	if (x2 > x(n)) then
		if (jerr2.eq.0) then
			write (6,*) 'TRAPZ2: X2 outside valid range:',x2,x(1),x(n),n
			jerr2=1
		endif
		ierr=2
		return
!		stop
	elseif (x2 == x(n)) then
		i2 = n
	else
!		Find last element in array X below X2
		do i=i1,n
		if (x(i) > x2) then
			i2=i-1
!			Find by interpolation values of Y2=y(X2)
			Y2 = ( (x(i2+1)-x2)*y(i2) + (x2-x(i2))*y(i2+1) ) / (x(i2+1)-x(i2))
!			Add area from x(i2) to X2
			trapz2=trapz2+(x2-x(i2))*(Y2+y(i2))/2.
			goto 2
		endif
		enddo
	endif

!	Compute area from x(i1) to x(i2)
2	do i=i1+1,i2
	trapz2=trapz2+(x(i)-x(i-1))*(y(i)+y(i-1))/2.
	enddo 
!	write (6,*) 'TRAPZ2: using = ',i1+1,i2,x(i1),x(i2),i2-i1,' points'
	return
	end
