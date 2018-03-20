	SUBROUTINE GAUSSIAN_V_DISP(X,Y,N,SIGMA)

!	Applies a gaussian filter to the sed in (X,Y) to reproduce
!	the effect of the velocity dispersion of stars in a galaxy
!	Enter sigma in Km/s

!	Array declaration
	real x(n),y(n),z(50000),u(10000),g(10000)

!	If sigma = 0, return
	if (sigma.le.0.) return

!	Speed of light in Km/s
	c=3.00E5

!	Number of sigmas to deviate from v=0
	m=6

!	Copy array y to array z
	if (n.gt.50000) then
		write (6,*) 'Too many data points in GAUSSIAN_V_DISP:',n
		stop
	endif
	do i=1,n
	z(i)=y(i)
	enddo

	do i=1,n
	xmax=c*x(i)/(c-m*sigma)
	m2 = min(1 + ilocate(x,n,xmax),n)
	m1 = max(2*i-m2,1)
	mt = m2-m1+1
	if (mt > 10000) then
		write (6,*) 'Too many data points in GAUSSIAN_V_DISP:',m1,m2,mt
		stop
	elseif (m1 > 1 .and. x(m2) <= 8750.) then
		k=0
		do j=m2,m1,-1
		v=(x(i)/x(j)-1.)*c
		k=k+1
		u(k)=v
		g(k)=z(j)*gauss(v,0.,sigma)
		enddo
		y(i)=trapz1(u,g,k)
	endif
	enddo
	return
	end
