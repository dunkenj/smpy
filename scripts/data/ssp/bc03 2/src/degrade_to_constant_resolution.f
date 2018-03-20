	SUBROUTINE DEGRADE_TO_CONSTANT_RESOLUTION(X,Y,N,S,S0,YS)

c	Degrades high resolution sed in (x,y) to a lower resolution
c	(of constant sigma S in wavelength) by applying a gaussian broadening function.

c	S0 = resolution (sigma) in Angstrom of input sed.

c	Returns smoothed spectrum in array (x,ys)

c	Array declaration
	real x(n),y(n),ys(n),s0(n),u(10000),g(10000)

c	Number of sigmas to perform integration
	m=6

c	Compute resulting sigma
	do i=1,n
	sigma=s
	if (sigma.le.s0(i)) then
c		write (6,*) 'sigma < sigma_model',sigma,s0(i)
c		stop
		ys(i)=y(i)
	else
		sigma=sqrt(sigma**2-s0(i)**2)
		xmax=x(i)+m*sigma
		call locate(x,n,xmax,i1)
		m2=i1+1
		m1=2*i-m2
c		write (81,*) i,x(i),sigma,xmax,m1,m2
		if (m1.lt.1) m1=1
		if (m2.gt.n) m2=n
		if (m2-m1+1.lt.10) then
			ys(i)=y(i)
		else
			k=0
			do j=m1,m2
			k=k+1
			u(k)=x(j)
			g(k)=y(j)*gauss(x(j),x(i),sigma)
			enddo
			ys(i)=trapz1(u,g,k)
		endif
	endif
	enddo
	return
	end
