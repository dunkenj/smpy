	SUBROUTINE DEGRADE_RESOLUTION(X,Y,N,FWHMS,XS,YS,NS)

c	Degrades STELIB sed in (x,y) to lower resolution sed by applying a
c	gaussian broadening function. The FWHM of the gaussian is supposed
c	to be constant (= FWHM) as a function of wavelength.

c	Returns smoothed spectrum in array (x,ys)

c	Array declaration
	real x(n),y(n),xs(n),ys(n),u(1000),g(1000)

c	Check if fwhms = 0
	if (fwhms.le.0) then
		ns=n
		do i=1,ns
		xs(i)=x(i)
		ys(i)=y(i)
		enddo
		return
	endif

c	Number of sigmas to perform integration
	m=6

c	FWHM of STELIB seds = 3 A (Leborgne et al. 2002)
	fwhm_stelib = 3.

c	Perform smoothing
	fwhm2 = fwhms**2 - fwhm_stelib**2
	fwhm  = sqrt(fwhm2)

c	Note: FWHM = 2*sqrt[2*ln(2)]*sigma = 2.3548*sigma
	sigma=fwhm/2.3548
	do i=1,n
	xmax=x(i)+m*sigma
	call locate(x,n,xmax,i1)
	m2=i1+1
	m1=2*i-m2
	if (m1.lt.1) m1=1
	if (m2.gt.n) m2=n
	k=0
	do j=m1,m2
	k=k+1
	u(k)=x(j)
	g(k)=y(j)*gauss(x(j),x(i),sigma)
	enddo
	ys(i)=trapz1(u,g,k)
	xs(i)=x(i)
	enddo

c	Suppress points at both ends of smoothed sed
	k=fwhms
c	Not really
	k=0
	ns=0
	do i=1+k,n-k
	ns=ns+1
	xs(ns)=xs(i)
	ys(ns)=ys(i)
	enddo

	return
	end
