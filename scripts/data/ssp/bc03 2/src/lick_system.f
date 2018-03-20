	SUBROUTINE LICK_SYSTEM(X,Y,N,YS)

c	Transforms sed in (x,y) to Lick/IDS system by applying a gaussian
c	broadening function to artificially simulate the Lick/IDS spectral
c	resolution. The FWHM of the gaussian as a function of wavelength is
c	given in Table 8 of Worthey and Ottaviani (1997, ApJS, 111, 377)

c	Returns smoothed spectrum in array (x,ys)

c	Array declaration
	real x(n),y(n),ys(n),u(10000),g(10000)

c	Number of sigmas to perform integration
	m=6

c	FWHM of STELIB seds = 3 A (Leborgne et al. 2002)
	fwhm_stelib = 3.

c	Perform smoothing
	do i=1,n
	fwhm2 = fwhm_ids(x(i))**2 - fwhm_stelib**2
	fwhm  = sqrt(fwhm2)

c       Note: FWHM = 2*sqrt[2*ln(2)]*sigma = 2.3548*sigma
	sigma=fwhm/2.3548

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
	enddo
	return
	end
