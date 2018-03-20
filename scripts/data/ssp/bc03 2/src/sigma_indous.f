	FUNCTION SIGMA_INDOUS(X)

c	Returns resolution (sigma) of IndoUS library at wavelength X in A.
c	Fit provided by the P. Coelho as FWHM vs lambda.

c	compute FWHM = d at wavelength X
	if (x.lt.3600.) then
c		d = 1.6
c		Updated 26/5/06 by Paula
		d = 1.57
	else
c		d = 2.15 - 1.54E-4*x
c		Updated 26/5/06 by Paula
		d = 1.92 - 9.84E-05*x
	endif
	c=2.*sqrt(2.*alog(2.))
	sigma_indous=d/c
	return
	end

	FUNCTION FWHM_INDOUS(X)

c	Returns resolution (FWHM) of IndoUS library at wavelength X in A.
c	Fit provided by the P. Coelho as FWHM vs lambda.

c	compute FWHM = d at wavelength X
	if (x.lt.3600.) then
c		d = 1.6
c		Updated 26/5/06 by Paula
		d = 1.57
	else
c		d = 2.15 - 1.54E-4*x
c		Updated 26/5/06 by Paula
		d = 1.92 - 9.84E-05*x
	endif
	fwhm_indous=d
	return
	end
