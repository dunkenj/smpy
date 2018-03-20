	FUNCTION FWHM_IDS(X)

c	Returns FWHM at wavelength X(A) according to Table 8 of
c	Worthey and Ottaviani (1997, ApJS, 111, 377)

c	Data points from Table 8
	parameter (nw=5)
	real w(nw),f(nw),linear
	data w/4000.,4400.,4900.,5400.,6000./
	data f/ 11.5,  9.2,  8.4,  8.4,  9.8/

c	Interpolate at wavelength X
	icall=0
	if (x.le.w(1)) then
		fwhm_ids=f(1)
	elseif (x.ge.w(nw)) then
		fwhm_ids=f(5)
	else
		icall=0
		fwhm_ids=linear(x,w,f,nw,icall)
	endif
	return
	end
