	FUNCTION GAUSS(X,X0,SIGMA)

c	Returns value of gaussian function normalized to area = 1

c	Normalization constant
	data a/0./
	if (a.eq.0) then
		pi=4.*atan(1.)
		a=1./sqrt(2.*pi)
	endif
        c=a/sigma
        gauss=c*exp(-(x-x0)**2/sigma**2/2)
!	write (998,*) x,x0,a,sigma,c
        return
        end
