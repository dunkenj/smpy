        FUNCTION TRAPZ3(X,Y,N,X1,X2,IERR)

c	Especially written to compute spectral indices in the Lick System

c       Finds area below SED y(x) from x=X1 to x =X2
c       X1 and X2 are not necessarily values of x(i)
c       but must fulfill X1>=x(1) and X2<=x(n)

c	Routine TRAPZ2   is used for high resolution seds (1 A sampling)
c	Routine TRAPZ_LR is used for low  resolution seds (> 1 A sampling)

c       Array declaration
        real x(n),y(n)

c	Check for wavelength scale
c	The following wavelength scales are in use in BC03 models
c	Low Res pure LBC library:         1221 wavelength points
c	Low res extended Pickles library: 6505 wavelength points
c	High Res pure STELIB library:     6700 wavelength points
c	High Res STELIB + LBC library:    6900 wavelength points
c	High Res STELIB+IUE+LBC library: 11288 wavelength points

c	Total number of points are not fed to this routine
c	if (n.lt.6600) then
c	Determine wavelength step
	dx=x(2)-x(1)
	if (dx.gt.5.) then
c		Low resolution sed
c		write (6,*) 'Using trapz_lr:',x(1),x(2),dx
		trapz3=trapz_lr(x,y,n,x1,x2,ierr)
	else
c		High resolution sed
c		write (6,*) 'Using trapz2:',x(1),x(2),dx
		trapz3=trapz2(x,y,n,x1,x2,ierr)
	endif
	return
	end
