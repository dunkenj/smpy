	FUNCTION GAL_MASS(IO,T,SFR)

c	Computes the galaxy mass up to time t as the integral of the SFR

c	Array declaration
	real*8 glast,trapz8,tg(20000),sf(20000)
	data ic,ig,tg,sf/0,0,40000*0./,glast/0./

c	Zero arrays on first call
	if (ic.eq.0) then
		do i=1,10000
		tg(i)=0.
		sf(i)=0.
		enddo
		gal_mass=0.
		ic=1
		ig=0
	endif

c	Check if SSP (IO=0)
	if (io.eq.0) then
c		SSP, 1 Mo of stars formed at t = 0
		gal_mass=1.
	else
c		Store SFR and compute integral
		ig=ig+1
		tg(ig)=t
		sf(ig)=sfr
		if (sfr > 0.) then
			glast =trapz8(tg,sf,ig)
			if (ig == 1) glast = sfr*t
		endif
		gal_mass=glast
	endif
!	if (t.gt.19.9E9) ic=0
	if (t >= 20.0E9) ic=0
	return
	end

	REAL*8 FUNCTION TRAPZ8(X,Y,N)
	REAL*8 X(N),Y(N)
	TRAPZ8=0.
	IF (N.LE.1) RETURN
	DO J=2,N
	TRAPZ8= TRAPZ8 + ABS(X(J)-X(J-1))*(Y(J)+Y(J-1))/2.
	ENDDO
	RETURN
	END
