	FUNCTION DISMOD(H0,Q0,Z)

c	Returns cosmological distance modulus

c	H0 = Ho in km/sec/Mpc
c	Q0 = qo
c	DL = Luminosity distance in Mpc

	dismod = 5*alog10(dl(h0,q0,z)*1.E6/10.)

	return
	end
