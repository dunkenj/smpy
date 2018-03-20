	real*4 function zx(tx, h, q, lamb)

c	Returns the value of Z = ZX (redshift) corresponding to a given
c	light travel time TX (measured in Gyr).

c	H = Ho in km/sec/Mpc
c	Q = qo (if problems with qo = 0, try 0.0001)

	implicit none
	real*4 ltt, z(47), age, tx, h, q, hh0, dz, t, omega0, omegainv, lamb, zi
	integer j
	external t
	data z /0.,.001,.002,.004,.006,.008,.01,.02,.04,.06,.08,.1,.2,.3,
     &      .4,.5,.545,.6,.7,.8,.9,.945,1.,1.2,1.4,1.6,1.8,2.,3.,4.,5.,6.,
     &      7.,8.,9.,10.,12.,14.,16.,18.,20.,30.,40.,60.,80.,100.,1000./

c     data z / 0., .001, .002, .004, .006, .008, .01, .02, .04, .06, .08
c    &, .1, .2, .3, .4, .5, .545, .6, .7, .8, .9, .945, 1., 1.2, 1.4,
c    &1.6, 1.8, 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18.
c    &, 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80.
c    &, 85., 90., 95., 100., 500., 1000. /

c	Define light travel time
	ltt(q,zi) = age - t(h,q,zi,lamb)

c	Compute omega
	if (lamb .ne. 0.) then
		omega0 = (2. * (q + 1.)) / 3.
		omegainv = 1. / omega0
	endif

c	Check for zero value
	zx = 0.
	if (tx .eq. 0.) return

c	Check for maximum age of universe
	age = t(h,q,0.,lamb)
	zx = -2.
	if (tx .ge. age) return

c	Express Ho in billion years ** (-1)
	hh0 = h * 0.001022

c	Check for q=0.5 case
	if (q .eq. 0.5) then
		zx = ((1. - (((3. * hh0) * tx) / 2.)) ** (- (2. / 3.))) - 1.
		return
	endif

c	General case
	do j = 1, 47
	if (tx .le. ltt(q,z(j))) goto 2
	enddo
2	if (j > 1) then
		zx = z(j - 1)
	else
		zx = z(1)
	endif
	zi=zx

c	Iterate to find best value
	do j = 1, 100000
	if (q .ge. 0.) then
c		zero cosmological constant
		dz = ((hh0 * ((tx - t(h,q,0.,lamb)) + t(h,q,zx,lamb))) *
     *		     ((1. + zx) ** 2)) * ((1. + ((2. * q) * zx)) ** 0.5)
	elseif (q .lt. 0.) then
c		non-zero cosmological constant
		dz = ((hh0 * ((tx - t(h,q,0.,lamb)) + t(h,q,zx,lamb))) *
     *		     (1. + zx)) * sqrt(omega0)
		dz = dz * sqrt((((1. + zx) ** 3.) + omegainv) - 1.)
	endif
c	if (abs(dz / zx) .lt. 0.00001) return
c	if (abs(dz / zx) .lt. 0.00005) return
c	if (abs(dz / zx) .lt. 0.00010)  return
	if (abs(dz / zx) .lt. 0.00003) return
	zx = zx + dz
	enddo
	return
	end
