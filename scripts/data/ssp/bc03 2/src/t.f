      real*4 function t(h, q, z, lamb)
c...... Returns age of universe at redshift z
      implicit none
      real*4 h, q, z, a, b, c, hh0, cosh, x, lamb
c... for lambda
      real*4 aa, bb, epsr, s0, s, funq, omega0
      integer npts
      logical success
      external funq
      common /cosm/ omega0

c	H = Ho in km/sec/Mpc
c	Q = qo  (if problems with qo = 0, try 0.0001)

      a(q,z) = (sqrt(1. + ((2. * q) * z)) / (1. - (2. * q))) / (1. + z)
      b(q) = q / (abs((2. * q) - 1.) ** 1.5)
      c(q,z) = ((1. - (q * (1. - z))) / q) / (1. + z)

      cosh(x) = alog(x + sqrt((x ** 2) - 1.))
      hh0 = h * 0.001022
c     in (billion years)**(-1)

      if (lamb .ne. 0.0) then
         	omega0 = (2. * (q + 1.)) / 3.
         	aa = 0.
         	bb = 1. / (1. + z)
	        success=.false.
		s0=1.e-10
		npts=0
		do while (.not.success)
		  npts=npts+1
		  callmidpnt(funq,aa,bb,s,npts)
		  epsr=abs(s-s0)/s0
		  if (epsr.lt.1.e-4) then
			success=.true.
		  else
			s0=s
		  endif
		enddo
		t=s
      else if (q .eq. 0.0) then
      	t = 1. / (1. + z)
      else if (q .lt. 0.5) then
      	t = a(q,z) - (b(q) * cosh(c(q,z)))
      else if (q .eq. 0.5) then
      	t = (2. / 3.) / ((1. + z) ** 1.5)
      else
      	t = a(q,z) + (b(q) * cos(c(q,z)))
      end if

      t = t / hh0

      return
      end
c=======================================================================
      real*4 function funq(x)
c.... For non-zero cosmological constant
      real*4 x, omega0, omegainv
      common /cosm/ omega0
      omegainv = 1. / omega0
      funq = sqrt(omegainv) / (x*sqrt((omegainv-1.)+(1./(x**3.))))
      return
      end
