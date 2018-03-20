      real*4 function dl(h,q,z)
c.....
c	  Computes luminosity distance corresponding to a redshift z.
c	  Uses Mattig formulae for qo both 0 and non 0
c	  Revised January 1991 to implement cosmolgical constant
c	Ho in km/sec/Mpc
c	******DL is in Mpc******
c.....
      implicit none
c     include 'cosmo2.dec'
      real*4 h, q, z, d1, d2
c......
      real*4 aa, bb, epsr, s, s0, funl
      real*4 dd1, dd2, omega0
      logical success
      integer npts
      external funl
      common /cosm/ omega0
c.....
      if (z.le.0.) then
c	10 pc
	dl=1.e-5
      	return
      endif

      if (q .eq. 0) then
      dl = ((3.e5 * z) * (1 + (z / 2.))) / h
      else if (q .gt. 0.) then
      d1 = (q * z) + ((q - 1.) * (sqrt(1. + ((2. * q) * z)) - 1.))
      d2 = ((h * q) * q) / 3.e5
      dl = d1 / d2
      else if (q .lt. 0.) then
         omega0 = (2. * (q + 1.)) / 3.
         aa = 1.
         bb = 1. + z
         success=.false.
         s0=1.e-10
         npts=0
         do while (.not.success)
             npts=npts+1
             call midpnt(funl,aa,bb,s,npts)
             epsr=abs(s-s0)/s0
             if (epsr.lt.1.e-4) then
                success=.true.
             else
                s0=s
             endif
         enddo
         dd1=s
	 dd2 = (3.e5 * (1. + z)) / (h * sqrt(omega0))
	 dl = dd1 * dd2
      end if
c.....
      return
      end
c=======================================================================
      real*4 function funl(x)
c.... For non-zero cosmological constant
      real*4 x, omega0, omegainv
      common /cosm/ omega0
      omegainv = 1. / omega0
      funl = 1. / sqrt(((x ** 3.) + omegainv) - 1.)
      return
      end
c=======================================================================
