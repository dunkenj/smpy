	FUNCTION GC_D4000_FIT(X,Y,N,TEFF,GLOG,FEH,FL,FR)

c	Computes 4000 A break amplitude according to the algorithm written
c	by Javier Gorgas and Nicolas Cardiel (GC), based on fitting functions.

c	Needs Teff, log g, and [Fe/H] for star of sed (x,y)

c	Returns value of fl = F(4000-) and fr = F(4000+)
c	D4000 = fr/fl

	real x(n),y(n),w(2000),z(2000)
!	save i1,i2
	save i3,i4
	data icall/0/

c	Find wavelength points
	gc_d4000_fit=0.
c	if (icall.eq.0) then
c		To allow different x scales in same program, always locate
c		icall=1

c		call locate (x,n,4050.,i1)
c		call locate (x,n,4250.,i2)
c		i2=i2+1

		call locate (x,n,3750.,i3)
		call locate (x,n,3950.,i4)
		i4=i4+1

c	endif

c	Transform to Fnu units
c	i=0
c	do j=i1,i2
c	i=i+1
c	w(i)=x(j)
c	z(i)=y(j)*x(j)**2
c	enddo
c	fr=trapz1(w,z,i)/(x(i2)-x(i1))

c	Compute F(4000-)
	i=0
	do j=i3,i4
	i=i+1
	w(i)=x(j)
	z(i)=y(j)*x(j)**2
	enddo
	fl=trapz1(w,z,i)/(x(i4)-x(i3))

c	Compute break amplitude (= findex) according to Gorgas and Cardiel
	call eff_d4000(teff,glog,feh,findex,eindex,iflag)

c	Compute F(4000+)
	fr=fl*findex
	gc_d4000_fit=findex

	return
	end

	FUNCTION GC_D4000_OBS(X,Y,N)

c	Measures 4000 A break from sed according to J. Gorgas and N. Cardiel

	parameter (jid=35)
	real x(n),y(n),w(11000),z(11000)
	save i1,i2,i3,i4
c	Common to store fluxes
	integer iim(jid)
	real ffb(jid),ffr(jid),ffc(jid),ffl(jid)
        common /fluxes/ iim,ffb,ffr,ffc,ffl
	data icall/0/

c	Find wavelength points
	gc_d4000_obs=0.
	if (icall.eq.0) then
c		To allow different x scales in same program, always locate
c		icall=1

		call locate (x,n,4050.,i1)
		call locate (x,n,4250.,i2)
		i2=i2+1

		call locate (x,n,3750.,i3)
		call locate (x,n,3950.,i4)
		i4=i4+1

	endif

c	Transform to Fnu units
c	Compute F(4000+)
	i=0
	do j=i1,i2
	i=i+1
	w(i)=x(j)
	z(i)=y(j)*x(j)**2
	enddo
	fr=trapz1(w,z,i)/(x(i2)-x(i1))

c	Compute F(4000-)
	i=0
	do j=i3,i4
	i=i+1
	w(i)=x(j)
	z(i)=y(j)*x(j)**2
	enddo
	fl=trapz1(w,z,i)/(x(i4)-x(i3))

c	Compute break amplitude
	gc_d4000_obs=fr/fl

c	Store fluxes for future use
	iim(30)=-1
	ffb(30)=fl
	ffr(30)=fr
	ffc(30)=gc_d4000_obs
	ffl(30)=gc_d4000_obs

	return
	end

c-----------------------------------------------------------------------------
c Version 02/12/98
c-----------------------------------------------------------------------------
c Program to evaluate the D4000 index using the empirical fitting functions
c from Gorgas, Cardiel, Pedraz, Gonzalez (1998, A&A Suppl., in preparation)
c
c INPUT:
c       t = effective temperature, in K                      (REAL)
c       g = logarithm (base 10) of the surface gravity       (REAL)
c       z = metallicity ([Fe/H])                             (REAL)
c  (A value higher or equal to 99 in g or z in input means that this 
c    parameter is unknown. In same cases, an index value can still be 
c    computed. T<=0 means that temperature is unknown) 
c  
c
c OUTPUT:
c       findex = index value                                 (REAL)
c       eindex = error in the index value                    (REAL)
c       iflag = indicates sucess of the functions:           (INTEGER)
c              0 -> OK
c              1 -> extrapolation to high temperatures
c              2 -> extrapolation to low temperatures
c              3 -> extrapolation in g (uncertain value)
c             -1 -> Error (no Teff)
c             -2 -> Error (z needed)
c             -3 -> Error (g needed)
c             -5 -> Error

      subroutine eff_d4000(t,g,z,findex,eindex,iflag)
      implicit none
      integer iflag
      real t,g,z,findex,eindex
      
      integer i,j
      real thet,x,findex0,eindex0
      real ctcg,ctcd,ctd,errt,errw
      double precision theta,thetaz,theta2,theta2z,xh(25)
      double precision ccg(25,25),ccd(25,25),cw(25,25),ch(25,25),cg(25)
      double precision cd(25,25),dumt,xh2(25)
      double precision fcg(25),fcd(25),fw(25),fh(25),fg(25),fd(25)
      double precision srcg,srcd,vh,srw,srh,dx,srg,xfg,srfg,srd,fw1
      double precision sri,fi(25),ci(25,25)
      logical nog,noz

c-----------------------------------------------------------------------
c Coefficients and Variance-covariance matrices:

c For hot stars
      srh=0.4615183238D-01
      fh(1)=0.6548051838D+00
      fh(2)=0.1340106583D+01
      ch(1,1)=0.2313283739D+00
      ch(1,2)=-0.4761853888D+00
      ch(2,2)=0.1142256111D+01

c For warm stars
      srw=0.3422608269D-01
      fw(1)=0.1822983515D+01
      fw(2)=-0.5068388040D+00
      cw(1,1)=0.1663359277D+02
      cw(1,2)=-0.2414081677D+02
      cw(2,2)=0.3519089299D+02
      errw=0.159526188D-01
      fw1= 1.4428544D+00
c For hot supergiants
      sri=0.7879815732D-01
      fi(1)=0.8613232907D+00
      fi(2)=0.1848690037D+01
      ci(1,1)=0.1731426778D+00
      ci(1,2)=-0.8470102471D+00
      ci(2,2)=0.7135431404D+01

c For cool giants:
      ctcg=0.9
      srcg=0.1309941549D+00
      fcg(1)=-0.5664970674D+01
      fcg(2)=0.9278637462D+01
      fcg(4)=-0.3273398709D+01
      fcg(5)=0.7322437356D+01
      fcg(7)=-0.3080110968D+01
      fcg(16)=-0.3694125492D+01
      data ((ccg(i,j),j=1,6),i=1,6)/0.1424049628D+03,-0.2516045719D+03,
     c 0.8616301985D+02,-0.1533506540D+03,0.1102731351D+03,
     c 0.6779636477D+02,0.D0,0.4460123942D+03,-0.1534799668D+03,
     c 0.2742738641D+03,-0.1961177449D+03,-0.1217490895D+03,0.D0,0.D0,
     c 0.3013629553D+03,-0.5639211958D+03,0.6795204379D+02,
     c 0.2608696279D+03,0.D0,0.D0,0.D0,0.1060365074D+04,
     c -0.1219393397D+03,-0.4926576419D+03,0.D0,0.D0,0.D0,0.D0,
     c 0.8652160543D+02,0.5435617680D+02,0.D0,0.D0,0.D0,0.D0,0.D0,
     c 0.2298069839D+03/

c For cool dwarfs:
      ctcd=0.9
      srcd=0.1647942839D+00
      fcd(1)=-0.8153624600D+01
      fcd(2)=0.1344994016D+02
      fcd(4)=-0.1706301430D+02
      fcd(5)=0.3811728640D+02
      fcd(7)=-0.4768742248D+01
      fcd(16)=-0.2068889239D+02
      data ((ccd(i,j),j=1,6),i=1,6)/0.1482840358D+03,-0.3158204285D+03,
     c 0.1389514149D+02,-0.2538831761D+02,0.1668009227D+03,
     c 0.1078900852D+02,0.D0,0.6745000058D+03,-0.3087237636D+02,
     c 0.5813276884D+02,-0.3571841966D+03,-0.2581987859D+02,0.D0,0.D0,
     c 0.1950227161D+04,-0.4230766254D+04,0.1667832891D+02,
     c 0.2278759127D+04,0.D0,0.D0,0.D0,0.9188955680D+04,
     c -0.3219517176D+02,-0.4955197523D+04,0.D0,0.D0,0.D0,0.D0,
     c 0.1896456546D+03,0.1480236474D+02,0.D0,0.D0,0.D0,0.D0,0.D0,
     c 0.2675388494D+04/

c For cold giants
      srg=0.2107104366D+00
      fg(1)=0.9524702185D+01
      fg(2)=-0.4093622777D+01
      cg(1)=0.1861588935D+01
      cg(2)=0.2986724620D+01
      xfg=0.1301180D+01
      srfg=0.1310919D+00

c For cold dwarfs
      ctd=1.15
      srd=0.1868498444D+00
      fd(1)=-0.2908460064D+01
      fd(2)= 0.6534652942D+01
      fd(7)=-0.2975879712D+01
      data ((cd(i,j),j=1,3),i=1,3)/0.2490535985D+03,-0.3573457710D+03,
     c 0.1258060822D+03,0.D0,0.5155536330D+03,-0.1824972759D+03,0.D0,
     c 0.D0,0.6497403497D+02/

c-----------------------------------------------------------------------

      nog=.false.
      noz=.false.
      if(g.ge.99.) nog=.true.
      if(z.ge.99.) noz=.true.

      iflag=0
      findex=0.
      eindex=0.
      if(t.le.0.) then
         iflag=-1
         return
      end if

c Program works in theta (=5040/Teff)
      theta=5040.D0/dble(t)
      thet=real(theta)
c Now, it prepares coefficients for error evaluation
      thetaz=theta*dble(z)
      theta2=theta*theta
      theta2z=theta2*dble(z)
      xh(1)=1.D0
      xh(2)=theta
      xh(3)=dble(z)
      xh(4)=thetaz
      xh(5)=theta2
      xh(6)=theta2z
      dumt=0.8D+00
      xh2(1)=1.D0
      xh2(2)=dumt
      xh2(3)=dble(z)
      xh2(4)=dumt*dble(z)
      xh2(5)=dumt*dumt
      xh2(6)=dumt*dumt*dble(z)
 
c Hot supergiants
      if(thet.le.0.7) then
         if(g.lt.(3.456933-1.221595*thet)) then
            if(t.lt.0.19) iflag=1
            findex=real(fi(1)+fi(2)*theta*theta*theta)
            xh(2)=theta*theta*theta
            call comperr(2,xh,ci,vh)
            eindex=real(sri*dsqrt(vh))
            return
         end if
      end if

c Theta=(0.127,0.6326)
c Not known dependence on g or z
c Even is theta is lower, the program extrapolates.
      if(thet.le.0.6326) then
         if(t.lt.0.127) iflag=1
         findex=real(fh(1)+fh(2)*theta)
         call comperr(2,xh,ch,vh)
         eindex=real(srh*dsqrt(vh))
         return
      end if
      
c Theta=(0.6326,0.75)
      if(thet.le.0.75) then
         findex=real(fw(1)+fw(2)*theta)
         call comperr(2,xh,cw,vh)
         eindex=real(srw*dsqrt(vh))
         return
      end if

      if(nog) then
         if(thet.le.0.8) then
            findex=real(fw(1)+fw(2)*theta)
            call comperr(2,xh,cw,vh)
            eindex=real(srw*dsqrt(vh))
            return
         else
            iflag=-3
            return
         end if
      end if

c Giants (logg < 3.5)
      if(g.le.3.5) then
c Theta=(0.75,1.3)
         if(thet.lt.real(xfg)) then
c If there is no value of z,g it cannot compute the index, except for theta 
c below 0.8
            if(noz) then
               if(thet.le.0.8) then
                  findex=real(fw(1)+fw(2)*theta)
                  call comperr(2,xh,cw,vh)
                  eindex=real(srw*dsqrt(vh))
                  return
               else
                  iflag=-2
                  return
               end if
            end if
            if(thet.ge.0.8) then
               findex=ctcg+real(dexp(fcg(1)+fcg(2)*theta+fcg(4)*dble(z)+
     c           fcg(5)*thetaz+fcg(7)*theta2+fcg(16)*theta2z))
               call comperr(6,xh,ccg,vh)
               errt=real(srcg*dsqrt(vh))
c This is the error in the mean index for these atm. par. If we wanted the 
c expected error for one measurement we should compute:
c               errt=sigres*sqrt(1+vh)
               eindex=errt*(findex-ctcg)
c If theta is between 0.75 and 0.8 we interpolate between this value and the 
c constant value of warm stars at theta=0.75 (0.
            else
               dumt=0.8D+00
               findex=ctcg+real(dexp(fcg(1)+fcg(2)*dumt+fcg(4)*dble(z)+
     c           fcg(5)*dumt*dble(z)+fcg(7)*dumt*dumt+fcg(16)*dumt*dumt*
     c           dble(z)))
               call comperr(6,xh2,ccg,vh)
               errt=real(srcg*dsqrt(vh))
               eindex=errt*(findex-ctcg)
               x=(thet-0.75)/0.05
               findex=real(fw1)*(1.-x)+findex*x
               eindex=errw*(1.-x)+eindex*x
            end if
         else
c Theta=(1.3,1.72) for giants
            findex=real(fg(1)+fg(2)*theta)
            dx=theta-xfg
	    eindex=real(dsqrt(srg*srg*(dx**2)*cg(1)+
     c        srfg*srfg*(1.D0-dx*cg(2))**2))
            if(g.gt.1.5) iflag=3
            if(thet.gt.1.73) iflag=2
         end if
         if(g.le.3.) then
            return
         else
            findex0=findex
            eindex0=eindex
         end if
      end if

c Dwarfs (logg > 3.)
      if(g.gt.3.) then
c Theta=(0.75,1.07)
c Cambiar lo siguiente:
         if(thet.le.1.0757) then
c If there is no value of z it cannot compute the index, except for theta 
c below 0.8
            if(noz) then
               if(thet.le.0.8) then
                  findex=real(fw(1)+fw(2)*theta)
                  call comperr(2,xh,cw,vh)
                  eindex=real(srw*dsqrt(vh))
                  return
               else
                  iflag=-2
                  return
               end if
            end if
            if(thet.ge.0.8) then
               findex=ctcd+real(dexp(fcd(1)+fcd(2)*theta+fcd(4)*dble(z)+
     c           fcd(5)*thetaz+fcd(7)*theta2+fcd(16)*theta2z))
               call comperr(6,xh,ccd,vh)
               errt=real(srcd*dsqrt(vh))
               eindex=errt*(findex-ctcd)
c If theta is between 0.75 and 0.8 we interpolate between this value and the 
c constant value (for 0.6-0.75) 
            else
               dumt=0.8D+00
               findex=ctcd+real(dexp(fcd(1)+fcd(2)*dumt+fcd(4)*dble(z)+
     c           fcd(5)*dumt*dble(z)+fcd(7)*dumt*dumt+fcd(16)*dumt*dumt*
     c           dble(z)))
               call comperr(6,xh2,ccd,vh)
               errt=real(srcd*dsqrt(vh))
               eindex=errt*(findex-ctcd)
               x=(thet-0.75)/0.05
               findex=real(fw1)*(1.-x)+findex*x
               eindex=errw*(1.-x)+eindex*x
            end if
         else
c Theta=(1.08,1.84) for dwarfs
            findex=ctd+real(dexp(fd(1)+fd(2)*theta+fd(7)*theta2))
            xh(3)=xh(5)
            call comperr(3,xh,cd,vh)
            errt=real(srd*dsqrt(vh))
            eindex=errt*(findex-ctd)
            if(g.lt.4.) iflag=3
            if(thet.gt.1.84) iflag=2
         end if
         if(g.ge.3.5) then
            return
         else
c Here we interpolate between this value and the one computed for giants 
c (this applies when log g is between 3.0 and 3.5):
            x=(g-3.)/0.5
            findex=findex0*(1.-x)+findex*x
            eindex=eindex0*(1.-x)+eindex*x
            return
         end if
      end if
c Error: No value was computed
      iflag=-5
      return

      end
c
c
c
      subroutine comperr(n,x,c,v)
      implicit none

      integer n
      double precision x(25),v,c(25,25)

      integer i,j
!     integer n1,n2,n3,n4
      double precision x1(25,25),x2(25,25),temp(25,25),fin(25,25)

c First it reconstructs the whole v-c matrix
      
      do i=2,n
         do j=1,i-1
            c(i,j)=c(j,i)
         end do
      end do

      do j=1,n
         x1(1,j)=x(j)
         x2(j,1)=x(j)
      end do

      call multmatrix(x1,c,temp,1,n,n,n)

      call multmatrix(temp,x2,fin,1,n,n,1)
      v=fin(1,1)

      return
      end

c Subroutine to multiply matrices
c Input: a -> matrix n1 x n2
c        b -> matrix n3 x n4
c n2 must be equal to n3
c Input : n1,n2,n3,n4 (integers, dimensions of matrices)
c Output: c -> matrix n1 x n4
c In the main program matrices must be defined as double precision with 
c   dimensions 25x25
 
      subroutine multmatrix(a,b,c,n1,n2,n3,n4)
      implicit none

      integer n1,n2,n3,n4
      integer i,j,k
      double precision a(25,25),b(25,25),c(25,25)

      if(n2.ne.n3) then
         write(6,*) 'Error in matrix dimensions'
         stop
      end if
      
      do i=1,n1
         do j=1,n4
            c(i,j)=0.D0
            do k=1,n2
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
         end do
      end do

      return
      end
