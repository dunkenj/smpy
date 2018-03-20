      FUNCTION CHABIMF(ALO,AUP,M1,M2,STARMASS)

c     Returns number of stars in range from M1 to M2 for Chabrier (2001) IMF
c     Returns the mass in stars from M1 to M2
c     ALO and AUP are the lower and upper cutoffs of the IMF, which
c     is normalized to 1 Mo over this range

      implicit real (A-H,O-Z)
      real*8 m1,m2,starmass

      common /fmlog/xml,sigmac,anormlog
      common /fmlin/anormlin

*******************************************************
*                 - Adjustable parameters -
*
*    Lower and upper masses for the entire mass
*    integration of the IMF: alo, aup
*    Index "alpha=1+x" in the power-law part > 1 Msol
*           dn/dm propto m**(-alpha)
      alpha=2.3
*******************************************************

      xml=-1.1
      sigmac=0.69
      call chabint(alo,aup,alpha)    ! calculates the normalization
      al1=1.-alpha
      al2=2.-alpha
      anorm=anormlin/log(10.)
      v1log=anormlog*exp(-xml**2/2./sigmac**2)
c      print*,'normlog, normlin',anormlog,anormlin
c      print*,'norm @ 1 Msol',v1log,anormlin
      if(v1log.ne.anormlin) then
                            print*,'normalization pb in routine chabimf!'
                            stop
      endif

	if (m2.le.m1) then
		chabimf=0.
		starmass=0.
		return
	endif

      IF(m1.LE.1.0) THEN
            if(m2.LE.1.0) then
                vnum=chabnum(m1,m2)     ! number-density
                vmas=chabmas(m1,m2)     ! mass-density
            else
                alim=1.0
                vnum=chabnum(m1,dble(alim))+anorm/al1*(m2**al1-1.0)
                vmas=chabmas(m1,dble(alim))+anorm/al2*(m2**al2-1.0)
            endif
      ELSE
                vnum=anorm/al1*(m2**al1-m1**al1)
                vmas=anorm/al2*(m2**al2-m1**al2)
      ENDIF

      chabimf=vnum
      starmass=vmas

      return
      end

      SUBROUTINE CHABINT(alo,aup,alpha)
**************************************************************************
**    chabmas calculates the normalization such that
**        [int m.(dn/dm)dm = 1]  from m=alo to m=aup
**************************************************************************
      implicit real (A-H,O-Z)
      common /fmlog/xml,sigmac,anormlog
      common /fmlin/anormlin

      pi=4.*atan(1.)
      cte2=exp(-xml**2/sigmac**2/2.)
      al2=2.-alpha
      clin=cte2/al2*(aup**al2-1.)
      if(aup.LE.1.) clin=0.

      sigmacp=sigmac*log(10.)
      xmlp=xml*log(10.)
      cte1=exp(xmlp+sigmacp**2/2.)
      x1=-(xmlp+sigmacp**2)/sqrt(2.)/sigmacp
      xlo=(log10(alo)*log(10.)-(xmlp+sigmacp**2))/sqrt(2.)/sigmacp
      clog=cte1*sqrt(2.*pi)/2.*sigmacp
     !       *  (erf(x1) - erf(xlo))

      anorm=1./(clog+clin)
      anormlog=anorm*log(10.)     ! cte de la partie lognormal
      anormlin=anormlog*cte2         ! cte de la partie loi de puissance

      return
      end

      FUNCTION CHABNUM(xmin,xmax)
**************************************************************************
**                         *    IMF Chabrier   *
**    chabnum returns [int (dn/dm)dm]  from m=xmin to m=xmax
*************************************************************************
      implicit real (A-H,O-Z)
	real*8 xmin,xmax
c	real*8 y1,y2,de
	include 'chabimf.dec'
      common /fmlog/xml,sigmac,anormlog

c	check mass range
	if (xmax.le.xmin) then
		chabnum=0.
		return
	else
		pi=4.*atan(1.)
c		Modified by GBA, April 19, 2003
c		chabnum=anormlog * sqrt(2.*pi)/2.*sigmac
c    !		*( erf((log10(xmax)-xml)/sqrt(2.)/sigmac )
c    !		-  erf((log10(xmin)-xml)/sqrt(2.)/sigmac ) )
		y1 =  (dlog10(xmax*1.D00)-xml)/sqrt(2.)/sigmac
		y2 =  (dlog10(xmin*1.D00)-xml)/sqrt(2.)/sigmac
		de =  erf(y1) - erf(y2)
		chabnum=anormlog * sqrt(2.*pi)/2.*sigmac * de
		return
	endif
	end

      FUNCTION CHABMAS(xmin,xmax)
**************************************************************************
**                         *    IMF Chabrier   *
**    chabmas returns [int m.(dn/dm)dm]  from m=xmin to m=xmax
*************************************************************************
      implicit real (A-H,O-Z)
	real*8 xmin,xmax
c	real*8 y1,y2,de
	include 'chabimf.dec'
      common /fmlog/xml,sigmac,anormlog

c	Check mass range
	if (xmax.le.xmin) then
		chabmas=0.
		return
	else
		pi=4.*atan(1.)
		sigmacp=log(10.)*sigmac
		xmlp=xml*log(10.)
		anormp=anormlog/log(10.)
		cte=exp(xmlp+sigmacp**2/2.)
c		Modified by GBA, April 19, 2003
c		x1=(log10(xmin)*log(10.)-(xmlp+sigmacp**2))/sqrt(2.)/sigmacp
c		x2=(log10(xmax)*log(10.)-(xmlp+sigmacp**2))/sqrt(2.)/sigmacp
c		chabmas=anormp*cte*sqrt(2.*pi)/2.*sigmacp * (erf(x2) - erf(x1))
		y1=(dlog10(xmin*1.D00)*dlog(10.D00)-(xmlp+sigmacp**2))/sqrt(2.)/sigmacp
		y2=(dlog10(xmax*1.D00)*dlog(10.D00)-(xmlp+sigmacp**2))/sqrt(2.)/sigmacp
		de=erf(y2) - erf(y1)
		chabmas=anormp * cte * sqrt(2.*pi)/2.*sigmacp * de
		return
	endif
	end
