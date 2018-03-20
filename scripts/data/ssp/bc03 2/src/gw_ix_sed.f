	FUNCTION GW_IX_SED(J,X,F,N,IGW)

c	Computes value of Guy Worthey index J directly from sed
c	Follows procedure outlined by Trager at el. (1998, ApJS 116, 1)
c	See their eq. 1-3, pag. 5.

c	See comments at the end of this file.

c	J = Index identification
c	X = wavelength scale
c	F = flux scale (sed)
c	N = dimension of X, F
c	GW_IX = index measured directly from sed (no formula)

c	Array declaration
	parameter (jid=35)
	character*40 idxlbl
	real x(n),f(n)
	real w(6),a(20),b(20),c(3),d(3),u(10000),z(10000)
        common /gwparam/ im,w,a,b,c,d,idxlbl
c	Common to store fluxes
	integer iim(jid)
	real ffb(jid),ffr(jid),ffc(jid),ffl(jid)
	common /fluxes/ iim,ffb,ffr,ffc,ffl

c	Equation of straight line joining (wb,fb) and (wr,fr)
	fc(v) = ( (wr-v)*fb + (v-wb)*fr ) / (wr-wb)

	if (j.eq.31) then
c		Return B4000VN index
		gw_ix_sed = b4000vn(x,f,n)
	elseif (j.eq.30) then
c		Returns Gorgas and Cardiel D4000 index
		gw_ix_sed = gc_d4000_obs(x,f,n)
	elseif (j.eq.35) then
c		Returns Brodie and Hanes HK index
		gw_ix_sed = bh_ix_sed(6,x,f,n)
	else
c		Find wavelength points everytime to allow calls to this routine in
c		the same program using different wavelength arrays.
c		Get Parameters
		if (igw.eq.0) then
			call gw_etal_dat(j)
		else
			call gw_thesis_dat(j)
		endif

c       	Compute mean height inside blue pseudo continuum band
        	fb=trapz3(x,f,n,w(1),w(2),ierr)/(w(2)-w(1))
		if (ierr.ne.0) then
			gw_ix_sed = 999.
			return
		endif
		wb=(w(1)+w(2))/2.
c		write (6,*)  j,wb,fb
c		write (66,*) j,wb,fb

c       	Compute mean height inside red pseudo continuum band
       		fr=trapz3(x,f,n,w(3),w(4),ierr)/(w(4)-w(3))
		if (ierr.ne.0) then
			gw_ix_sed = 888.
			return
		endif
		wr=(w(3)+w(4))/2.
c		write (6,*)  j,wr,fr
c		write (66,*) j,wr,fr

c		Compute equivalent width according to eq. (2) of
c		Trager at el. (1998, ApJS 116, 1)
c       	Locate central band in array x: w(5) and w(6)
        	call locate (x,n,w(5),i1)
        	call locate (x,n,w(6),i2)
c		Fill in auxiliary array from i1-5 to i2+5 with ratio f/fc
		nu=0
		do i=i1-5,i2+5
		nu=nu+1
		u(nu)=x(i)
		z(nu)=f(i)/fc(x(i))
c 		write (6,*)  j,i,nu,u(nu),fc(x(i)),z(nu)
c		write (66,*) j,i,nu,u(nu),fc(x(i)),z(nu)
		enddo
c       	Compute integral of ratio f(i)/fc(x(i)) from w(5) to w(6)
       		fm=trapz3(u,z,nu,w(5),w(6),ierr)
		if (ierr.ne.0) then
			gw_ix_sed = 777.
			return
		endif

c		Compute index
		if (im.eq.0) then
			gw_ix_sed = w(6) - w(5) - fm
		else
			fm=fm / (w(6) - w(5))
			gw_ix_sed = -2.5*alog10(fm)
			fm=fm * (w(6) - w(5))
		endif

c		Store flux values:
		iim(j)=im
		ffb(j)=fb
		ffr(j)=fr
c		Compute central pseudo continuum flux fc = gw_fc
c		according to G Worthey definition
		wc=(w(5)+w(6))/2.
		cc=(wc-wb)/(wr-wb)
		ffc(j)=(1.-cc)*fb + cc*fr
		ffl(j)=trapz3(x,f,n,w(5),w(6),ierr)
c		Implemented by SC 1/10/03: the above appears to
c		correspond to the old Burstein et al.definition
c	        Eqs (2) and (3) of Worthey et al. (1994) do perform
c		the integration. Store equivalent continuum level
c		write (*,*) j,im,ffc(j),ffl(j)/fm
		ffc(j)=ffl(j)/fm
	endif
	return
	end

c	From gustavo@MPA-Garching.MPG.DE Wed Oct  2 16:20:24 2002
c	Date: Thu, 3 Oct 2002 10:12:54 +0200
c	From: Gustavo Bruzual <gustavo@MPA-Garching.MPG.DE>
c	To: bruzual@cida.ve
c	Subject: [charlot@MPA-Garching.MPG.DE: 2 updates]
c	
c	Forwarded message from Stephane Charlot <charlot@MPA-Garching.MPG.DE>
c	
c	Date: Wed, 2 Oct 2002 02:32:22 +0200
c	From: Stephane Charlot <charlot@MPA-Garching.MPG.DE>
c	To: bruzual@cida.ve
c	Subject: 2 updates
c	User-Agent: Mutt/1.2i
c	
c	Hi Gustavo!
c	You must be arriving just about now!
c	I am working on this last analysis about line indices. I realized
c	the fluxes you stored are from the old Burstein definition, while
c	Worthey et al. 1994 (eq. 2-3) do include an integral. Of course
c	we cannot store the integrls, but I stored something which is
c	closer to it than what you had, for the continuum. I attach
c	the updated gw_ix_sed.f routine.
c	I also attach the updated vel_disp.f program, in which I also
c	store these fluxes in .8lsindx_sed_fluxes
c	
c	cheers!!
c	stephane
