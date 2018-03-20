	SUBROUTINE GW_INDICES(T,FI,FC,LUN,IGW)

c	Computes spectral indices from the fluxes FI and FC using fitting functions
c	Worthey and Ottaviani indices added June 2002, indices 22 to 25.
c	Writes down the answers

c	Array declarations
	character*40 idxlbl
	parameter (jid=26)
	real gwix(jid),w(6),a(20),b(20),c(3),d(3)
	real fi(jid),fc(jid),dw(jid)
	integer ix(jid)
	common /gwparam/ im,w,a,b,c,d,idxlbl
	data igwcall/0/,dw/jid*0./,ix/jid*0/

	if (t.eq.0.) return

c	Fill arrays
	if (igwcall.eq.0) then
		do j=1,25
		if (igw.eq.0) then
			call gw_etal_dat(j)
		else
			call gw_thesis_dat(j)
		endif
		dw(j)=w(6)-w(5)
		ix(j)=im
		enddo
		igwcall=1
	endif

c	Compute indices (G Worthey)
	do j=1,25
	if (ix(j).eq.0) then
		gwix(j)=dw(j)*(1.-fi(j)/fc(j))
	else
		
		gwix(j)=-2.5*alog10(fi(j)/fc(j))
	endif
c	if (j.eq.1) write (49,*) tl,fi(j),fc(j),fi(j)/fc(j)
	enddo

c	Compute 4000 A break, according to Gorgas and Cardiel procedure
c	d4000=fc(jid)/fi(jid)
	gwix(jid)=fc(jid)/fi(jid)

	tl=alog10(t)
!	write (lun+6,100) tl,(gwix(j),j= 1, 13)
	write (lun+6,100) tl,(gwix(j),j= 1, 21)
!	write (lun+7,101) tl,(gwix(j),j=14,jid)
	write (lun+7,101) tl,(gwix(j),j=22,jid)
100	format (f10.6,21f9.4)
!101	format (f10.6,8f8.3,5f12.3)
101	format (f10.6,1x,5f11.4)
	return
	end
