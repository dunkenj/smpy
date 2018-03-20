	SUBROUTINE CONVOLVE_TX(Y,AGE,BOL,STR,EVF,SNR,PNR,BH,SN,WD,RM,GS,GL,DR,DP)

!	Modified by G. Bruzual on 01-Aug-2013
!	In this version of the routine the AGE in yr is entered in the arguments
!	The values of IC is found below

c	Written by G. Bruzual on 25-March-1999

c	Uses different approach from previous versions. The SFR is
c	used as given (from t=0 to t=age) and the SSP spectra are
c	interpolated at required age. This allows for arbitrarily
c	short bursts to be described correctly.

c	Computes sed at age tb(k) by performing convolution integral
c	of the sed for a SSP and the chosen SFR using trapezoidal rule

c	Convolution integral according to trapezoidal rule

c	BOL = convolved bolometric flux at age tb(k)
c	STR = convolved total mass in stars at age tb(k)
c	EVF = convolved evolutionary flux at age tb(k)
c	SNR = convolved Super Nova rate at age tb(k)
c	PNR = convolved PN birth rate at age tb(k)
c	BH  = convolved total number of BH at age tb(k)
c	SN  = convolved total number of NS at age tb(k)
c	WD  = convolved total number of WD at age tb(k)
c	RM  = convolved total mass in Remnants at age tb(k)

c	Array declaration
	include 'jb.dec'
	include 'csp.dec'
	real y(imw),taux(1000000)
	real fxu(0:jts),fxb(0:jts),fxi(0:jts),fxk(0:jts)
	real ufwa,bfwa,ifwa,kfwa,mwa,ufwla,bfwla,ifwla,kfwla,mwla
	common /w_ages/ fxu,fxb,fxi,fxk,ufwa,bfwa,ifwa,kfwa,mwa,ufwla,bfwla,ifwla,kfwla,mwla

c	Modified by G. Bruzual on 24/9/93 to include common /totgas/ where
c	the total mass in gas (new + recycled) vs time is stored. This
c	quantity is used in models with SFR which recycles gas.
	include 'recycle.dec'

c	Clear buffers
	bol=0.
	str=0.
	evf=0.
	snr=0.
	pnr=0.
	bh=0.
	sn=0.
	wd=0.
	rm=0.
	gs=0.
	gl=0.
	ufwa=0.
	bfwa=0.
	ifwa=0.
	kfwa=0.
	mwa=0.
	duwa=0.
	dbwa=0.
	diwa=0.
	dkwa=0.
	dmwa=0.
	do j=1,inw
	y(j)=0.
	enddo

c	Build time scale for accurate integration
	kt=0
!	age=tb(ic)		! age entered as a calling argument
!	Find corresponding value of ic
	ic=0
	do while (tb(ic) < age)
	ic=ic+1
	enddo
!	Preform convolution integral at t=age
	do k=0,ic
	kt=kt+1
	taux(kt)=max(0.,age-tb(k))
	kt=kt+1
	taux(kt)=tb(k)
	enddo
c	Add more time steps under option 7 (user sfr entered in table)
	if (io.eq.7) then
		do i=1,isfrc
		if (age.ge.time(i)) then
			kt=kt+1
			taux(kt)=age-time(i)
		endif
		kt=kt+1
		taux(kt)=time(i)
		enddo
	endif
c       Add more steps if tcut.lt.20.E9
	if (age.gt.tcut.and.tcut.lt.20.E9) then
		do k=0,ic
		if (tb(k).le.tcut) last=k
		enddo
		tlast=tb(last)
		tnext=tb(last+1)
        	dt=0.0005*tcut
        	do i=1,100000
	        tlast=tlast+dt
		if (tlast.gt.tnext) goto 2
		kt=kt+1
		taux(kt)=age-tlast
		kt=kt+1
		taux(kt)=tlast
		enddo
	endif
2	call sort(kt,taux)
c	Suppress duplicate entries below tb(ic)
	kn=1
	do k=2,kt
	if (taux(k).le.tb(ic).and.taux(k).gt.taux(k-1)) then
		kn=kn+1
		taux(kn)=taux(k)
	endif
	enddo
	if (kn.eq.1) return
c	Suppress duplicate entries
	call usort(kn,taux)
c	Write in reverse order
	call xreverse(taux,kn)

c	Compute weights to assign to each sed and compute convolution integral
c	do k=0,ic
	do k=1,kn
c	sr=sfr(tb(k))
	sr=sfr(max(0.,age-taux(k)))
	if (sr.gt.0.) then
		if (k.eq.1) then
c			dt=tb(k+1)-tb(k)
			dt=taux(1)-taux(2)
		elseif (k.eq.kn) then
c			dt=tb(ic)-tb(ic-1)
			dt=taux(kn-1)-taux(kn)
		else
c			dt=tb(k+1)-tb(k-1)
			dt=taux(k-1)-taux(k+1)
		endif
		w=sr*dt/2.
		if (io.eq.0) w=1.
c		tx=age-tb(k)
		tx=taux(k)
		do j1=0,ic
		if (tx.eq.tb(j1)) then
c			Used sed and other quantities at tx = age-tb(k)
			do m=1,inw
			y(m)=y(m)+w*fl(m,j1)
			enddo
			bol=bol+w*bflx(j1)
			str=str+w*strm(j1)
			evf=evf+w*evfl(j1)
			snr=snr+w*snbr(j1)
			pnr=pnr+w*pnbr(j1)
			sn =sn +w*sntn(j1)
			bh =bh +w*bhtn(j1)
			wd =wd +w*wdtn(j1)
			rm =rm +w*rmtm(j1)
			gs =gs +w*gasms(j1)
			gl =gl +w*galms(j1)
			ufwa = ufwa + w*fxu(j1)*tb(j1)
			bfwa = bfwa + w*fxb(j1)*tb(j1)
			ifwa = ifwa + w*fxi(j1)*tb(j1)
			kfwa = kfwa + w*fxk(j1)*tb(j1)
			mwa  = mwa  + w*tb(j1)
			duwa = duwa + w*fxu(j1)
			dbwa = dbwa + w*fxb(j1)
			diwa = diwa + w*fxi(j1)
			dkwa = dkwa + w*fxk(j1)
			dmwa = dmwa + w
			goto 1
		endif
		if (tx.ge.tb(j1).and.tx.lt.tb(j1+1)) then
			j2=j1+1
			if (tb(j1).gt.0.) then
				a=alog10(tx/tb(j1))/alog10(tb(j2)/tb(j1))
			else
				a=tx/tb(j2)
			endif
			b=1.-a
			wa=w*a
			wb=w*b
c			Interpolate sed and other quantities at tx = age-tb(k)
			do m=1,inw
			y(m)=y(m)+wb*fl(m,j1)+wa*fl(m,j2)
			enddo
			bol=bol+wb*bflx(j1)+wa*bflx(j2)
			str=str+wb*strm(j1)+wa*strm(j2)
			evf=evf+wb*evfl(j1)+wa*evfl(j2)
			snr=snr+wb*snbr(j1)+wa*snbr(j2)
			pnr=pnr+wb*pnbr(j1)+wa*pnbr(j2)
			sn =sn +wb*sntn(j1)+wa*sntn(j2)
			bh =bh +wb*bhtn(j1)+wa*bhtn(j2)
			wd =wd +wb*wdtn(j1)+wa*wdtn(j2)
			rm =rm +wb*rmtm(j1)+wa*rmtm(j2)
			gs =gs +wb*gasms(j1)+wa*gasms(j2)
			gl =gl +wb*galms(j1)+wa*galms(j2)
			ufwa = ufwa + wb*fxu(j1)*tb(j1) + wa*fxu(j2)*tb(j2)
			bfwa = bfwa + wb*fxb(j1)*tb(j1) + wa*fxb(j2)*tb(j2)
			ifwa = ifwa + wb*fxi(j1)*tb(j1) + wa*fxi(j2)*tb(j2)
			kfwa = kfwa + wb*fxk(j1)*tb(j1) + wa*fxk(j2)*tb(j2)
			mwa  = mwa  + wb*tb(j1)         + wa*tb(j2)
			duwa = duwa + wb*fxu(j1)        + wa*fxu(j2)        
			dbwa = dbwa + wb*fxb(j1)        + wa*fxb(j2)        
			diwa = diwa + wb*fxi(j1)        + wa*fxi(j2)        
			dkwa = dkwa + wb*fxk(j1)        + wa*fxk(j2)        
			dmwa = dmwa + wb                + wa
			goto 1
		endif
		enddo
		write (6,*) 'Problem interpolating tx,j1 = ',tx,j1
		stop
	endif
1	continue
	enddo

c	Compute flux weighted age returned in common
	ufwa = ufwa/duwa
	bfwa = bfwa/dbwa
	ifwa = ifwa/diwa
	kfwa = kfwa/dkwa
	mwa  = mwa /dmwa

c	Store amount of processed gas
	if (io.eq.1) then
		str=amin1(1.,str)
c		unprocessed gas so far
		ugas=exp(-age/tau)
c		processed gas = gas formed into stars - mass in stars - remnants
		prgas = 1. - ugas - str - rm
		if (prgas.lt.0.) prgas=0.
c		lgas=ic+1
		lgas=lgas+1
		so(lgas)=prgas
		to(lgas)=age
	endif
	return
	end
