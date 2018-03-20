	SUBROUTINE CONVOLVE_NEW(Y,IC,BOL,STR,EVF,SNR,PNR,BH,SN,WD,RM,GS,GL,MR,DP)

!	Written by G. Bruzual on 25-March-1999

!	Uses different approach from previous versions. The SFR is
!	used as given (from t=0 to t=age) and the SSP spectra are
!	interpolated at required age. This allows for arbitrarily
!	short bursts to be described correctly.

!	Computes sed at age tb(k) by performing convolution integral
!	of the sed for a SSP and the chosen SFR using trapezoidal rule

!	Convolution integral according to trapezoidal rule

!	BOL = convolved bolometric flux at age tb(k)
!	STR = convolved total mass in stars at age tb(k)
!	EVF = convolved evolutionary flux at age tb(k)
!	SNR = convolved Super Nova rate at age tb(k)
!	PNR = convolved PN birth rate at age tb(k)
!	BH  = convolved total number of BH at age tb(k)
!	SN  = convolved total number of NS at age tb(k)
!	WD  = convolved total number of WD at age tb(k)
!	RM  = convolved total mass in Remnants at age tb(k)
!	MLR = mass loss rate
!	DPR = dust production rate

!	Array declaration
	include 'jb.dec'
	include 'csp.dec'
	real y(imw),taux(100000)
	real fxu(0:jts),fxg(0:jts),fxr(0:jts),fxi(0:jts),fxz(0:jts),fxk(0:jts)
	real ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwa,mwla
	common /w_ages/ fxu,fxg,fxr,fxi,fxz,fxk,ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwa,mwla

!	Modified by G. Bruzual on 24/9/93 to include common /totgas/ where
!	the total mass in gas (new + recycled) vs time is stored. This
!	quantity is used in models with SFR which recycles gas.
	include 'recycle.dec'

!	Clear buffers
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
	mr=0.
	dp=0.
	ufwa=0.
	gfwa=0.
	rfwa=0.
	ifwa=0.
	zfwa=0.
	kfwa=0.
	ufwla=0.
	gfwla=0.
	rfwla=0.
	ifwla=0.
	zfwla=0.
	kfwla=0.
	mwa=0.
	mwla=0.
	duwa=0.
	dgwa=0.
	drwa=0.
	diwa=0.
	dzwa=0.
	dkwa=0.
	dmwa=0.
	do j=1,inw
	y(j)=0.
	enddo

!	Build time scale for accurate integration
	kt=0
	age=tb(ic)
	do k=0,ic
	kt=kt+1
	taux(kt)=age-tb(k)
	kt=kt+1
	taux(kt)=tb(k)
	enddo
!	Add more time steps under option 7 (user sfr entered in table)
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
!       Add more steps if tcut.lt.20.E9
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
!	Suppress duplicate entries below tb(ic)
	kn=1
	do k=2,kt
	if (taux(k).le.tb(ic).and.taux(k).gt.taux(k-1)) then
		kn=kn+1
		taux(kn)=taux(k)
	endif
	enddo
	if (kn.eq.1) return
!	Suppress duplicate entries
	call usort(kn,taux)
!	Write in reverse order
	call xreverse(taux,kn)

!	Compute weights to assign to each sed and compute convolution integral
!	do k=0,ic
	do k=1,kn
!	sr=sfr(tb(k))
	sr=sfr(age-taux(k))
	if (sr.gt.0.) then
		if (k.eq.1) then
!			dt=tb(k+1)-tb(k)
			dt=taux(1)-taux(2)
		elseif (k.eq.kn) then
!			dt=tb(ic)-tb(ic-1)
			dt=taux(kn-1)-taux(kn)
		else
!			dt=tb(k+1)-tb(k-1)
			dt=taux(k-1)-taux(k+1)
		endif
		w=sr*dt/2.
		if (io.eq.0) w=1.
!		tx=age-tb(k)
		tx=taux(k)
		do j1=0,ic
		if (tx.eq.tb(j1)) then
!			Used sed and other quantities at tx = age-tb(k)
			do m=1,inw
			y(m)=y(m)+w*fl(m,j1)
			enddo
			bol = bol + w*bflx(j1)
			str = str + w*strm(j1)
			evf = evf + w*evfl(j1)
			snr = snr + w*snbr(j1)
			pnr = pnr + w*pnbr(j1)
			sn  = sn  + w*sntn(j1)
			bh  = bh  + w*bhtn(j1)
			wd  = wd  + w*wdtn(j1)
			rm  = rm  + w*rmtm(j1)
			gs  = gs  + w*gasms(j1)
			gl  = gl  + w*galms(j1)
			mr  = mr  + w*rml(j1)
			dp  = dp  + w*rdp(j1)
			ufwa  = ufwa  + w*fxu(j1)*tb(j1)
			gfwa  = gfwa  + w*fxg(j1)*tb(j1)
			rfwa  = rfwa  + w*fxr(j1)*tb(j1)
			ifwa  = ifwa  + w*fxi(j1)*tb(j1)
			zfwa  = zfwa  + w*fxz(j1)*tb(j1)
			kfwa  = kfwa  + w*fxk(j1)*tb(j1)
			ufwla = ufwla + w*fxu(j1)*tlog(tb(j1))
			gfwla = gfwla + w*fxg(j1)*tlog(tb(j1))
			rfwla = rfwla + w*fxr(j1)*tlog(tb(j1))
			ifwla = ifwla + w*fxi(j1)*tlog(tb(j1))
			zfwla = zfwla + w*fxz(j1)*tlog(tb(j1))
			kfwla = kfwla + w*fxk(j1)*tlog(tb(j1))
			mwa   = mwa   + w*tb(j1)
			mwla  = mwla  + w*tlog(tb(j1))
			duwa  = duwa  + w*fxu(j1)
			dgwa  = dgwa  + w*fxg(j1)
			drwa  = drwa  + w*fxr(j1)
			diwa  = diwa  + w*fxi(j1)
			dzwa  = dzwa  + w*fxz(j1)
			dkwa  = dkwa  + w*fxk(j1)
			dmwa  = dmwa  + w
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
!			Interpolate sed and other quantities at tx = age-tb(k)
			do m=1,inw
			y(m)=y(m)+wb*fl(m,j1)+wa*fl(m,j2)
			enddo
			bol = bol + wb*bflx(j1)  + wa*bflx(j2)
			str = str + wb*strm(j1)  + wa*strm(j2)
			evf = evf + wb*evfl(j1)  + wa*evfl(j2)
			snr = snr + wb*snbr(j1)  + wa*snbr(j2)
			pnr = pnr + wb*pnbr(j1)  + wa*pnbr(j2)
			sn  = sn  + wb*sntn(j1)  + wa*sntn(j2)
			bh  = bh  + wb*bhtn(j1)  + wa*bhtn(j2)
			wd  = wd  + wb*wdtn(j1)  + wa*wdtn(j2)
			rm  = rm  + wb*rmtm(j1)  + wa*rmtm(j2)
			gs  = gs  + wb*gasms(j1) + wa*gasms(j2)
			gl  = gl  + wb*galms(j1) + wa*galms(j2)
			mr  = mr  + wb*rml(j1)   + wa*rml(j2)
			dp  = dp  + wb*rdp(j1)   + wa*rdp(j2)
			ufwa  = ufwa  + wb*fxu(j1)*tb(j1) + wa*fxu(j2)*tb(j2)
			gfwa  = gfwa  + wb*fxg(j1)*tb(j1) + wa*fxg(j2)*tb(j2)
			rfwa  = rfwa  + wb*fxr(j1)*tb(j1) + wa*fxr(j2)*tb(j2)
			ifwa  = ifwa  + wb*fxi(j1)*tb(j1) + wa*fxi(j2)*tb(j2)
			zfwa  = zfwa  + wb*fxz(j1)*tb(j1) + wa*fxz(j2)*tb(j2)
			kfwa  = kfwa  + wb*fxk(j1)*tb(j1) + wa*fxk(j2)*tb(j2)
			ufwla = ufwla + wb*fxu(j1)*tlog(tb(j1)) + wa*fxu(j2)*tlog(tb(j2))
			gfwla = gfwla + wb*fxg(j1)*tlog(tb(j1)) + wa*fxg(j2)*tlog(tb(j2))
			rfwla = rfwla + wb*fxr(j1)*tlog(tb(j1)) + wa*fxr(j2)*tlog(tb(j2))
			ifwla = ifwla + wb*fxi(j1)*tlog(tb(j1)) + wa*fxi(j2)*tlog(tb(j2))
			zfwla = zfwla + wb*fxz(j1)*tlog(tb(j1)) + wa*fxz(j2)*tlog(tb(j2))
			kfwla = kfwla + wb*fxk(j1)*tlog(tb(j1)) + wa*fxk(j2)*tlog(tb(j2))
			mwa   = mwa   + wb*tb(j1)               + wa*tb(j2)
			mwla  = mwla  + wb*tlog(tb(j1))         + wa*tlog(tb(j2))
			duwa  = duwa  + wb*fxu(j1) + wa*fxu(j2)        
			dgwa  = dgwa  + wb*fxg(j1) + wa*fxg(j2)        
			drwa  = drwa  + wb*fxr(j1) + wa*fxr(j2)        
			diwa  = diwa  + wb*fxi(j1) + wa*fxi(j2)        
			dzwa  = dzwa  + wb*fxz(j1) + wa*fxz(j2)        
			dkwa  = dkwa  + wb*fxk(j1) + wa*fxk(j2)        
			dmwa  = dmwa  + wb + wa
			goto 1
		endif
		enddo
		write (6,*) 'Problem interpolating tx,j1 = ',tx,j1
		stop
	endif
1	continue
	enddo

!	Compute flux weighted age and log age returned in common
	ufwa  = ufwa/duwa
	gfwa  = gfwa/dgwa
	rfwa  = rfwa/drwa
	ifwa  = ifwa/diwa
	zfwa  = zfwa/dzwa
	kfwa  = kfwa/dkwa
	mwa   = mwa /dmwa
	ufwla = ufwla/duwa
	gfwla = gfwla/dgwa
	rfwla = rfwla/drwa
	ifwla = ifwla/diwa
	zfwla = zfwla/dzwa
	kfwla = kfwla/dkwa
	mwla  = mwla /dmwa

!	Store amount of processed gas
	if (io.eq.1) then
		str=amin1(1.,str)
!		unprocessed gas so far
		ugas=exp(-age/tau)
!		processed gas = gas formed into stars - mass in stars - remnants
		prgas = 1. - ugas - str - rm
		if (prgas.lt.0.) prgas=0.
!		lgas=ic+1
		lgas=lgas+1
		so(lgas)=prgas
		to(lgas)=age
	endif
	return
	end

	FUNCTION TLOG(T)

!	avoids log of 0

	if (t <= 0.) then
		tlog = 0.
	else
		tlog = alog10(t)
	endif
	return
	end
