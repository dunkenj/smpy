	SUBROUTINE EXPAND_TIME_STEPS
 
c	Expands time scale as iequired for model with tcut < 20 Gyr
 
	include 'jb.dec'
	include 'csp.dec'
	real t(0:jts)
	real fxu(0:jts),fxg(0:jts),fxr(0:jts),fxi(0:jts),fxz(0:jts),fxk(0:jts)
	real ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwa,mwla
	common /w_ages/ fxu,fxg,fxr,fxi,fxz,fxk,ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwa,mwla
 
c	Build enlarged time scale
	k=-1
	do i=0,nsteps-1
	if (tb(i).le.tcut) then
		k=k+1
		last=i
		t(k)=tb(i)
		tlast=t(k)
		next=i+1
		tnext=tb(i+1)
	endif
	enddo
	m1=last
 
c	Add desired steps
	dt=0.002E9
	dt=0.002*tcut
	if (dt.lt.2.E6) dt=2.E6
	dt=dt/4.
	if (tlast < tcut) then
		k=k+1
		t(k)=tcut
		tlast=tcut
	endif
	do i=1,100000
	tlast=tlast+dt
	if (tlast.lt.tnext) then
		k=k+1
		t(k)=tlast
	endif
	if (i.eq.20) dt=dt*4.
	enddo
	k1=k+1
 
c	Add remaining steps
	do i=next,nsteps-1
	k=k+1
	t(k)=tb(i)
	enddo
 
c	Shift sed''s and other quantities for i>=next
	ks=nsteps-next
	do n=1,ks
	knew=k-n+1
	kold=nsteps-n
	do i=1,inw
	fl(i,knew)  = fl(i,kold)
	enddo
	bflx(knew)  = bflx(kold)
	strm(knew)  = strm(kold)
	evfl(knew)  = evfl(kold)
	snbr(knew)  = snbr(kold)
	pnbr(knew)  = pnbr(kold)
	bhtn(knew)  = bhtn(kold)
	sntn(knew)  = sntn(kold)
	wdtn(knew)  = wdtn(kold)
	rmtm(knew)  = rmtm(kold)
	gasms(knew) = gasms(kold)
	galms(knew) = galms(kold)
	fxu(knew)   = fxu(kold)
	fxg(knew)   = fxg(kold)
	fxr(knew)   = fxr(kold)
	fxi(knew)   = fxi(kold)
	fxz(knew)   = fxz(kold)
	fxk(knew)   = fxk(kold)
	enddo
	if (ks == 0) knew = k
	m2=knew

c	Interpolate SSP sed''s and other quantities at new time steps
	do n=m1+1,m2-1
	a=(t(m2)-t(n))/(t(m2)-t(m1))
	b=1.-a
	do i=1,inw
	fl(i,n)  = a*fl(i,m1)  + b*fl(i,m2)
	enddo
	bflx(n)  = a*bflx(m1)  + b*bflx(m2)
	strm(n)  = a*strm(m1)  + b*strm(m2)
	evfl(n)  = a*evfl(m1)  + b*evfl(m2)
	snbr(n)  = a*snbr(m1)  + b*snbr(m2)
	pnbr(n)  = a*pnbr(m1)  + b*pnbr(m2)
	bhtn(n)  = a*bhtn(m1)  + b*bhtn(m2)
	sntn(n)  = a*sntn(m1)  + b*sntn(m2)
	wdtn(n)  = a*wdtn(m1)  + b*wdtn(m2)
	rmtm(n)  = a*rmtm(m1)  + b*rmtm(m2)
	gasms(n) = a*gasms(m1) + b*gasms(m2)
	galms(n) = a*galms(m1) + b*galms(m2)
	fxu(n)   = a*fxu(m1)   + b*fxu(m2)
	fxg(n)   = a*fxg(m1)   + b*fxg(m2)
	fxr(n)   = a*fxr(m1)   + b*fxr(m2)
	fxi(n)   = a*fxi(m1)   + b*fxi(m2)
	fxz(n)   = a*fxz(m1)   + b*fxz(m2)
	fxk(n)   = a*fxk(m1)   + b*fxk(m2)
	enddo
 
c	Redefine time scale
	nsteps=k+1
	do i=0,nsteps-1
	tb(i)=t(i)
	enddo
 
	return
	end
 
	SUBROUTINE ADD_TIME_STEPS
 
c	Includes in basic time scale the singular time stpes in usr_sfr
 
	include 'jb.dec'
	include 'csp.dec'
	real t(0:nsfrp)
	real fxu(0:jts),fxg(0:jts),fxr(0:jts),fxi(0:jts),fxz(0:jts),fxk(0:jts)
	real ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwa,mwla
	common /w_ages/ fxu,fxg,fxr,fxi,fxz,fxk,ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwa,mwla
 
c	Verify if any point has to be added
	if (jadd.eq.0) return
	write (6,*) jadd,' time steps to add. Computing spectra at these new time steps.'
	
c	Add desired steps
	do j=1,jadd
	k=-1
	do i=0,nsteps-1
	k=k+1
	t(k)=tb(i)
	if (tb(i).lt.tadd(j).and.tadd(j).lt.tb(i+1)) then
		m1=i
                next=i+1
		ks=nsteps-next
 
c		Add remaining steps
		k=k+1
		t(k)=tadd(j)
        	do n=next,nsteps-1
		k=k+1
        	t(k)=tb(n)
        	enddo
c		if (k.ne.nsteps+1) then
c			write (6,*) 'wrong number of time steps',k,nsteps
c			stop
c		endif
 
c		Shift sed''s and other quantities for i>=next
		do n=1,ks
		knew=k-n+1
		kold=nsteps-n
		do l=1,inw
		fl(l,knew)  = fl(l,kold)
		enddo
		bflx(knew)  = bflx(kold)
		strm(knew)  = strm(kold)
		evfl(knew)  = evfl(kold)
		snbr(knew)  = snbr(kold)
		pnbr(knew)  = pnbr(kold)
		bhtn(knew)  = bhtn(kold)
		sntn(knew)  = sntn(kold)
		wdtn(knew)  = wdtn(kold)
		rmtm(knew)  = rmtm(kold)
		gasms(knew) = gasms(kold)
		galms(knew) = galms(kold)
		fxu(knew)   = fxu(kold)
		fxg(knew)   = fxg(kold)
		fxr(knew)   = fxr(kold)
		fxi(knew)   = fxi(kold)
		fxz(knew)   = fxz(kold)
		fxk(knew)   = fxk(kold)
		enddo
		m2=knew
 
c		Interpolate SSP sed''s and other quantities at new time steps
		n=m1+1
		a=(t(m2)-t(n))/(t(m2)-t(m1))
		b=1.-a
		do l=1,inw
		fl(l,n)  = a*fl(l,m1)  + b*fl(l,m2)
		enddo
		bflx(n)  = a*bflx(m1)  + b*bflx(m2)
		strm(n)  = a*strm(m1)  + b*strm(m2)
		evfl(n)  = a*evfl(m1)  + b*evfl(m2)
		snbr(n)  = a*snbr(m1)  + b*snbr(m2)
		pnbr(n)  = a*pnbr(m1)  + b*pnbr(m2)
		bhtn(n)  = a*bhtn(m1)  + b*bhtn(m2)
		sntn(n)  = a*sntn(m1)  + b*sntn(m2)
		wdtn(n)  = a*wdtn(m1)  + b*wdtn(m2)
		rmtm(n)  = a*rmtm(m1)  + b*rmtm(m2)
		gasms(n) = a*gasms(m1) + b*gasms(m2)
		galms(n) = a*galms(m1) + b*galms(m2)
		fxu(n)   = a*fxu(m1)   + b*fxu(m2)
		fxg(n)   = a*fxg(m1)   + b*fxg(m2)
		fxr(n)   = a*fxr(m1)   + b*fxr(m2)
		fxi(n)   = a*fxi(m1)   + b*fxi(m2)
		fxz(n)   = a*fxz(m1)   + b*fxz(m2)
		fxk(n)   = a*fxk(m1)   + b*fxk(m2)
 
c		Redefine time scale
		nsteps=nsteps+1
		do l=0,nsteps-1
		tb(l)=t(l)
		enddo
	endif
	enddo
	enddo
	return
	end
