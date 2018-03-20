	SUBROUTINE SFR_0_B(IDEF,ZOBS)

!	Selects parameters for SFR

	include 'jb.dec'
	include 'csp.dec'
	character ans*1,aux*128,bux*128

!	Ask for selection
	io=0
	tcut=20.E9
	isingle=0.
	if (idef.ne.0) then
		write (6,*)
		write (6,'(x,a)') 'Choose SFR: 0 = SSP (Delta Burst = zero length burst)'
		write (6,'(x,a)') '            1 = Exponential (enter Tau)'
		write (6,'(x,a)') '           -1 = Exponential (enter mu_SFR parameter)'
		write (6,'(x,a)') '            2 = Single Burst of finite length'
		write (6,'(x,a)') '            3 = Constant'
		write (6,'(x,a)') '            4 = Delayed'
		write (6,'(x,a)') '            5 = Linearly decreasing'
		write (6,'(x,a)') '            6 = Read SFR(t) from ASCII file'
		write (6,'(x,a)') '            7 = Single or double exponential + burst (after Chen et al.)'
20		write (6,'(x,a,$)') '   Choice = '
		read (5,'(i10)',err=20,end=6) io
		if (io.eq.7) io=9
		if (io.eq.6) io=7
		if (io.eq.5) io=8
		if (io.eq.4) io=6
	endif
	if (io.eq.1) then
!		Exponential SFR (enter tau)
1		write (6,'(x,a,$)') 'Exponential with e-folding time TAU (Gyr) = '
		read (5,'(f10.0)',err=1,end=6) tau
		tmu=1.-exp(-1./tau)
		write (6,'(x,a,f10.4,5x,a,f10.4)') 'mu_SFR =',tmu,'tau =',tau
		tau=1.e9*tau
		write (6,'(x,a,$)') 'Recycle gas ejected by stars Y/[N] = '
		read (5,'(a)',end=6) ans
		if (ans.eq.'Y'.or.ans.eq.'y') then
			call get_epsilon(io)
		else
			io=4
		endif
	elseif (io.eq.-1) then
!		Exponential SFR (enter mu)
2		write (6,'(x,a,$)') 'Exponential with mu_SFR parameter = '
		read (5,'(f10.0)',err=2,end=6) tmu
		tau=-1.e9/alog(1.-tmu)
		write (6,'(x,a,f10.4,5x,a,f10.4)') 'mu_SFR =',tmu,'tau =',tau*1.E-9
		write (6,'(x,a,$)') 'Recycle gas ejected by stars Y/[N] = '
		read (5,'(a)',end=6) ans
		if (ans.eq.'Y'.or.ans.eq.'y') then
			io=1
			call get_epsilon(io)
		else
			io=4
		endif
	elseif (io.eq.2) then
!		Long Burst SFR (enter duration of burst)
3		write (6,'(x,a,$)') 'Duration of burst (Gyr) = '
		read (5,'(f10.0)',err=3,end=6) tau
		tau=1.e9*tau
		tcut=tau
	elseif (io.eq.3) then
!		Constant SFR (enter constant value)
		tau=1.
5		write (6,'(x,a,$)') 'Enter SFR in Mo/yr [1] = '
		read (5,'(f10.0)',err=5,end=6) str
		if (str.gt.0.) tau=str
	elseif (io.eq.6) then
!		Delayed SFR (enter time at which maximum SFR occurs)
 		write (6,'(x,a  )') 'Delayed SFR as defined in Bruzual (1983)'
7		write (6,'(x,a,$)') 'Maximum in SFR at time TAU (Gyr) = '
		read (5,'(f10.0)',err=7,end=6) tau
		tau=1.e9*tau
	elseif (io.eq.7) then
!		Function SFR(t) read 2 column from ASCII file
 		write (6,'(x,a)') 'SFR(t) read from disk file.'
 		write (6,'(x,a)') 't(yr), SFR(t) in Mo/yr read from 2 column ascii file.'
 		write (6,'(x,a)') 'Linear interpolation between listed points.'
 		write (6,'(x,a)') 'Skip lines starting with #'
8 		write (6,'(x,a,$)') 'SFR in file name = '
		tau=usrsfr(-1.)
		if (tau.lt.0) goto 8
	elseif (io.eq.8) then
!		Linear SFR (enter time at which SFR  = 0)
 		write (6,'(x,a  )') 'Linearly decreasing SFR'
9		write (6,'(x,a,$)') 'SFR = 0 at time TAU (Gyr) = '
		read (5,'(f10.0)',err=9,end=6) tau
		tau=1.e9*tau
	elseif (io.eq.9) then
		write (6,*)
10		write (6,'(x,a,$)') 'For t <= Tj, use exponential SFR with e-folding time TAU_1 (Gyr)   = '
		read (5,'(f20.0)',err=10,end=6) tau1
		tau1 = 1.e9*tau1
		con1 = 1./tau1
11		write (6,'(x,a,$)') 'For t >  Tj, use exponential SFR with e-folding time TAU_2 (Gyr)   = '
		read (5,'(f20.0)',err=11,end=6) tau2
		if (tau2 /= 0.) then
			tau2 = 1.e9*tau2
!12			write (6,'(x,a,$)') '                  Join continuously both SFR''s at time T_J (Gyr)   = '
12			write (6,'(x,a,$)') '        Join continuously both SFR''s at look back time T_J (Gyr)   = '
			read (5,'(f20.0)',err=12,end=6) tauj
			tauj = 1.e9*tauj
		endif
13		write (6,'(x,a,$)') '   Look back time for this galaxy from z = 0, a.k.a. Tform (Gyr)   = '
		read (5,'(f20.0)',err=13,end=6) tauf
		tauf = 1.e9*tauf
		if (tau2 /= 0.) then
!		    Next instruction added at Alfredo and Gladis request (Jan-23-2015):
!			Ayer Gladis y yo encontramos una inconsistencia entre la SFR que calculas con el programa para
!			correr las SFH a un z dado y los parámetros de la SSAG. El problema parece ser la interpretación
!			del tiempo de truncado (tcut, cuando comienza a decaer más rápidamente la SFR). Ivan calcula estas
!			cosas en look-back-time, lo que quiere decir que si él dice que una galaxia que se formó hace 8Gaños
!			tiene un tiempo de truncado de 1Gaño, quiere decir que comenzó a truncarse hace 1Gaño a partir del
!			presente. Mientras en tus SFR aparecería el tiempo de truncado como 8-1Gaños. Afortunadamente esto solo
!			afecta a las galaxias con truncado, solo 6 (+1 a z>0) en la muestra y puedo continuar sin estas galaxias.
			tauj = tauf - tauj	! added 03-Jun-2015
			if (tauj/tau2 < 100.) then
				con2 = con1 * exp(-tauj/tau1) / exp(-tauj/tau2)
			else
				con2 = 0.
			endif
		else
			tau2 = tau1
			tauj = tauf
			con2 = con1
		endif

!		Mass formed up to t = tauf
!		if (tauf <= tauj) then
!			xm1 = 1. - exp(-tauf/tau1)
!			xm2 = 0.
!		else
!			xm1  = 1. - exp(-tauj/tau1)
!			xm2  = con2*tau2*(exp(-tauj/tau2) - exp(-tauf/tau2))
!		endif

!		Compute SFR for both exponentials plus burst
		isfrc = 0
		do i=0,nsteps-1
		if (tb(i) <= tauf) then
			isfrc = isfrc + 1
			time(isfrc) = tb(i)
		endif
		enddo
!		Ask for burst parameters
14		write (6,'(x,a,$)') 'Add burst of fractional amplitude A = mass in burst/subyacent mass = '
		read (5,'(f20.0)',err=14,end=6) ampl
		if (ampl > 0.) then
15			write (6,'(x,a,$)') '                                      Burst starts at time (Gyr)   = '
			read (5,'(f20.0)',err=15,end=6) taui
			taui = 1.e9*taui
16			write (6,'(x,a,$)') '                                         Duration of burst (Gyr)   = '
			read (5,'(f20.0)',err=16,end=6) tdur
			tdur = 1.e9*tdur
			isfrc = isfrc + 1
			time(isfrc) = 0.9999*taui
			isfrc = isfrc + 1
			time(isfrc) = 0.99999*taui
			if (taui+tdur < tauf) then
				isfrc = isfrc + 1
				time(isfrc) = 1.0001*(taui+tdur)
				isfrc = isfrc + 1
				time(isfrc) = 1.00001*(taui+tdur)
!			else
!				isfrc = isfrc + 1
!				time(isfrc) = 1.0001*tauf
!				isfrc = isfrc + 1
!				time(isfrc) = 1.00001*tauf
			endif
			nbr=25
			dtbr=tdur/nbr
			do i=1,nbr+1
			if (taui+(i-1)*dtbr < tauf) then
				isfrc = isfrc + 1
				time(isfrc) = taui+(i-1)*dtbr
			endif
			enddo
		else
			ampl = 0.
			tburst=0.
		endif
		ampr=ampl

!		Store t = tauf
		isfrc = isfrc + 1
		time(isfrc) = tauf			! tauf = age of model galaxy requested by user
!		Store t = tobs in time array.
		tobs=age_model(tauf,0.,znew,zform,tlt1,tlt2)	! tobs = age in SSP model closest to tauf at z = 0
		isfrc = isfrc + 1
		time(isfrc) = tobs
!		Store t = tzob in time array.
		tzob=age_model(tobs,zobs,znew,zfor2,tlt1,tlt2)! tzob = age in SSP model closest to tauf-LTT(0:zobs); znew = z corresponding to age tzob, may differ slightly from zobs
		isfrc = isfrc + 1
		time(isfrc) = tzob
		tsng = tzob

!		Fill in SFR array
		do i=1,isfrc
		if (time(i) <= tauj) then
			usr_sfr(i) = con1*exp(-time(i)/tau1)
		else
			usr_sfr(i) = con2*exp(-time(i)/tau2)
		endif
		enddo
!		Sort and suppress duplicate entries in SFR
		call sort2(isfrc,time,usr_sfr)
		j=1
		do i=2,isfrc
		if (time(i) > time(i-1)) then
			j=j+1
			time(j)=time(i)
			usr_sfr(j)=usr_sfr(i)
			if (time(j) == tauf) jsfrc=j
			if (time(j) == tobs) jsfro=j
			if (time(j) == tzob) jsfrz=j
		endif
		enddo
		isfrc=j
!		Compute mass of subyacent population at t = tform
		tsubya=trapz1(time,usr_sfr,jsfrc)
!		Scale burst amplitude to obtain required fraction
!		Compute mass of stars formed in burst
		if (ampr > 0.) then
			if (taui+tdur <= tauf) then
				ampl   = ampr*tsubya/tdur
				tburst = ampl*tdur
			else
				ampl   = ampr*tsubya/(tauf-taui)
				tburst = ampl*(tauf-taui)
			endif
!			Add burst population
			do i=1,isfrc
			if (time(i) >= taui .and. time(i) <= taui+tdur) then
				usr_sfr(i) = usr_sfr(i) + ampl
			endif
			enddo
		endif
!		Compute total mass in stars at t = tform
		tstelm=trapz1(time,usr_sfr,jsfrc)
!		Compute total mass in stars at t = tobs (at z = 0)
		tstelo=trapz1(time,usr_sfr,jsfro)
!		Compute total mass in stars at t = tzob (at z = zobs)
		tstelz=trapz1(time,usr_sfr,jsfrz)
!		Compute mass of stars formed during last Gyr
		tlast1=trapz2(time,usr_sfr,isfrc,tauf-1.e9,tauf,ierr)
!		Write ascii file with SFR
!		close (501)
		write (aux,'(a)') '# SFR parameters after Chen et al. (2012)'
		call s500('a',aux,0.)
		write (aux,'(a)') '#'
		call s500('a',aux,0.)
!		read  (501,*) zx
		call r500(1,'f',aux,zx)
		if (zx == 0.) then
			write (aux,101)   'Z/Zsun'
			call s500('a',aux,0.)
		else
			write (aux,101)   'Z/Zsun',zx
			call s500('a',aux,0.)
		endif
!		read  (501,*) w1
		call r500(2,'f',aux,w1)
!		read  (501,*) aux
		call r500(3,'a',bux,zap)
		if (w1 > 0.) then
			write (aux,104)   'SSP-1:','   '//bux(1:largo(bux)),'   weight =',w1
			call s500('a',aux,0.)
		endif
!		read  (501,*) w2
		call r500(4,'f',aux,w2)
!		read  (501,*) aux
		call r500(5,'a',bux,zap)
		if (w2 > 0.) then
			write (aux,104)   'SSP-2:','   '//bux(1:largo(bux)),'   weight =',w2
			call s500('a',aux,0.)
		endif
!		read  (501,*) vtau,vmud
		call r500(6,'f',aux,vtau)
		call r500(7,'f',aux,vmud)
!		close (501)
		write (aux,101)   'tau_v ',vtau
		call s500('a',aux,0.)
		write (aux,101)   'mu_d  ',vmud
		call s500('a',aux,0.)
		write (aux,'(a)') '#'
		call s500('a',aux,0.)
		write (aux,100)   'Tform ',tauf*1.E-9
		call s500('a',aux,0.)
		if (tau1 < 1.E11) then
			write (aux,100)   'tau1  ',tau1*1.E-9
		else
			write (aux,107)   'tau1  ',tau1*1.E-9
		endif
		call s500('a',aux,0.)
		if (tau2 < 1.E11) then
			write (aux,100)   'tau2  ',tau2*1.E-9
		else
			write (aux,107)   'tau2  ',tau2*1.E-9
		endif
		call s500('a',aux,0.)
		write (aux,100)   'tauj  ',tauj*1.E-9
		call s500('a',aux,0.)
		write (aux,'(a)') '#'
		call s500('a',aux,0.)
		write (aux,101)   'Burst amplitude',ampr
		call s500('a',aux,0.)
		write (aux,100)   'Burst starts at',taui*1.e-9
		call s500('a',aux,0.)
		write (aux,100)   'Burst lasts    ',tdur*1.e-9
		call s500('a',aux,0.)
		write (aux,'(a)') '#'
		call s500('a',aux,0.)
!		Report stellar mass
		write (6,*)
		write (6,*) '              Mass computed numerically form SFR:'
		write (6,*)
		write (6,*) '              Mass of stars in subyacent population = ',tsubya
		write (6,*) '              Mass of stars formed in burst         = ',tburst
		write (6,*) '              Mass of stars formed during last Gyr  = ',tlast1
		if (ampr > 0.) then
			write (6,*) '              Mass ratio burst/subyacent population = ',tburst/tsubya,' = burst amplitude. % difference with requested amplitude =',(tburst/tsubya -ampr)/ampr*100.
		else
			write (6,*) '              Mass ratio burst/subyacent population = ',tburst/tsubya,' = burst amplitude'
		endif
		write (6,*)
		write (6,*) '              Mass of stars in total population     = ',tstelm,'   at t = Tform       = ',tauf*1.e-09,' Gyr  (requested)'
		write (6,*) '              Mass of stars in total population     = ',tstelo,'   at t = Tform''      = ',tobs*1.e-09,' Gyr  (used)     <---'
		write (6,*) '              Mass of stars in total population     = ',tstelz,'   at t = Tform''-LTT  = ',tzob*1.e-09,' Gyr  (used)     <---'
		write (6,*) '                                        At z = zobs = ',zobs  ,'   light travel time  = ',tlt1,' Gyr  (requested)'
		write (6,*) '                                        At z =        ',znew  ,'                 LTT  = ',tlt2,' Gyr  (used)     <---'
		write (6,*) '                                              zform = ',zform
		write (6,*)
		write (aux,'(a)') '# Mass computed numerically from SFR in table below.'
		call s500('a',aux,0.)
		write (aux,'(a)') '#'
		call s500('a',aux,0.)
		write (aux,102) 'Mass of stars in subyacent population',tsubya
		call s500('a',aux,0.)
		write (aux,102) 'Mass of stars formed in burst        ',tburst
		call s500('a',aux,0.)
		write (aux,102) 'Mass of stars formed during last Gyr ',tlast1
		call s500('a',aux,0.)
		write (aux,'(a)') '#'
		call s500('a',aux,0.)
		if (ampr > 0.) then
			write (aux,103) 'Mass ratio burst/subyacent population',tburst/tsubya,' = burst amplitude. % difference with requested amplitude =',(tburst/tsubya -ampr)/ampr*100.
		else
			write (aux,103) 'Mass ratio burst/subyacent population',tburst/tsubya,' = burst amplitude.'
		endif
		call s500('a',aux,0.)
		write (aux,'(a)') '#'
		call s500('a',aux,0.)
		write (aux,102) 'Mass of stars in total population    ',tstelm,'   at t = Tform       = ',tauf*1.e-09,' Gyr  (requested)'
		call s500('a',aux,0.)
		write (aux,102) 'Mass of stars in total population    ',tstelo,'   at t = Tform''      = ',tobs*1.e-09,' Gyr  (used)     <---'
		call s500('a',aux,0.)
		write (aux,102) 'Mass of stars in total population    ',tstelz,'   at t = Tform''-LTT  = ',tzob*1.e-09,' Gyr  (used)     <---'
		call s500('a',aux,0.)
		write (aux,105) '                         At z = zobs ',zobs  ,'   light travel time  = ',tlt1,' Gyr  (requested)'
		call s500('a',aux,0.)
		write (aux,105) '                         At z        ',znew  ,'                 LTT  = ',tlt2,' Gyr  (used)     <---'
		call s500('a',aux,0.)
		write (aux,105) '                               zform ',zform
		call s500('a',aux,0.)
		write (aux,'(a)') '#'
		call s500('a',aux,0.)
		write (aux,106) 'summary:',tauf,tstelm,tobs,tstelo,zform,zobs,znew,tzob,tstelz
		call s500('a',aux,0.)
		write (aux,'(a)') '#'
		call s500('a',aux,0.)
!		Write SFR properly
		write (aux,'(a)') '#    t(yr)           SFR (Mo/yr)'
		call s500('a',aux,0.)
		isfrc=jsfrz
		tauf=tzob
		do j=1,isfrc
		write (aux,*) time(j),usr_sfr(j)
		call s500('a',aux,0.)
		enddo
!		close (500)
100		format ('# ',a,' = ',f10.5,' Gyr')
107		format ('# ',a,' = ',f10.1,' Gyr')
101		format ('# ',a,' = ',f10.5)
102		format ('# ',a,' = ',1pe12.5,' Mo',a,0pf10.5,a)
103		format ('# ',a,' = ',f8.5,a,f6.3)
104		format ('# ',a,' = ',2a,f10.5)
105		format ('# ',a,' = ',f8.5,7x,a,0pf10.5,a)
106		format ('# ',a,1p4e14.5,0p3f10.5,1p2e14.5)
!		Read SFR from ascii file fort.500 written above and continue as in option 7
		io = 7
		tau=usrsfr(-500.)
		isingle=1
	endif
	last=0
	if (io.eq.0.or.io.eq.2.or.io.eq.7) return
17	write(6,'(x,a,f4.1,a,$)')'Make SFR = 0 at time > TCUT [',tcut*1.e-9,' Gyr] = '
	read (5,'(f10.0)',err=17,end=6) ttcut
	if (ttcut.gt.0.) then
		tcut=1.e9*ttcut
		if (io.eq.1.or.io.eq.4) then
			agas=exp(-tcut/tau)
		else
			agas=rgas(tcut)
		endif
		write (6,'(x,a,f6.4,a)') 'Unprocessed gas at TCUT = ',agas,' Mo'
	endif
	return
6	stop
	end

	Function rgas(tcut)

!	Computes amount of unprocessed gas remaining at t=tcut, not including
!	gas ejected by dying stars.

	real tfe(500)
	dt=0.05E9
	tt=-dt
	do i=1,500
	tt=tt+dt
	if (tt.ge.tcut) goto 1
	tfe(i)=sfr(tt)
	enddo
1	rgas=1.-trapez(tfe,i-1,dt)
	if (rgas.lt.0.) rgas=0.
	return
	end

	SUBROUTINE GET_EPSILON(IO)

!	Gets fraction of ejected gas to be recycled in stars = epsilon

	include 'recycle.dec'

	write (6,'(1x,a)') 'Epsilon = fraction of ejected gas to be recycled in stars = '
	write (6,'(1x,a)') 'Values from 0.001 to 1 have been explored'
	write (6,'(1x,a)') 'Epsilon = 0.001 reproduces old_galaxev mu_SFR=0.50 model'
	write (6,'(1x,a)') 'Epsilon > 1 emulates gas infall'
	write (6,'(1x,a)') 'Epsilon < 0 emulates galactic wind'
1	write (6,'(1x,a,$)') 'epsilon = '
	read (5,'(f10.0)',err=1,end=2) epsilon
	if (epsilon.eq.0.) io=4
	lgas=0
	return
2	stop
	end

	FUNCTION AGE_MODEL(TX,ZOBS,ZNEW,ZFORM,TLT1,TLT2)

!	Returns last time step to be used in csp_galaxev integration for Chen et al. SFH model
!	Last time step is the closest value of tb(i) to TX

!	Variables
	include 'jb.dec'
	include 'csp.dec'

!	Check for zobs > 0; Cosmology: data h/70./,omega/0.30/,omega_lambda/0.70/
	h0=70.
!	Obtain cosmological constant and parameter q
	clambda=cosmol_c(h0,0.30,0.70,q)
	if (zobs > 0.) then
		tuni=t(h0,q,0.00,clambda)	! age of universe at z = 0
		tobs=t(h0,q,zobs,clambda)	! age of universe at z = zobs
		tlt1=(tuni-tobs)*1.E9		! light travel time from z = 0 to z=zobs
		tobs=tx-tlt1			! age of galaxy at z = zobs
	else
		tlt1=0.       			! light travel time from z = 0
		tobs=tx				! requested age of galaxy at z = 0
	endif
!	Find closest age to t = tobs in array tb
	nfinal = kfinal(tobs)
!	Age of model at i = nfinal
	age_model= tb(nfinal)
!	Store age of galaxy (used) at z = 0.
!	Redshift corresponding to LTT = tobs-tb(nfinal) in Gyr
	if (age_model < tx) then
		tlt1  = tlt1*1.E-9
		tlt2  = (tx-age_model)*1.E-9
		znew  = zx(tlt2,h0,q,clambda)
	else
		tlt2  = 0.
		znew  = 0.
		zform = zx(tx*1.E-9,h0,q,clambda)
		if (zform == -2.) zform = 99.99999
	endif
	return
	end

	INTEGER FUNCTION KFINAL(TX)

!	Returns i for time step tb(i) closest in value to TX
!	Integration in csp_galaxev is performed up to this value for Chen et al. SFH model

!	Variables
	include 'jb.dec'
	include 'csp.dec'

!	Find closes time step to t = TX
	kfinal=0
	do i=0,nsteps-1
	if (tb(i) <= tx) then
		kfinal = i
	endif
	enddo

!	Refine position of age tb(i) closest to tx
	if (abs(tb(kfinal+1)-tx) < abs(tb(kfinal)-tx)) kfinal=kfinal+1
	return
	end
