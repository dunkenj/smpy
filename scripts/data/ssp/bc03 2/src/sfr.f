	FUNCTION SFR(TP)

c	Returns SFR at time TP

	include 'SSP_13.dec'
	include 'recycle.dec'

	sfr=0.
	if (tp.gt.tcut) return
	if (tp.lt.0.) then
		write (6,'(x,a,1pe12.4)') 'SFR called with tp < 0',tp
		return
	endif

c	Check for tcut
	if (io.eq.0) then
		if (tp.eq.0.) sfr=1.
	elseif (io.eq.1) then
		if (tau.gt.0.) then
			sfr=(1. + epsilon*pgas(tp))*exp(-tp/tau)/tau
		elseif (tau.lt.0.) then
			sfr=(1./(1.-exp(1.)) + epsilon*pgas(tp))*exp(-tp/tau)/tau
		endif
	elseif (io.eq.2) then
		if (tp.le.tau) sfr=1./tau
	elseif (io.eq.3) then
		sfr=tau
	elseif (io.eq.4) then
		if (tau.gt.0.) then
			sfr = exp(-tp/tau)/tau
		elseif (tau.lt.0.) then
			sfr = exp(-tp/tau)/tau/(1.-exp(1.))
		endif

	elseif (io.eq.6) then
		sfr=tp*exp(-tp/tau)/tau**2
	elseif (io.eq.7) then
		sfr=usrsfr(tp)
	elseif (io.eq.8) then
		sfr=2./tau*(1.-tp/tau)
		if (sfr.lt.0.) sfr=0.
	elseif (io.eq.9) then
		if (ampl > 0.) then
!			If burst has been added to double exponential SFR, use
			sfr=usrsfr(tp)
		else
!			Use analytic double exponential, joined at t = tauj
			if (tp < = tauj) then
				sfr = exp(-tp/tau1)/tau1
			else
				sfr = con2*exp(-tp/tau2)
			endif
		endif
	endif
	return
	end
