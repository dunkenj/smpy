
	function pgas(tx)

c	Returns amount of processed gas at time tx

	include 'recycle.dec'

c	Interpolate array (so,to) if possible
	if (tx.le.to(1).or.lgas.lt.1) then
		pgas=0.
	elseif (tx.lt.to(lgas)) then
		call locate(to,lgas,tx,i)
		a=(tx-to(i))/(to(i+1)-to(i))
		pgas=(1.-a)*so(i)+a*so(i+1)
	else
		pgas=so(lgas)
	endif
	return
	end
