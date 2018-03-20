	SUBROUTINE EVOL_ISED(KL,Z,H0,Q0,CLAMBDA,Y,ICALL,LAST)

c	New version 4-JULY-2003, includes cosmological constant clambda

c...	Returns (in array Y) the rest frame s.e.d. observed at redshift Z
c...	A linear interpolation is performed between two adjacents s.e.d.s
c...	It is assumed that Z=0 corresponds to s.e.d. in file record KL.
c...	ICALL must be set to zero the first time the routine is called.

	include 'read_ised.dec'
	Real y(imw),zz(2000)
	Save zz,tg,ageuni
	last=0
	if (icall.eq.0) then
		icall=1
		ageuni=T(h0,q0,0.,clambda)
		tg=ta(kl)
		write (6,*)
		do k=1,kl
		call percent(6,k,kl,'EVOL_ISED')
		i1=kl-k+1
		t1=(tg-ta(i1))*1.e-9
		zz(i1)=ZX(t1,h0,q0,clambda)
		enddo
	endif
		
	do k=1,kl-1
	i1=kl-k+1
	if (zz(i1).le.z.and.z.lt.zz(i1-1)) goto 1
	enddo
	if (z.eq.zz(1)) then
		i1=1
		i2=1
		a1=1.
		a2=0.
		goto 3
	endif
c	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c	If z > zz(i1), use first sed in model (=> no evolution at this z)
c	Modified as suggested by Rogier Windhorst to allow high z
c	colors to be computed. This means that after the t=0
c	sed from galaxev is used, only the k-correction will be
c	computed (no more evolution possible because galaxy
c	did not exist at higher z''s)
c	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	last=1
	i1=1
	i2=1
	i1=2
	i2=2
	a1=1.
c	write (6,'(/x,a,f6.2,a,f6.2,a/)') '***   z =',z,' >',zz(1),', using t = 0 sed (no evolution)   ***'
	goto 2
1	i2=i1-1
	tx=ageuni-T(h0,q0,z,clambda)
	t1=(tg-ta(i1))*1.e-9
	t2=(tg-ta(i2))*1.e-9
	a1=(t2-tx)/(t2-t1)
2	a2=1.-a1
3	do i=1,iw
	y(i)=a1*f(i,i1)+a2*f(i,i2)
	enddo
	return
	end

	FUNCTION ZK(KZ,ZF)
c	Returns redshift used at step kz
	if (kz.le.21) then
		zk =  0.000 + (kz  -1)*0.005
	elseif (kz.le.97) then
		zk =  0.100 + (kz -21)*0.025
	elseif (kz.le.177) then
		zk =  2.000 + (kz -97)*0.100
	else
		zk = 10.000 + (kz-177)*0.020
	endif
	if (zk.gt.zf) zk=zf
	return
	end
