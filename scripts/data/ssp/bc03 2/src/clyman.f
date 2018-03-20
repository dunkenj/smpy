	FUNCTION CLYMAN(X,Y,N)

!	Compute number of Lyman continuum photons in sed Y(X).
!	Assumes X is in Angstroms and Y in ergs/sec/Angstroms (physical flux)

!	variables
	real x(n),y(n)
	real*4 c,h,const,trapz1,w(200000),f(200000)
	save c,h,const
	data icall/0/,const/0./,wly/912./

!	compute proportionality constant
	if (icall.eq.0) then
		c=2.997925E10
		h=6.6262E-27
		const=1.0E-8/h/c
		icall=1
	endif

!	Find Lyman limit in sed and log(number of photons)
	if (x(1).ge.wly) then
		clyman=0.
		return
	endif
	do i=1,n
	if (x(i).lt.wly) then
		w(i)=x(i)
		f(i)=w(i)*y(i)
	elseif (x(i).eq.wly) then
		w(i)=x(i)
		f(i)=w(i)*y(i)
		goto 1
	elseif (x(i).gt.wly) then
		w(i)=wly
		f(i)=y(i-1)+(w(i)-x(i-1))*(y(i)-y(i-1))/(x(i)-x(i-1))
		f(i)=w(i)*f(i)
		goto 1
	endif
	enddo
1	clyman=const*trapz1(w,f,i)
	return
	end
