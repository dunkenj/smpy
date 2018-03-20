	FUNCTION CHELIUM(X,Y,N,CHE2)

c	Compute number of Helium ionizing photons in sed Y(X).
c	Assumes X is in Angstroms and Y in ergs/sec/Angstroms (physical flux)

	real x(n),y(n),w(100000),f(100000)
	save c,h,const
	data icall/0/,const/0./

c	compute proportionality constant
	if (icall.eq.0) then
		c=2.997925E10
		h=6.6262E-27
		const=1.0E-8/h/c
		icall=1
	endif

c	Find Helium limit in sed and log(number of photons)
c	First electron, all photons shortward of whelium = 504.4 A
	whelium = 504.4
	if (x(1).ge.whelium) then
		chelium=0.
		return
	endif
	do i=1,n
	if (x(i).lt.whelium) then
		w(i)=x(i)
		f(i)=w(i)*y(i)
	elseif (x(i).eq.whelium) then
		w(i)=x(i)
		f(i)=w(i)*y(i)
		goto 1
	elseif (x(i).gt.whelium) then
		w(i)=whelium
		f(i)=y(i-1)+(w(i)-x(i-1))*(y(i)-y(i-1))/(x(i)-x(i-1))
		f(i)=w(i)*f(i)
		goto 1
	endif
	enddo
1	chelium=const*trapz1(w,f,i)

!	Second electron, all photons shortward of whelium = 228. A
	whelium = 228.
	if (x(1).ge.whelium) then
		che2=0.
		return
	endif
	do i=1,n
	if (x(i).lt.whelium) then
		w(i)=x(i)
		f(i)=w(i)*y(i)
	elseif (x(i).eq.whelium) then
		w(i)=x(i)
		f(i)=w(i)*y(i)
		goto 2
	elseif (x(i).gt.whelium) then
		w(i)=whelium
		f(i)=y(i-1)+(w(i)-x(i-1))*(y(i)-y(i-1))/(x(i)-x(i-1))
		f(i)=w(i)*f(i)
		goto 2
	endif
	enddo
2	che2=const*trapz1(w,f,i)
	return
	end
