	FUNCTION B912(X,Y,N)

c	Returns 912 A break amplitude for sed in (x,y).
c	The break is defined as the ratio of the average Fnu flux densities:
c	(800-900)/(1000-1100) (Bruzual (1983).

	real x(n),y(n),w(500),z(500)
	save a1,a2,a3,a4,i1,i2,i3,i4,x1,x2,x3,x4
	data icall/0/
	b912=0.
	if (x(1).ge.1100.) return
	if (icall.eq.0) then
		x1=800.
		call locate(x,n,x1,i1)
		a1=(x1-x(i1))/(x(i1+1)-x(i1))
		x2=900.
		call locate(x,n,x2,i2)
		a2=(x2-x(i2))/(x(i2+1)-x(i2))
		x3=1000.
		call locate(x,n,x3,i3)
		a3=(x3-x(i3))/(x(i3+1)-x(i3))
		x4=1100.
		call locate(x,n,x4,i4)
		a4=(x4-x(i4))/(x(i4+1)-x(i4))
		icall=1
	endif

c	Compute fluxes and transform to Fnu units
	i=1
	w(i)=x1
	z(i)=a1*y(i1+1)+(1.-a1)*y(i1)
	z(i)=z(i)*w(i)**2
	do j=i1+1,i2
	i=i+1
	w(i)=x(j)
	z(i)=y(j)*x(j)**2
	enddo
	i=i+1
	w(i)=x2
	z(i)=a2*y(i2+1)+(1.-a2)*y(i2)
	z(i)=z(i)*w(i)**2
	fl=trapz1(w,z,i)/(w(i)-w(1))
	i=1
	w(i)=x3
	z(i)=a3*y(i3+1)+(1.-a3)*y(i3)
	z(i)=z(i)*w(i)**2
	do j=i3+1,i4
	i=i+1
	w(i)=x(j)
	z(i)=y(j)*x(j)**2
	enddo
	i=i+1
	w(i)=x4
	z(i)=a4*y(i4+1)+(1.-a4)*y(i4)
	z(i)=z(i)*w(i)**2
	fr=trapz1(w,z,i)/(w(i)-w(1))
	if (fl.gt.0.) b912=fr/fl
	return
	end
