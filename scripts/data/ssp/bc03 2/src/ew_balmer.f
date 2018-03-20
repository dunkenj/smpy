	FUNCTION EW_BALMER(X,Y,N,H)

c	Computes equivalent width of the H_betha, H_gamma and H_delta lines
c	Computes the average of these 3 values

	real x(n),y(n),w(2000),z(2000),h(3)
	integer ip(3),im(3),io(3)
c	data icall/0/,ip/4,3,4/,im/4,10,4/,io/3*0/
	data icall/0/,ip/4,3,4/,im/4,10,4/,io/3*0/

c	Find position of Balmer lines
	ew_balmer=0.
	if (icall.eq.0) then
		icall=1
		do i=1,n
c		H_delta
		if (x(i).le.4100.) io(1)=i
c		H_gamma
		if (x(i).le.4340.) io(2)=i
c		H_betha
		if (x(i).ge.4860.) then
			io(3)=i
			goto 1
		endif
		enddo
1		if (n.gt.3500) then
c			High resolution mode
			call locate (x,n,4060.,im(1))
			im(1)=io(1)-im(1)
			call locate (x,n,4140.,ip(1))
			ip(1)=ip(1)-io(1)
			im(2)=16
			ip(2)=16
			call locate (x,n,4820.,im(3))
			im(3)=io(3)-im(3)
			call locate (x,n,4900.,ip(3))
			ip(3)=ip(3)-io(3)
		endif
	endif

c	Compute equivalent widths
	ew_balmer=0.
	do k=1,3
	i1=io(k)-im(k)
	i2=io(k)+ip(k)
	i=0
	do j=i1,i2
	a=(x(i2)-x(j))/(x(i2)-x(i1))
	fc=a*y(i1)+(1.-a)*y(i2)
	i=i+1
	w(i)=x(j)
	z(i)=1.-y(j)/fc
	enddo
	if (k.eq.2) then
c		Supress G-band (make profile symmetrical)
		if (n.lt.3500) then
c			G & S atlas
			w(1)=w(7)
			z(1)=z(14)
			w(2)=w(8)
			z(2)=z(13)
			do i=3,8
			w(i)=w(i+6)
			z(i)=z(i+6)
			enddo
			i=8
		else
c			High Resolution
			do l=1,12
			z(l)=z(33-l+1)
			enddo
		endif
	endif
	h(k)=trapz1(w,z,i)
	if (h(k).lt.0.) h(k)=0.
	ew_balmer=ew_balmer+h(k)
	enddo
	ew_balmer=ew_balmer/3.
	return
	end
