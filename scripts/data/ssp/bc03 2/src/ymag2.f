	FUNCTION YMAG2(WL,F,NI)

c	Computes Mg2 Index defined by Burstein, Faber, et al.
c	Written by G. Magris (May 1991)

	real wl(ni),f(ni),w0(6),wc(3),fl(200),wi(200),fi(3),wn(200),delta(3)
	real linear
	data nw,w0/3,4897.,4958.25,5156.,5197.25,5303.,5366.75/

	ymag2=0.

	wc(1)=(w0(1)+w0(2))/2.
	wc(3)=(w0(5)+w0(6))/2.
	wc(2)=(w0(3)+w0(4))/2.
	delta(1)=w0(2)-w0(1)
	delta(3)=w0(6)-w0(5)
	delta(2)=w0(4)-w0(3)

	np=21
	do j=1,3
	k=2*j-1
	del=(w0(k+1)-w0(k))/np
		do i=1,np+1
		l=i+(np+1)*(j-1)
		wi(l)=w0(k)+del*(i-1)
		enddo
	enddo

	do j=1,3
		i0=0
		do i=1,np+1
		l=i+(np+1)*(j-1)
		wn(i)=wi(l)
		fl(i)= linear(wi(l),wl,f,ni,i0)
		enddo
	fi(j)=trapz1(wn,fl,np+1)/delta(j)
	enddo

	fc=fi(1)+(fi(3)-fi(1))/(wc(3)-wc(1))*(wc(2)-wc(1))
	ymag2=-2.5*alog10(fi(2)/fc)

	return
	end
