	FUNCTION TRAPZ_LR(X,Y,N,X1,X2,IERR)

c	Finds area below y(x) from x=X1 to x =X2
c	X1 and X2 are not necessarily values of x(i)
c	but must fulfill X1>=x(1) and X2<=x(n)

c	This routine was especially written for coarse X arrays
c	where there are few points (or none) from X1 to X2 (like in Kurucz
c	model atmospheres). Interpolates Y(X) every unit in X.

c	Array declaration
	real x(n),y(n),xa(1000),ya(1000),linear

c       Check for values outside valid range
        TRAPZ_LR=0.
        ierr=0.
        if (n.le.1) return
        if (x1.lt.x(1)) then
                write (6,*) 'TRAPZ_LR: X1 outside valid range:',x1,x(1),1
                ierr=1
                return
c               stop
        endif
        if (x2.gt.x(n)) then
                write (6,*) 'TRAPZ_LR: X2 outside valid range:',x2,x(n),n
                ierr=2
                return
c               stop
        endif

c	Build auxiliary arrays (xa,ya) by linear interpolation
	j=0
	icall=0
	w1=x1-1.
	do i=1,n
	w1=w1+1
	if (w1.gt.x2) goto 1
	j=j+1
	xa(j)=w1
	ya(j)=linear(w1,x,y,n,icall)
	enddo
1	j=j+1
	xa(j)=x2
	ya(j)=linear(x2,x,y,n,icall)
	trapz_lr=trapz1(xa,ya,j)
	return
	end
