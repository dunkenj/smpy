	SUBROUTINE NREAD(X,NA,*,*)

c	This routine is useful when reading several parameters separated by
c	commas (not supported by current version of f77 in Sun 4/110).

c	Returns NA arguments in array X (read from the screen).
c	RETURN 1: Statement to execute if error while reading.
c	RETURN 2: Statement to execute if EOF detected.

	parameter (npar=24)
	character b*132
	real x(npar)

c	Clear buffer
	do i=1,npar
	x(i)=0.
	enddo

c	Read string b
	read (5,'(a)',end=2) b

c	Number of characters read into b
c	na=index(b,' ')-1
	na=largo(b)
	if (na.eq.0) return

c	Adds "," at the end of string b to guarantee extraction of last value
	b(na+1:na+1)=','

c	Extract numbers from b
	n=0
	l1=0
	do i=1,100
	l1=l1+n+1
	n=index(b(l1:),',')-1
	if (n.eq.0) then
		x(i)=0.
	elseif (n.gt.0) then
      		read(b(l1:l1+n-1),'(G132.0)',err=1) x(i)
	else
		na=i-1
		return
	endif
	enddo
1	return 1
2	return 2

	ENTRY QREAD(V,K,*,*)

c	Emulates q format of VMS/Fortran

c	V = required value
c	K = number of characters typed
c	RETURN 1: Statement to execute if error while reading.
c	RETURN 2: Statement to execute if EOF detected.

c	Read string b and assign value to v
	read (5,'(a)',end=4) b
c	k = index(b,' ') - 1
	k=largo(b)
	v=0.
	if (k.gt.0) then
      		read(b(:k),'(G132.0)',err=3) v
	endif
	return
3	return 1
4	return 2
	end
