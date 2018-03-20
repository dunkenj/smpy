	FUNCTION B4000(X,Y,N)

c	Returns 4000 A break amplitude for sed in (x,y)
c	The break is defined as the ratio of the average Fnu flux densities:
c	(4050-4250)/(3750-3950) (Bruzual (1983), Hamilton 1985 ApJ 297, 371).

	real x(n),y(n),w(10000),z(10000)
	data icall/0/
	save i1,i2,i3,i4,j1,j2
	b4000=0.

	IF (N.EQ.1206.or.N.EQ.1977) THEN

c		BC93 (G+S) wavelength scale
		if (icall.eq.0) then
			icall=1
			do i=1,n

c			if (x(i).eq.3770.) i1=i
c			if (x(i).eq.3800.) i2=i
c			if (x(i).eq.3890.) i3=i
c			if (x(i).eq.3920.) i4=i

			if (x(i).le.3761.) i1=i
			if (x(i).eq.3810.) i2=i
			if (x(i).eq.3850.) i3=i
			if (x(i).eq.3920.) i4=i

			if (x(i).eq.4050.) j1=i
			if (x(i).eq.4250.) then
				j2=i
				goto 1
			endif
			enddo
		endif

c		Transform to Fnu units
1		i=0
		do j=j1,j2
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fr=trapz1(w,z,i)/(x(j2)-x(j1))

		i=0
		do j=i1,i2
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fl1=trapz1(w,z,i)/(x(i2)-x(i1))
		
		i=0
		do j=i3,i4
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fl2=trapz1(w,z,i)/(x(i4)-x(i3))
	
		fl=(fl1+fl2)/2.
	
		b4000=fr/fl

	ELSE

c		General procedure
		if (icall.eq.0) then
c			To allow different wavelength scales in same code always locate
c			icall=1
			call locate (x,n,4050.,i1)
			call locate (x,n,4250.,i2)
			call locate (x,n,3750.,i3)
			call locate (x,n,3950.,i4)
			i2=i2+1
			i4=i4+1
		endif

c		Transform to Fnu units
		i=0
		do j=i1,i2
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fr=trapz1(w,z,i)/(x(i2)-x(i1))

		i=0
		do j=i3,i4
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fl=trapz1(w,z,i)/(x(i4)-x(i3))

		b4000=fr/fl

	ENDIF
	return
	end

	FUNCTION B4000VN(X,Y,N)

c	MPA definition, via S Charlot

c	Returns 4000 A break amplitude for sed in (x,y)
c	The break is defined as the ratio of the average Fnu flux densities:
c	(4050-4250)/(3750-3950) (Bruzual (1983), Hamilton 1985 ApJ 297, 371).

	parameter (jid=35)
	real x(n),y(n),w(10000),z(10000)
	save i1,i2,i3,i4,j1,j2
c	Common to store fluxes
	integer iim(jid)
	real ffb(jid),ffr(jid),ffc(jid),ffl(jid)
	common /fluxes/ iim,ffb,ffr,ffc,ffl
	data icall/0/
	b4000vn=0.

	IF (N.EQ.1206.or.N.EQ.1977) THEN

c		BC93 (G+S) wavelength scale
		if (icall.eq.0) then
			icall=1
			do i=1,n

c			if (x(i).eq.3770.) i1=i
c			if (x(i).eq.3800.) i2=i
c			if (x(i).eq.3890.) i3=i
c			if (x(i).eq.3920.) i4=i

			if (x(i).le.3761.) i1=i
			if (x(i).eq.3810.) i2=i
			if (x(i).eq.3850.) i3=i
			if (x(i).eq.3920.) i4=i

			if (x(i).eq.4050.) j1=i
			if (x(i).eq.4250.) then
				j2=i
				goto 1
			endif
			enddo
		endif

c		Transform to Fnu units
1		i=0
		do j=j1,j2
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fr=trapz1(w,z,i)/(x(j2)-x(j1))

		i=0
		do j=i1,i2
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fl1=trapz1(w,z,i)/(x(i2)-x(i1))
		
		i=0
		do j=i3,i4
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fl2=trapz1(w,z,i)/(x(i4)-x(i3))
	
		fl=(fl1+fl2)/2.
	
		b4000vn=fr/fl

	ELSE
c		General procedure
		if (icall.eq.0) then
c			To allow different wavelength scales in same code always locate
c			icall=1
			call locate (x,n,4000.,i1)
			call locate (x,n,4100.,i2)
			call locate (x,n,3850.,i3)
			call locate (x,n,3950.,i4)
			i2=i2+1
			i4=i4+1
		endif

c		Transform to Fnu units
		i=0
		do j=i1,i2
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fr=trapz1(w,z,i)/(x(i2)-x(i1))

		i=0
		do j=i3,i4
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fl=trapz1(w,z,i)/(x(i4)-x(i3))

		b4000vn=fr/fl

c		Store fluxes for future use
		iim(31)=-1
		ffb(31)=fl
		ffr(31)=fr
		ffc(31)=b4000vn
		ffl(31)=b4000vn

	ENDIF
	return
	end

	FUNCTION B4000_SDSS(X,Y,N)

c	Returns 4000 A break amplitude for sed in (x,y)
c	Uses definition by Stoughton et al. (Ap. J), for SDSS data

	real x(n),y(n),w(10000),z(10000)
	save i1,i2,i3,i4,j1,j2
	data icall/0/
	b4000_sdss=0.

	IF (N.EQ.1206.or.N.EQ.1977) THEN

c		BC93 (G+S) wavelength scale
		if (icall.eq.0) then
			icall=1
			do i=1,n

c			if (x(i).eq.3770.) i1=i
c			if (x(i).eq.3800.) i2=i
c			if (x(i).eq.3890.) i3=i
c			if (x(i).eq.3920.) i4=i

			if (x(i).le.3761.) i1=i
			if (x(i).eq.3810.) i2=i
			if (x(i).eq.3850.) i3=i
			if (x(i).eq.3920.) i4=i

			if (x(i).eq.4050.) j1=i
			if (x(i).eq.4250.) then
				j2=i
				goto 1
			endif
			enddo
		endif

c		Transform to Fnu units
1		i=0
		do j=j1,j2
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fr=trapz1(w,z,i)/(x(j2)-x(j1))

		i=0
		do j=i1,i2
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fl1=trapz1(w,z,i)/(x(i2)-x(i1))
		
		i=0
		do j=i3,i4
		i=i+1
		w(i)=x(j)
		z(i)=y(j)*x(j)**2
		enddo
		fl2=trapz1(w,z,i)/(x(i4)-x(i3))
	
		fl=(fl1+fl2)/2.
	
		b4000_sdss=fr/fl

	ELSE

c		General procedure
		if (icall.eq.0) then
c			To allow different wavelength scales in same code always locate
c			icall=1
			call locate (x,n,4051.,i1)
			call locate (x,n,4251.,i2)
			call locate (x,n,3751.,i3)
			call locate (x,n,3951.,i4)
			i2=i2+1
			i4=i4+1
		endif

c		Store in new array
		i=0
		do j=i1,i2
		i=i+1
		w(i)=x(j)

c		Transform to Fnu units
c		z(i)=y(j)*x(j)**2

		z(i)=y(j)
		enddo
		fr=trapz1(w,z,i)/(x(i2)-x(i1))

		i=0
		do j=i3,i4
		i=i+1
		w(i)=x(j)

c		Transform to Fnu units
c		z(i)=y(j)*x(j)**2

		z(i)=y(j)
		enddo
		fl=trapz1(w,z,i)/(x(i4)-x(i3))

		b4000_sdss=fr/fl

	ENDIF
	return
	end

	FUNCTION B4000VN3(X,Y,N)

c	Returns 4000 A break amplitude for sed in (x,y)

c	The break is defined as the ratio of the average Fnu flux densities:
c	(4050-4250)/(3750-3950) (Bruzual (1983), Hamilton 1985 ApJ 297, 371).

c	Narrow definition from MPA: the break is defined as the ratio of the
c	average Fnu flux densities: (4000-4100)/(3850-3950)

c	Array declaration
	real x(n),y(n),w(10000),z(10000)
	b4000vn=0.

c	Find length of sed to be transformed to Fnu
	call locate (x,n,3850.,i1)
	call locate (x,n,4100.,i2)

c	Transform to Fnu units
	i=0
	do j=i1-10,i2+10
	i=i+1
	w(i)=x(j)
	z(i)=y(j)*x(j)**2
	write (21,*) w(i),z(i),y(j)
	enddo

c	Compute average flux on the red side
	x1=4000.
	x2=4100.
	fr = trapz2(w,z,i,x1,x2,ierr)/(x2-x1)

c	Compute average flux on the blue side
	x1=3850.
	x2=3950.
	fb = trapz2(w,z,i,x1,x2,ierr)/(x2-x1)

c	Compute b4000
	b4000vn3=fr/fb
        write (6,*) fr,fb,b4000vn3
        write (21,*) fr,fb,b4000vn3

	return
	end
