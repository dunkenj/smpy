	FUNCTION LARGO(A)
!	Returns significant length of string a
	character*(*) a
!	largo=index(a,' ')-1
	largo=0
	do i=1,len(a)
	if (a(i:i).ne.' ') largo=i
	enddo
	return
	end
	FUNCTION MARGO(A)
!	Returns max(1,largo(a)) to avoid error in DEC fortran
	character*(*) a
	margo=max(1,largo(a))
	return
	end
	FUNCTION NARGO(A)
!	Returns position where first dot occurs
	character*(*) a
	nargo=largo(a)
	do i=1,nargo
	if (a(i:i).eq.'.') then
		nargo=i-1
		return
	endif
	enddo
	return
	end
	FUNCTION LASTD(A)
!	Returns position of last dot
	character*(*) a
	lastd=0
	lastd=len(a)
	do i=1,len(a)
	if (a(i:i) == '.') lastd=i
	enddo
	return
	end
