	FUNCTION LDOT(S)

!	Returns position of last dot "." in array S

!	Variables
	character*(*) s

!	Finds last dot in array s
	j=index(s,'.')
	ldot = j
	if (ldot>0) then
		do while (j>0)
		j=index(s(ldot+1:),'.')
		ldot=ldot+j
		enddo
	else
		s=trim(s) // '.'
		ldot = len(trim(s))
	endif
	return
	end

	SUBROUTINE CHEXT(S,E)

!	Changes extension to filename s

!	Variables
	character*(*) s,e

!	Finds last dot in array s
	s(ldot(s)+1:)=e
	return
	end

	SUBROUTINE CHAEXT(S,E,N)

!	Changes extension to filename s

!	Variables
	character*(*) s,e

!	Use function chext
	call chext(s,e)
	return
	end

!	SUBROUTINE OLD_CHAEXT(S,E,N)

!	Changes extension (suffix) to filename s
!	SUN/FORTRAN version. G. Bruzual. 08-NOV-1988

!	character*(*) s,e
!	l=len(s)
!	k=len(e)
!	ib=largo(s)	
!	id=0
!	if (n >= 0) then
!!		Finds last dot in array s
!		do i=ib-1,1,-1
!		if (s(i:i).eq.'.') then
!			id=i
!			if (s(i+1:i+1).eq.'/') id=ib
!			goto 1
!		endif
!		enddo
!	endif
!1	if (id.eq.0) id=ib+1
!	s(id:id)='.'
!	if (id+k.gt.l) k=l-id
!	do i=1,k
!	s(id+i:id+i)=e(i:i)
!	enddo
!	n=id+k
!	do i=id+k+1,l
!	s(i:i)=' '
!	enddo
!	return
!	end
