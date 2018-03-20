	SUBROUTINE S500(T,A,X)

!	Stores in array S the variable A for later use

!	Variables
	character*(*) A,T*1,s(1000)*128
	common /save500/ k,s
	data k/0/

!	Store variable in array s
	if (t == 'c') then
!		Clean buffer
		do i=1,1000
		s(i)=' '
		enddo
	elseif (t == 'k') then
!		Reset pointer
		k = nint(x)
	elseif (t == 'f') then
!		Store floating point x
		k=k+1
		write (s(k),*) x
	elseif (t == 'a') then
!		Store alpha numeric a
		k=k+1
		write (s(k),'(a)') a(:largo(a))
	elseif (t == 'r') then
!		Store SFR in file name.sfr
		s(k+1)=a
		call chext(a,'sfr')
		open (500,file=a)
		do i=8,k
		write (500,'(a)') s(i)(:largo(s(i)))
		enddo
		close (500)
		a=s(k+1)
	endif
	return
	end

	SUBROUTINE R500(J,T,A,X)

!	Reads floating point or alpha numeric previously stored in array S

!	Variables
	character T*1,a*128,s(1000)*128
	common /save500/ k,s

!	Check range
	if (j > k) then
		j=-k
		return
	endif

!	Read variable
	if (t == 'f') then
!		Read floating point
		read (s(j),*) x
	elseif (t == 'a') then
!		Read alpha numeric
		read (s(j),'(a)') a
	endif
	return
	end
