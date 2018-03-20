	SUBROUTINE USORT(N,RA)

c	Sorts array RA(N) and suppreses duplicate entries

	DIMENSION RA(N)

	if (n.le.1) return
	call sort(n,ra)
	j=1
	do i=2,n
	if (ra(i).gt.ra(i-1)) then
		j=j+1
		ra(j)=ra(i)
	endif
	enddo
	n=j
	return
	end

	SUBROUTINE UISORT(N,RA)

c	Sorts array RA(N) and suppreses duplicate entries

	INTEGER RA(N)

	if (n.le.1) return
	call isort(n,ra)
	j=1
	do i=2,n
	if (ra(i).gt.ra(i-1)) then
		j=j+1
		ra(j)=ra(i)
	endif
	enddo
	n=j
	return
	end

	SUBROUTINE REL_SORT(N,RA,EPS)

c	Sorts array RA(N) and suppreses entries that differ in less than eps %

	DIMENSION RA(N)

	if (n.le.1) return
	call sort(n,ra)
	j=2
	do i=3,n
	if ( (ra(i)-ra(i-1))/ra(i-1).gt.eps) then
		j=j+1
		ra(j)=ra(i)
	endif
	enddo
	n=j
	return
	end
