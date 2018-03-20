	SUBROUTINE INTRP_ISED(NAME,AGE,X,Y,NW)

!	Interpolates .ised file at any arbitrary age entered by user

!	Variables
	include 'jb.dec'
	include 'csp.dec'
	character*1024 name,save
	real x(imw),y(imw)
	data iread/0/,save/''/

!	Open and read input file, if required
	if (name /= save) then
		close (777)
		open (777,file=name,form='unformatted',status='old',err=1)
!		Read basic parameters from input file
		read (777) nsteps,(tb(i),i=1,nsteps)	! ml,mu,iseg,(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),totm,totn,avs,jo,tauo,id,tcut,ttt,ttt,ttt,id,id,igw,stelib
!		Read sed from file
		read (777) nw,(x(i),i=1,nw)
		do n=1,nsteps
		read (777) nw,(fl(i,n),i=1,nw)
		enddo
!		write (6,'(2a)') ' Reading file ',name(1:largo(name))
!		write (6,*) nw,' wavelength points per record'
!		write (6,*) nsteps,' time steps'
		close (777)
		save=name
	endif

!	Interpolate at required age
	if (age <= tb(1)) then
		do l=1,nw
		y(l) = fl(l,1)
		enddo
	elseif (age >= tb(nsteps)) then
		do l=1,nw
		y(l) = fl(l,nsteps)
		enddo
	else
		call locate(tb,nsteps,age,i)
		if (tb(i-1) > 0.) then
			a = alog10(age/tb(i-1))/alog10(tb(i)/tb(i-1))
		else
			a = age/tb(i)
		endif
		b=1.-a
		if (a == 0.) then
			do l=1,nw
			y(l) = fl(l,i-1)
			enddo
		elseif (a == 1.) then
			do l=1,nw
			y(l) = fl(l,i)
			enddo
		else
			do l=1,nw
			y(l) = b*fl(l,i-1) + a*fl(l,i)
			enddo
		endif
!		write (6,*) b,i-1,tb(i-1),a,i,tb(i)
	endif
	return
1	write (6,'(x,5a)') 'File ',name(1:largo(name)),' not found. Stopping code'
	stop
	end
