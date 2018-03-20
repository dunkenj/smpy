	SUBROUTINE READ_ISED(file,tg,k0,ierr)

c	Reads galaxy s.e.d.s = output file from program GALAXEV.
c	Returns K0 = record number at which galaxy age = TG.

c	The array and common declaration are in file 'read_ised.dec'
	include 'read_ised.dec'
	include 'SSP_4.dec'
	include 'SSP_13.dec'
	logical stelib
	character*(*) file

c	Open and read file *.ised
	open (2,file=file,form='unformatted',status='old',err=1)
c	Read time scale and model parameters
	read (2,err=1) ks,(ta(i),i=1,ks),ml,mu,iseg,
     &	(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),
     &	totm,totn,avs,jo,tauo,id,tcut,ttt,ttt,ttt,id,id,igw,stelib
c	Read wavelength scale
	read (2) iw,(w(i),i=1,iw)
c	Read model sed
	do k=1,ks
	read (2) iw,(f(i,k),i=1,iw)
	enddo

c	Determine records bracketing age tg
	do k=1,ks-1
	if (tg.ge.ta(k).and.tg.lt.ta(k+1)) then
		a=(tg-ta(k))/(ta(k+1)-ta(k))
		b=1.-a
		k0=k+1
		do i=1,iw
		f(i,k0)=a*f(i,k)+b*f(i,k0)
		enddo
		ta(k0)=tg
!		write (6,100) k0,ta(k0)
!100		format (' Exit READ_ISED. Last record: ',i3,' Age = ',1pe12.4)
		close (2)
		ierr=0
		return
	endif
	enddo
1	write (6,'(x,3a)') 'File ',file(1:largo(file)),' not found. Please try again'
	ierr=1
	return
	end
