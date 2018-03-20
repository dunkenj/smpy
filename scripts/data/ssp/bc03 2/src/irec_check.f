	function irec_check(lun,a)

c	Get record length in *.ised file

	character*(*) a

c	Open and read file
	close (lun)
	open (lun,file=a,form='unformatted',status='old',err=3)
	read (lun)
	read (lun)
	read (lun,err=1) inl,(fl,i=1,inl),inx
!	write (6,*) 'irec_check: ',inl,inx
	irec_check=inx
	goto 2
1	irec_check=0
2	close (lun)
	return
3	write (6,'(2a)') 'File not found: ',a(1:largo(a))
	stop
	end
