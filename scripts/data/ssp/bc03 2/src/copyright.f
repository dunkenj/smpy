	SUBROUTINE COPYRIGHT(JUN)
	data icount/0/
	lun=abs(jun)
	write (lun,*)
	write (lun,'(x,a)') 'Galaxy Spectral Evolution Library (GALAXEV)'
	write (lun,'(x,a)') 'UNIX Version (C) 1995-2017 - G. Bruzual and S. Charlot - All Rights Reserved'
	write (lun,*)
	if (jun.lt.0) return
	icount=icount+1
	write (lun,'(x,a,i3)') 'Computing model No.',icount
	write (lun,*)
	return
	end
