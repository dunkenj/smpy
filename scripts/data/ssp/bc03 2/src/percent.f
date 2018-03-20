	subroutine percent(lun,i,n,label)

!	Reports percentage done

	character*(*) label,aux*2048
	common /kpercent/ ipcall,iplast,iphead,jc
	data ipcall/0/,iplast/0/,iphead/0/,jc/0/

	if (ipcall.eq.0) then
	  ipcall=1
	  iplast=0
	  iphead=1
	  jc=0
!	  aux = label(1:lastd(label)) // ':'
	  l = max(index(label,'.ised')-1,index(label,'.wfpc')-1)
	  if (l <= 0) l = largo(label)
	  aux = label(1:l) // ':'
	  write (lun,'(/x,a)') aux(1:l+1)
	  write (lun,'(a,i10,a)') ' % done: ...10...20...30...40...50...60...70...80...90..100 (',n,' steps)'
	  write (lun,'(a,$)')  	  '         '
	endif
	xi=i
	xn=n
	ip=nint(xi/xn*50)
	if (ip.gt.iplast) then
	  iplast=ip
	  write (lun,'(a,$)') '.'
	  call flush(lun)
!	  jc=jc+1
!	  if (jc.eq.(jc/10)*10) write (lun,*)
	endif
	if (i.ge.n) then
		ipcall=0
		iphead=0
		write (lun,*)
		write (lun,*)
	endif
	return
	end

	subroutine vpercent(lun,i,n,label)

!	Reports percentage done in vertical mode

	character*(*) label,aux*24
	common /kpercent/ ipcall,iplast,iphead,jc
!	data ipcall/0/,iplast/0/,iphead/0/,jc/0/

	if (ipcall.eq.0) then
	  ipcall=1
	  iplast=0
	  iphead=1
	  jc=0
	  aux = label(1:largo(label)) // ':'
	  write (lun,'(a)') aux(1:largo(aux))
	  write (lun,'(a,i5,a)') ' % done: ...  0  (',n,' steps left)'
!	  write (lun,'(a,i5,a)') ' % done: ... 10  (',n,' steps)'
!	  write (lun,'(a,$)') '         '
	endif
	xi=i
	xn=n
	ip=nint(xi/xn*50)
	if (ip.gt.iplast) then
	  iplast=ip
!	  write (lun,'(a,$)') '.'
	  jc=jc+1
	  if (jc.eq.(jc/5)*5) write (lun,'(a,i3,a,i4)') '         ...',2*jc,'    ',n-i
	  call flush(lun)
	endif
	if (i.ge.n) then
		ipcall=0
		iphead=0
		jc=0
!		write (lun,*)
!		write (lun,*)
	endif
	return
	end

	subroutine evpercent(lun,i,n,label)

!	Reports percentage done in vertical mode (en español)

	character*(*) label,aux*24
	common /kpercent/ ipcall,iplast,iphead,jc
!	data ipcall/0/,iplast/0/,iphead/0/,jc/0/

	if (ipcall.eq.0) then
	  ipcall=1
	  iplast=0
	  iphead=1
	  jc=0
	  aux = label(1:largo(label)) // ':'
	  write (lun,'(a)') aux(1:largo(aux))
	  write (lun,'(a,i5,a)') ' % hecho: ...  0  (',n,' por hacer)'
!	  write (lun,'(a,i5,a)') ' % hecho: ... 10  (',n,' steps)'
!	  write (lun,'(a,$)') '         '
	endif
	xi=i
	xn=n
	ip=nint(xi/xn*50)
	if (ip.gt.iplast) then
	  iplast=ip
!	  write (lun,'(a,$)') '.'
	  jc=jc+1
	  if (jc.eq.(jc/5)*5) write (lun,'(a,i3,a,i4)') '          ...',2*jc,'    ',n-i
	  call flush(6)
	endif
	if (i.ge.n) then
		ipcall=0
		iphead=0
		jc=0
!		write (lun,*)
!		write (lun,*)
	endif
	return
	end
