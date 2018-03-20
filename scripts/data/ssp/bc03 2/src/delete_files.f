	subroutine delete_files(name,iw,iop)

c	Writes command file to delete some files produced by program bc97_galaxev

c	Array declaration
	character*(*) name

c	Open command file
	close (91)
	open (91,file='bc03.rm',status='unknown',access='append')

c	Check type of file
	if (abs(iw).eq.1238.or.abs(iw).eq.6528) then
c		low resolution mode
		ires=0
	elseif (iw.eq.9000) then
c		various resolution modes called by downgrade_resolution
		ires=-2
	elseif (iw.lt.0) then
c		high resolution mode called by sub_mix_stelib
		ires=-1
	elseif (iw.eq.9848.or.iw.eq.11288.or.iw.eq.20305) then
c		high resolution mode called by sub_mix_stelib_pickles
		ires=-3
	elseif (iw.eq.4367.or.iw.eq.6700.or.iw.eq.15011) then
c		high resolution mode (hrs) called by any program
		ires=+1
	else
c		high resolution mode in general
		ires=+2
	endif

c	Write delete file command
	if (ires.eq.0) then
c		Keep only *.?color, *.1ABmag, and *.ised files in _lr_ mode
		write (91,'(a)') '# Delete _lr_ lsindx files.'
		call chaext(name,'?lsindx_ffn',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'?lsindx_sed',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'?lsindx_sed_lick_system',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'8lsindx_sed_fluxes',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
	elseif (ires.eq.-2) then
c		Keep only *.1color, *.2color, *.1ABmag, *.?lsindx_sed, and *.ised files
		write (91,'(a)') '# Delete *.?color files.'
		call chaext(name,'3color',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'4color',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		write (91,'(a)') '# Delete *.?lsindx_ffn files and similar.'
		call chaext(name,'?lsindx_ffn',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'?lsindx_sed_lick_system',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'8lsindx_sed_fluxes',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
	elseif (ires.eq.-1) then
c		Delete or rename *_hrs_* files after usage by sub_mix_stelib
		write (91,'(a)') '# Delete or rename _hrs_ color and lsindx files.'
		call chaext(name,'?color',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'1ABmag',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'?lsindx_ffn',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'?lsindx_sed',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'ised',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'6lsindx_sed_lick_system',nn)
		j=index(name,'_hrs_')
		write (91,'(5a)') 'mv -f ',name(1:largo(name)),' ',name(1:j+2),name(j+4:largo(name))
		call chaext(name,'7lsindx_sed_lick_system',nn)
		write (91,'(5a)') 'mv -f ',name(1:largo(name)),' ',name(1:j+2),name(j+4:largo(name))
		call chaext(name,'8lsindx_sed_fluxes',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
c		write (91,'(3a)') 'rm -f ',name(1:j+2),name(j+4:largo(name))
	elseif (ires.eq.-3) then
c		Delete excess files produced by sub_mix_stelib_pickles
		write (91,'(a)') '# Delete excess files produced by sub_mix_stelib_pickles'
		call chaext(name,'3color',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'4color',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'5color',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'?lsindx_ffn',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'?lsindx_sed',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'?lsindx_sed_lick_system',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'8lsindx_sed_fluxes',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		write (91,'(2a)') 'rm -f bc2003_lr_p?_*'
	elseif (ires.eq.+1) then
c		Delete or rename *_hrs_* color files after created by bc97_galaxev
		write (91,'(a)') '# Delete _hrs_ color and 8lsindx files.'
		call chaext(name,'?color',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'1ABmag',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
		call chaext(name,'8lsindx_sed_fluxes',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
	elseif (ires.eq.+2) then
		write (91,'(a)') '# Delete 8lsindx file.'
		call chaext(name,'8lsindx_sed_fluxes',nn)
c		write (91,'(2a)') 'rm -f ',name(1:largo(name))
	endif
	if (iop.ne.0) then
c		Delete *.5color file for non SSP's
		write (91,'(a)') '# Delete 5color file for non SSP''s.'
		call chaext(name,'5color',nn)
		write (91,'(2a)') 'rm -f ',name(1:largo(name))
	endif
	return
	end
