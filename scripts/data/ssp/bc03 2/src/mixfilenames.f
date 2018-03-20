	SUBROUTINE MIXFILENAMES(B,A)

!	Builds script to perform correct file assignments for program mix_stelib

!	Variables
	character*(*) a,b

!	Find first and third occurrence of character '_'
	if (index(b,'BaSeL').gt.0) then
		if (index(a,'BaSeL') == 0) then
			a = a(:largo(a)) // '_lr_BaSeL_ssp'
		endif
		i=index(a,'_lr')
		j=index(a(i+1:),'_')
		k=index(a(i+j+1:),'_')
		l=largo(a)
		open (120,file='mix_stelib.tmp')
!		Write script to use program mix_stelib
		write (120,'(5a)') 'setenv name1  ',a(:i),'hrs_stelib' ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name2  ',a(:i),'lr_BaSeL'   ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name3  ',a(:i),'hr_stelib'  ,a(i+j+k:l),'.ised'
		write (120,'(1a)') '#'
		write (120,'(5a)') 'setenv name4  ',a(:i),'lrs_hngsl'  ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name5  ',a(:i),'lr_BaSeL'   ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name6  ',a(:i),'lr_xhngsl'  ,a(i+j+k:l),'.ised'
		write (120,'(1a)') '#'
		write (120,'(5a)') 'setenv name7  ',a(:i),'hrs_indous' ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name8  ',a(:i),'lr_BaSeL'   ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name9  ',a(:i),'hr_xindous' ,a(i+j+k:l),'.ised'
		write (120,'(1a)') '#'
		write (120,'(5a)') 'setenv name10 ',a(:i),'hrx_miles'  ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name11 ',a(:i),'hr_stelib'  ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name12 ',a(:i),'hr_xmiless' ,a(i+j+k:l),'.ised'
		write (120,'(1a)') '#'
		write (120,'(5a)') 'setenv name13 ',a(:i),'hrx_miles'  ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name14 ',a(:i),'hr_xindous' ,a(i+j+k:l),'.ised'
		write (120,'(5a)') 'setenv name15 ',a(:i),'hr_xmilesi' ,a(i+j+k:l),'.ised'
		write (120,'(1a)') '#'
!		write (120,'(5a)') 'setenv name16 ',a(:i),'hrl_miles'  ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name17 ',a(:i),'hr_stelib'  ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name18 ',a(:i),'hr_xmilesls',a(i+j+k:l),'.ised'
!		write (120,'(1a)') '#'
!		write (120,'(5a)') 'setenv name19 ',a(:i),'hrt_miles'  ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name20 ',a(:i),'hr_stelib'  ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name21 ',a(:i),'hr_xmilests',a(i+j+k:l),'.ised'
!		write (120,'(1a)') '#'
!		write (120,'(5a)') 'setenv name22 ',a(:i),'hrw_miles'  ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name23 ',a(:i),'hr_stelib'  ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name24 ',a(:i),'hr_xmilesws',a(i+j+k:l),'.ised'
!		write (120,'(1a)') '#'
!		write (120,'(5a)') 'setenv name25 ',a(:i),'hrl_miles'  ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name26 ',a(:i),'hr_xindous' ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name27 ',a(:i),'hr_xmilesli',a(i+j+k:l),'.ised'
!		write (120,'(1a)') '#'
!		write (120,'(5a)') 'setenv name28 ',a(:i),'hrt_miles'  ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name29 ',a(:i),'hr_xindous' ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name30 ',a(:i),'hr_xmilesti',a(i+j+k:l),'.ised'
!		write (120,'(1a)') '#'
!		write (120,'(5a)') 'setenv name31 ',a(:i),'hrw_miles'  ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name32 ',a(:i),'hr_xindous' ,a(i+j+k:l),'.ised'
!		write (120,'(5a)') 'setenv name33 ',a(:i),'hr_xmileswi',a(i+j+k:l),'.ised'
!		write (120,'(1a)') '#'
		close (120)
	endif
	return
	end
