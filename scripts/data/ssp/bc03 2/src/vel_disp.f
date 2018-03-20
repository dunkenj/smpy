        PROGRAM VEL_DISP

c	Applies gaussian velocity dispersion to all sed''s in *.ised file.
c	Computes line strength indices after smoothing sed.

c	Array declaration
        parameter (jb=2,jid=35)
	include 'csp.dec'
        character name*256,aux*9,library*10
	logical stelib
        real t(jts),x(imw),y(imw),w(imw),z(imw),gwx(jid),sf(0:jts)
!	real snr(0:jts),pnr(0:jts),bh(0:jts),sn(0:jts),wd(0:jts),rm(0:jts)
!	real bol(0:jts),str(0:jts),sf(0:jts),evf(0:jts)
c       Common to store fluxes used to compute indices
        integer iim(jid)
        real ffb(jid),ffr(jid),ffc(jid),ffl(jid)
        common /fluxes/ iim,ffb,ffr,ffc,ffl

c	Ask for file name
1       call copyright(6)
        write (6,'(x,a,$)') 'BC_GALAXEV SSP sed in file = '
	read (5,'(a)',end=10) name
2	write (6,'(x,a,$)') 'Sigma for velocity dispersion (km/s) = '
	read (5,'(f10.0)',err=2,end=10) sigma
	if (sigma.ge.100.) then
	  write (aux,'(a,i3)') '_sigma',nint(sigma)
	else
	  write (aux,'(a,i2)') '_sigma0',nint(sigma)
	endif
	call chaext(name,'ised',nm)
	open (1,file=name,form='unformatted',status='old')
c	Read basic parameters from SSP input file
	read (1) ks,(t(i),i=1,ks),ml,mu,iseg,
     &  (xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),
     &  totm,totn,avs,io,tauo,id,tcut,ttt,ttt,ttt,id,id,igw,stelib
	igw=0
	stelib=.false.
	read (1) n,(x(i),i=1,n)

c	Identify stellar library used in model
	imod=1
	if (n == 1238) then
!		BaSeL 3.1 _lr_ model. Nothing to do
		library = 'lr_BaSeL'
		write (6,*)
		write (6,*) 'You have entered a low resolution model *.ised file: ',library
		write (6,*) 'No velocity dispersion will be applied'
		close (1)
		goto 1
	elseif (n == 6700) then
!		Pure Stelib model
		library = 'hrs_stelib'
		imod=0
	elseif (n == 6917) then
!		Extended Stelib model. Smooth hr segment.
		library = 'hr_stelib'
		v1=3322.
		v2=9300.
	elseif (n == 4300) then
!		Pure Miles model
		library = 'hrs_miles'
		imod=0
	elseif (n == 12511) then
!		Extended Miles model. Smooth hr segment.
		library = 'hr_xmiles'
		v1=3550.
		v2=7400.
		v1=x(1)
		v2=9300.	! sed's have been extended with stelib
	elseif (n == 15011) then
		library = 'hrs_indous'
		imod=0
	elseif (n == 21243) then
		library = 'hr_indous'
!		Extended IndoUS model. Smooth hr segment.
		v1=3465.
		v2=9465.
	else
		write (6,*)
		write (6,*) 'You have entered an unknown type of model with',n,' points per record'
		write (6,*) 'Wavelength range:',x(1),x(n)
		close(1)
		goto 1
	endif

!	Proceed according to user choice
	if (imod == 1) then
!		_hr_ model. Smooth _hrs_ segment.
		write (6,*)
		write (6,*) 'You have entered a mixed low and high resolution model *.ised file: ',library
		write (6,*) 'Velocity dispersion will be applied to the high resolution segment'
		write (6,*) 'Low resolution part will be left untouched'
		do i=1,n
		if (x(i).le.v1) i1=i
		if (x(i).le.v2) i2=i
		enddo
	else
!		_hrs_ model. Smooth full sed.
		write (6,*)
		write (6,*) 'You have entered a pure high resolution model *.ised file: ',library
		write (6,*) 'Velocity dispersion will be applied to the full spectrum'
		write (6,*) 'Please ignore a few points at each end of the smoothed sed'
		v1=x(1)
		v2=x(n)
		i1=1
		i2=n
	endif
!	write (6,*) n,x(1),x(n)
!	write (6,*) i1,x(i1)
!	write (6,*) i2,x(i2)
	write (6,*)
	write (6,*) 'Spectral indices past',v2,' A cannot be computed'

c	Open output files
	i=index(name,'.')-1
	name=name(1:i) // aux(1:largo(aux))
	call chaext(name,'6lsindx_sed',nm)
	open (8,file=name,form='formatted',status='unknown')
	do i=1,27
	write (8,'(a,f10.0)') '# sigma = ',sigma
	enddo
	call chaext(name,'7lsindx_sed',nm)
	open (9,file=name,form='formatted',status='unknown')
	do i=1,27
	write (9,'(a,f10.0)') '# sigma = ',sigma
	enddo
	call chaext(name,'8lsindx_sed_fluxes',nm)
	open (10,file=name,form='formatted',status='unknown')
	do i=1,27
	write (10,'(a,f10.0)') '# sigma = ',sigma
	enddo
        write (8,801)
        write (8,802)
        write (8,803)
        write (8,804)
801     format ('# Index_No.     1:      2:      3:      4:      5:      6:      7:      8:      9:     ',
     *          '10:     11:     12:     13:')
802     format ('# log-age     CN_1    CN_2  Ca4227   G4300  Fe4383  Ca4455  Fe4531  Fe4668  H\\beta  ',
     *          'Fe5015    Mg_1    Mg_2    Mg-b')
803	format ('#   (yr)     (mag)   (mag)   (\\AA)   (\\AA)   (\\AA)   (\\AA)   (\\AA)   (\\AA)   ',
     *          '(\\AA)   (\\AA)   (mag)   (mag)   (\\AA)')
804	format ('#    (1)      (2)     (3)     (4)     (5)     (6)     (7)     (8)     (9)     (10)    ',
     *          '(11)    (12)    (13)    (14)')

        write (9,901)
        write (9,902)
        write (9,903)
        write (9,904)
901	format ('# Index_No.    14:     15:     16:     17:     18:     19:     20:     21:',
     *          '       WO-1:       WO-2:       WO-3:       WO-4:       GC-1:       4000A',
     *          '     DTT-Ca1     DTT-Ca2     DTT-Ca3     DTT-MgI       DM-04       DM-04       DM-04         BH:')
902	format ('# log-age   Fe5270  Fe5335  Fe5406  Fe5709  Fe5782    Na-D   TiO_1   TiO_2',
     *          '   H\\delta_A   H\\gamma_A   H\\delta_F   H\\gamma_F     D(4000)       B4_VN',
     *          '    CaII8498    CaII8542    CaII8662     MgI8807     H8_3889     H9_3835    H10_3798       BH-HK')
903	format ('#   (yr)     (\\AA)   (\\AA)   (\\AA)   (\\AA)   (\\AA)   (\\AA)   (mag)   (mag)',
     *          '       (\\AA)       (\\AA)       (\\AA)       (\\AA)          .           .',
     *          '        (\\AA)       (\\AA)       (\\AA)       (\\AA)       (\\AA)       (\\AA)       (\\AA)       (\\AA)')
904	format ('#    (1)      (2)     (3)     (4)     (5)     (6)     (7)     (8)     (9)',
     *          '         (10)        (11)        (12)        (13)        (14)        (15)',
     *          '        (16)        (17)        (18)        (19)        (20)        (21)        (22)        (23)')
	write (10,100)
100	format ('# log-age   Nx  Im   Flux_Blue   Flux_Red    Flux_Ctrl   Flux_Line       Index')

c	Open binary file to write new sed
	call chaext(name,'ised',nm)
	open (11,file=name,form='unformatted',status='unknown')
	write (11) ks,(t(i),i=1,ks),ml,mu,iseg,
     &  (xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),
     &  totm,totn,avs,io,tauo,id,tcut,ttt,ttt,ttt,id,id,igw,stelib
	write (11) n,(x(i),i=1,n)
	read   (1) n,(y(i),i=1,n),inx,(y(i),i=n+1,n+inx)
	write (11) n,(y(i),i=1,n),inx,(y(i),i=n+1,n+inx)

c	Smooth each sed in file and compute indices.
	do k=2,ks
	read (1) n,(y(i),i=1,n),inx,(y(i),i=n+1,n+inx)
	if (imod.eq.0) then
c		Smooth entire sed
		call gaussian_v_disp(x,y,n,sigma)
		write (11) n,(y(i),i=1,n),inx,(y(i),i=n+1,n+inx)
c		Compute Worthey indices plus other indices directly from smoothed sed
        	do j=1,jid
        	gwx(j)=gw_ix_sed(j,x,y,n,0)
        	enddo
	elseif (imod.eq.1) then
c		Smooth central segment
		j=0
		do i=i1,i2
		j=j+1
		w(j)=x(i)
		z(j)=y(i)
		enddo
		call gaussian_v_disp(w,z,j,sigma)
c		Ignore 10 smoothed points at each extreme to avoid border effects
		j2=j-10
		i3=i2-10
		write (11) n,(y(i),i=1,i1+9),(z(i),i=11,j2),(y(i),i=i3+1,n),inx,(y(i),i=n+1,n+inx)
c		Compute Worthey indices plus other indices directly from smoothed sed
        	do i=1,jid
        	gwx(i)=gw_ix_sed(i,w,z,j,0)
        	enddo
	endif

c       Write results for gwx(j)
	tl=alog10(t(k))
	write (8,105) tl,(gwx(j),j= 1,13)
        write (9,106) tl,(gwx(j),j=14,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,jid)
        write (10,107) (tl,j,iim(j),ffb(j),ffr(j),ffc(j),ffl(j),gwx(j),j=1,jid)
105     format (f10.6,13f8.3)
106     format (f10.6,8f8.3,14f12.3)
107	format (f10.6,i4,i4,1p4e12.3,0pf12.3)
        enddo

c	Add remaining records
	read (1,end=4) m,(bflx(i),i=1,m)
	write (11) m,(bflx(i),i=1,m)
	read (1,end=4) m,(strm(i),i=1,m)
	write (11) m,(strm(i),i=1,m)
	read (1,end=4) m,(sf(i),i=1,m)
	write (11) m,(sf(i),i=1,m)
	read (1,end=4) m,(evfl(i),i=1,m)
	write (11) m,(evfl(i),i=1,m)
	read (1,end=4) m,(snbr(i),i=1,m)
	write (11) m,(snbr(i),i=1,m)
	read (1,end=4) m,(pnbr(i),i=1,m)
	write (11) m,(pnbr(i),i=1,m)
	read (1,end=4) m,(bhtn(i),i=1,m)
	write (11) m,(bhtn(i),i=1,m)
	read (1,end=4) m,(sntn(i),i=1,m)
	write (11) m,(sntn(i),i=1,m)
	read (1,end=4) m,(wdtn(i),i=1,m)
	write (11) m,(wdtn(i),i=1,m)
	read (1,end=4) m,(rmtm(i),i=1,m)
	write (11) m,(rmtm(i),i=1,m)
	read (1,end=4) m,(toff(i),i=1,m)
	write (11) m,(toff(i),i=1,m)
	read (1,end=4) m,(bolms(i),i=1,m)
	write (11) m,(bolms(i),i=1,m)
	read (1,end=4) m,(gasms(i),i=1,m)
	write (11) m,(gasms(i),i=1,m)
	read (1,end=4) m,(galms(i),i=1,m)
	write (11) m,(galms(i),i=1,m)
4	close (1)
	close (8)
	close (9)
	close (10)
	close (11)
c	Delete unwanted files
	call delete_files(name,n,1)
	goto 1
10	end
