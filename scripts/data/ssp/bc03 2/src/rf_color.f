	SUBROUTINE RF_COLOR(io,t,x,y,inw,lun,bolflux,strmass,sf,evflux,snbr,pnbr,bh,sn,wd,rm,xm,bolms,gasmass,galmass,tmlr,tdpr)

!	Array declarations
	parameter (nc=100,lb=2*nc,jid=35,its=400)
	character genfile*96,envfile*96,rfcolorfile*96,h(0:200)*24
!	Negative indices allow for mandatory colors
	integer p,n1(-4:nc),n2(-4:nc),no(lb),ly(6),kerr(nc)
	real zp(-4:nc),col(-4:nc),phot(0:200),f(200),w(200),s(200)
	real umag,bmag,rmag,jmag,kmag,bsun,vsun,ksun,blr1,vlr1,klr1,blr,vlr,klr
	real x(inw),y(inw),fx(lb),balm(3),gwx(jid),ab(0:20)
	real m_pn,m_hm,m_bh,m_ns,m_wd
	save no,ly,zp,sl,mb,nv,mc,kb,nb,v0
	include 'stelib.dec'
	include 'filter.dec'
	common /kpercent/ ipcall,iplast,iphead,jc
!	Common to store fluxes used to compute indices
	integer iim(jid)
	real ffb(jid),ffr(jid),ffc(jid),ffl(jid)
	common /fluxes/ iim,ffb,ffr,ffc,ffl
	common /massrem/ m_pn,m_hm,m_bh,m_ns,m_wd
	real*8 mabo(its),mbel(its),mcut,rmup(its),rmlo(its)
	common /topimf/ mabo,mbel,rmup,rmlo,mcut,ic
	data icall/0/,ly/6*0/,last/0/
!	Define mandatory colors
	data n1/ 15,15,15,12,14, nc*0/
	data n2/125,57,84,13,15, nc*0/
	data kerr/nc*0/,jfits/0/
	common /f1001/ mab,ab

!	Check for zero flux sed
	if (t == 0. .or. bolflux <= 0.) return
	izero=0
	do i=1,inw
	if (y(i).gt.0) then
		izero=izero+1
	endif
	enddo
	if (izero.eq.0) then
!		write (6,*) 'Exiting rf_color at t =',t,' because of zero flux'
!		return
!		To allow using .?color files with a single number of lines
		do i=1,inw
		y(i) = 1.E-33
		enddo
	endif

!	Compute galaxy mass and mass in gas (except for add_bursts io = 5)
!	if (io.ne.10) then
!		galmass=gal_mass(io,t,sf)
!	endif

!	Compute Worthey indices directly from sed plus other indices
!	Store fluxes used to compute spectral indices in binary file
	do j=1,jid
	gwx(j)=gw_ix_sed(j,x,y,inw,0)
	enddo
!	Write results for gwx(j)
	tl=alog10(t)
	write (lun+8,105)  tl,(gwx(j),j= 1,21)
	write (lun+9,106)  tl,(gwx(j),j=22,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,jid)
	write (lun+10,110) (tl,j,iim(j),ffb(j),ffr(j),ffc(j),ffl(j),gwx(j),j=1,jid)
105     format (f10.6,21f9.4)
106     format (f10.6,1x,6f9.4,3f11.4,4f10.4,f10.4)
110	format (f10.6,i4,i4,1p4e12.3,0pf12.4)
!106    format (f10.6,8f8.3,14f12.3)
!	write (lun+9,106) tl,(gwx(j),j=14,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,jid)
!	write (lun+8,105) tl,(gwx(j),j= 1,13)

!	Check for stelib model
	if (stelib) then
		if (iphead.eq.0) then
!			Report what we are doing
			write (6,*)
!			write (6,'(x,a,$)') 'Computing SEDs and Indices'
			write (6,'(x,a/ )') 'Computing magnitudes, colours and indices'
			iphead=1
		endif
	endif
	if (stelib.or.inw.gt.6600) then
!		Compute indices in Lick system only for stelib or hr csp_galaxev models
		ism=0
		do j=1,jid
		gwx(j)=gw_ix_sed_lick_system(j,x,y,inw,0,ism)
		enddo
!		Write results for gwx(j)
		write (lun+11,105) tl,(gwx(j),j= 1,21)
		write (lun+12,106) tl,(gwx(j),j=22,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,jid)
!		write (lun+11,105) tl,(gwx(j),j= 1,13)
!		write (lun+12,106) tl,(gwx(j),j=14,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,34)
!		Do not compute colors for _hrs_ models
		if ((inw.eq.6700.or.inw.eq.15011).and.io.gt.0) return
	endif

!	Check if number of points has changed
	if (inw.ne.last) then
!		write (6,*) 'Resetting:',inw,last
!		ireset=1
		last=inw
	endif

	if (icall.eq.0) then
		icall=1

!		Format of file modified 07/01/2003:
!		To avoid recompiling this routine the arrays n1,n2 of
!		up to nc elements each are now read from a file.
!		Get file name from environment variable RF_COLORS_ARRAYS
		envfile='RF_COLORS_ARRAYS'
		call getenv(envfile,rfcolorfile)
		close (1)
		open (1,file=rfcolorfile,status='old',form='formatted',err=2)
!		write (6,*) 'List of filters in file: ',rfcolorfile(1:largo(rfcolorfile))
!		Read one line of header
		read (1,'(a)') genfile
!		write (6,*) genfile(1:largo(genfile))
		kb=0
!		Read filter pairs
		do i=1,nc
		read (1,'(2i4)',err=10) n1(i),n2(i)
		kb=kb+1
		fx(kb)=n1(i)
		kb=kb+1
		fx(kb)=n2(i)
		enddo
10		mc=i-1
		close (1)
!		Add mandatory filters
		do i=-4,0
		kb=kb+1
		fx(kb)=n1(i)
		kb=kb+1
		fx(kb)=n2(i)
		enddo
!		Sort filters in numerical order
		call sort(kb,fx)
!		Find independent filters in arrays n1 and n2, and store in no
		nb=1
		no(1)=fx(1)
		do i=2,kb
		if (fx(i).gt.fx(i-1)) then
			nb=nb+1
			no(nb)=nint(fx(i))
		endif
		enddo
!		Write filter ID
		call filterid(fid)
!		write (6,'(x,a)') 'Selected filters:'
!		do i=1,nb
!		write (6,'(i4,'': '',i4,3x,a)') i,no(i),fid(no(i))
!		enddo
!		write (6,'(i3,a)') mc,' colors selected:'
!		write (6,'(32i4)') (n1(i),i=1,mc)
!		write (6,'(32i4)') (n2(i),i=1,mc)

!		Log of solar luminosity
		sl=33.+alog10(3.826)

!		Read filter file and compute zero points
!		write (6,*)
!		write (6,*) 'Computing Zero Points:'
!		Compute V magnitude Vega zero point
		v0=vega_0p_n(15)
		do i=-4,mc
		zp(i)=zerop_n(n1(i),n2(i))

!		Fill arrays with filter numbers
		do j=1,nb
		if (n1(i).eq.no(j)) n1(i)=j
		if (n2(i).eq.no(j)) n2(i)=j
		if (no(j).eq.14) mb=j
		if (no(j).eq.15) nv=j
		enddo
		enddo

!		Find in array x the position of points used to define the
!		continuum at Lyman alpha
		do i=1,inw
		if (x(i).le.1120.) ly(1)=i
		if (x(i).le.1140.) ly(2)=i
		if (x(i).le.1160.) ly(3)=i
		if (x(i).le.1280.) ly(4)=i
		if (x(i).le.1300.) ly(5)=i
		if (x(i).ge.1320.) then
			ly(6)=i
			goto 1
		endif

		enddo
	endif
1	continue

!	Report what we are doing
	if (iphead.eq.0) then
		write (6,*)
!		write (6,'(x,a,$)') 'Computing SEDs and Colors. '
		write (6,'(x,a,$)') 'Computing magnitudes, colours and indices'
		iphead=1
	endif

!	Compute colors for SDSS
!	ireset=1
	call sdss_color(t,x,y,inw,lun,bolflux)

!	Compute flux through each of nb filters
	do i=1,nb
	if (kerr(i).eq.0) then
		fx(i)=filter_n(no(i),x,y,inw,0.,kerr(i))
	else
		fx(i)=0.
	endif
	enddo

!	Compute colors in Vega system
	do i=-4,mc
	if (fx(n1(i)).gt.0..and.fx(n2(i)).gt.0.) then
		col(i)=zp(i)-2.5*alog10(fx(n1(i))/fx(n2(i)))
!		write (6,*) i,n1(i),kerr(n1(i)),n2(i),kerr(n2(i))
	else
		col(i)=-99.99
	endif
	enddo

!	Compute bolometric magnitude
	bolmag=4.75-2.5*alog10(bolflux)
	bolpms=bolflux-bolms
	if (bolms.gt.0.) bolrat=bolpms/bolms

!	Compute V magnitude for a 1 Mo galaxy
!	It is -27.5 magnitudes brighter for a 1E11 Mo galaxy
	vmag=v0-2.5*alog10(fx(nv))

!	Compute U, B, R, K, and J2MASS magnitudes
	bmag=vmag+col( 0)
	umag=bmag+col(-1)
	rmag=vmag-col(-2)
	kmag=vmag-col(-3)
	jmag=vmag-col(-4)

!	Compute mass-to-visual-light ratio in solar units SUPERSEDED (Feb. 27, 2004).
!	Using a G2 V sed and the filters number 14 and 15
!	in the filter file, one derives:
!		fblue(sun) = 0.138Lo
!		fvis (sun) = 0.113Lo
!	This numbers apply only to the filters in this filter library.
!	Express blue and visual flux of model galaxy (also measured in
!	Lo) in units of the blue and visual flux of the sun:
!	fblu=fx(mb)/0.138
!	fvis=fx(nv)/0.113
!	fblu=filter(14,x,y,inw,0.)/0.138
!	fvis=filter(15,x,y,inw,0.)/0.113
!	Total mass in galaxy  = 1 Mo
!	Compute mass-to-visual-light ratio
!	blr=1./fblu
!	vlr=1./fvis
!	Compute stellar-mass-to--light ratios
!	blr=strmass/fblu
!	vlr=strmass/fvis

!	The solar absolute magnitudes for U,B,V,R,I,J,H,K were calibrated against the values
!	of Binney and Merrifield 1998, Galactic Astronomy, Table 2.1 (page 53), assuming
!	Bessell filters, and the offsets used to calibrate the entire set of filters.
!	Some values need checking - particularly those using UV filters (FOCA and Galex).
!	Taken from: http://mips.as.arizona.edu/~cnaw/sun.html

!		Filter	 B&M	 here	 difference
!		U	 5.61	 5.55	 0.06
!		B	 5.48	 5.45	 0.03
!		V	 4.83	 4.80	 0.03
!		R	 4.42	 4.46	 -0.04
!		I	 4.08	 4.11	 -0.03
!		J	 3.64	 3.67	 -0.02
!		H	 3.32	 3.33	 0.01
!		K	 3.28	 3.29	 0.01


!	Compute mass-to-visual-light ratio in solar units. Improved definition (Feb. 27, 2004).
!	Use solar absolute V and B magnitudes and (B-V)sun = 0.65
	vsun=4.80
	bsun=5.45
	ksun=3.29
!	M/L using the total mass in stars (M*), i.e. no remmnants
	blr1=strmass*10.**(0.4*(bmag-bsun))
	vlr1=strmass*10.**(0.4*(vmag-vsun))
	klr1=strmass*10.**(0.4*(kmag-ksun))
!	Modified Jan. 2011 to add mass of remmnants (rm = mBH + mNS + mWD)
	blr=(strmass+rm)*10.**(0.4*(bmag-bsun))
	vlr=(strmass+rm)*10.**(0.4*(vmag-vsun))
	klr=(strmass+rm)*10.**(0.4*(kmag-ksun))

!	Number of Lyman Continuum photons = Cly (log)
!	Flux in Lyman alpha from recombination theory
!	E(Lalpha) = 4.78E-13 * 33.1 * Nuv
!	log E = log(Nuv) -10.8 = cly -10.8
!	Number of Lyman continuum photons
	phly=clyman(x,y,inw)
	if (phly.gt.0.) then
		cly=sl+alog10(phly)
		fa=cly-10.8
	else
		cly=0.
		fa=0.
	endif

!	Number of Helium ionizing photons
	phe=chelium(x,y,inw,phe2)
	if (phe.gt.0.) then
		che=sl+alog10(phe)
	else
		che=0.
	endif
	if (phe2.gt.0.) then
		che2=sl+alog10(phe2)
	else
		che2=0.
	endif

!	Stellar continuum at Lyman alpha
	if (x(1).ge.1320.) then
		scly=0.
	else
		scly=(y(ly(1))+y(ly(2))+y(ly(3))+y(ly(4))+y(ly(5))+y(ly(6)))/6.
	endif
	if (scly.gt.0.) then
		fc=sl+alog10(scly)
	else
		fc=0.
	endif

!	Ly alpha equivalent width assuming that the continuum is the stellar continuum
	if (fa.gt.0..and.fc.gt.0.) then
		ew=10.**(fa-fc)
		ew2=phly/scly/10.**(10.8)
	else
		ew=0.
		ew2=0.
	endif

!	Compute Mg2 index
	ymg2=ymag2(x,y,inw)

!	Compute 912 A break
	b9=b912(x,y,inw)

!	Compute 4000 A break
	b4=b4000(x,y,inw)

!	Compute narrow version of D4000
	b4_n=b4000vn(x,y,inw)

!	Compute SDSS version of D4000
	b4_s=b4000_sdss(x,y,inw)

!	Compute equivalent width of Balmer lines (Hgamma, Hdelta, Hbeta)
	ewbl=ew_balmer(x,y,inw,balm)

!	Compute bolometric magnitude
	bolmag=4.75-2.5*alog10(bolflux)

!	Compute specific flux, snbr, pnbr
	evf=evflux/bolflux
!	write (6,*) t,evflux,evf,xm
	snr=snbr/bolflux
	pnr=pnbr/bolflux

!	SFR/year
!	sf=sfr(t)

!	Compute quantities requested by C. Popescu
!	call popescu(tl,sf,x,y,inw)

!	Compute flux from 1500 to 2800A
!	it=0
!	do i=1,inw
!	if (x(i).ge.1500.0.and.x(i).le.2800.0) then
!		it=it+1
!		xx(it)=x(i)
!		yy(it)=y(i)
!	endif
!	enddo
!	fuv=trapz1(x,y,it)
!	write (9,*) tl, sf,fuv

!	Write results
	write (lun+3 ,102) tl,b4,b4_n,b4_s,b9,cly,che,che2,bolmag,bolflux,snr,bh,sn,pnr,wd,rm
	write (lun+4 ,103) tl,bolmag,bmag,vmag,kmag,strmass,rm,gasmass,galmass,sf,strmass+rm,blr,vlr,klr,blr1,vlr1,klr1
	write (lun+5 ,104) tl,bolmag,evflux,evf,xm,bolrat
	write (lun+1 ,101) tl,bolmag,umag,bmag,vmag,kmag,(col(i),i=1,9)
	write (lun+2 ,101) tl,rmag,jmag,kmag,            (col(i),i=10,20)
	write (lun+14,111) tl,     kmag,                 (col(i),i=21,31),tmlr,tdpr
	write (lun+84,108) tl,vmag,kmag,                 (col(i),i=32,42),cly,blr,vlr,klr
	write (lun+85,108) tl,vmag,kmag,                 (col(i),i=43,53),cly,blr,vlr,klr
	write (lun+86,108) tl,vmag,kmag,                 (col(i),i=54,64),cly,blr,vlr,klr
	write (lun+87,108) tl,vmag,                      (col(i),i=65,80)
	write (lun+88,108) tl,vmag,kmag,                 (col(i),i=81,92),cly,blr,vlr,klr
!	write (lun+17,171) tl,mabo(ic),mbel(ic),mabo(ic)+mbel(ic),mcut,rmup(ic),rmlo(ic)    ! Esto escribe los datos requeridos por Alba (opcion kdist = 5)
!171	format (f10.6,1p6e12.3)
101	format (f10.6,14f10.4,1pe13.4)
102	format (f10.6,3f10.4,f11.4,2x,3f10.4,2x,f10.4,1p7e12.4)
103	format (f10.6,4f10.4,1p6e13.4,1x,6e12.4)
104	format (f10.6,f10.4,1p3e12.4,0pf10.4)
108	format (f10.6,18f10.4)
111	format (f10.6,12f10.4,1p2e13.4)

!	Organize photometric magnitudes in order to be listed in fits file
!	Write temporary file fort.1000 with photometric magnitudes to be inserted in fits file
!	return					 ! comment this statement to get file fort.1000 written

!	Johnson UBVRIJKL filters
	p =   0 ; phot(p) = bolmag               ;                h(p) = 'Mbol'          !  Mbol
	p = p+1 ; phot(p) = umag                 ; f(p) =  12   ; h(p) = 'U_Johnson'     !  U
	p = p+1 ; phot(p) = bmag                 ; f(p) =  14   ; h(p) = 'B_Johnson'     !  B
	p = p+1 ; phot(p) = vmag                 ; f(p) =  15   ; h(p) = 'V_Johnson'     !  V
	p = p+1 ; phot(p) = vmag - col(76)       ; f(p) =  32   ; h(p) = 'R_Johnson'     !  R = V - (V-R)
	p = p+1 ; phot(p) = vmag - col(77)       ; f(p) =  33   ; h(p) = 'I_Johnson'     !  I = V - (V-I)
	p = p+1 ; phot(p) = vmag - col(78)       ; f(p) =  34   ; h(p) = 'J_Johnson'     !  J = V - (V-J)
	p = p+1 ; phot(p) = vmag - col(79)       ; f(p) =  35   ; h(p) = 'K_Johnson'     !  K = V - (V-K)
	p = p+1 ; phot(p) = vmag - col(80)       ; f(p) =  36   ; h(p) = 'L_Johnson'     !  L = V - (V-L)

!	Cousins RI filters
	p = p+1 ; phot(p) = rmag                 ; f(p) =  84   ; h(p) = 'R_Cousins'
	p = p+1 ; phot(p) = vmag - col(11)       ; f(p) =  85   ; h(p) = 'I_Cousins'     !  I = V - (V-I)

!	Palomar JHK filters
	p = p+1 ; phot(p) = vmag - col(12)       ; f(p) =  55   ; h(p) = 'J_Palomar'     !  J = V - (V-J)
	p = p+1 ; phot(p) = phot(p-1) - col(15)	 ; f(p) =  56   ; h(p) = 'H_Palomar'     !  H = J - (J-H)
	p = p+1 ; phot(p) = kmag                 ; f(p) =  57   ; h(p) = 'K_Palomar'

!	Cowie K' filter
	p = p+1 ; phot(p) = vmag - col(17) 	 ; f(p) =  88   ; h(p) = 'Kprime_Cowie'  !   K' = V - (V-K')

!	2Mass JHKs filters
	p = p+1 ; phot(p) = vmag - col(18) 	 ; f(p) = 125   ; h(p) = 'J_2Mass'       !  J2M = V - (V-J2M)    2M = 2Mass
	p = p+1 ; phot(p) = vmag - col(19) 	 ; f(p) = 126   ; h(p) = 'H_2Mass'       !  H2M = V - (V-H2M)    2M = 2Mass
	p = p+1 ; phot(p) = vmag - col(20) 	 ; f(p) = 127   ; h(p) = 'Ks_2Mass'      ! Ks2M = V - (V-Ks2M)   2M = 2Mass

!	IRAC + IRAS + MIPS IR filters
	p = p+1 ; phot(p) = kmag      - col(21)	 ; f(p) = 128   ; h(p) = 'I3p6_IRAC'  ! I3.5 = K    - (K-I3.5)      IRAC
	p = p+1 ; phot(p) = phot(p-1) - col(22)	 ; f(p) = 129   ; h(p) = 'I4p5_IRAC'  ! I4.5 = I3.5 - (I3.5-I4.5)   IRAC
	p = p+1 ; phot(p) = phot(p-1) - col(23)	 ; f(p) = 130   ; h(p) = 'I5p7_IRAC'  ! I5.7 = I4.5 - (I4.5-I5.7)   IRAC
	p = p+1 ; phot(p) = phot(p-1) - col(24)	 ; f(p) = 131   ; h(p) = 'I7p9_IRAC'  ! I7.9 = I5.7 - (I5.7-I7.9)   IRAC
	p = p+1 ; phot(p) = phot(p-1) - col(25)	 ; f(p) =  71   ; h(p) = 'I12_IRAS'   ! I12  = I7.9 - (I7.9 - I12)  IRAS
	p = p+1 ; phot(p) = phot(p-1) - col(26)	 ; f(p) =  72   ; h(p) = 'I25_IRAS'   ! I25  = I12  - (I12  - I25)  IRAS
	p = p+1 ; phot(p) = phot(p-1) - col(27)	 ; f(p) =  73   ; h(p) = 'I60_IRAS'   ! I60  = I25  - (I25  - I60)  IRAS
	p = p+1 ; phot(p) = phot(p-1) - col(28)	 ; f(p) =  74   ; h(p) = 'I100_IRAS'  ! I100 = I60  - (I60 - I100)  IRAS
	p = p+3 ; phot(p) = phot(p-3) - col(31)	 ; f(p) = 134   ; h(p) = 'M160_MIPS'  ! M160 = I100 - (I100 - M160) MIPS
	p = p-1 ; phot(p) = phot(p+1) + col(30)	 ; f(p) = 133   ; h(p) = 'M70_MIPS'   ! M70  = M160 + (M70  - M160) MIPS
	p = p-1 ; phot(p) = phot(p+1) + col(29)	 ; f(p) = 132   ; h(p) = 'M24_MIPS'   ! M24  = M70  + (M24  - M70 ) MIPS

!	HST ACS wide filters
	p = p+3 ; phot(p) = vmag - col(54)	 ; f(p) = 219   ; h(p) = 'ACS_WFC_F220w' ! ACS WFC wide filter F220w = V - (V-F220w)
	p = p+1 ; phot(p) = vmag - col(55)	 ; f(p) = 220   ; h(p) = 'ACS_WFC_F250w' ! ACS WFC wide filter F250w = V - (V-F250w)
	p = p+1 ; phot(p) = vmag - col(56)	 ; f(p) = 221   ; h(p) = 'ACS_WFC_F330w' ! ACS WFC wide filter F330w = V - (V-F330w)
	p = p+1 ; phot(p) = vmag - col(57)	 ; f(p) = 222   ; h(p) = 'ACS_WFC_F410w' ! ACS WFC wide filter F410w = V - (V-F410w)
	p = p+1 ; phot(p) = vmag - col(58)	 ; f(p) = 223   ; h(p) = 'ACS_WFC_F435w' ! ACS WFC wide filter F435w = V - (V-F435w)
	p = p+1 ; phot(p) = vmag - col(59)	 ; f(p) = 224   ; h(p) = 'ACS_WFC_F475w' ! ACS WFC wide filter F475w = V - (V-F475w)
	p = p+1 ; phot(p) = vmag - col(60)	 ; f(p) = 225   ; h(p) = 'ACS_WFC_F555w' ! ACS WFC wide filter F555w = V - (V-F555w)
	p = p+1 ; phot(p) = vmag - col(61)	 ; f(p) = 226   ; h(p) = 'ACS_WFC_F606w' ! ACS WFC wide filter F606w = V - (V-F606w)
	p = p+1 ; phot(p) = vmag - col(62)	 ; f(p) = 227   ; h(p) = 'ACS_WFC_F625w' ! ACS WFC wide filter F625w = V - (V-F625w)
	p = p+1 ; phot(p) = vmag - col(63)	 ; f(p) = 228   ; h(p) = 'ACS_WFC_F775w' ! ACS WFC wide filter F775w = V - (V-F775w)
	p = p+1 ; phot(p) = vmag - col(64)	 ; f(p) = 229   ; h(p) = 'ACS_WFC_F814w' ! ACS WFC wide filter F814w = V - (V-F814w)

!	HST UVIS1 filters
	p = p+1 ; phot(p) = vmag - col(81)	 ; f(p) = 242   ; h(p) = 'UVIS1_f225w' ! UVIS1 f225w = V - (V - f225w.UVIS1)
	p = p+1 ; phot(p) = vmag - col(82)	 ; f(p) = 243   ; h(p) = 'UVIS1_f275w' ! UVIS1 f275w = V - (V - f275w.UVIS1)
	p = p+1 ; phot(p) = vmag - col(83)	 ; f(p) = 244   ; h(p) = 'UVIS1_f336w' ! UVIS1 f336w = V - (V - f336w.UVIS1)
	p = p+1 ; phot(p) = vmag - col(84)	 ; f(p) = 245   ; h(p) = 'UVIS1_f438w' ! UVIS1 f438w = V - (V - f438w.UVIS1)
	p = p+1 ; phot(p) = vmag - col(85)	 ; f(p) = 246   ; h(p) = 'UVIS1_f547m' ! UVIS1 f547m = V - (V - f547m.UVIS1)
	p = p+1 ; phot(p) = vmag - col(86)	 ; f(p) = 247   ; h(p) = 'UVIS1_f555w' ! UVIS1 f555w = V - (V - f555w.UVIS1)
	p = p+1 ; phot(p) = vmag - col(87)	 ; f(p) = 248   ; h(p) = 'UVIS1_f606w' ! UVIS1 f606w = V - (V - f606w.UVIS1)
	p = p+1 ; phot(p) = vmag - col(88)	 ; f(p) = 249   ; h(p) = 'UVIS1_f625w' ! UVIS1 f625w = V - (V - f625w.UVIS1)
	p = p+1 ; phot(p) = vmag - col(89)	 ; f(p) = 250   ; h(p) = 'UVIS1_f656n' ! UVIS1 f656n = V - (V - f656n.UVIS1)
	p = p+1 ; phot(p) = vmag - col(90)	 ; f(p) = 251   ; h(p) = 'UVIS1_f657n' ! UVIS1 f657n = V - (V - f657n.UVIS1)
	p = p+1 ; phot(p) = vmag - col(91)	 ; f(p) = 252   ; h(p) = 'UVIS1_f658n' ! UVIS1 f658n = V - (V - f658n.UVIS1)
	p = p+1 ; phot(p) = vmag - col(92)	 ; f(p) = 253   ; h(p) = 'UVIS1_f814w' ! UVIS1 f814w = V - (V - f814w.UVIS1)

!	AB mags
	p = p+1 ; phot(p) = ab (0)        	 ; f(p) = 120   ; h(p) = 'u_SDSS_AB'    ! SDSS u AB mag
	p = p+1 ; phot(p) = ab (1)        	 ; f(p) = 121   ; h(p) = 'g_SDSS_AB'    ! SDSS g AB mag
	p = p+1 ; phot(p) = ab (2)        	 ; f(p) = 122   ; h(p) = 'r_SDSS_AB'    ! SDSS r AB mag
	p = p+1 ; phot(p) = ab (3)        	 ; f(p) = 123   ; h(p) = 'i_SDSS_AB'    ! SDSS i AB mag
	p = p+1 ; phot(p) = ab (4)        	 ; f(p) = 124   ; h(p) = 'z_SDSS_AB'    ! SDSS z AB mag
	p = p+1 ; phot(p) = ab (5)        	 ; f(p) = 237   ; h(p) = 'u_CFHT_MC_AB' ! CFHT MC u AB mag
	p = p+1 ; phot(p) = ab (6)        	 ; f(p) = 238   ; h(p) = 'g_CFHT_MC_AB' ! CFHT MC g AB mag
	p = p+1 ; phot(p) = ab (7)        	 ; f(p) = 239   ; h(p) = 'r_CFHT_MC_AB' ! CFHT MC r AB mag
	p = p+1 ; phot(p) = ab (8)        	 ; f(p) = 255   ; h(p) = 'i_CFHT_MC_AB' ! CFHT MC i AB mag
	p = p+1 ; phot(p) = ab (9)        	 ; f(p) = 115   ; h(p) = 'y_CFHT_MC_AB' ! CFHT MC y AB mag
	p = p+1 ; phot(p) = ab (10)        	 ; f(p) = 241   ; h(p) = 'z_CFHT_MC_AB' ! CFHT MC z AB mag
	p = p+1 ; phot(p) = ab (11)        	 ; f(p) = 256   ; h(p) = 'Ks_CFHT_MC_AB'! CFHT MC Ks AB mag
	p = p+1 ; phot(p) = ab (12)        	 ; f(p) = 139   ; h(p) = 'FUV_GALEX_AB'  ! GALEX 1500 AB mag
	p = p+1 ; phot(p) = ab (13)        	 ; f(p) = 140   ; h(p) = 'NUV_GALEX_AB'  ! GALEX 2300 AB mag

!	Write temporary file fort.1000 with photometric magnitudes to be inserted in fits file
	if (jfits == 0) then
		jfits = 1
!		Compute effective wavelength for each of the p filtes
		do i=1,p
		s(i) = i
		area = areafilt(int(f(i)),w(i))
		enddo
!		Sort in increasing order of effective wavelength
		call sort2(p,w,s)
!		  do i=1,p
!		  j = s(i)
!		  write (1006,*) i,w(i),'   ',h(j)
!		  write (1006,*) i+2,'   ',h(j)
!		  enddo
!		  stop
!		  write (1000, '(a)') '# logage_yr   Mbol    U       B       V       R       I       J       K       L       Rc      Ic      Jp      Hp      Kp      Kpr     J2M     H2M    Ks2M    I3p6    I4p5    I5p7    I7p9    I12     I25     I60     I100     M24     M70    M160  A_F220w A_F250w A_F330w A_F410w A_F435w A_F475w A_F555w A_F606w A_F625w A_F775w A_F814w U_f225w U_f275w U_f336w U_f438w U_f547m U_f555w U_f606w U_f625w U_f656n U_f657n U_f658n U_f814w'
!		  write (1001, '(a)') '    u       g       r       i       z      u_MC    g_MC    r_MC    i_MC    y_MC    z_MC    Ks_MC   FUV     NUV'
!		write (1000, '(201(a,2x))') '# logage_yr',(trim(h(i)),i=0,p)
		write (1000, '(201(a,2x))') '# logage_yr',trim(h(0)),(trim(h(int(s(i)))),i=1,p)
		write (6,*) p+1,' magnitudes stored in temporary file fort.1000'
	endif
!	write (1000,109) tl,(phot(i),i=0,p)
	write (1000,109) tl,phot(0),(phot(int(s(i))),i=1,p)
109	format (f10.6,90f8.4)

	return
2	write (6,*) 'Error opening file: ',rfcolorfile(1:largo(rfcolorfile))
	stop
	end
