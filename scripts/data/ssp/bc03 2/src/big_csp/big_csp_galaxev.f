	PROGRAM BIG_CSP_GALAXEV

!	Computes spectral energy distribution as a function of time for a
!	Composite Stellar Population (CSP) by performing a convolution integral
!	of the spectral energy distributions for a Single Stellar Population SSP
!	with the chosen Star Formation Rate (SFR).

!	Modified in Toulouse: 6/8/97 to include both convolution of sed and spectral indices

!	Modified version of csp_galaxev to allow for large SFH files entered by user.
!	GBA, August 1st, 2013

!	Array declarations
	include 'jb.dec'
	include 'csp.dec'
	character ans,name*256,save*256,aux*512,atlas
	logical stelib
	real w(imw),h(imw)
	real tauv,mu_d,tv,fext,fi(40),fc(40),ta(0:jts)
	real snr(0:jts),pnr(0:jts), bh(0:jts), sn(0:jts), wd(0:jts), rm(0:jts)
	real bol(0:jts),str(0:jts), sf(0:jts),evf(0:jts),gas(0:jts),gal(0:jts),xml(0:jts),xdp(0:jts)
	real fxu(0:jts),fxg(0:jts),fxr(0:jts),fxi(0:jts),fxz(0:jts),fxk(0:jts),tofr(0:jts),bolmr(0:jts)
	real ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwa,mwla
	common /w_ages/ fxu,fxg,fxr,fxi,fxz,fxk,ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwa,mwla
	data ihrd/0/,atlas/' '/

!	Check if correct filter file is in use.
	j=ifilter()
	if (j.eq.0) then
		write (6,'(2a)') char(7),'Please assign correct filter file'
		write (6,'(2a)') char(7),'Use command: stdfilt'
		stop
	endif

!	Ask for file name
	save='_'
1	call copyright(6)
2	l=margo(save)
	write (6,'(x,3a,$)') 'BC_GALAXEV SSP sed in file [',save(1:l),'] = '
	read (5,'(a)',end=10) name
	if (largo(name).eq.0) then
		name=save
	else
!		Check if interpolation in metallicity of SSP is requested
		call interpolate_ssp(name)
	endif
	call s500('k','',5.)
	call chaext(name,'ised',mm)
	iread=0
	if (name.ne.save) iread=1
!	lr=index(name,'_lr_')
	open (1,file=name,form='unformatted',status='old',err=3)
	close (1)
	jdef=jjdef(name)

!       Ask for dust attenuation parameters
        write (6,'(/x,a,$)')   'Include attenuation by dust? Y/[N] '
        read (5,'(a)',end=10) ans
        if (ans.ne.'y'.and.ans.ne.'Y') then
           tauv=0.
	   mu_d=0.
        else
13         write (6,'(/x,a)')     '[using simple 2-component model of Charlot & Fall (2000)]'
           write (6,'(x,a,$)') 'Enter total effective attenuation optical depth: tau_V [1.0] = '
!          read (5,'(f10.0)',err=13,end=10) tauv
           read (5,'(a)',end=10) aux
	   if (largo(aux).eq.0) then
                tauv=1.
           else
                read (aux,*,err=13) tauv
           endif
14         write (6,'(x,a,$)') 'Enter fraction of tau_V arising from the ambient ISM: mu [0.3] = '
!          read (5,'(f10.0)',err=14,end=10) mu_d
	   read (5,'(a)',end=10) aux
           if (largo(aux).eq.0) then
                mu_d=0.3
           else
                read (aux,*,err=14) mu_d
           endif
	   write (6,'(x,a,f7.4)') '...using tau_V = ',tauv
	   write (6,'(x,a,f7.4)') '            mu = ',mu_d
        endif
!	write (501,*) tauv,mu_d
	call s500('f','',tauv)
	call s500('f','',mu_d)

!	Open input sed file (read again, tb modified in case of truncated SFR)
!	if (iread.ne.0) then
	   save=name
	   close (1)

!	   Check length of records (needed in Linux)
	   open (1,file=name,form='unformatted',status='old',err=3)
	   read (1)
	   read (1)
	   read (1,err=55) inl,(h(i),i=1,inl),inx
!	   write (6,*) inl,inx
	   imx=1
	   goto 56
55	   imx=0
56	   close (1)
	   open (1,file=name,form='unformatted',status='old',err=3)
!	   Read basic parameters from SSP input file
	   read (1) nsteps,(tb(i),i=0,nsteps-1),ml,mu,iseg,
     &	   (xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),
     &	   totm,totn,avs,jo,tauo,id,tcut,ttt,ttt,ttt,id,id,igw,stelib

	   if (jo.ne.0) then
	      write (6,'(x,a$)')'File does not contain an SSP. Proceed Y/[N] ? '
	      read (5,'(a)',end=10) ans
	      if (ans.ne.'y'.and.ans.ne.'Y') goto 1
	   endif
!	   Read sed from SSP file
	   write (6,'(2a)') ' Reading file ',name
	   read (1) inl,(w(i),i=1,inl)
	   do n=0,nsteps-1
	   if (imx.gt.0) then
		read (1) inl,(fl(i,n),i=1,inl),inx,(fl(i,n),i=inl+1,inl+inx)
	   else
		read (1) inl,(fl(i,n),i=1,inl)
		inx=0
	   endif
!          Attenuate by dust if requested (use Charlot & Fall, 2000, ApJ, 539, 718)
           if (tauv.gt.0.) then
              if (tb(n).le.1.e7) then
                 tv=tauv
              else
                 tv=mu_d*tauv
              endif
              do i=1,inl
                fl(i,n)=fl(i,n)*fext(w(i),tv)
              enddo
           endif
	   enddo
	   inw=inl+inx
	   write (6,*) inl,' wavelength points per record'
	   write (6,*) inx,' spectral index points per record'
	   write (6,*) inw,' total points per record'
	   write (6,*) nsteps,' time steps'
	   read (1,end=4) nsteps,(bflx(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(strm(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(sf(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(evfl(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(snbr(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(pnbr(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(bhtn(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(sntn(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(wdtn(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(rmtm(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(tofr(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(bolmr(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(gasms(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(galms(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(rml(i),i=0,nsteps-1)
	   read (1,end=4) nsteps,(rdp(i),i=0,nsteps-1)
	   do i=0,nsteps-1 ; toff(i)  = tofr(i) ; enddo
	   do i=0,nsteps-1 ; bolms(i) = bolmr(i) ; enddo
4	   close (1)

!	   Store flux in the u, g, r, i, z, and K bands for each SSP age in order
!	   to compute the flux weighted age of the composite population.
!		filter = 120 SDSS Camera u Response Function, airmass = 1.3 (June 2001)
!		filter = 121 SDSS Camera g Response Function, airmass = 1.3 (June 2001)
!		filter = 122 SDSS Camera r Response Function, airmass = 1.3 (June 2001)
!		filter = 123 SDSS Camera i Response Function, airmass = 1.3 (June 2001)
!		filter = 124 SDSS Camera z Response Function, airmass = 1.3 (June 2001)
!		filter =  57 IR K filter + Palomar 200 IR detectors + atmosphere.57
5	   write (6,'(/x,a,$)') 'Compute flux weighted age in the galaxy rest frame at z [0] = '
	   read (5,'(f12.0)',err=5,end=10) z
!	   read (5,*,err=5,end=10) z
	   do n=0,nsteps-1
	   do i=1,inl
	   h(i)=fl(i,n)
	   enddo
!	   In the observer frame use the value of z entered above
!	   zu = z
!	   In the galaxy restframe use z = 0
	   zu = 0.
	   fxu(n) = filter_n(120,w,h,inl,zu,kerr)
	   fxg(n) = filter_n(121,w,h,inl,zu,kerr)
	   fxr(n) = filter_n(122,w,h,inl,zu,kerr)
	   fxi(n) = filter_n(123,w,h,inl,zu,kerr)
	   fxz(n) = filter_n(124,w,h,inl,zu,kerr)
	   fxk(n) = filter_n( 57,w,h,inl,zu,kerr)
	   enddo

!	endif

!	Ask for SFR
	call sfr_0_b(1,z)

!	Define time steps to compute new sed, for instance:
!	tbeg=2.5E9 ; tend= 3.5e9 ; tstp=5.0e6 ; msteps=1+(tend-tbeg)/tstp ; do n=0,msteps-1 ; ta(n)=tbeg+n*tstp ; enddo
!	tbeg=2.5E9 ; tend= 3.5e9 ; tstp=1.0e6 ; msteps=1+(tend-tbeg)/tstp ; do n=0,msteps-1 ; ta(n)=tbeg+n*tstp ; enddo
	tbeg=0.0E9 ; tend=14.0e9 ; tstp=1.0e6 ; msteps=1+(tend-tbeg)/tstp ; do n=0,msteps-1 ; ta(n)=tbeg+n*tstp ; enddo
	write (6,*) msteps,' time steps in new grid, from:',ta(0),' to',ta(msteps-1)

!	Expand time scale if required
!	if ((io.ne.2.and.tcut.lt.20.E9).or.(io.eq.2.and.tcut.gt.2.E9)) call expand_time_steps
	if ((io.ne.2.and.tcut.lt.20.E9).or.(io.eq.2                 )) call expand_time_steps
	if (io.eq.7) call add_time_steps

!	Ask for output file name. Open files.
!	Write time scale, IMF, and wavelength scale in output file
	lun=-1
	ioption=99
	call name_sed(lun,jun,kdist,ihrd,ioption,jdef,ldef,ml,mu,name,atlas)

!	Compute convolution integral, rest frame colors, and write results for CSP
	if (isingle == 0) then
!		Use for standard csp models
		nstart = 0
		nfinal = msteps-1
		write (lun) msteps,(ta(i),i=0,msteps-1),ml,mu,iseg,(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),totm,totn,avs,io,tau,id,tau,tau,1.,1.,id,id,igw,stelib
	else
!		Use for Chen et al. SFH model
		call s500('r',name,0.)
		nfinal=kfinal(tsng)
		nstart = nfinal
		write (lun) 1,tb(nfinal),ml,mu,iseg,(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),totm,totn,avs,io,tau,id,tau,tau,1.,1.,id,id,igw,stelib
	endif
	write (lun) inl,(w(i),i=1,inl)
	do n=nstart,nfinal
	age=ta(n)
	if (io.gt.0) then
!		Perform convolution with chosen SFR
		call convolve_tx(h,age,bol(n),str(n),evf(n),snr(n),pnr(n),bh(n),sn(n),wd(n),rm(n),gas(n),gal(n),xml(n),xdp(n))
		call file_w_ages(io,name,zu,age,w,h,inl,ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,mwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwla)
	else
!		Interpolate in log t input SSP model at this age
		if (age <= tb(0)) then
			i=1
			a=0.
		elseif (age >= tb(nsteps-1)) then
			i=nsteps
			a=0.
		else
			call locate(tb,nsteps-1,age,i)
			if (tb(i-1) > 0.) then
				a = alog10(age/tb(i-1))/alog10(tb(i)/tb(i-1))
			else
				a = age/tb(i)
			endif
		endif
		b=1.-a
		do l=1,inw
		h(l)     = b*fl(l,i-1)   + a*fl(l,i)
		enddo
		bol(n)   = b*bflx(i-1)   + a*bflx(i)
		str(n)   = b*strm(i-1)   + a*strm(i)
		evf(n)   = b*evfl(i-1)   + a*evfl(i)
		snr(n)   = b*snbr(i-1)   + a*snbr(i)
		pnr(n)   = b*pnbr(i-1)   + a*pnbr(i)
		bh (n)   = b*bhtn(i-1)   + a*bhtn(i)
		sn (n)   = b*sntn(i-1)   + a*sntn(i)
		wd (n)   = b*wdtn(i-1)   + a*wdtn(i)
		rm (n)   = b*rmtm(i-1)   + a*rmtm(i)
		toff(n)  = b*tofr(i-1)   + a*tofr(i)
		bolms(n) = b*bolmr(i-1)  + a*bolmr(i)
		gas(n)   = b*gasms(i-1)  + a*gasms(i)
		gal(n)   = b*galms(i-1)  + a*galms(i)
		xml(n)   = b*rml(i-1)    + a*rml(i)
		xdp(n)   = b*rdp(i-1)    + a*rdp(i)
!		write (6,*) age,b,i-1,tb(i-1),a,i,tb(i)
!		read (5,*)
	endif
!	Compute galaxy mass
	do i=1,n
!	sf(n)=sfr(tb(n))
	sf(n)=sfr(ta(n))	! just for new time scale in big_csp_galaxev
	enddo
!	galmass=gal_mass(io,tb(n),sf(n))
	galmass=gal_mass(io,ta(n),sf(n))

!	Store sed. Report standard colors and Guy Worthey Spectral indices
	if (isingle == 0 .or. n == nfinal) then
		write (lun) inl,(h(i),i=1,inl),inx,(h(i),i=inl+1,inl+inx)
		call rf_color(io,age,w,h,inl,lun,bol(n),str(n),sf(n),evf(n),snr(n),pnr(n),bh(n),sn(n),wd(n),rm(n),toff(n),bolms(n),gas(n),gal(n),xml(n),xdp(n))
!		Guy Worthey Spectral indices
		inx2=inx/2
		do i=1,inx2
		fi(i)=h(inl+i)
		fc(i)=h(inl+inx2+i)
		enddo
		igw=0
!		call gw_indices(tb(n),fi,fc,lun,igw)
		call gw_indices(ta(n),fi,fc,lun,igw)
	endif

!	Report percent done
	if (n.gt.1) call percent(6,n,nfinal,'BIG_CSP_GALAXEV ' // name(1:nargo(name)))
	enddo

!	Add time behaviour of various quantities
	if (isingle == 0) then
		write (lun) msteps,(  bol(i),i=0,msteps-1)
		write (lun) msteps,(  str(i),i=0,msteps-1)
		write (lun) msteps,(   sf(i),i=0,msteps-1)
		write (lun) msteps,(  evf(i),i=0,msteps-1)
		write (lun) msteps,(  snr(i),i=0,msteps-1)
		write (lun) msteps,(  pnr(i),i=0,msteps-1)
		write (lun) msteps,(   bh(i),i=0,msteps-1)
		write (lun) msteps,(   sn(i),i=0,msteps-1)
		write (lun) msteps,(   wd(i),i=0,msteps-1)
		write (lun) msteps,(   rm(i),i=0,msteps-1)
		write (lun) msteps,( toff(i),i=0,msteps-1)
		write (lun) msteps,(bolms(i),i=0,msteps-1)
		write (lun) msteps,(  gas(i),i=0,msteps-1)
		write (lun) msteps,(  gal(i),i=0,msteps-1)
		write (lun) msteps,(  xml(i),i=0,msteps-1)
		write (lun) msteps,(  xdp(i),i=0,msteps-1)
	endif

!	Write command file to delete unwanted files
	call delete_files(name,inl,1)
	call file_w_ages(-1,name,zu,age,w,h,inl,ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,mwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwla)
	goto 1

3	write (6,'(x,5a)') char(7),'File ',name(1:largo(name)),' not found',char(7)
	goto 2
10	end

        REAL FUNCTION FEXT(X,TV)
        real x,tv,tau

        fext=1.

        if (tv.eq.0.) return

        tau=tv*( (5500./x)**0.7 )
        fext=exp(-tau)

        return
        end
