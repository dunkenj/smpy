	SUBROUTINE FILE_HEADER(IUN,NAME,LISTED)

c	Writes header information in file iun

	include 'SSP_4.dec'
	include 'SSP_13.dec'
	include 'cosmo.dec'
	character name*(*),und,buffer*96,blank*32
	real*8 cux
	common /vazdekis/ jvaz,xmu
	data und/'-'/,blank/'                                '/

c	Date file
	call dattim(iun,0,'L')

c	Clear file name
	ldot=index(name,'.')-1
	buffer(1:ldot)=name(1:ldot)
	buffer(ldot+1:64)=blank(ldot+1:32)

c	Check if model is for Vazdekis et al. (1996) IMF
	if (jvaz == 0) then
		i = index(buffer,'_v')
		if (i > 0 .and. buffer(i+3:i+3) == 'p' .and. buffer(i+6:i+6) == '_') jvaz =1
	endif

c	Write top lines of header
	write (iun,200) (und,i=1,120)
200	format ('#',10x,'I',120a,'I')
	write (iun,201)
201	format ('#',10x,'I',120x,'I')
	write (iun,202) buffer(1:42)
c202	format ('#',10x,'I      BC_GALAXEV  ---  MODEL PARAMETERS:   Generic file name for this model = ',a,10X,'I')
202	format ('#',10x,'I      BC_GALAXEV  ---  MODEL PARAMETERS:   Generic file name for this model = ',a,    'I')
	write (iun,201)
	write (iun,203) id
203	format ('#',10x,'I      TRACKS: ',a,26X,'I')
	write (iun,201)

c	Write information about SFR:
	if (io.eq.0) then
		write (iun,204)
204		format ('#',10x,'I      S.F.R.: SSP = Zero Length Burst at t = 0',74X,'I')
		write (iun,201)
	elseif (io.eq.1) then
		t9=tau*1.e-9
		tmu=1.-exp(-1./t9)
		write (iun,205) tmu,t9
205		format ('#',10x,'I      S.F.R.: Exponential with MU9 = ',F5.3,'/Gyr',6X,'TAU = ',F7.3,' Gyr',
     *                  ' (includes processed gas recycling)',16X,'I')
		write (iun,201)
	elseif (io.eq.2) then
		write (iun,206) tau
206		format ('#',10x,'I      S.F.R.: Finite Burst of Duration = ',1pe9.3,' yr',67X,'I')
		write (iun,201)
	elseif (io.eq.3) then
		write (iun,207) tau
207		format ('#',10x,'I      S.F.R.: Constant = ',1pe10.3,' Mo/yr',79X,'I')
		write (iun,201)
	elseif (io.eq.4) then
		t9=tau*1.e-9
		tmu=1.-exp(-1./t9)
		write (iun,208) tmu,t9
208		format ('#',10x,'I      S.F.R.: Exponential with MU9 = ',F5.3,'/Gyr',6X,'TAU = ',F7.3,' Gyr',
     *                  ' (does not include processed gas recycling)', 8X,'I')
		write (iun,201)
	elseif (io.eq.5) then
		ldot=index(name1,' ')-1
		buffer(1:ldot)=name1(1:ldot)
		buffer(ldot+1:32)=blank(ldot+1:32)
		write (iun,209) t0(1)*1.E-9,s(1),buffer(1:32)
209		format ('#',10x,'I      S.F.R.: 2 Bursts: Burst 1 at ',F6.3,'/Gyr, Strength = ',F8.5,', in file = ',a,10x,'I')
		ldot=index(name2,' ')-1
		buffer(1:ldot)=name2(1:ldot)
		buffer(ldot+1:32)=blank(ldot+1:32)
		write (iun,210) t0(2)*1.E-9,s(2),buffer(1:32)
210		format ('#',10x,'I                        Burst 2 at ',F6.3,'/Gyr, Strength = ',F8.5,', in file = ',a,10x,'I')
	elseif (io.eq.9 .or. isingle > 0) then
		write (iun,2091) tau1*1.E-9
		write (iun,2092) tau2*1.E-9
		write (iun,2093) tauj*1.E-9
		write (iun,2094) tauf*1.E-9
		write (iun,2095) ampr
		write (iun,2096) taui*1.E-9
		write (iun,2097) tdur*1.E-9
2091		format ('#',10x,'I      S.F.R.: Double Exponential with TAU_1 = ',F7.3,' Gyr (after Chen et al., does not include processed gas recycling) I')
2092		format ('#',10x,'I                                      TAU_2 = ',F7.3,' Gyr                                                               I')
2093		format ('#',10x,'I               which join smoothly at TAU_J = ',f7.3,' Gyr                                                               I')
2094		format ('#',10x,'I                    Look-back-time = Tform  = ',f7.3,' Gyr                                                               I')
2095		format ('#',10x,'I                            Burst amplitude = ',f7.3,'     ( = ratio of mass in burst / subyacent mass at t = Tform)     I')
2096		format ('#',10x,'I                            Burst starts at = ',f7.3,' Gyr                                                               I')
2097		format ('#',10x,'I                             Burst duration = ',f7.3,' Gyr                                                               I')
	else
		write (iun,201)
		write (iun,201)
	endif
	write (iun,201)

c	Write information about IMF
	if (iseg.eq.0) then
	   jseg=2
	   write (iun,221)
221	   format ('#',10x,'I      I.M.F.: Lognormal + power law   ',82X,'I')
	   write (iun,222)
222	   format ('#',10x,'I                                                                      mass in      number',31x,'I')
	   write (iun,223)
223	   format ('#',10x,'I                                                  from m     to m     segment    of stars      c',24x,'I')
	   aux=1.
           bux=chabimf(ml,mu,dble(ml),dble(aux),cux)
           write (iun,224) ml,aux,cux,bux
           bux=chabimf(ml,mu,dble(aux),dble(mu),cux)
           write (iun,225) aux,mu,cux,bux
224        format ('#',10x,'I',40x,'lognormal',f7.2,f9.2,2f12.4,31x,'I')
225        format ('#',10x,'I',40x,'  x=1.3  ',f7.2,f9.2,2f12.4,31x,'I')
	elseif (jvaz > 0) then
	   jseg=3
	   write (iun,321)
321	   format ('#',10x,'I      I.M.F.: Vazdekis et al. (1996) bimodal IMF',72x,'I')
	   write (iun,322)
322	   format ('#',10x,'I                                                                       total       number',31x,'I')
	   write (iun,323)
323	   format ('#',10x,'I                                              x   from m     to m       mass     of stars      c',24x,'I')
	   write (iun,324)
324	   format ('#',10x,'I                                           flat     0.10     0.20',55x,'I')
	   write (iun,325)
325	   format ('#',10x,'I                                         spline     0.20     0.60',55x,'I')
           write (iun,214) xx(1),0.60,um(1),baux(1),cn(1),cc(1)
	else
	   jseg=iseg
	   write (iun,211) iseg
211	   format ('#',10x,'I      I.M.F.: Power law (',I2,' segments):',82X,'I')
	   write (iun,212)
212	   format ('#',10x,'I                                                                      mass in      number',31x,'I')
	   write (iun,213)
213	   format ('#',10x,'I                                              x   from m     to m     segment    of stars      c',24x,'I')
           do i=1,iseg
           write (iun,214) xx(i),lm(i),um(i),baux(i),cn(i),cc(i)
214        format ('#',10x,'I',40x,f7.2,f9.2,f9.2,2f12.4,1pe12.4,19x,'I')
           enddo
	endif
	write (iun,201)
        write (iun,215) totm,totn,avs
215     format ('#',10x,'I',59x,'totals',2f12.4,f8.4,' Mo/star',15x,'I')
	do i=1,6-jseg
	write (iun,201)
	enddo
	if (listed.eq.0) then
		write (iun,216)
216		format ('#',10x,'I      LISTED: Rest frame properties',85X,'I')
		write (iun,201)
	else
                write (iun,217) h,omega,omega_lambda,q,clambda,tu,ttg,zf
217		format ('#',10x,'I      LISTED: Observer frame properties for the following cosmology:',52x,'I'/,
     *                  '#',10x,'I',14x,'Ho=',f4.0,2x,'Omega=',f4.2,2x,'Omega_lambda=',f4.2,2x,'qo=',f5.2,2x,
     *                  'Lambda=',1pE9.2,2x,'tu=',0pf5.2,' Gyr',2x,'tg=',f5.2,' Gyr',2x,'zf=',f5.2,2x,'I')
	endif
	write (iun,201)
	write (iun,218)
218	format ('#',10x,'I',30X,'(C) 1995-2015 G. Bruzual A. & S. Charlot - All Rights Reserved',28x,'I')
	write (iun,200) (UND,I=1,120)
	write (iun,301)
301	format ('#')
	return
	end
