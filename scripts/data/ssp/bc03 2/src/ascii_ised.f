	PROGRAM ASCII_ISED

c	Writes *.ised file in ASCII format

c	Array declarations
	logical stelib
	character id*80,id2*80,id3*80,name*256

c	Maximum number of wavelength points
	PARAMETER (imw=100000)
	real w(imw),h(imw),tb(0:500),f(100)

c	Array declaration for IMF data
	parameter (imf=10)
	real ml,mu,xx(imf),lm(imf),um(imf),baux(imf),cc(imf),cn(imf)

c	Ask for file name if not enterd in command line
	j=iargc()
        if (j.eq.0) then
		call copyright(6)
		write (6,'(x,a,$)') 'BC model *.ised file = '
		read (5,'(a)',end=10) name
	else
		call getarg(1,name)
		call chaext(name,'ised',nn)
		write (6,*)
		write (6,*) 'Input file:  ',name(1:largo(name))
	endif

c	Read basic parameters from input file
	open (1,file=name,form='unformatted',status='old',err=3)
	read  (1)       nsteps,(tb(i),i=1,nsteps),
     &		        ml,mu,iseg,(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),
     &			totm,totn,avs,jo,tau,id,tau1,tau2,tau3,tau4,id2,id3,iop,stelib
	write (6,*)     nsteps,' time steps'

c	Read wavelength scale from input file and write to ASCII file
	read  (1)       inw,(w(i),i=1,inw)
	write (6,'(i8,a,f10.2,a,f10.2,a)')     inw,' wavelength points per record, from',w(1),' to',w(inw),' A'

c	Write basic parameters to output file
	name=name(1:largo(name)) // '_ASCII'
	open  (2,file=name,form='formatted',status='unknown')
	write (2,*) nsteps,(tb(i),i=1,nsteps)
	write (2,*) ml,mu,iseg,(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg)
     	write (2,*) totm,totn,avs,jo,tau,   tau1,tau2,tau3,tau4,        iop,stelib
	write (2,'(a)') id
	write (2,'(a)') id2
	write (2,'(a)') id3

c	Write wavelength scale to binary file
	write (2,*) inw,(w(i),i=1,inw)

c	Read sed''s from input file and write to ASCII file
c	Include fluxes for fitting functions
	do n=1,nsteps
	read  (1,err=1)   ik,(h(i),i=1,ik),ix,(f(i),i=1,ix)
1	write (2,*)       ik,(h(i),i=1,ik),ix,(f(i),i=1,ix)
	enddo
	write (6,*)     ik,' spectral fluxes per record,', ix,' fitting function fluxes per record'
c	Copy 12 records after the sed''s.
	do n=1,12
	read  (1,end=4)   ik,(h(i),i=1,ik)
	write (2,*)       ik,(h(i),i=1,ik)
	enddo
	write (6,*)  12,' extra records,',ik,' points each.'
4	write (6,*) 'Output file:  ',name(1:largo(name))
	write (6,*)

c	Close files
	close (1)
	close (2)
	stop
c	goto 1
3	write (6,'(3a)') 'File ',name(1:largo(name)),' not found'
c	goto 1
10	end
