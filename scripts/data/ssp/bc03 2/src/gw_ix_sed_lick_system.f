	FUNCTION GW_IX_SED_LICK_SYSTEM(J,X,F,N,IGW,ISM)

c	Compute Lick index J value in the Lick system.
c	Stelib sed (x,f) is transformed previously to the
c	Lick system if ISM=0 and then the index is computed

c	array declaration
	real x(n),f(n),shift(25)
!	real y2(16000)
!	save y2

c       Transform input sed to Lick/IDS system and store seds
c	UPDATED NOV 8, 2008
c	Input sed is assumed to be already degraded to Lick IDS spectral resolution
c	if (ism.eq.0) then
c        	call lick_system(x,f,n,y2)
c		ism=1
c		do i=1,n
c		y2(i)=f(i)
c		enddo
c	endif

c	Define median shift obtained from STELIB(smoothed) - LICK
c	i.e., to transform to Lick System, subtract shift to value
c	obtained from smoothed stelib sed

c	November 2002 version (STELIB wavelength scale uncorrected for radial velocity)
c	I used routine trapz2 for integration of line indices
c	shift( 1) = -0.0101
c	shift( 2) = -0.0068
c	shift( 3) = -0.0704
c	shift( 4) = -0.0048
c	shift( 5) =  0.3649
c	shift( 6) = -0.1808
c	shift( 7) = -0.2267
c	shift( 8) =  0.3792
c	shift( 9) =  0.1065
c	shift(10) =  0.2799
c	shift(11) = -0.0210
c	shift(12) = -0.0188
c	shift(13) = -0.0634
c	shift(14) =  0.0424
c	shift(15) = -0.0864
c	shift(16) =  0.0762
c	shift(17) =  0.0236
c	shift(18) =  0.0097
c	shift(19) = -0.0393
c	shift(20) =  0.0032
c	shift(21) =  0.0037
c	shift(22) =  0.9586
c	shift(23) = -0.8760
c	shift(24) =  0.1282
c	shift(25) = -0.2813

c	April 2003 version (STELIB wavelength scale corrected for radial velocity)
c	Changed to routine trapz3 for integration of line indices
c	sigma = FWHM/2.0000 (wrong)
c	shift( 1) = -0.0106
c	shift( 2) = -0.0079
c	shift( 3) = -0.0891
c	shift( 4) = -0.0135
c	shift( 5) =  0.3401
c	shift( 6) = -0.1877
c	shift( 7) = -0.2389
c	shift( 8) =  0.3413
c	shift( 9) =  0.1110
c	shift(10) =  0.3286
c	shift(11) = -0.0208
c	shift(12) = -0.0188
c	shift(13) = -0.0416
c	shift(14) =  0.0765
c	shift(15) = -0.0217
c	shift(16) =  0.1154
c	shift(17) =  0.0068
c	shift(18) =  0.0080
c	shift(19) = -0.0034
c	shift(20) =  0.0028
c	shift(21) =  0.0035
c	shift(22) =  0.9580
c	shift(23) = -0.8710
c	shift(24) =  0.1758
c	shift(25) = -0.3096

c	August 2003 version (STELIB wavelength scale corrected for radial velocity)
c	Changed to routine trapz3 for integration of line indices
c	sigma = FWHM/2.3548
	shift( 1) = -0.0107
	shift( 2) = -0.0012
	shift( 3) = -0.0187
	shift( 4) =  0.0543
	shift( 5) =  0.5247
	shift( 6) = -0.1342
	shift( 7) = -0.1204
	shift( 8) =  0.4260
	shift( 9) =  0.1269
	shift(10) =  0.4657
	shift(11) = -0.0197
	shift(12) = -0.0179
	shift(13) =  0.0264
	shift(14) =  0.1740
	shift(15) =  0.0720
	shift(16) =  0.2040
	shift(17) =  0.0267
	shift(18) =  0.0420
	shift(19) =  0.0117
	shift(20) =  0.0029
	shift(21) =  0.0038
	shift(22) =  0.8330
	shift(23) = -0.8895
	shift(24) =  0.2019
	shift(25) = -0.2922

c	Compute index J using lick system (smoothed) sed
c	Apply median shift for the index
	if (j.le.25) then
c		gw_ix_sed_lick_system=gw_ix_sed(j,x,y2,n,igw)-shift(j)
		gw_ix_sed_lick_system=gw_ix_sed(j,x,f,n,igw)-shift(j)
	else
c		Lick system transformation only for indices 1 to 25
		gw_ix_sed_lick_system=gw_ix_sed(j,x,f,n,igw)
	endif
	return
	end
