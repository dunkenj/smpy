	FUNCTION BH_IX_SED(J,X,F,N)

c	Computes value of Brodie and Hanes (1986) index J directly from sed
c	Follows procedure outlined by Brodie and Hanes (1986, ApJ 300, 258)
c	See their eqs. 1-2, and Table 2 in pag. 260.

c	J = Index identification
c	X = wavelength scale
c	F = flux scale (sed)
c	N = dimension of X, F
c	BH_IX_SED = index measured directly from sed.

c	Array declaration
	character*40 idxlbl
	real x(n),f(n),w(4)
	parameter (jid=35)
	integer iim(jid)
	real ffb(jid),ffr(jid),ffc(jid),ffl(jid)
	common /fluxes/ iim,ffb,ffr,ffc,ffl

c	Define wavelength range for selected index
	if (j.eq.1) then
		idxlbl = 'Delta (UV Blanketing)'
		w(1) = 3800.
		w(2) = 4000.
		w(3) = 4000.
		w(4) = 4200.
c		Mean flux through red window
       		fr=trapz3(x,f,n,w(3),w(4),ierr)
c		Mean flux through blue window
       		fb=trapz3(x,f,n,w(1),w(2),ierr)
c		Compute index 1
		bh_ix_sed = 2.5 * alog10(fr/fb)
		return
	elseif (j.eq.2) then
		idxlbl = 'H10 (Balmer line)'
		w(1) = 3785.
		w(2) = 3795.
		w(3) = 3810.
		w(4) = 3825.
	elseif (j.eq.3) then
		idxlbl = 'H9 (Balmer line)'
		w(1) = 3810.
		w(2) = 3825.
		w(3) = 3850.
		w(4) = 3865.
	elseif (j.eq.4) then
		idxlbl = 'Hpsi (Balmer line)'
		w(1) = 3865.
		w(2) = 3885.
		w(3) = 3910.
		w(4) = 3925.
	elseif (j.eq.5) then
		idxlbl = 'CN (UV cyanogen)'
		w(1) = 3785.
		w(2) = 3810.
		w(3) = 3910.
		w(4) = 3925.
	elseif (j.eq.6) then
		idxlbl = 'HK index (Ca II H and K lines)'
		w(1) = 3910.
		w(2) = 3925.
		w(3) = 3995.
		w(4) = 4015.
	elseif (j.eq.7) then
		idxlbl = 'Hdelta (Balmer line)'
		w(1) = 4075.
		w(2) = 4090.
		w(3) = 4125.
		w(4) = 4140.
	elseif (j.eq.8) then
		idxlbl = 'Ca (Ca I line)'
		w(1) = 4200.
		w(2) = 4215.
		w(3) = 4245.
		w(4) = 4260.
	elseif (j.eq.9) then
		idxlbl = 'G (CH G band)'
		w(1) = 4275.
		w(2) = 4285.
		w(3) = 4315.
		w(4) = 4325.
	elseif (j.eq.10) then
		idxlbl = 'Hgamma (Balmer line)'
		w(1) = 4315.
		w(2) = 4325.
		w(3) = 4360.
		w(4) = 4375.
	elseif (j.eq.11) then
		idxlbl = 'Hbeta (Balmer line)'
		w(1) = 4800.
		w(2) = 4830.
		w(3) = 4890.
		w(4) = 4920.
	elseif (j.eq.12) then
		idxlbl = 'Mg (Mg b triplet)'
		w(1) = 5125.
		w(2) = 5150.
		w(3) = 5195.
		w(4) = 5220.
	elseif (j.eq.13) then
		idxlbl = 'MH (Mg b + MgH)'
		w(1) = 4740.
		w(2) = 4940.
		w(3) = 5350.
		w(4) = 5550.
	elseif (j.eq.14) then
		idxlbl = 'FC (Fe I + Ca I)'
		w(1) = 5225.
		w(2) = 5250.
		w(3) = 5280.
		w(4) = 5305.
	elseif (j.eq.15) then
		idxlbl = 'Na (Na D lines)'
		w(1) = 5385.
		w(2) = 5865.
		w(3) = 5920.
		w(4) = 5950.
	elseif (j.eq.16) then
		idxlbl = 'FTC (Fe I + TiO + Ca I)'
		w(1) = 5950.
		w(2) = 6100.
		w(3) = 6330.
		w(4) = 6480.
	elseif (j.eq.17) then
		idxlbl = 'Halpha (Balmer line)'
		w(1) = 6480.
		w(2) = 6615.
		w(3) = 6575.
		w(4) = 6610.
	else
		write (6,*) 'Undefined BH index'
		bh_ix_sed = 999.99
	endif

c	Compute index according to BH definition
c	Mean flux through blue window
       	fb=trapz3(x,f,n,w(1),w(2),ierr)/(w(2)-w(1))
c	Mean flux through red window
       	fr=trapz3(x,f,n,w(3),w(4),ierr)/(w(4)-w(3))
c	Mean flux through central window
       	fc=trapz3(x,f,n,w(2),w(3),ierr)/(w(3)-w(2))
c	Compute index
	bh_ix_sed = 1.0 - 2.0*fc/(fb+fr)

c	Store values for HK index (Ca II H and K lines)
	if (j.eq.6) then
		iim(35)=1
		ffb(35)=fb
		ffr(35)=fr
		ffc(35)=fc
c		fl is not defined for this index
c		ffl(35)=fl
		ffl(35)=0.
	endif
	return
	end
