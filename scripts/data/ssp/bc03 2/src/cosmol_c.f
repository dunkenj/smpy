	real*4 function cosmol_c(h,omega,omega_lambda,q)

c	Returns cosmological constant = cosmol_c and parameter q

c       Omega is entered by the user
c       omega=1.-omega_lambda

c       Cosmological constant
        cosmol_c=omega_lambda/(3.*h**2)

c       Compute q=q0
        if (omega_lambda.eq.0.) then
                q=omega/2.
        else
                q=(3.*omega/2.) - 1.
        endif
	return
	end
