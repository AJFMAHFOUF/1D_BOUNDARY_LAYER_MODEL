subroutine pbl_height(ustar,fluxh,ths,qvs,tha_nlev,za_nlev,zi,phi_m,phi_h,wstar,gamma_cg)
!----------------------------------------------------------------------------
!
! Computation of surface boundary layer stability functions according to
! Monin-Obukhov similarity theory and planetary boundary layer height
! according to Deardorff (1974) equation and used in Mahrer and Pielke (1975)
!
!                                 Jean-Francois MAHFOUF (01/22)
! 
!----------------------------------------------------------------------------
 use setup, only : dt, f
 use const, only : karman, grav
 implicit none
 real, intent(in)    :: ustar, fluxh, ths, qvs, tha_nlev, za_nlev
 real, intent(inout) :: zi
 real, intent(out)   :: phi_m, phi_h, gamma_cg, wstar
 real, parameter     :: dthadzp = 8.5E-3  ! potential temperature gradient above inversion
 real                :: lmo, zeta
!
!  Monin-Obukhov length   
!   
  lmo = -ustar**3/(grav*karman/(ths*(1.0 + 0.608*qvs))*fluxh)
!
  zeta = za_nlev/lmo
!
!  Stability functions (Dyer, 1974) - unstable regime
!    
  phi_m = (1.0 - 16.0*zeta)**(-0.25)
  phi_h = (1.0 - 16.0*zeta)**(-0.50)   
!  
!  Evolution of planetary boundary layer depth  (Deardorff, 1974)
!   
  wstar = (grav/tha_nlev*fluxh*zi)**(1./3.)
  zi = zi + dt*(1.8*(wstar**3 + 1.1*ustar**3 - 3.3*ustar**2*abs(f)*zi) &
 &     /(grav*zi*zi/tha_nlev*dthadzp + 9.0*wstar**2 + 7.2*ustar*ustar))   
!
!  Countergradient for potential temperature flux (Therry Lacarrère 1983)
!   
   gamma_cg = 5.0*fluxh/(wstar*zi) 
   
 return
end subroutine pbl_height

