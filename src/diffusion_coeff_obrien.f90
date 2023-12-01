subroutine diffusion_coeff_obrien(nlev,zi,ustar,phi_m,phi_h,theta,u,v,z,km,kh)
!-------------------------------------------------------------------------
!
! Computation of the eddy diffusion coefficients from O'Brien (1970) 
!
!                                 Jean-Francois MAHFOUF (01/22))
! 
!--------------------------------------------------------------------------
 use const, only : karman
 implicit none         
 integer,                 intent(in)  :: nlev
 real,                    intent(in)  :: zi, ustar, phi_m, phi_h 
 real, dimension(nlev),   intent(in)  :: theta, u, v
 real, dimension(0:nlev), intent(in)  :: z
 real, dimension(nlev),   intent(out) :: km,kh
 real             :: z_nk, zki, zkm, dkzm, zkh, dkzh, dths
!
 integer          :: jk
!
 z_nk = z(nlev-1)
! 
 zki = 1.0E-4
 km(:) = 1.0E-4
 kh(:) = 1.0E-4
! 
 zkm = karman*ustar*z_nk/phi_m
 dkzm =karman*ustar/phi_m
!
 zkh = karman*ustar*z_nk/phi_h
 dkzh =karman*ustar/phi_h
! 
 do jk=1,nlev
   km(jk) = zki + ((zi - z(jk-1))/(zi - z_nk))**2 &
             &*(zkm - zki + (z(jk-1) - z_nk)*(dkzm + 2.*(zkm - zki)/(zi - z_nk)))    
   km(jk) = 3.*km(jk)          
   if (z(jk) > zi) km(jk) = zki

   kh(jk) = zki + ((zi - z(jk-1))/(zi - z_nk))**2 &
             &*(zkh - zki + (z(jk-1) - z_nk)*(dkzh + 2.*(zkh - zki)/(zi - z_nk)))
   kh(jk) = 3.*kh(jk)
   if (z(jk) > zi) kh(jk) = zki
 enddo 
 return
end subroutine diffusion_coeff_obrien
