subroutine diffusion_coeff_louis(nlev,theta,qv,u,v,z,km,kh)
!-------------------------------------------------------------------------
!
! Computation of the eddy diffusion coefficients
! from Louis et al. (1981) with modifications described in Coiffier (2011)
!
!                                 Jean-Francois MAHFOUF (01/22)
! 
!--------------------------------------------------------------------------
 use const, only : karman, grav
 use surf1, only : z0, z0h
 implicit none         
 integer,                 intent(in)  :: nlev
 real, dimension(nlev),   intent(in)  :: theta, qv, u, v
 real, dimension(0:nlev), intent(in)  :: z
 real, dimension(nlev),   intent(out) :: km, kh
 real                                 :: zm, lm, lh, dz, du, dv, thetam, dtheta, fm, fh, cm, ch
 real, dimension(nlev)                :: shear, rib, dthdz
!
 real, parameter  :: cons_b = 5.0, cons_c = 5.0, cons_d = 5.0
 real, parameter  :: lambda_m = 150.0, lambda_h = lambda_m*sqrt(1.5*cons_d)
 real, parameter  :: beta = 0.05, Zmax = 1500.0
 integer          :: jk
!
 cm = cons_c
 ch = cons_c
! 
 km(:) = 0.0
 kh(:) = 0.0
 rib(:) = 0.0
 shear(:) = 0.0
 dthdz(:) = 0.0
 do jk=2,nlev
   zm = 0.5*(z(jk) + z(jk-1))
   lm = karman*(zm + z0)/(1.0 + karman*(zm + z0)/lambda_m)*(beta + (1.0 - beta)/(1.0 + ((zm + z0)/Zmax)**2))
   lh = karman*(zm + z0h)/(1.0 + karman*(zm + z0h)/lambda_h)*(beta + (1.0 - beta)/(1.0 + ((zm + z0h)/Zmax)**2)) 
   dz = 1.0/(z(jk) - z(jk-1))
   du = u(jk) - u(jk-1)
   dv = v(jk) - v(jk-1)
   shear(jk) = max(1E-8,sqrt((du*dz)**2 + (dv*dz)**2))
   thetam = 0.5*(theta(jk)*(1.0 + 0.608*qv(jk)) + theta(jk-1)*(1.0 + 0.608*qv(jk-1)))
   dtheta = (theta(jk) - theta(jk-1))
   dthdz(jk) = dtheta*dz
   rib(jk) = grav/thetam*dtheta*dz/shear(jk)**2
   if (rib(jk) > 100.) rib(jk) = 100.0
   if (rib(jk) >  0.0) then
     fm = 1.0/(1.0 + 2.0*cons_b*rib(jk)/sqrt(1.0 + cons_d*rib(jk)))
     fh = 1.0/(1.0 + 3.0*cons_b*rib(jk)*sqrt(1.0 + cons_d*rib(jk)))
   else
     fm = 1.0 - 2.0*cons_b*rib(jk)/(1.0 + 3.0*cons_b*cm*sqrt(abs(rib(jk))/27.0*(lm/(zm + z0))**2))
     fh = 1.0 - 2.0*cons_b*rib(jk)/(1.0 + 3.0*cons_b*ch*sqrt(abs(rib(jk))/27.0*(lh/(zm + z0h))**2))
   endif  
   km(jk) = lm*lm*shear(jk)*fm
   kh(jk) = lm*lh*shear(jk)*fh
 enddo 
 km(1) = km(2)
 kh(1) = kh(2)
end subroutine diffusion_coeff_louis
