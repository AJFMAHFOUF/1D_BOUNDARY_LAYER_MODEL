subroutine vertical_diffusion(nlev,x0,z,rho,kdiff,flux,x1,gamma_cg,icg)
!-------------------------------------------------------------------------
!
! Solve vertical diffusion equation by inversion of tridiagonal matrix
!
!                                 Jean-Francois MAHFOUF (01/22)
! 
!-------------------------------------------------------------------------
 use setup, only : dt
 implicit none         
 integer,                 intent(in)  :: nlev
 real,                    intent(in)  :: flux, gamma_cg
 real, dimension(nlev),   intent(in)  :: x0, rho, kdiff
 integer,                 intent(in)  :: icg ! countergradient term for theta only
 real, dimension(nlev+1), intent(in)  :: z
 real, dimension(nlev),   intent(out) :: x1

 real, dimension(nlev)   :: a, b, c, d, rhoh   
 real, dimension(nlev)   :: ai, bi, ci 
 real, dimension(nlev)   :: ae, be, ce     
 real, dimension(nlev+1) :: zm      
 integer                 :: jk
 real, parameter         :: beta = 1.0 ! impliciteness factor
!
 do jk=2,nlev
  rhoh(jk) = 0.5*(rho(jk) + rho(jk-1))
 enddo
 rhoh(1) = rhoh(2) 
!
 do jk=2,nlev
   zm(jk) = 0.5*(z(jk) + z(jk-1))
 enddo  
 zm(1) = z(1)
 zm(nlev+1) = z(nlev+1)
!
 a(1) = 0.0
 ae(1) = 0.0
 ai(1) = 0.0
 do jk=2,nlev
   a(jk) = 2.0*dt*kdiff(jk)/((zm(jk) - zm(jk+1))*(zm(jk+1) - zm(jk-1)))
   ae(jk) = (beta - 1.0)*a(jk)
   ai(jk) = beta*a(jk)
 enddo 
! 
 c(nlev) = 0.0
 ce(nlev) = 0.0
 ci(nlev) = 0.0 
 do jk=1,nlev-1
   c(jk) = 2.0*dt*kdiff(jk+1)/((zm(jk) - zm(jk+1))*(zm(jk+2) - zm(jk)))
   ce(jk) = (beta - 1.0)*c(jk)
   ci(jk) = beta*c(jk)
 enddo
! 
 do jk=1,nlev
   b(jk) = 1.0 - a(jk) - c(jk)
   be(jk) = 1.0 - ae(jk) - ce(jk)
   bi(jk) = 1.0 - ai(jk) - ci(jk)
 enddo
! 
 do jk=2,nlev-1
   d(jk) = ae(jk)*x0(jk-1) + be(jk)*x0(jk) + ce(jk)*x0(jk+1) - &
         & float(icg)*dt*(kdiff(jk) - kdiff(jk+1))*gamma_cg/(zm(jk) - zm(jk+1))
 enddo
 d(nlev) = ae(nlev)*x0(nlev-1) + be(nlev)*x0(nlev) + dt*flux/((zm(nlev) - zm(nlev+1))) 
 d(1) = be(1)*x0(1) + ce(1)*x0(2) 
!
 call inv_mat_tri(nlev,a,b,c,d,x1)
!   
 return
end subroutine vertical_diffusion
