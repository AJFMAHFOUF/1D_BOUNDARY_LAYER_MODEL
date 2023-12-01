subroutine soilt_vertical_diffusion(nlevs,x0,zs,rhocs,Lambda_s,tsk,x1)
!-------------------------------------------------------------------------
!
! Solve soil heat diffusion equation by inversion of tridiagonal matrix
!
!                                 Jean-Francois MAHFOUF (04/22)
! 
!-------------------------------------------------------------------------
 use setup, only : dt
 use surf1, only : Lambda_sk
 implicit none         
 integer,                  intent(in)  :: nlevs
 real,                     intent(in)  :: tsk
 real, dimension(nlevs),   intent(in)  :: x0, rhocs, Lambda_s
 real, dimension(nlevs),   intent(in)  :: zs
 real, dimension(nlevs),   intent(out) :: x1

 real, dimension(nlevs)   :: a0, b0, c0, d0, Lambda_sm  
 real, dimension(nlevs)   :: ai, bi, ci 
 real, dimension(nlevs)   :: ae, be, ce  
 real, dimension(nlevs+1) :: zsm     
 integer                  :: jk
 real, parameter          :: zbeta = 1.5 ! impliciteness factor 
!
! Definition of half level grid where fluxes are estimated
!
 zsm(1) = 0.0 ! first half level at ground surface
 do jk=2,nlevs
   zsm(jk) = 0.5*(zs(jk) + zs(jk-1))
 enddo
! last flux level extrapolated below last model level
 zsm(nlevs+1) =  2.0*zs(nlevs) - zsm(nlevs) 
!
! Mean thermal conductivity (at flux levels)
!
 do jk=2,nlevs
  Lambda_sm(jk) = 0.5*(Lambda_s(jk) + Lambda_s(jk-1))
 enddo
 Lambda_sm(1) = Lambda_s(1)
!
! Fill the elements of the tridiagonal matrix
!
 a0(1) = 0.0
 ae(1) = 0.0
 ai(1) = 0.0
 do jk=2,nlevs
   a0(jk) = -dt*Lambda_sm(jk)/(rhocs(jk)*(zs(jk) - zs(jk-1))*(zsm(jk+1) - zsm(jk)))
   ae(jk) = (zbeta - 1.0)*a0(jk)
   ai(jk) = zbeta*a0(jk)
 enddo 
! 
 c0(nlevs) = 0.0
 ce(nlevs) = 0.0
 ci(nlevs) = 0.0
 do jk=1,nlevs-1
   c0(jk) = -dt*Lambda_sm(jk+1)/(rhocs(jk)*(zs(jk+1) - zs(jk))*(zsm(jk+1) - zsm(jk)))
   ce(jk) = (zbeta - 1.0)*c0(jk)
   ci(jk) = zbeta*c0(jk)
 enddo
! 
 do jk=1,nlevs
   b0(jk) = 1.0 - a0(jk) - c0(jk)
   be(jk) = 1.0 - ae(jk) - ce(jk)
   bi(jk) = 1.0 - ai(jk) - ci(jk)
 enddo
!
 do jk=2,nlevs-1
   d0(jk) = ae(jk)*x0(jk-1) + be(jk)*x0(jk) + ce(jk)*x0(jk+1)
 enddo 
 d0(nlevs) = ae(nlevs)*x0(nlevs-1) + be(nlevs)*x0(nlevs)
!
! Surface boundary condition with skin temperature (similar to ECMWF) 
!
 b0(1) = 1.0 - a0(1) - c0(1) + dt*Lambda_sk/(rhocs(1)*(zsm(2) - zsm(1)))
! 
 d0(1) = be(1)*x0(1) + ce(1)*x0(2) + dt*Lambda_sk*tsk/(rhocs(1)*(zsm(2) - zsm(1))) 
!
 call inv_mat_tri(nlevs,a0,b0,c0,d0,x1)
!   
 return
end subroutine soilt_vertical_diffusion
