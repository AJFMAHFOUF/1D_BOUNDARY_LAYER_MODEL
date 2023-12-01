subroutine soilw_vertical_diffusion(nlevs,x0,zs,K_w,D_w,EmP,Root_ext,x1)
!-------------------------------------------------------------------------
!
! Solve soil water diffusion equation by inversion of tridiagonal matrix
!
!                                 Jean-Francois MAHFOUF (04/22)
! 
!-------------------------------------------------------------------------
 use setup, only : dt
 use soil, only : wsat, Kwsat, b, Psis
 implicit none         
 integer,                  intent(in)  :: nlevs
 real,                     intent(in)  :: EmP, Root_ext
 real, dimension(nlevs),   intent(in)  :: x0, K_w, D_w
 real, dimension(nlevs),   intent(in)  :: zs
 real, dimension(nlevs),   intent(out) :: x1

 real, dimension(nlevs)   :: a0, b0, c0,  d0  
 real, dimension(nlevs)   :: ai, bi, ci 
 real, dimension(nlevs)   :: ae, be, ce  
 real, dimension(nlevs)   :: K_wm, D_wm, S
 real, dimension(nlevs+1) :: zsm 
 integer                  :: jk
 real, parameter          :: zbeta = 1.5 ! impliciteness factor
 real                     :: infiltr_max, EmP_max, runoff_s, Sn
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
!  Empirical root extraction profile
! 
 Sn = 0.0
 S(:) = 0.0
 do jk=2,nlevs
   Sn = Sn + exp(-2.5*(zs(jk) - zs(2)))*D_w(jk)
   S(jk) = Root_ext*exp(-2.5*(zs(jk) - zs(2)))*D_w(jk)
 enddo  
 S(:) = S(:)/Sn
!
! Mean hydraulic properties (at flux levels)
!
 do jk=2,nlevs
   K_wm(jk) = 0.5*(K_w(jk) + K_w(jk-1))
   D_wm(jk) = 0.5*(D_w(jk) + D_w(jk-1))
!  K_wm(jk) = sqrt(K_w(jk)*K_w(jk-1))
!  D_wm(jk) = sqrt(D_w(jk)*D_w(jk-1))
   K_wm(jk) = amax1(K_w(jk),K_w(jk-1))
   D_wm(jk) = amax1(D_w(jk),D_w(jk-1))
 enddo
 D_wm(1) = D_w(1)
 K_wm(1) = K_w(1)
!
! Fill the elements of the tridiagonal matrix
!
 a0(1) = 0.0
 ae(1) = 0.0
 ai(1) = 0.0
 do jk=2,nlevs
   a0(jk) = -dt*D_wm(jk)/((zs(jk) - zs(jk-1))*(zsm(jk+1) - zsm(jk)))
   ae(jk) = (zbeta - 1.0)*a0(jk)
   ai(jk) = zbeta*a0(jk)
 enddo 
! 
 c0(nlevs) = 0.0
 ce(nlevs) = 0.0
 ci(nlevs) = 0.0 
 do jk=1,nlevs-1
   c0(jk) = -dt*D_wm(jk+1)/((zs(jk+1) - zs(jk))*(zsm(jk+1) - zsm(jk)))
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
   d0(jk) = ae(jk)*x0(jk-1) + be(jk)*x0(jk) + ce(jk)*x0(jk+1) - dt*(K_wm(jk+1) - K_wm(jk))/(zsm(jk+1) - zsm(jk)) - dt*S(jk)
 enddo
 d0(nlevs) = ae(nlevs)*x0(nlevs-1) + be(nlevs)*x0(nlevs)  ! + dt*K_wm(nlevs)/(zsm(nlevs+1) - zsm(nlevs))
! 
! Surface boundary condition for water flux (precipitation minus evaporation)
!
! Maximum soil infiltration proposed by Mahrt and Pan (1984)
! Useless since it leads to too high values => threshold never reached
! 
 infiltr_max = (-2.0*b*Kwsat*Psis/wsat*(wsat - x0(1))/(zsm(2) - zsm(1)) + Kwsat)
 EmP_max = amin1(EmP,infiltr_max)
 runoff_s = (EmP - EmP_max)*dt*1000.0 ! runoff expressed in mm 
 d0(1) = be(1)*x0(1) + ce(1)*x0(2) + dt*(EmP_max - K_wm(2))/((zsm(2) - zsm(1))) 
!
 call inv_mat_tri(nlevs,a0,b0,c0,d0,x1)
!
! Superficial runoff in case of saturation of first soil layer
! 
 x1 = amax1(0.05,x1)
 if (x1(1) > wsat) then
   runoff_s = runoff_s + (x1(1) - wsat)*zsm(2)*1000.
   x1(1) = wsat
 endif
!
 return
end subroutine soilw_vertical_diffusion
