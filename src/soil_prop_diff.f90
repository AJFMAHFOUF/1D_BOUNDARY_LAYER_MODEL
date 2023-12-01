subroutine soil_prop_diff(nlevs,wsoil,K_w,D_w,rhocs,Lambda_s)
!--------------------------------------------------------
!
! Estimation of soil thermal and hydraulic properties
! for the resolution of soil diffusion equations
!
!--------------------------------------------------------
 use soil
 implicit none
 integer, intent(in)                  :: nlevs
 real, dimension (nlevs), intent(in)  :: wsoil  
 real, dimension (nlevs), intent(out) :: K_w, D_w, rhocs, Lambda_s 
 integer :: j
 real    :: zke, zlambda_sat
! 
 do j=1,nlevs
   K_w(j) = Kwsat*(wsoil(j)/wsat)**(2.0*b + 3.0)
   D_w(j) = -b*Kwsat*Psis/wsat*(wsoil(j)/wsat)**(b + 2.0)
   rhocs(j) = rhoc_soil*(1.0 - wsat) + 4.186E6*wsoil(j)
   if (wsoil(j)/wsat > 0.1) then
     zke = alog10(wsoil(j)/wsat) + 1.0
   else
     zke = 0.0
   endif
   zlambda_sat = lsoil**(1.0 - wsat)*lwater**(wsoil(j))
   Lambda_s(j) = zke*(zlambda_sat - Lambda_dry) + Lambda_dry       
 enddo 
!
 return
end subroutine soil_prop_diff
