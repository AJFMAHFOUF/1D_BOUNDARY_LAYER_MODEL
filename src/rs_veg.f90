subroutine rs_veg(w2,ps,qa,ta,rg,rs) 
!-------------------------------------------------------------
!
! Computation of the canopy surface resistance RS
! according to Noilhan and Planton (1989) and
! Noilhan and Mahfouf (1996) except for the
! saturation water vapor deficit dependency 
!
!                                Jean-Francois MAHFOUF (11/06)   
!-------------------------------------------------------------
 use surf1
 use soil
 implicit none 
 interface
  real function qsat(p,t)
   implicit none
   real, intent(in)  :: p,t
  end function qsat
 end interface
 real, intent(in)  :: w2, ps, qa, ta, rg
 real, intent(out) :: rs
 real              :: zz, zf, zf1, zf2, zf3, zf4
 zf  = 1.1 * rg / (Rgl * LAI)
 zf1 = (1. + zf)/(zf + Rsmin/5000.)
 zf2 = (w2 - wwilt)/(wfc - wwilt)
 zf2 = min(1.,max(1.0E-4,zf2))
 zf3 = 1.0 + gamma*sqrt(max(0.,qsat(ps,ta) - qa))
 zz  = max(0.,qsat(ps,ta) - qa)
 zf3 = 1./(1.0 - gamma*zz)
 if (zf3 < 0.0) then
   zf3 = 5000.0 ! saturation deficit too large
   print *,' *** WARNING SATURATION DEFICIT TOO LARGE ***',zz
   stop
 endif
 zf4 = 1. - 0.0016*(298.0 - ta)**2
 rs = (Rsmin/LAI)*zf1*zf3/(zf2*zf4)
 return
end subroutine rs_veg
