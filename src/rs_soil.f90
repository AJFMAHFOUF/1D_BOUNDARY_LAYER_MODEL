subroutine rs_soil(wg,veg,rsoil)
!--------------------------------------------------------------------------
!
! Computation of a soil surface resistance according
! to ECMWF formulation
! Replaces the Hu formulation from ISBA (to avoid dew flux pb)
! Minimum threshold modified (wwilt -> veg*wwilt) - Albergel et al. (2012)
!
!                                            Jean-Francois MAHFOUF (11/06)
!                                            Modified (10/21)          
!--------------------------------------------------------------------------
 use soil
 implicit none  
 real, intent(in)  :: wg, veg
 real, intent(out) :: rsoil                   
 real              :: zf2   
!
 zf2 = (wg - veg*wwilt)/(wfc - veg*wwilt)
 zf2 = min(1.0,max(1.0E-4,zf2))
 rsoil = 100./zf2
 return
end subroutine rs_soil
