subroutine drag_coeff_z0h(zref,tha,ts,qa,qs,ua,va,ra,ustar)
!-------------------------------------------------------------------------
!
! Computation of the surface aerodynamic resistance
! from Louis et al. (1981) stability functions generalized
! by Mascart et al. (1995) when z0m /= ZOh
!
!                                 Jean-Francois MAHFOUF (11/06)
!
! Modification (JFM 03/07) : virtual temperature in Rib 
!              (JFM 06/07) : inclusion of momentum coefficients to get u*
!--------------------------------------------------------------------------
 use const
 use surf1
 implicit none         
 real, intent(in)  :: zref, tha, qa, ts, qs, ua, va
 real, intent(out) :: ra, ustar
 real              :: zcd, zum, zrib, zcorh, zmu, zchs, zcms, &
&                     zcorm, zpm, zph,                        &
&                     zch, zcm, zratio, ztvs, ztva
 real, parameter   :: zcons1=10.0, zcons2=5.0 ! zcons1 = 10 - zcons2 = 1 in Viterbo et al. (1999)
!
 zmu  = log(z0/z0h)
 zchs = 3.2165 + 4.3431*zmu + 0.5360*zmu*zmu - 0.0781*zmu*zmu*zmu
 zph  = 0.5802 - 0.1571*zmu + 0.0327*zmu*zmu - 0.0026*zmu*zmu*zmu
 zcms = 6.8741 + 2.6933*zmu - 0.3601*zmu*zmu + 0.0154*zmu*zmu*zmu
 zpm  = 0.5233 - 0.0815*zmu + 0.0135*zmu*zmu - 0.0010*zmu*zmu*zmu 
 zcd = (karman/log((zref + z0)/z0))**2
 zratio = log((zref + z0)/z0)/log((zref + z0)/z0h)
 zch  = 15.*zchs*zcd*((zref + z0)/z0h)**zph*zratio
 zcm  = 10.*zcms*zcd*((zref + z0)/z0h)**zpm
 zum  = max(0.01,sqrt(ua*ua+va*va))
 ztva = tha*(1. + 0.608*qa)
 ztvs =  ts*(1. + 0.608*qs)
 zrib = 2.*grav*zref*(ztva - ztvs)/((ztvs + ztva)*zum*zum)
 if (zrib > 0.) then 
   zcorh = (1./(1. + 1.5*zcons1*zrib*sqrt(1. + zcons2*zrib)))*zratio
   zcorm = (1./(1. + zcons1*zrib/sqrt(1. + zcons2*zrib)))
 else
   zcorh = (1. -  1.5*zcons1*zrib/(1. + zch*sqrt(abs(zrib))))*zratio
   zcorm =  1. -  zcons1*zrib/(1. + zcm*sqrt(abs(zrib)))
 endif
 ra = 1./(zcd * zum * zcorh)
 ustar = sqrt(zcd*zcorm)*zum
 return
end subroutine drag_coeff_z0h
