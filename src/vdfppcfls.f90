subroutine vdfppcfls(zref,ps,ta,ts,ua,va,qa,qs,u10,v10,t2m,q2m,rh2m)
!--------------------------------------------------------------------
!
! Interpolation of atmospheric U, V, T and q at observation level
! formulation from Geleyn (1988) using Louis et al. (1981)
! stability functions of the surface boundary layer generalized
! by Mascart et al. (1995) when ZOM and Z0H are different
!
!
!                                     Jean-Francois MAHFOUF (11/06)
!--------------------------------------------------------------------
 use const
 use surf1
 implicit none
 interface
  real function qsat(p,t)
   implicit none
   real, intent(in)  :: p,t
  end function qsat
  real function esat(t)
   implicit none
   real, intent(in)  :: t
  end function esat
 end interface
 real, intent(in)  :: zref, ps, ta, ts, ua, va, qa, qs
 real, intent(out) :: u10, v10, t2m, q2m, rh2m
 real              :: zrib, z1dz0m, z1dz0h, zxlnm, zxlnh, zcdnm, zcdnh, &
&                     zmu, zchs, zph, zcms, zpm, zcm, zch, zum, zcfm, zcfh, &
&                     zbn, zbnh, zbd, zbh, zru, zrs, zlogu, zlogs, zcoru, zcors, &
&                     z10uiv, z2siv, zcpt2m, z2m, z10m, ztvs, ztva, zeps, zev
!
 z2m  =  2.0
 z10m = 10.0
 zeps = 0.622
!
 z1dz0m = (zref + z0)/z0
 z1dz0h = (zref + z0)/z0h
 zxlnm  = log(z1dz0m)
 zxlnh  = log(z1dz0h)
 zcdnm  = (karman/zxlnm)**2
 zcdnh  = (karman/zxlnh)**2   
 zmu  = log(z0/z0h)
 zchs = 3.2165 + 4.3431*zmu + 0.5360*zmu*zmu - 0.0781*zmu*zmu*zmu
 zph  = 0.5802 - 0.1571*zmu + 0.0327*zmu*zmu - 0.0026*zmu*zmu*zmu
 zcms = 6.8741 + 2.6933*zmu - 0.3601*zmu*zmu + 0.0154*zmu*zmu*zmu
 zpm  = 0.5233 - 0.0815*zmu + 0.0135*zmu*zmu - 0.0010*zmu*zmu*zmu
 zcm  = zcms*z1dz0m**zpm
 zch  = zchs*z1dz0h**zph
 zum = max(0.01,sqrt(ua*ua+va*va))
 ztva = (ta + grav*zref/Cp)*(1. + 0.608*qa)
 ztvs =  ts * (1. + 0.608*qs)
 zrib = 2.*zref*grav*(ztva - ztvs)/((ztvs+ztvs)*zum*zum)
 if (zrib > 0.0) then
   zcfm = zcdnm /(1. + 10.*zrib/sqrt(1. + 5.*zrib))
   zcfh = zcdnh /(1. + 15.*zrib*sqrt(1. + 5.*zrib))
 else                     
   zcfm = zcdnm*(1. - 10.*zrib/(1. + 10.*zcdnm*zcm*sqrt(abs(zrib))))
   zcfh = zcdnh*(1. - 15.*zrib/(1. + 15.*zcdnh*zch*sqrt(abs(zrib))))
 endif   
 zbn   = karman/sqrt(zcdnm)
 zbnh  = karman*sqrt(zcdnm)/zcdnh
 zbd   = karman/sqrt(zcfm)
 zbh   = karman*sqrt(zcfm)/zcfh
 zru   = z10m/zref
 zrs   = z2m/zref
 zlogu = log(1. + zru*(exp(zbn ) - 1.))
 zlogs = log(1. + zrs*(exp(zbnh) - 1.))
 if (zrib > 0.0) then
   zcoru= zru*(zbn-zbd)
   zcors = zrs*(zbnh-zbh)
 else                  
   zcoru = log(1. + zru*(exp(max(0.,zbn -zbd))-1.))
   zcors = log(1. + zrs*(exp(max(0.,zbnh-zbh))-1.))
 endif   
 z10uiv = max(0.,min(1.,(zlogu-zcoru)/zbd))
 z2siv  = max(0.,min(1.,(zlogs-zcors)/zbh))
 u10=ua*z10uiv
 v10=va*z10uiv
 zcpt2m = Cp*ts + (Cp*ta + grav*zref - Cp*ts)*z2siv
 t2m = (zcpt2m - grav*2.0)/Cp
 q2m = qs + (qa - qs)*z2siv
 !zev = q2m/qsat(ps,t2m)
 zev = ps*q2m/((zeps + (1. - zeps)*q2m)*esat(t2m))
 rh2m = max(0.,min(1.0,zev))
 return
end subroutine vdfppcfls
