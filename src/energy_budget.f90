subroutine energy_budget(rho,ts,t2,tsk,tsoil1,ra,rs,rsoil,hu,delta,ct,ps,tha,qa,rg,rl,wg,tsn,t2n,tskn)
!-------------------------------------------------------------------------
!
! Solve (implicitely) the force-restore equations for Ts and T2
! Modified to account for Hu formulation of bare soil evaporation 
! and evaporation from interception reservoir
!
!
!                                         Jean-Francois MAHFOUF (11/06)
!                                                               (10/21)
!-------------------------------------------------------------------------
 use surf1
 use const
 use soil
 use setup, only : dt
 implicit none
 interface
  real function qsat(p,t)
   implicit none
   real, intent(in)  :: p,t
  end function qsat
  real function dqsat(p,t)
   implicit none
   real, intent(in)  :: p,t
  end function dqsat
 end interface
 real, intent(in) :: rho, ts, t2, tsk, tsoil1, ra, rs, ps, tha, qa, rg, rl, wg, delta, ct
 real, intent(inout) :: rsoil, hu
 real, intent(out)   :: tsn, t2n, tskn
 real :: zbeta1, zbeta2, zqs, zdqsdt, za, zb, zc, ztau2, zhu, zra
 real :: zres, zresp
 integer :: niter
! 
 ztau2 = tau/(2.0*rpi)
!
 zdqsdt = dqsat(ps,ts)
 zqs = qsat(ps,ts)                
! 
 if (wg < wfc) then
   zhu = 0.5*(1.0 - cos(wg/wfc*rpi))
 else
   zhu = 1.0
 endif
 if (zqs < qa) then ! Dew flux at potential rate
  zhu = 1.0
  rsoil = 0.0
 endif
 if (zhu*zqs < qa .and. zqs > qa) then ! Avoid dew flux on very dry and warm soils
  zhu = qa/zqs
 endif
 hu = zhu
 zbeta1 = (1. - delta)*veg/(rs + ra) + delta*veg/ra + (1. - veg)*zhu/(rsoil + ra)
 zbeta2 = (1. - delta)*veg/(rs + ra) + delta*veg/ra + (1. - veg)/(rsoil + ra)
!
! Solve implicitly the surface energy balance (for ts)
!
 za  = -ct*(4.*emis*Stefan*ts**3 + rho*(Cp/ra + Lv*zbeta1*zdqsdt))
 zb  =  ct*(3.*emis*Stefan*ts**3 + rho*Lv*zbeta1*zdqsdt)
 zc  =  ct*((1.-alpha)*rg + emis*rl + rho*(Cp*tha/ra + Lv*(zbeta2*qa - zbeta1*zqs)))
 tsn = ((1. + zb*dt)*ts + dt*t2/ztau2 + zc*dt)/(1. - za*dt + dt/ztau2) 
!
! Solve surface energy budget for skin temperature (tsk) with Newton iteration loop
!
 tskn = tsk
 do niter = 1,2
   zres  = Lambda_sk*tskn + emis*Stefan*tskn**4 + rho*Cp*tskn/ra + rho*Lv*zbeta1*qsat(ps,tskn)
   zres  = zres - (Lambda_sk*t2 + (1.-alpha)*rg + emis*rl + rho*Cp*tha/ra + rho*Lv*zbeta2*qa)
   zresp = Lambda_sk   + 4.0*emis*Stefan*tskn**3 + rho*Cp/ra   + rho*Lv*zbeta1*dqsat(ps,tskn)
   tskn   = tskn - zres/zresp
 enddo  
!
!
! Solve explicitly the surface energy balance (for ts)
! 
! tsn = ts + ct*dt*((1.-alpha)*rg + emis*rl - emis*Stefan*ts**4 - rho*Cp*(ts - tha)/ra - &
!                  &                rho*Lv*(zbeta1*zqs - zbeta2*qa)) - dt*(ts - t2)/ztau2
 
! tsn = 2*tsn - 1*ts ! over implicit
!
! Evolution of mean surface temperature (t2)
!
 t2n = (t2 + tsn*dt/tau)/(1. + dt/tau)
 return
end subroutine energy_budget
