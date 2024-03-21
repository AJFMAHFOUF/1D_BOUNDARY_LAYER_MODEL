subroutine water_budget(wg,w2,wr,pr,leg,lev,letr,wgn,w2n,wrn,ro,EmP) 
!---------------------------------------------------------------------------------
!
! Solve (explicitely) the force-restore equations for ws, w2 and wr
!
! Modified to account the evolution of the interception reservoir
! and the corresponding evaporation fluxes and runoff components
! Addition of the net water flux at the surface for the multi-level soil scheme
!
!
!                                 Jean-Francois MAHFOUF (11/06)
!                                                       (10/21)
!                                                       (04/22)
!----------------------------------------------------------------------------------
 use soil
 use const
 use surf1
 use setup, only : dt
 implicit none
 real,  intent(in)  :: pr, leg, lev, letr
 real,  intent(in)  :: wg, w2, wr
 real,  intent(out) :: wgn, w2n, wrn, ro, EmP
 real               :: runoff1, runoff2, runoff3, runoff4
!
! Evolution of the interception reservoir
!
 wrn = wr + dt*(veg*pr - lev/Lv) 
!
 if (wrn > Wrmax) then
  runoff4 = (wrn - Wrmax)/dt
  wrn = Wrmax
 else 
  runoff4 = 0.
 endif     
 if (wrn < 0.) then
  runoff4 = (0.0 - wr)/dt ! conserve water by getting the missing amount in the soil
  wrn = 0.0
 endif    
!
! Evolution of the surface volumetric water content
!
 wgn= (wg +  dt*(C1*((1.-veg)*pr + runoff4 - leg/Lv)/rhow + C2*wgeq/tau))/(1. + C2*dt/tau)
!
 wgn = max(wl,wgn)
 if (wgn > wsat) then
  runoff1 = (wgn - wsat)*d1*rhow
  wgn = wsat
 else
  runoff1 = 0.
 endif 
!
! Evolution of mean soil moisture content
!
 w2n = w2 + dt*((1.-veg)*pr + runoff4 - leg/Lv - letr/Lv)/(d2*rhow) - dt*C3/tau*max(0.,w2 - wfc)  
!
 runoff3 = dt*C3/tau*max(0.,w2 - wfc)*d2*rhow
! 
 w2n = max(wl,w2n)
 if (w2n > wsat) then
  runoff2 = (w2n - wsat)*d2*rhow
  w2n = wsat
 else
  runoff2 = 0.
 endif  
 ro = runoff1 + runoff2 + runoff3 
 !
 ! Net water flux at the surface for the multi-layer soil scheme
 !
 EmP = ((1.-veg)*pr + runoff4 - leg/Lv)/rhow ! - letr/Lv
end subroutine water_budget 
