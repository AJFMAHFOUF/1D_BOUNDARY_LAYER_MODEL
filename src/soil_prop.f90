subroutine soil_prop(wg,w2,ts)
!--------------------------------------------------------
!
! Estimation of soil thermal and hydraulic properties
! for the resolution of force-restore equations
!
!--------------------------------------------------------
 use soil
 implicit none
 real, intent(in) :: wg, w2, ts        
 real             :: zeta, zwmax, zc1max, zsigma2, zx 
!
 zeta = (-1.815E-2*ts + 6.41)*wwilt+(6.5E-3*ts - 1.4)
 zwmax = zeta*wwilt
 zc1max = (1.19*wwilt - 5.09)*ts*0.01 + (1.464*wwilt + 17.86)
 zsigma2 = -zwmax**2/(2.0*log(0.01/zc1max))
 if (wg < wwilt) then 
   C1 = zc1max*exp(-(wg - zwmax)**2/(2.0*zsigma2))/0.01
 else
   C1 = C1sat*(wsat/wg)**(0.5*b + 1.0)
 endif
 C2 = C2ref*(w2/(wsat - w2 + wl))
 zx = w2/wsat
 wgeq = wsat*(zx - a*(zx**p*(1.0 - zx**(8.*p))))
 Cg = Cgsat*zx**(-b/(2.*log(10.))) 
!
 return
end subroutine soil_prop
