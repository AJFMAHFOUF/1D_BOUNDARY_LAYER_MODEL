subroutine fluxes(rho,tsn,t2n,tha,qs,qa,rg,rl,ra,rs,rsoil,delta,hu,ct,h,lev,leg,letr,le,rn,g)
 use const
 use surf1
 implicit none
 real, intent(in)  :: rho, tsn, tha, t2n, qa, qs, rg, rl, ra, rs, rsoil, delta, hu, ct
 real, intent(out) :: h, lev, leg, letr, le, rn, g
 h    = rho*Cp*(tsn - tha)/ra
 lev  = rho*Lv*delta*veg*(qs - qa)/ra
 letr = rho*Lv*(1. - delta)*veg*(qs - qa)/(ra + rs)
 leg  = rho*Lv*(1. - veg)*(hu*qs - qa)/(rsoil + ra)
 le   = lev + letr + leg
 rn   = -emis*Stefan*tsn**4 + (1.- alpha)*rg + emis*rl 
 !g  = rn - h - le
 ! g  = 2.0*rpi/(tau*ct)*(tsn - t2n) ! soil heat flux instead of energy balance residual
 g  = Lambda_sk*(tsn - t2n)
 return
end subroutine fluxes
