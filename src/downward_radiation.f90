subroutine downward_radiation(ta,qa,ps,clf,pmu,rg,rl)
!-------------------------------------------------------------------------
!
! Downward shortwave and longwave radiation fluxes at the surface
!
!
!                                         Jean-Francois MAHFOUF (01/22)
!                                                             
!-------------------------------------------------------------------------
 use const, only : S0, Stefan
 use setup, only : Trans
 implicit none
 real, intent(in)  :: ta, qa, ps, clf, pmu
 real, intent(out) :: rg, rl
 real              :: zea, zemis1, zemis2, zemis3
! 
 if (pmu > 0) then
   rg = S0*pmu*(Trans)**(1./pmu)
 else
   rg = 0.0
 endif
!
 zea = qa*ps/(0.378*qa + 0.622)
 !
 ! Effective emissivity
 !
 zemis1 = (clf + (1.0 - clf)*(9.36E-6*ta**2))          ! Swinbank (1963)
 zemis2 = (clf + (1.0 - clf)*(1.24*(zea/ta)**(1./7.))) ! Brustaert (1975)
 zemis3 = (clf + (1.0 - clf)*(0.67*(1670*qa)**(0.08))) ! Staley & Jurica (1972)
 !
 rl = zemis1*Stefan*ta**4
!
 return
end subroutine downward_radiation
