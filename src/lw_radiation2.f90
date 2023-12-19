subroutine lw_radiation2(nlev,ts,qvs,ps,tha,qva,pa,dtdt)
!-------------------------------------------------------------------------
!
! Compute comprehensive longwave radiative cooling 
!
!                                 Jean-Francois MAHFOUF (02/22)
! 
!-------------------------------------------------------------------------
 use const, only : Stefan, grav, Cp, Rd, p00
 implicit none    
 interface
  real function emis1_co2(u)
   implicit none
   real, intent(in)  :: u
  end function emis1_co2
  real function emis1_h2o(u)
   implicit none
   real, intent(in)  :: u
  end function emis1_h2o
  real function emis2_h2o(u)
   implicit none
   real, intent(in)  :: u
  end function emis2_h2o
  real function emis3_h2o(u)
   implicit none
   real, intent(in)  :: u
  end function emis3_h2o
  real function emis4_h2o(u)
   implicit none
   real, intent(in)  :: u
  end function emis4_h2o
 end interface     
 integer,                 intent(in)  :: nlev
 real,                    intent(in)  :: ts, qvs, ps
 real, dimension(nlev),   intent(in)  :: tha, qva, pa
 real, dimension(nlev),   intent(out) :: dtdt

 real, dimension(nlev)      :: ta, tah, qvah   
 real, dimension(nlev+1)    :: zpa, zta, pah, flux_u, flux_d      
 real                       :: emis_totm, emis_totp, invcpdp, zp0, zt0, &
                            &  total_column_h2o, total_column_co2, delta_flux, &
                            &  zscale_h2o, zscale_co2, du_h2o, du_co2
 integer                    :: jk, jk1, jk2
 logical                    :: l_wvcont
!
!
 l_wvcont=.true.
!
! Reference values for scaling optical path for H2O and CO2
!
 zp0 = 1.0E5
 zt0 = 273.0
!
 do jk=1,nlev
   ta(jk) = tha(jk)*(pa(jk)/p00)**(Rd/Cp)
 enddo  
! 
 do jk=1,nlev
   zpa(jk) = pa(jk)
   zta(jk) = ta(jk)
 enddo
 zpa(nlev+1) = ps
 zta(nlev+1) = ts
!
! Define half level values for p, T and qv
!   
 do jk=1,nlev
  if (jk /= nlev) then
    tah(jk)  = 0.5*(ta(jk) + ta(jk+1))
    qvah(jk) = 0.5*(qva(jk) + qva(jk+1))
    pah(jk)  = 0.5*(pa(jk) + pa(jk+1))
  else
    tah(jk)  = 0.5*(ta(jk) + ts)
    qvah(jk) = 0.5*(qva(jk) + qvs)
    pah(jk)  = 0.5*(pa(jk) + ps)
  endif 
 enddo
!
! Downward longwave flux
! 
 flux_d(1) = 0.0 
 do jk1=1,nlev
   total_column_h2o = 0.0
   total_column_co2 = 0.0
   emis_totm = 0.0
   delta_flux = 0.0
   do jk2=jk1,1,-1 
     zscale_h2o = (pah(jk2)/zp0)**1.20*(zt0/tah(jk2))**0.5
     zscale_co2 = (pah(jk2)/zp0)**0.75*(zt0/tah(jk2))**0.0
     du_h2o = 0.1/grav*qvah(jk2)*zscale_h2o*(zpa(jk2+1) - zpa(jk2))
     total_column_h2o = total_column_h2o + du_h2o
     du_co2 = 0.00612*zscale_co2*(zpa(jk2+1) - zpa(jk2))
     total_column_co2 = total_column_co2 + du_co2
     emis_totp = emis1_co2(total_column_co2) + emis3_h2o(total_column_h2o)
     delta_flux = delta_flux + Stefan*tah(jk2)**4*(emis_totp - emis_totm)
     emis_totm = emis_totp
   enddo
   flux_d(jk1+1) = delta_flux ! + Stefan*ta(1)**4*emis_totm
 enddo  
!
! Updward longwave flux
! 
 do jk1=nlev,1,-1
   total_column_h2o = 0.0
   total_column_co2 = 0.0
   emis_totm = 0.0
   delta_flux = 0.0   
   do jk2=jk1,nlev
     zscale_h2o = (pah(jk2)/zp0)**1.20*(zt0/tah(jk2))**0.5
     zscale_co2 = (pah(jk2)/zp0)**0.75*(zt0/tah(jk2))**0.0
     du_h2o = 0.1/grav*qvah(jk2)*zscale_h2o*(zpa(jk2+1) - zpa(jk2))
     total_column_h2o = total_column_h2o + du_h2o
     du_co2 = 0.00612*zscale_co2*(zpa(jk2+1) - zpa(jk2))
     total_column_co2 = total_column_co2 + du_co2
     emis_totp = emis1_co2(total_column_co2) + emis3_h2o(total_column_h2o)
     delta_flux = delta_flux + Stefan*tah(jk2)**4*(emis_totp - emis_totm)
     emis_totm = emis_totp
   enddo
   flux_u(jk1) = delta_flux + Stefan*ts**4*(1.0 - emis_totm)
 enddo 
 flux_u(nlev+1) = Stefan*ts**4 ! not used for radiative tendencies   
! 
!  Radiative cooling at each level
!  
 do jk=2,nlev
   invcpdp = grav/(Cp*(pah(jk) - pah(jk-1)))
   dtdt(jk) = invcpdp*(flux_u(jk) - flux_u(jk-1) + flux_d(jk-1) - flux_d(jk))    
!
!  Empirical correction fo water vapour continuum (Savijarvi, 1990)
!   
   if (l_wvcont) then
     dtdt(jk) = dtdt(jk) - (0.5*1.E-3*(1.E3*qva(jk))**3 + 0.1)/86400.0  
   endif
 enddo
! 
!   Empirical extrapolation for the upper most level
!
 dtdt(1)  = dtdt(2) + (pa(1) - pa(2))/(pa(2) - pa(3))*(dtdt(2) - dtdt(3)) 
!
 return
end subroutine lw_radiation2
