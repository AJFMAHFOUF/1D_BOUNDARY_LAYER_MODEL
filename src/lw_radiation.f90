subroutine lw_radiation(nlev,ts,qvs,ps,tha,qva,pa,rho,z,dtdt)
!-------------------------------------------------------------------------
!
! Compute simplified longwave radiative cooling (Sasamori, 1972)
!
!                                 Jean-Francois MAHFOUF (01/22)
!                                 modified by JFM (12/23)
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
 end interface     
 integer,                 intent(in)  :: nlev
 real,                    intent(in)  :: ts, qvs, ps
 real, dimension(nlev),   intent(in)  :: tha, qva, pa, rho
 real, dimension(nlev+1), intent(in)  :: z
 real, dimension(nlev),   intent(out) :: dtdt

 real, dimension(nlev)   :: rhoh, ta, tah, qvah, dtdth
 real, dimension(nlev)   :: uh2o_d, uco2_d, uh2o_u, uco2_u      
 real, dimension(nlev+1) :: zpa      
 real                    :: zrhos, emis_up1, emis_up2, emis_dn1, emis_dn2, invcpdz
 real                    :: zt_top, zp0, zt0, zscale_h2o, zscale_co2 
 integer                 :: jk, jk1, jk2
 logical                 :: l_cts, l_wvcont, l_simtend
!
! Options for computing radiative tendencies
! 
 l_cts = .false.
 l_wvcont = .true.
 l_simtend = .false.
!
! Reference values for scaling optical path for H2O and CO2
!
 zp0 = 1.0E5
 zt0 = 273.0
! 
 do jk=1,nlev-1
  rhoh(jk) = 0.5*(rho(jk) + rho(jk+1))
 enddo
 zrhos = ps/(Rd*ts*(1.0 + 0.608*qvs))
 rhoh(nlev) = 0.5*(rho(nlev) + zrhos) 
!
 do jk=1,nlev
   ta(jk) = tha(jk)*(pa(jk)/p00)**(Rd/Cp)
 enddo  
! 
 do jk=1,nlev
   zpa(jk) = pa(jk)
 enddo
 zpa(nlev+1) = ps
!   
 do jk=1,nlev
  if (jk /= nlev) then
    tah(jk)  = 0.5*(ta(jk) + ta(jk+1))
    qvah(jk) = 0.5*(qva(jk) + qva(jk+1))
  else
    tah(jk)  = 0.5*(ta(jk) + ts)
    qvah(jk) = 0.5*(qva(jk) + qvs)
  endif 
 enddo
! 
! Temperature at model top - set to zero with a full profile (not only PBL)
!
 zt_top = tah(1)
! 
! Path length for water vapour and carbon dioxide
!
! a) from a given level jk1 down to the surface
!
 uh2o_d(:) = 0.0
 uco2_d(:) = 0.0
 do jk1=1,nlev
   do jk2=jk1,nlev 
     zscale_h2o = (zpa(jk2)/zp0)**0.85*(zt0/tah(jk2))**0.5
     zscale_co2 = (zpa(jk2)/zp0)**0.75*(zt0/tah(jk2))**0.0
     uh2o_d(jk1) = uh2o_d(jk1) + 0.1/grav*qvah(jk2)*zscale_h2o*(zpa(jk2+1) - zpa(jk2))
     uco2_d(jk1) = uco2_d(jk1) + 0.00612*zscale_co2*(zpa(jk2+1) - zpa(jk2))
   enddo  
 enddo 
!
! b) from a given level jk1 up to model top  
! 
 uh2o_u(:) = 0.0
 uco2_u(:) = 0.0
 do jk1=nlev,1,-1
   do jk2=jk1,1,-1 
     zscale_h2o = (zpa(jk2)/zp0)**0.85*(zt0/tah(jk2))**0.5
     zscale_co2 = (zpa(jk2)/zp0)**0.75*(zt0/tah(jk2))**0.0
     uh2o_u(jk1) = uh2o_u(jk1) + 0.1/grav*qvah(jk2)*zscale_h2o*(zpa(jk2+1) - zpa(jk2))
     uco2_u(jk1) = uco2_u(jk1) + 0.00612*zscale_co2*(zpa(jk2+1) - zpa(jk2))
   enddo  
   uh2o_u(jk1) = max(1.E-4,uh2o_u(jk1))
 enddo  
!
! Effective emissivities and longwave radiative cooling rate (at half levels)
! 
 do jk=1,nlev   
   emis_up2 = emis3_h2o(uh2o_u(jk)) + emis1_co2(uco2_u(jk))
   if (jk .ne. 1) then    
     emis_up1 = emis3_h2o(uh2o_u(jk-1)) + emis1_co2(uco2_u(jk-1))
   else
     emis_up1 = 1.0*emis_up2
   endif 
! 
   emis_dn1 = emis3_h2o(uh2o_d(jk)) + emis1_co2(uco2_d(jk))
   if (jk .ne. nlev) then    
     emis_dn2 = emis3_h2o(uh2o_d(jk+1)) + emis1_co2(uco2_d(jk+1))
   else
     emis_dn2 = 1.0*emis_dn1
   endif      
!   
   !invcpdz = -Stefan/(rhoh(jk)*Cp*(z(jk+1) - z(jk)))
! Modified formulation to prevent from too low density values at high levels   
   invcpdz = grav*Stefan/(Cp*(zpa(jk+1) - zpa(jk)))
   dtdth(jk) = invcpdz*((tah(jk)**4 - ts**4)*(emis_dn2 - emis_dn1) + &
            &        (zt_top**4 - tah(jk)**4)*(emis_up2 - emis_up1))          
!
!  Cooling to space approximation
!             
   if (l_cts) then
     dtdth(jk) = invcpdz*tah(jk)**4*(emis_up2 - emis_up1)         
   endif
!
!  Simplified tendencies proposed by Pielke (1984)
!     
   if (l_simtend) then
     dtdth(jk) = -(0.017*(tah(jk)-273.15) + 1.8)/86400.0  
   endif
!
!  Empirical correction fo water vapour continuum (Savijarvi, 1990)
!   
   if (l_wvcont) then
     dtdth(jk) = dtdth(jk) - (1.E-3*(1.E3*qvah(jk))**3 + 0.1)/86400.0  
   endif
!    
 enddo
! 
 do jk=2,nlev
   dtdt(jk) = 0.5*(dtdth(jk) + dtdth(jk-1))
 enddo
 dtdt(1)  = 0.5*dtdth(1)
!
! write (*,*) '----------------------------------------------------------------'
! do jk=1,nlev
!  write(*,*) 'heating rate',jk,pa(jk)/100.,dtdt(jk)*86400.0,-(0.017*(ta(jk)-273.15) + 1.8)
! enddo
! stop  
 return
end subroutine lw_radiation
