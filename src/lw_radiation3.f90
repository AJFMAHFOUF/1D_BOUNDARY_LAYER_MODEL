subroutine lw_radiation3(nlev,ts,qvs,ps,tha,qva,pa,dtdt)
!-------------------------------------------------------------------------------------
!
! Compute simplified longwave radiative cooling (Sasamori, 1972)
!
!                                 Jean-Francois MAHFOUF (01/22)
!                                 modified by JFM (12/23)
!                                 modified by JFM (08/26) - tendencies at full levels  
! 
!--------------------------------------------------------------------------------------
 use const, only : Stefan, grav, Cp, Rd, p00
 implicit none    
 interface
  real function emis1_co2(u)
   implicit none
   real, intent(in)  :: u
  end function emis1_co2
  real function emis2_co2(u)
   implicit none
   real, intent(in)  :: u
  end function emis2_co2
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

 real, dimension(nlev)   :: uh2o_d, uco2_d, uh2o_u, uco2_u, ta     
 real, dimension(nlev+1) :: zpah      
 real                    :: emis_up1, emis_up2, emis_dn1, emis_dn2, invcpdp
 real                    :: zt_top, zp0, zt0, zscale_h2o, zscale_co2, zdelta_emis, zt_surf
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
! Empirical correction of emissivity at model top 
! 
 zdelta_emis = 0.0 ! 3.0E-5 for ICRCCM
!
! Mean atmospheric variables at half-levels
!
 do jk=1,nlev
   ta(jk) = tha(jk)*(pa(jk)/p00)**(Rd/Cp)
 enddo  
!   
! Define pressure at half levels - extrapolation at model top
!
 do jk=2,nlev
  zpah(jk) = exp(0.5*(log(pa(jk)) + log(pa(jk-1))))
 enddo
 zpah(1) = 2.0*pa(1) - zpah(2)
 zpah(nlev+1) = ps
! 
! Temperature at model top - set to zero with a full profile (not only PBL)
!
 zt_top = ta(1)
! 
! Effective surface temperature for radiation
! 
 zt_surf = (0.5*(ta(nlev)**4 + ts**4))**(0.25)
! 
! Path length for water vapour and carbon dioxide
!
! a) from a given level jk1 down to the surface
!
 uh2o_d(:) = 0.0
 uco2_d(:) = 0.0
 do jk1=1,nlev
   do jk2=jk1,nlev 
     zscale_h2o = (pa(jk2)/zp0)**1.20*(zt0/ta(jk2))**0.5
     zscale_co2 = (pa(jk2)/zp0)**0.75*(zt0/ta(jk2))**0.0
     uh2o_d(jk1) = uh2o_d(jk1) + 0.1/grav*qva(jk2)*zscale_h2o*(zpah(jk2+1) - zpah(jk2))
     uco2_d(jk1) = uco2_d(jk1) + 0.00612*zscale_co2*(zpah(jk2+1) - zpah(jk2))
   enddo  
 enddo 
!
! b) from a given level jk1 up to model top  
! 
 uh2o_u(:) = 0.0
 uco2_u(:) = 0.0
 do jk1=nlev,1,-1
   do jk2=jk1,1,-1 
     zscale_h2o = (pa(jk2)/zp0)**1.20*(zt0/ta(jk2))**0.5
     zscale_co2 = (pa(jk2)/zp0)**0.75*(zt0/ta(jk2))**0.0
     uh2o_u(jk1) = uh2o_u(jk1) + 0.1/grav*qva(jk2)*zscale_h2o*(zpah(jk2+1) - zpah(jk2))
     uco2_u(jk1) = uco2_u(jk1) + 0.00612*zscale_co2*(zpah(jk2+1) - zpah(jk2))
   enddo  
 enddo  
!
! Effective emissivities and longwave radiative cooling rate (at half levels)
! 
 do jk=1,nlev   
   emis_up2 = emis3_h2o(uh2o_u(jk)) + emis1_co2(uco2_u(jk))
   if (jk .ne. 1) then    
     emis_up1 = emis3_h2o(uh2o_u(jk-1)) + emis1_co2(uco2_u(jk-1))
   else
     emis_up1 = 1.0*emis_up2 - zdelta_emis
   endif 
! 
   emis_dn1 = emis3_h2o(uh2o_d(jk)) + emis1_co2(uco2_d(jk))
   if (jk .ne. nlev) then    
     emis_dn2 = emis3_h2o(uh2o_d(jk+1)) + emis1_co2(uco2_d(jk+1))
   else
     emis_dn2 = 1.0*emis_dn1
   endif      
!   
   invcpdp = grav*Stefan/(Cp*(zpah(jk+1) - zpah(jk)))
   dtdt(jk) = invcpdp*((ta(jk)**4 - zt_surf**4)*(emis_dn2 - emis_dn1) + &
            &      (zt_top**4 - ta(jk)**4)*(emis_up2 - emis_up1))          
!
!  Cooling to space approximation
!             
   if (l_cts) then
     dtdt(jk) = -invcpdp*ta(jk)**4*(emis_up2 - emis_up1)         
   endif
!
!  Simplified tendencies proposed by Pielke (1984)
!     
   if (l_simtend) then
     dtdt(jk) = -(0.017*(ta(jk)-273.15) + 1.8)/86400.0  
   endif
!
!  Empirical correction for water vapour continuum (Savijarvi, 1990)
!   
   if (l_wvcont) then
     dtdt(jk) = dtdt(jk) - (1.E-3*(1.E3*qva(jk))**3 + 0.1)/86400.0  
   endif
!    
 enddo
!
! write (*,*) '----------------------------------------------------------------'
! do jk=1,nlev
!  write(*,*) 'heating rate',jk,pa(jk)/100.,dtdt(jk)*86400.0,-(0.017*(ta(jk)-273.15) + 1.8)
! enddo
 return
end subroutine lw_radiation3
