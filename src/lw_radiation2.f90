subroutine lw_radiation2(nlev,ts,qvs,ps,tha,qva,pa,rho,z,dtdt)
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
  real function emis2_h2o(u)
   implicit none
   real, intent(in)  :: u
  end function emis2_h2o
  real function emis1_h2o(u)
   implicit none
   real, intent(in)  :: u
  end function emis1_h2o
 end interface     
 integer,                 intent(in)  :: nlev
 real,                    intent(in)  :: ts, qvs, ps
 real, dimension(nlev),   intent(in)  :: tha, qva, pa, rho
 real, dimension(nlev+1), intent(in)  :: z
 real, dimension(nlev),   intent(out) :: dtdt

 real, dimension(nlev)      :: rhoh, ta, tah, qvah, dtdth, flux_u, flux_d
 real, dimension(nlev,nlev) :: uh2o, uco2      
 real, dimension(nlev+1)    :: zpa      
 real                       :: zrhos, emis_1, emis_2, invcpdz 
 integer                    :: jk, jk1, jk2
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
! Path length for water vapour and carbon dioxide (between two levels)
!
 uh2o(:,:) = 0.0
 uco2(:,:) = 0.0
 do jk1=1,nlev
   do jk2=jk1+1,nlev 
     uh2o(jk1,jk2) = uh2o(jk1,jk2-1) + 0.1/grav*qvah(jk2)*(zpa(jk2+1) - zpa(jk2))
     uco2(jk1,jk2) = uco2(jk1,jk2-1) + 0.00612*(zpa(jk2+1) - zpa(jk2))
   enddo  
 enddo 
 do jk1=1,nlev
   do jk2=1,jk1
     uh2o(jk1,jk2) = uh2o(jk2,jk1) 
     uco2(jk1,jk2) = uco2(jk2,jk1) 
   enddo  
 enddo 
!
! Downward longwave flux
! 
 do jk1=1,nlev
   flux_d(jk1) = 0.0
   do jk2=2,jk1
     emis_1 = emis1_co2(uco2(jk1,jk2))   + emis1_h2o(uh2o(jk1,jk2))
     emis_2 = emis1_co2(uco2(jk1,jk2-1)) + emis1_h2o(uh2o(jk1,jk2-1))
     flux_d(jk1) = flux_d(jk1) + Stefan*tah(jk2-1)**4*(emis_2 - emis_1) 
   enddo
   flux_d(jk1) = flux_d(jk1) + Stefan*ta(1)**4*(1.0 - emis1_co2(uco2(jk1,1)) - emis1_h2o(uh2o(jk1,1)))
 enddo  
!
! Updward longwave flux
! 
 do jk1=nlev,2,-1
   flux_u(jk1) = 0.0
   do jk2=nlev,jk1,-1
     emis_1 = emis1_co2(uco2(jk1,jk2))   + emis1_h2o(uh2o(jk1,jk2))
     emis_2 = emis1_co2(uco2(jk1,jk2-1)) + emis1_h2o(uh2o(jk1,jk2-1))
     flux_u(jk1) = flux_u(jk1) + Stefan*tah(jk2-1)**4*(emis_1 - emis_2) 
   enddo
   flux_u(jk1) = flux_u(jk1) + Stefan*ts**4*(1.0 - emis1_co2(uco2(jk1,nlev)) - emis1_h2o(uh2o(jk1,nlev)))
 enddo    
! 
!  Radiative cooling at each level
!  
 do jk=1,nlev-1
   invcpdz = 1.0/(rhoh(jk)*Cp*(z(jk+1) - z(jk)))
   dtdth(jk) = -invcpdz*(flux_u(jk+1) - flux_u(jk) + flux_d(jk) - flux_d(jk+1))    
 enddo
! 
 do jk=2,nlev-1
   dtdt(jk) = 0.5*(dtdth(jk) + dtdth(jk-1))
 enddo
 dtdt(1)  = dtdth(2)
 dtdt(nlev) = dtdth(nlev-1)
!
 !write (*,*) '----------------------------------------------------------------'
 !do jk1=1,nlev
 !  do jk2=1,nlev
    ! write(*,*) jk1,jk2,uh2o(jk1,jk2),uco2(jk1,jk2)
  !enddo  
 !enddo  
 stop
 return
end subroutine lw_radiation2
