!----------------------------------------
 function emis1_co2(u)
!****************************************
! CO2 flux emissivity (Kondratyev, 1969)
!
! input : u in cm
! ouput : emis1_co2 
!
! Jean-Francois Mahfouf (15/01/2022)
!
!****************************************
 implicit none
 real :: emis1_co2 
 real, intent(in) :: u
 if (u > 0.0) then
  emis1_co2 = 0.185*(1.0 - exp(-0.3919*u**(0.4)))
 else
  emis1_co2 = 0.0
 endif
 end function emis1_co2
!---------------------------------------
 function emis1_h2o(u)
!****************************************
! H2O flux emissivity (Staley & Jurica, 1972, fit @ 20 °C)
!
! input : u in cm 
! ouput : emis1_h2o
!
! Jean-Francois Mahfouf (15/01/2022)
!
!*****************************************
 implicit none
 real :: emis1_h2o 
 real, intent(in) :: u   
 real :: logu
 if (u > 0.0) then
   logu = log10(u)                          
   emis1_h2o = 0.61746 + 0.21334*logu + 0.018539*logu*logu
 else
   emis1_h2o = 0.0
 endif
 end function emis1_h2o
!------------------------------------------------------------
 function emis2_h2o(u)
!****************************************
! H2O flux emissivity (Kuhn, 1963 ; Jacobs et al, 1974)
!
! input : u in cm 
! ouput : emis1_h2o
!
! Jean-Francois Mahfouf (15/01/2022)
!
!*****************************************
 implicit none
 real :: emis2_h2o 
 real, intent(in) :: u   
 real :: logu
 if (u > 0.0) then
  logu = log10(u) 
  if (logu < -4.0) then
    emis2_h2o = 0.11288*log10(1 + 12.63*u)
  else if (logu < -3.0) then
    emis2_h2o = 0.104*logu + 0.440
  else if (logu < -1.5) then
    emis2_h2o = 0.121*logu + 0.491
  else if (logu < -1.0) then
    emis2_h2o = 0.146*logu + 0.527
  else if (logu <  0.0) then
    emis2_h2o = 0.161*logu + 0.542
  else
    emis2_h2o = 0.136*logu + 0.542
  endif
 else
  emis2_h2o = 0.0
 endif  
 end function emis2_h2o

