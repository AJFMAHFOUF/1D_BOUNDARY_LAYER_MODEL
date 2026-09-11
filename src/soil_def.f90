subroutine soil_def(clay,sand)    
!--------------------------------------------------------
!
! Initialisation of soil thermal and hydraulic properties
! (part 1)
!
!--------------------------------------------------------
 use soil
 implicit none
 real, intent(in) :: clay, sand 
 real             :: lquartz, l0, rhod
!
 lquartz = 7.7 ! Quartz thermal conductivity
 lwater = 0.57 ! Water thermal conductivity
! 
 wsat = (-108.*sand + 494.305)*1.E-3
 wwilt = 37.1342E-3*sqrt(clay*100.)
 wfc = 89.0467E-3*(clay*100.)**(0.3496)
 b = 13.7*clay + 3.501
 Cgsat = (-1.557*sand - 1.441*clay + 4.7021)*1.0E-6 
 C1sat = 5.58*clay + 0.8488
 C2ref = 13.815*(clay*100.)**(-0.954) 
 C3 = 5.327*(clay*100.)**(-1.043)
 a = 0.73242*(clay*100.)**(-0.539)
 p = 13.4*clay + 3.4
 wl = 1.E-5
! 
! Kwsat = 1.E-6*10**(-8.967 + 6.83*(clay*100)**0.2)  ! formula of unknown origin
!
! Fit of Clapp - Hornberger (1978) values - bounded by upper clay fraction   
 Kwsat = 10**(8.95*min(0.7,clay)**2 - 9.17*min(0.7,clay) - 3.6)
! Kwsat = 10**(2.425*sand - 6.171)
 Psis = -10**(-0.883*sand - 0.143)  ! Decharme et al. (2011) very bad correlation coefficient (modified)
! 
 rhoc_soil = (1.398 - 0.532*clay)*1.E6
 if (sand < 0.2) then
   l0 = 3.0
 else
   l0 = 2.0
 endif  
 lsoil = lquartz**sand*l0**(1.0 - sand)
 rhod = (1.0 - wsat)*2700.0
 Lambda_dry = (0.135*rhod + 64.7)/(2700.0 - 0.947*rhod)  
 
! print *,'*** Soil properties ***'
! print *,'Kwsat=',Kwsat,' Psis=',Psis,' b=',b
! print *,'Lambda_dry=',Lambda_dry,' wsat=',wsat,' wfc=',wfc,' wwilt=',wwilt
!
 return
end subroutine soil_def
