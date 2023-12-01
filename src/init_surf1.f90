subroutine init_surf1
!---------------------------------------------------------------
!
! Definition of constant surface properties 
! Addition of maximum storage of interception reservoir WRMAX
!
!                                Jean-Francois MAHFOUF (11/06)
!                                                      (10/21)
!---------------------------------------------------------------
 use surf1
 implicit none
 z0 = 0.10           ! surface roughness length
 z0h = 0.01          ! surface roughness length for heat
 emis = 0.97         ! surface emissivity
 alpha = 0.2         ! surface albedo
 Rsmin = 40.         ! minimum stomatal resistance
 LAI = 1.            ! leaf area index
 d1 = 0.01           ! depth of surface soil layer
 d2 = 1.00           ! depth of deep soil layer
 gamma = 20.         ! dependency of RS with saturation vapor deficit
 Rgl = 100.          ! dependency of RS with solar radiation
 Cv = 2.E-5          ! vegetation thermal coefficient 
 veg = 0.85          ! vegetation fractionnal cover
 Wrmax = 0.2*veg*LAI ! maximum capacity of interception reservoir
 Lambda_sk = 7.0     ! Thermal properties of skin layer
 return
end subroutine init_surf1 
