subroutine main
 use setup 
 use const, only : Rd, Cp, Lv, p00, rpi
 use surf1
 use soil
 implicit none
 interface
  real function qsat(p,t)
   implicit none
   real, intent(in)  :: p,t
  end function qsat
 end interface
 integer, parameter :: nlev=80        ! number of atmospheric vertical levels
 integer, parameter :: nlevs=14       ! number of soil vertical levels
 real, parameter    :: t_length = 48. ! duration of integration (hours)
 real, parameter    :: t_freq = 1.0   ! time frequency of geostrophic forcing (hours)
 integer, parameter :: ngeos = t_length/t_freq + 1 ! number of geostrophic time slots
 logical, parameter :: isba_fr = .false. ! choice of the soil scheme 
 real, dimension (nlev)         :: pa, tha, qva, ua, va
 real, dimension (nlev+1)       :: za, ug, vg
 real, dimension (nlev)         :: km, kh, rho
 real, dimension (nlev)         :: than, qvan, uan, van 
 real, dimension (nlev)         :: dtdt_lw, dudt_cor, dvdt_cor
 real, dimension (nlev+1,ngeos) :: ug1, vg1
 real, dimension (nlevs)        :: tsoil, tsoiln, wsoil, wsoiln, zsoil 
 real, dimension (nlevs)        :: K_w, D_w, rhocs, Lambda_s
 real :: ps, zs, ths, qvs, us, vs, rtime
 real :: fluxh, fluxq, fluxu, fluxv, ustar, zi, phi_m, phi_h
 real :: clf, zta, zqa, zua, zva, zps, pmu, rg, rl 
 real :: clay, sand
 real :: wg, w2, wr, ts, t2, tsk, delta 
 real :: wgn, w2n, wrn, tsn, t2n, tskn, pr, ro
 real :: rsoil, rs, ra, hu, ct, ustar2
 real :: h, lev, leg, letr, le, rn, g, EmP
 real :: z1, z2, qs
 real :: u10, v10, t2m, q2m, rh2m, zslope
 real :: swig, swi2, delta_ts, delta_t2
 real :: wg_in, w2_in, wsoil_m, Root_ext
 integer :: iloop, i, jk, ii, jj, jjp1
!
! Define soil vertical grid
!
 data zsoil/0.025,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.70,0.80,1.00/ 
!  
 namelist/surfprop/clay,sand,z0,z0h,emis,alpha,Rsmin,LAI,d2,veg
 namelist/simsetup/expid,lat,lon,Trans,f,dt,rtime0,iday0,lgeowind,llwrad
 namelist/soilinit/swig,swi2,delta_ts,delta_t2

 open (action='read',file='namelist',unit=8)
!
! Open input files
!
 open (unit=20,file='../data_in/Wangara_day33_0900_IC_L80.dat',status='old')
 open (unit=30,file='../data_in/Wangara_day33-34_UG_VG_L80.dat',status='old')
!
! Read initial atmospheric conditions
!
 do jk=1,nlev
   read (20,*) za(jk),pa(jk),tha(jk),qva(jk),ua(jk),va(jk)
 enddo
 read (20,*) za(nlev+1),ps,ths,qvs,us,vs  
! 
 if (ps == pa(nlev)) ps = pa(nlev) + 200.0
! 
 close (unit=20)
!
! Read hourly geostrophic wind profiles
! 
 do jj=1,49
  do jk=1,nlev+1
    read (30,*) ug1(jk,jj),vg1(jk,jj)
  enddo
 enddo   
! 
 close (unit=30)
! 
! Define initial surface and soil properties  
!
 call init_surf1
!
 clay = 0.38 ; sand = 0.50 
!
 read (unit=8,nml=surfprop)
!  
 call soil_def(clay,sand)
! 
! Initial soil conditions for the ISBA-2L scheme
! 
 swig = 0.5
 swi2 = 0.5
 delta_ts = 0.0
 delta_t2 = 0.0
! 
!  Read namelists to modify simulation set-up and initial conditions
!
 read (unit=8,nml=simsetup)
 read (unit=8,nml=soilinit)
! 
 wg = wwilt + swig*(wfc - wwilt)
 w2 = wwilt + swi2*(wfc - wwilt)
 ts = ths + delta_ts
 t2 = ths + delta_t2
 tsk = ts
 tsoil(:) = t2
 wsoil(:) = w2
 wsoil(1) = wg
 wr = 0.0 
! 
 iloop = int(t_length*3600.0/dt) ! number of time steps 
!
! Define surface relative humidity (bare soil) for first time step
! 
 if (isba_fr) then
   wg_in = wg
 else
   wg_in = wsoil(1)
 endif
 if (wg_in < wfc) then
   hu = 0.5*(1.0 - cos(wg_in/wfc*rpi))
 else
   hu = 1.0
 endif 
!
! Define initial depth of planetary boundary layer
!
 ustar = 0.25
 zi = 0.25*ustar/abs(f)
 zi = 100.0
! 
! Surface precipitation
!
 pr = 0.0
!
! Print experimental set-up
!
 write (*,113) '***********************************************************************'
 write (*,113) '*                                                                     *'
 write (*,113) '*   SINGLE COLUMN MODEL - PLANETARY BOUNDARY LAYER/SOIL/VEGETATION    *'
 write (*,113) '*   --------------------------------------------------------------    *'
 write (*,113) '*                                                                     *'
 write (*,114) '*                 EXPERIMENTAL SET-UP  : ',expid,'                       *'
 write (*,113) '*                 ==============================                      *'
 write (*,113) '*                                                                     *'
 write (*,100) '* LATITUDE =',lat,' LONGITUDE =',lon,' DAY OF YEAR =',iday0,'      *'
 write (*,101) '* INITIAL TIME (SEC) =',rtime0,' SIMULATION LENGTH (HOURS) =',t_length,'*' 
 write (*,102) '* MODEL TIME STEP (SEC) =',dt,' GEOSTROPHIC FORCING =',lgeowind,'       *'
 write (*,103) '* LONGWAVE RADIATION =',llwrad,'                                       *'
 write (*,104) '* CLAY FRACTION =',clay*100.,' SAND FRACTION =',sand*100.,'             *'
 write (*,105) '* WSAT =',wsat,' WFC = ',wfc,' WWILT =',wwilt,'                *'
 write (*,106) '* MOMENTUM ROUGHNESS LENGTH =',z0,' HEAT ROUGHNESS LENGTH =',z0h,' *'
 write (*,107) '* SURFACE EMISSIVITY =',emis,' SURFACE ALBEDO =',alpha,'        *'
 write (*,108) '* MINIMUM STOMATAL RESISTANCE =',Rsmin,' LEAF AREA INDEX =',LAI,'   *'
 write (*,109) '* VEGETATION FRACTION =',veg,' SOIL DEPTH (CM) =',d2*100.0,'            *'
 write (*,110) '* SUPERFICIAL SWI =',swig*100.0,' DEEP SWI =',swi2*100.,'                    *'
 write (*,111) '* SUPERFICIAL MOISTURE CONTENT =',wg,' DEEP MOISTURE CONTENT =',w2,'  *'
 write (*,112) '* SURFACE TEMPERATURE =',ts,' DEEP SOIL TEMPERATURE =',t2,'   *'
 write (*,113) '*                                                                     *'
 write (*,113) '***********************************************************************'
 100 format(1X,A12,1X,F7.2,A12,1X,F7.2,A14,1X,I3,6X,A7)
 101 format(1X,A22,1X,F6.0,A28,1X,F4.0,6X,A3)
 102 format(1X,A25,1X,F5.0,A22,1X,L3,6X,A8)
 103 format(1X,A22,L3,6X,A40)
 104 format(1X,A17,1X,F4.0,A16,1X,F4.0,6X,A22)
 105 format(1X,A8,1X,F5.3,A7,1X,F5.3,A8,1X,F5.3,6X,A24)
 106 format(1X,A29,1X,F4.2,A24,1X,F5.3,5X,A2)
 107 format(1X,A22,1X,F3.2,A17,1X,F3.2,6X,A18)
 108 format(1X,A31,1X,F5.1,A18,1X,F4.1,6X,A5)
 109 format(1X,A23,1X,F4.2,A18,1X,F5.0,6X,A13)
 110 format(1X,A19,1X,F4.0,A11,1X,F4.0,6X,A25)
 111 format(1X,A32,1X,F5.3,A24,1X,F5.3,A3)
 112 format(1X,A23,1X,F6.2,A24,1X,F6.2,6X,A4)
 113 format(A72)
 114 format(1X,A41,A6,A24)
!
!  Open output files
!
 open (unit=100,file='../data_out/atmospheric_profiles_exp_'//expid//'.dat') 
 open (unit=200,file='../data_out/soil_profiles_exp_'//expid//'.dat') 
 open (unit=300,file='../data_out/surface_fluxes_and_variables_exp_'//expid//'.dat')
 open (unit=400,file='../data_out/radiative_heating_rate_exp_'//expid//'.dat')  
!
 ii = 0 
! ----------------------------------------------------------
! Start integration of column model
!-----------------------------------------------------------
 do i = 0, iloop
   rtime = rtime0 + dt*i
!
! Time interpolation of geostrophic wind
!  
   jj = int(i*dt/(t_freq*3600.)) + 1
   zslope = mod(i*dt,t_freq*3600.)/(t_freq*3600.0)
   jjp1 = min(jj+1,ngeos)
   do jk=1,nlev+1
     ug(jk) = ug1(jk,jj) + zslope*(ug1(jk,jjp1) - ug1(jk,jj))
     vg(jk) = vg1(jk,jj) + zslope*(vg1(jk,jjp1) - vg1(jk,jj))
   enddo    
!
! Compute tendencies from Coriolis and gradient pressure forces
!   
  call coriolis(nlev,ua,va,ug,vg,dudt_cor,dvdt_cor)
  
  if (lgeowind) then
    ua(:) = ua(:) + dt*dudt_cor(:)
    va(:) = va(:) + dt*dvdt_cor(:)
  endif
!   
! Define fraction of wet vegetation   
!   
   if (wr > 0.0 .and. Wrmax /= 0.0) then
     delta = (wr/Wrmax)**(2./3.)
   else
     delta = 0.0
   endif
!
!  Compute soil properties
!   
   call soil_prop(wg,w2,ts)
!   
   if (isba_fr) then
     wg_in = wg
   else
     wg_in = wsoil(1)
   endif 
!      
   call rs_soil(wg_in,veg,rsoil)  
!
!  Multi-layer soil thermal and hydraulic properties
!   
   call soil_prop_diff(nlevs,wsoil,K_w,D_w,rhocs,Lambda_s)
!
!  Soil/vegetation thermal inertia
!    
   ct = 1./(veg/Cv + (1. - veg)/Cg)
!
!  Compute solar zenith angle 
!        
   call solar_angle(rtime,pmu)
!   
   clf = 0.0 ! cloud fraction set to zero
   zta = tha(nlev)*(pa(nlev)/p00)**(Rd/Cp)
   zqa = qva(nlev)
   zua = ua(nlev)
   zva = va(nlev)
   zps = ps   
   zs  = za(nlev+1)
!
!  Air density profile
!   
   rho(:) = pa(:)/(Rd*tha(:)*(pa(:)/p00)**(Rd/Cp)*(1.0 + 0.608*qva(:)))
! 
!  Downward longwave and shorwave radiation fluxes
!   
   call downward_radiation(zta,zqa,zps,clf,pmu,rg,rl)
!
!  Radiative cooling temperature tendency
!   
   call lw_radiation3(nlev,tsk,qvs,zps,tha,qva,pa,rho,za,dtdt_lw)   
!   
   if (llwrad) then
     tha(:) = tha(:) + dt*(p00/pa(:))**(Rd/Cp)*dtdt_lw(:)
   endif  
!   
!  Canopy resistance and aerodynamic resistances 
!
   wsoil_m = 0.5*wsoil(1)*zsoil(1)
   do jk=2,nlevs-1
     wsoil_m = wsoil_m + 0.5*wsoil(jk)*(zsoil(jk+1) - zsoil(jk-1))
   enddo
   wsoil_m = wsoil_m + wsoil(nlevs)*(zsoil(nlevs) - zsoil(nlevs-1))
   wsoil_m = wsoil_m*(1.5*zsoil(nlevs) - 0.5*zsoil(nlevs-1))
   if (isba_fr) then
     w2_in = w2
   else
     w2_in = wsoil_m 
   endif  
   call rs_veg(w2_in,zps,zqa,zta,rg,rs)    
   call drag_coeff_z0h(za(nlev),tha(nlev),tsk,zqa,qvs,zua,zva,ra,ustar)
!
!  Surface fluxes (to be provided by land surface scheme) 
!   
   qs = qsat(zps,tsk)
   call fluxes(rho(nlev),tsk,tsoil(1),tha(nlev),qs,zqa,rg,rl,ra,rs,rsoil,delta,hu,ct, &
             & h,lev,leg,letr,le,rn,g)    
!
!   fluxh = (0.4*rg)/Cp 
!   fluxq = (0.1*rg)/Lv
   fluxh = h/Cp
   fluxq = le/Lv
   fluxu = -ustar*ustar*zua/sqrt(zua*zua + zva*zva)
   fluxv = -ustar*ustar*zva/sqrt(zua*zua + zva*zva)
!   
!  Compute exchange coefficients differently according to surface stability
!
   if (h > 0.0) then 
!
!  Height of the PBL and stability functions in the constant flux layer 
!
     call pbl_height(ustar,fluxh,tsk,qvs,tha(nlev),za(nlev),zi,phi_m,phi_h)   
!
!  Eddy diffusivity exchange coefficients (0'Brien 1970)
!   
     call diffusion_coeff_obrien(nlev,zi,ustar,phi_m,phi_h,tha,ua,va,za,km,kh)
!
   else
!   
!  Eddy diffusivity exchange ceofficients (Louis et al., 1981) 
!      
      call diffusion_coeff_louis(nlev,tha,qva,ua,va,za,km,kh)
!      
      zi = 20. ! set PBL height at arbitrary small value
!      
   endif         
!   
!  Implicit vertical diffusion equation (tridiagonal solver)
!
   call vertical_diffusion(nlev,tha,za,rho,kh,fluxh,than)
   call vertical_diffusion(nlev,qva,za,rho,kh,fluxq,qvan)
   call vertical_diffusion(nlev,ua,za,rho,km,fluxu,uan)
   call vertical_diffusion(nlev,va,za,rho,km,fluxv,van)
!
!  Solve surface energy budget - evolution of soil temperatures 
!   
   if (isba_fr) then
     wg_in = wg
   else
     wg_in = wsoil(1)
   endif 
   call energy_budget(rho(nlev),ts,t2,tsk,tsoil(1),ra,rs,rsoil,hu,delta,ct, &
                    & zps,tha(nlev),zqa,rg,rl,wg_in,tsn,t2n,tskn)
!
!  Solve multi-layer soil heat transfers (diffusion equation)      
!
   call soilt_vertical_diffusion(nlevs,tsoil,zsoil,rhocs,Lambda_s,tskn,tsoiln)                     
!
!  Solve surface water budget - evolution of soil moisture contents + interception reservoir  
!   
   call water_budget(wg,w2,wr,pr,leg,lev,letr,wgn,w2n,wrn,ro,EmP) 
!
!  Transpiration flux for soil root extraction (expressed in m/s)
! 
   Root_ext = 1.E-3*letr/Lv 
!
!  Solve multi-layer soil water transfers (diffusion equation)   
!   
   call soilw_vertical_diffusion(nlevs,wsoil,zsoil,K_w,D_w,EmP,Root_ext,wsoiln)
!
!  Write results in file
!
   if (amod(dt*i,3.*3600.) == 0.0) then 
     do jk=1,nlev
      write (100,*) za(jk),tha(jk),qva(jk),ua(jk),va(jk),km(jk),kh(jk)
     enddo
     do jk=1,nlevs
      write (200,*) zsoil(jk),tsoil(jk),wsoil(jk),Lambda_s(jk),D_w(jk)
     enddo
   ii = ii + 1   
   endif  
!
!  Swapp over time steps
!    
   tha(:) = than(:)
   qva(:) = qvan(:)
   ua(:)  = uan(:)
   va(:)  = van(:)
!   
   ts  = tsn
   tsk = tskn
   t2  = t2n
   wg  = wgn
   w2  = w2n
   wr  = wrn
   tsoil(:) = tsoiln(:)
   wsoil(:) = wsoiln(:)
!
!  Surface specific humidity
!
  z1 = ra/(ra + ra) ; z2 = ra/(ra + rsoil)
  qs = qsat(zps,tsk)
  qvs = qva(nlev) + veg*delta*(qs - qva(nlev)) + veg*(1. - delta)*(qs -qva(nlev))*z1 + &
      & (1. -veg)*(hu*qs - qva(nlev))*z2   
!
!  Vertical interpolation to get screen-level parameters: wind at 10m and T, q, RH at 2m
!   
   call vdfppcfls(za(nlev),zps,zta,tsk,ua(nlev),va(nlev),qva(nlev),qvs,u10,v10,t2m,q2m,rh2m)
!
!  Write outputs in files for graphics
!   
   write (300,200) i,rg,rl,rn,h,le,leg,lev,letr,g,tsk,tsoil(1),wg,w2,t2m,wsoil(1)
   200 format(1X,I5,15(1X,F10.6))
!   
   do jk=1,nlev
     write (400,*) za(jk),dtdt_lw(jk)*86400.0
   enddo  
!   
!   write (*,*) ts,qvs*1000.0,t2m,q2m,u10,ua(nlev)
!   
 enddo
!
! Close output files
!
 close (unit=100) 
 close (unit=200)
 close (unit=300)
 close (unit=400)
! 
 return
end subroutine main
