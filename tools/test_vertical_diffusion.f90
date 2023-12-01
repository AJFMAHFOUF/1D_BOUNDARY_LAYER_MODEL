program test
integer, parameter :: nlev = 10.
real, dimension (nlev) :: theta, t, qv, u, v, z, km, kh
real, dimension (nlev) :: thetam, um, vm, zm
data z /100.,200.,400.,500.,800.,1000.,1300.,1500.,1800.,2000./
data t /280.,278.,276.,274.,270.,268.,266.,262.,260.,258./
qv = 0.0
data u /2.0,3.0,5.0,7.0,10.0,11.0,9.0,7.0,9.0,12.0/
v = u
call init_surf1
do jk =1, nlev
 theta(jk) = t(jk) + 9.81/1000.0*z(jk)
enddo 
do jk = 1,nlev
  thetam(jk) = theta(nlev - jk + 1)
  um (jk) = u(nlev - jk + 1)
  zm(jk) = z(nlev - jk + 1)
enddo  
vm = um
call diffusion_coeff(nlev,thetam,qv,um,vm,zm,km,kh)
do jk = 1,nlev
 print *,'km kh',zm(jk),thetam(jk),um(jk),km(jk),kh(jk)
enddo 
stop
end program test
