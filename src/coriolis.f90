subroutine coriolis(nlev,ua,va,ug,vg,dudt,dvdt)
!-------------------------------------------------------------------------
!
! Compute wind tendencies from Coriolis and pressure gradient forces
!
!                                 Jean-Francois MAHFOUF (01/22)
! 
!-------------------------------------------------------------------------
 use setup, only : dt, f
 implicit none       
 integer,                   intent(in)  :: nlev
 real, dimension(nlev),     intent(in)  :: ua, va
 real, dimension(nlev+1),   intent(in)  :: ug, vg
 real, dimension(nlev),     intent(out) :: dudt, dvdt

 real                    :: zdtf, zdtf2, zopdtf2
 integer                 :: jk
 logical, parameter      :: noimp=.false.
!
 zdtf = dt*f
 zdtf2 = zdtf*zdtf
!
! Control implicit solution of equations 
!  
 if (noimp) zdtf2 = 0.0 
 zopdtf2 = 1.0 + zdtf2
! 
 do jk=1,nlev
  dudt(jk) =  (zdtf*(va(jk) - vg(jk)) - zdtf2*(ua(jk) - ug(jk)))/(zopdtf2*dt)
  dvdt(jk) = -(zdtf*(ua(jk) - ug(jk)) + zdtf2*(va(jk) - vg(jk)))/(zopdtf2*dt)
 enddo
!
 return
end subroutine coriolis
