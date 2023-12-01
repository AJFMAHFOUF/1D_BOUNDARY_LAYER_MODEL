subroutine solar_angle(rtime,pmu)
!-----------------------------------------------------------------------
!   Computation of solar zenith angle (simplified formula)
!   
!
!                                    Jean-François Mahfouf (04/01/2022)
!-----------------------------------------------------------------------
 use const, only : rpi, tau
 use setup, only : lat, lon, rtime0, iday0
 implicit none
 real,    intent(in)  :: rtime      ! time in seconds since beginning of forecast
 real,    intent(out) :: pmu        ! cosinus of zenith angle
 real :: rdecli, rwsovr, iday, rtime2, rsidec, rcodec, pgemu, pgelam, rcovsr, rsivsr
!
 rtime2 = mod(rtime + rtime0,tau)
 iday   = iday0 + int((rtime + rtime0)/tau)
 rdecli = 0.409*cos(2.0*rpi*(iday - 172.0)/365.25) ! solar declinaison 
 rwsovr = 2.0*rpi*rtime2/tau
 rsidec = sin(rdecli)
 rcodec = cos(rdecli)
 pgemu  = sin(lat*rpi/180.)
 pgelam = lon*rpi/180.
 rcovsr = cos(rwsovr)
 rsivsr = sin(rwsovr)
! 
 pmu = max(rsidec*pgemu &
     &   - rcodec*rcovsr*sqrt(1.- pgemu**2)*cos(pgelam) &
     &   + rcodec*rsivsr*sqrt(1.- pgemu**2)*sin(pgelam), 0.)
 if (pmu > 0.0) then 
   pmu = sqrt(1224.*pmu*pmu + 1.)/35. ! Magnification factor
 endif  
 return
end subroutine solar_angle
