program vertical_grid
 real, parameter :: a1=0.0376*1.044, a2=0.417/3.0, a3=0.01*2.
 real, dimension (0:80) :: Z, P, T, TH, TH2, Q, U, V
 real, dimension (31) :: z1, p1, t1, th1, q1, u1, v1
 real, parameter :: z0 = 0.01, g=9.81, Rd=287.0, Cp=1005.0, p00=1.0E5  
 real, dimension(18) :: ug0, vg0
 real, dimension(6)  :: dug1,dug2,dvg1,dvg2 
 real, dimension(0:48) :: zug0, zvg0, zdug1, zdvg1, zdug2, zdvg2
 real, dimension(0:80,0:48) :: ug, vg
! data ug0/-5.5,-5.7,-6.2,-6.3,-6.0,-6.8,-8.0,-7.2,-7.7,-8.0,-9.0,-8.8,-8.0,-8.4,-10.0,-8.7,-9.2,0.0/
 data ug0/-5.5,-5.7,-6.2,-6.3,-6.0,-6.8,-8.0,-7.2,-7.7,-8.0,-9.0,-8.6,-8.0,-8.4,-10.0,-8.7,-9.2,0.0/
! data vg0/-0.5,1.2,0.9,-0.7,-2.0,-3.0,-3.5,-5.0,-4.7,-3.7,-4.5,-6.1,-5.9,-8.4,-10.0,-9.9,-9.2,0.0/ 
 data vg0/-0.5,1.2,0.9,-0.7,-2.0,-3.0,-3.5,-5.0,-4.7,-3.7,-4.5,-6.1,-5.9,-8.4,-10.0,-9.6,-9.2,0.0/ 
 data dug1/2.98,2.81,1.87,1.91,1.70,0.0/
 data dvg1/-0.04,-0.67,0.59,1.30,1.21,0.0/
 data dug2/1.49,1.32,1.04,-0.52,2.91,0.0/
 data dvg2/0.26,0.45,0.98,1.19,0.04,0.0/
 open (unit=10,file='../data_in/wangara_day33_0900_uvs.dat',status='old')
 open (unit=25,file='../data_in/Wangara_day33_0900_IC_L80.dat')
 open (unit=30,file='../data_in/Wangara_day33-34_UG_VG_L80.dat')
!
! Vertical grid according to Yamada and Mellor (1975) with 
! modified constants (a1,a2,a3) 
!
! first model level at 10 m
!
 Z(1) = 10.0
!
 do k=2,80
 zn = Z(k-1)
   do 
     zn1 = zn - (a1*zn + a2*log(zn/a3) - float(k))/(a1 + a2/zn)
     if (abs(zn1 - zn) < 1.E-3) then
       exit
     else
       zn = zn1
     endif
   enddo
 Z(k) = zn
 enddo
!
! read profiles on the observations grid
!
 do j=1,31
   read (10,*) z1(j), p1(j), t1(j), th1(j), q1(j), u1(j), v1(j)
 enddo 
!
! Perform linear interpolation on model grid
!
 do k=1,80
   do j=1,30
     if (z1(j+1) > Z(k) .and. z1(j) <= Z(k)) then
        !print*,'altitudes',Z(k),z1(j+1),z1(j)
        alpha = (Z(k) - z1(j))/(z1(j+1) -z1(j)) 
        T(k)  = t1(j)  + alpha*(t1(j+1) - t1(j))
        TH(k) = th1(j) + alpha*(th1(j+1) - th1(j))
        Q(k)  = q1(j)  + alpha*(q1(j+1) - q1(j))
        U(k)  = u1(j)  + alpha*(u1(j+1) - u1(j))
        V(k)  = v1(j)  + alpha*(v1(j+1) - v1(j))
     endif  
   enddo
!
! Modified vertical interpolation in surface boundary layer
!
   if (Z(k) < z1(2)) then
     print *,Z(k),z1(2)
     U(k)  = u1(2)*log(Z(k)/z0)/log(z1(2)/z0)
     V(k)  = v1(2)*log(Z(k)/z0)/log(z1(2)/z0)
     TH(k) = th1(1) + (th1(2)-th1(1))*log(Z(k)/z0)/log(z1(2)/z0)
     Q(k)  = q1(1) + (q1(2) - q1(1))*log(Z(k)/z0)/log(z1(2)/z0)
   endif 
 enddo
!
! define pressure at full levels for altitude grid
!
 P(1) = p1(1)
 do k=2,80
   Tvm = 0.5*(T(k)*(1.0 + 0.608*Q(k)) + T(k-1)*(1.0 + 0.608*Q(k-1))) 
   P(k)=P(k-1)*exp(-g/(Rd*Tvm)*(Z(k) - Z(k-1)))
 enddo
! 
 Z(0) = 0
 P(0) = p1(1)
 T(0) = t1(1)
 TH(0) = t1(1)
 TH2(0) = t1(1)*(p00/p1(1))**(Rd/Cp)
 Q(0) = q1(1)
 U(0) = 0.0
 V(0) = 0.0
!
! define potential temperature
!
 do k=1,80
   TH2(k) = T(k)*(p00/P(k))**(Rd/Cp)
 enddo
!
!  Write initial profiles on model grid
!
! do k=0,80
!   print *,k,Z(k),P(k),T(k),U(k),V(k),TH(k),TH2(k),Q(k)
!   kk = 80 - k
!   write (25,*) Z(kk),P(kk),TH2(kk),Q(kk),U(kk),V(kk)
! enddo
!
! Temporal interpolation of geostrophic wind (hourly values)
!
 do i=0,48
  j = int(i/3) + 1
  im = 3*(j-1) 
  ip = 3*j
  jj = int(i/12) + 1
  iim = 12*(jj-1) 
  iip = 12*jj
  zvg0(i)  = vg0(j) + float(i-im)*(vg0(j+1) - vg0(j))/float(ip-im)
  zug0(i)  = ug0(j) + float(i-im)*(ug0(j+1) - ug0(j))/float(ip-im)
  zdug1(i) = dug1(jj) + float(i - iim)*(dug1(jj+1) - dug1(jj))/float(iip - iim)
  zdug2(i) = dug2(jj) + float(i - iim)*(dug2(jj+1) - dug2(jj))/float(iip - iim)
  zdvg1(i) = dvg1(jj) + float(i - iim)*(dvg1(jj+1) - dvg1(jj))/float(iip - iim)
  zdvg2(i) = dvg2(jj) + float(i - iim)*(dvg2(jj+1) - dvg2(jj))/float(iip - iim)
 enddo 
!
! Vertical interpolation (parabolic profiles)
! 
 do i=0,48
  do k=80,0,-1
    ug(k,i)= (zdug2(i) -  zdug1(i))/2.0E6*Z(k)*Z(k) + (3.*zdug1(i) - zdug2(i))/2.0E3*Z(k) + zug0(i)
    vg(k,i)= (zdvg2(i) -  zdvg1(i))/2.0E6*Z(k)*Z(k) + (3.*zdvg1(i) - zdvg2(i))/2.0E3*Z(k) + zvg0(i)
    write (30,*) ug(k,i),vg(k,i)
  enddo
 enddo
! do k=80,0,-1
!   write(*,fmt=100) (ug(k,i),i=25,48)
! enddo
! 100 format(25(1X,F6.2))  
!  
 stop
end program vertical_grid
