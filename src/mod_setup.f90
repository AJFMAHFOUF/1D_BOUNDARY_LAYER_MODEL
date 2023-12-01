module setup
character (6) :: expid        
real :: lat = -34.5        , &
&       lon = 144.93       , &
&       Trans =  0.9       , &
&       f  =  -0.826E-4    , &
&       dt = 600.0          , &
&       rtime0 = 23.*3600.    
integer :: iday0 = 227 
logical :: lgeowind = .true., &
           llwrad =.true.
end module setup
