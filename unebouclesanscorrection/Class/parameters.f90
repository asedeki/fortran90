module parameters
  !64 bits integer, parameter :: wp = selected_real_kind(33,4931)
  !integer, parameter :: wp = selected_real_kind(15,307)!32 bits
  integer, parameter :: wp = selected_real_kind(10,50)
  REAL(wp), PARAMETER :: PI=3.141592653589793238462643383279502884197_wp
  REAL(wp), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_wp
  !integer,parameter  :: N_patche=32
  
contains
  subroutine timestamp_file ( unite)
    implicit none
    integer,intent(in):: unite
    character ( len = 40 ) string

    call timestring ( string )
    write ( unite, '(a)' ) trim ( string )
    
    return
  end subroutine timestamp_file
  subroutine timestring ( string )
    implicit none

    character ( len = 8 ) ampm
    integer d
    character ( len = 8 ) date
    integer h
    integer m
    integer mm
    character ( len = 9 ), parameter, dimension(12) :: month = (/ &
         '##January  ', '##February ', '##March    ', '##April    ', &
         '##May      ', '##June     ', '##July     ', '##August   ', &
         '##September', '##October  ', '##November ', '##December ' /)
    integer n
    integer s
    character ( len = * ) string
    character ( len = 10 ) time
    integer values(8)
    integer y
    character ( len = 5 ) zone

    call date_and_time ( date, time, zone, values )
    
    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)
    
    if ( h < 12 ) then
       ampm = 'AM'
    else if ( h == 12 ) then
       if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
       else
          ampm = 'PM'
       end if
    else
       h = h - 12
       if ( h < 12 ) then
          ampm = 'PM'
       else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
             ampm = 'Midnight'
          else
             ampm = 'AM'
          end if
       end if
    end if
    
    write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
         trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )
    
    return
  end subroutine timestring
  integer function  wrap (i,N_patche)
    !integer :: wrap
    integer , intent(in) :: i
    integer , intent(in) :: N_patche 
    integer:: i_sup, i_inf
    i_sup = N_patche/2-1
    i_inf = -N_patche/2
    
    wrap=i
    
    if (wrap < i_inf )then
       wrap=wrap+N_patche
    elseif(wrap > i_sup )then
       wrap=wrap-N_patche
    end if
    
  end function wrap
  integer function  wrap2(i,N_patche)
    !integer :: wrap
    integer , intent(in) :: i
    integer , intent(in) :: N_patche 
    integer:: i_sup, i_inf,itmp
    i_sup = N_patche/2-1
    i_inf = -N_patche/2
    itmp=i
    do while(itmp < i_inf .or. itmp > i_sup)
       if (itmp < i_inf )itmp=itmp+N_patche
       if (itmp > i_sup )itmp=itmp-N_patche
    end do
    wrap2=itmp
  end function wrap2
  subroutine nickel_calcul(r,arg1,arg2)
    implicit none
    real(wp),intent(out) ::  r
    real(wp),intent(in)  ::  arg1,arg2
    real(wp) :: arg,c,D,ftanh
    D= 1.0_wp+arg2/arg1
    r = (tanh(arg1)+tanh(arg2))/D 
  end subroutine nickel_calcul
  subroutine nickelThanh(r,arg1,arg2)
    implicit none
    real(wp),intent(out) ::  r
    real(wp),intent(in)  ::  arg1,arg2
    real(wp) :: arg,c,D,ftanh
    D= 1.0_wp+arg2/arg1
    if(abs(D) >= 0.5 ) then 
       r = (tanh(arg1)+tanh(arg2))/D 
    else
       if (abs(arg1) > 24 ) then
          r = 0.0_wp
       else
          arg = arg1+arg2
          if (abs(arg) > 0.1_wp) then
             r = sinh(arg)/cosh(arg1)
             r = r/cosh(arg2)
             r = r/D
          else
             c = arg*arg
             r = 1.0_wp + 1.666667D-1*c + 8.333333D-3*c*c
             r = r/cosh(arg1)
             r = r*arg1
             r = r/cosh(arg2)
          end if
       end if
    end if
  end subroutine nickelThanh
  function theta(x)
    implicit none
    real(wp),intent(in)::x
    real(wp):: theta
    if (x > 0)then
       theta=1_wp
    elseif( x < 0) then
       theta=0.0_wp
    else
       theta=0.5_wp
    end if
  end function theta
end module parameters
