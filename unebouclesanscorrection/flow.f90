program test
  use parameters
  use class_system
  use class_interaction
  use class_bubles
  use class_quasi1D
  implicit none

  character(20),parameter:: file="data.tpl"
  character(2)::ecrire
  type(quasi1D):: q
  real(kind=wp):: li=2,lf=50,Temp,x,dT,TempC
  integer::i,j
  open (unit=1,file="Tc.dat",access="append",status="unknown")
  q=quasi1D(file,ecrire)

  if(q%sys%in%find_TC)then
     call flow_Temp
  else
     call flow_one_Temp
  end if
  !call f1D
  call q%End
contains
  subroutine flow_Temp
    implicit none
    
    dT=1000.0
    Temp=q%Temperature
    do while(dT>=Temp)
       dT=dT/10.0
    end do
    
    do i=1,100
       li=0.0;lf=50.0
       
       if(.not.q%flow(li,lf))then
          call q%WriteOutputs()
       else
          
          !if(dT/q%Temperature<=1E-2)then
          if(abs(dT-0.01)<=1E-6)then
             write(1,'(f9.6,2x,e13.4)')q%t_perp2,q%Temperature
             exit
          else
             TempC=Temp
             Temp=Temp+dT
             dT=dT/10.0
          end if
       end if
       q=quasi1D(file)
       q%Temperature=Temp
       !if(dT<0.0001)exit
       if(abs(dT-q%Temperature)<=1E-6)dT=dT/10.0
       if(abs(q%Temperature-dT-TempC)<=1E-6)dT=dT/10.0
       q%Temperature=q%Temperature-dT
       Temp=q%Temperature
    end do
  end subroutine flow_Temp

  subroutine flow_one_Temp
    
    li=0.0;lf=50.0
    if(.not.q%flow(li,lf))then
       !print*,sum(q%interaction%g(:)%g2)!
       call q%WriteOutputs()
    end if
  end subroutine flow_one_Temp
  subroutine f1D
    integer::i
    
    li=0.0;lf=1
    do i=1,10
       q=quasi1D(file)
       if(.not.q%flow(li,lf))then
          print*,sum(q%interaction%g(:)%g2)!
          call q%WriteOutputs()
       end if
       li=0
       lf=lf+4
    end do
  end subroutine f1D
  
end program test
    
