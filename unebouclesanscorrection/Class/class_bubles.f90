module class_bubles
  use parameters
  use ode
  use class_system
  use class_indice
  implicit none
  type, public :: buble
     real(kind = wp),allocatable::IP(:,:),IC(:,:)
     real(kind = wp)::vec,t_perp_ini=0.0,t_perp2_ini=0.0
     real(kind = wp)::E0=0.0,Temperature=0.0
     integer::N_patche
   contains
     procedure :: create
     procedure :: destruct
     procedure::values
     procedure::values_tp20
     procedure::calcul  
     procedure::Eperp
  end type buble
  interface buble
     module procedure new_buble
  end interface buble
contains
  function new_buble(sys)
    type(system),intent(in)::sys 
    type(buble)::new_buble
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    new_buble%vec=2*pi/real(sys%in%N_patche)
    new_buble%N_patche=sys%in%N_patche
    new_buble%E0=sys%in%E0
    call new_buble%create()
  end function new_buble
  subroutine create(this)
    class(buble),intent(inout) ::this
    if(.not.allocated(this%IP))then
       allocate(this%IP(-this%N_patche/2:this%N_patche/2,-this%N_patche/2:this%N_patche/2))
       allocate(this%IC(-this%N_patche/2:this%N_patche/2,-this%N_patche/2:this%N_patche/2))
    end if
    this%IP=0.0_wp
    this%IC=0.0_wp
  end subroutine create
  subroutine destruct(this)
    class(buble),intent(inout) ::this
    deallocate(this%IP,this%IC)
  end subroutine destruct
  subroutine values(this,x,Temp,tpi,tp2i)
    class(buble),intent(inout) ::this
    real(kind = wp), intent(in) :: x,Temp,tpi,tp2i
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = wp) E
    integer :: kp,qp,k1,sigma,i_inf,i_sup
    character(6)::typ(2)
    
    this%t_perp_ini=tpi
    this%t_perp2_ini=tp2i
    this%Temperature=Temp
    typ(1)="C";typ(2)="P"
    E=this%E0
    this%E0=this%E0*exp(-x)
    this%IC=0.0_wp
    this%IP=0.0_wp
    
    i_inf=-this%N_patche/2
    i_sup=this%N_patche/2!-1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(qp,kp) 
!$OMP DO
    do qp=i_inf,i_sup
       do kp = i_inf,0
          this%IC(kp,qp)=this%calcul(typ(1),kp,qp)
          this%IP(kp,qp)=this%calcul(typ(2),kp,qp)
             
       end do
    end do
    
!$OMP END DO 
!$OMP END PARALLEL     
    this%E0=E
    this%IC(1:i_sup,i_inf:i_sup)=this%IC(-1:i_inf:-1,i_sup:i_inf:-1)
    this%IP(1:i_sup,i_inf:i_sup)=this%IP(-1:i_inf:-1,i_sup:i_inf:-1)

    !print*,x,32*this%IP(0,0,0),32*this%IC(0,0,0)

  end subroutine values
  subroutine values_tp20(this,x,Temp,tpi,tp2i)
    class(buble),intent(inout) ::this
    real(kind = wp), intent(in) :: x,Temp,tpi,tp2i
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = wp) E
    integer :: kp,qp,k1,sigma,i_inf,i_sup
    character(6)::typ(2)
    
    this%t_perp_ini=tpi
    this%t_perp2_ini=tp2i
    this%Temperature=Temp
    typ(1)="C";typ(2)="P"
    E=this%E0
    this%E0=this%E0*exp(-x)
    this%IC=0.0_wp
    this%IP=0.0_wp
    
    i_inf=-this%N_patche/2
    i_sup=this%N_patche/2!-1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(kp,qp) 
!$OMP DO
    do qp = i_inf,i_sup
       do kp=i_inf,0
          this%IC(kp,qp)=this%calcul(typ(1),kp,qp)
       end do
    end do
    
!$OMP END DO 
!$OMP END PARALLEL     
    
    this%E0=E
    this%IC(1:i_sup,i_inf:i_sup)=this%IC(-1:i_inf:-1,i_sup:i_inf:-1)
    
    this%IP(i_inf:i_sup,i_inf:0)=this%IC(i_sup:i_inf:-1,0:i_sup)
    this%IP(i_inf:i_sup,0:i_sup)=this%IC(i_sup:i_inf:-1,i_inf:0)
    !print*,x,32*this%IP(0,0,0),32*this%IC(0,0,0)

  end subroutine values_tp20
  function calcul(this,typ,kp,qp)
    class(buble),intent(inout) ::this
    character(*)::typ
    integer,intent(in) ::kp,qp
    real(kind = wp) :: calcul
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = wp) :: err=1E-7,abserr=1E-7
    integer,parameter::neq=1
    real ( kind = wp ) work(100+21*neq)
    integer :: iwork(5)
    integer:: iflag
    real (kind = wp) :: k_perp, q_perp
    real(kind=wp) :: k_perp_start, k_perp_end,y(neq),ytmp(neq)
    

    
    k_perp =real(kp); q_perp =real(qp)
    k_perp_start = k_perp-0.5_wp ; k_perp_end = k_perp+0.5_wp
    y=0.0_wp
    iflag = 1
100 call ode_rkf (DI, neq, y ,  k_perp_start ,k_perp_end, err,abserr, iflag, work, iwork )
    if(iflag==4)goto 100
    !call DI(k_perp,ytmp,y) 
    if (abs(iflag) == 2)then 
       calcul= y(1)/real(this%N_patche)
    else
       print*,"Error in integration : Class Buble"
       call error(iflag)
       stop
    end if
  contains
    subroutine DI(k_perp,y1,yp) 
      
      real ( kind = wp ),intent(out) :: yp(neq)
      real(kind = wp), intent(in) ::y1(neq), k_perp
      
      real(kind = wp) :: A,arg1,arg2,res
      integer :: i,mu
      A=this%Eperp(k_perp,q_perp,typ)
      
      mu =1
      yp = 0.0_wp
      arg1 = 0.5_wp*this%E0/this%Temperature
      do i=1,2
         arg2 = 0.5_wp*this%E0/this%Temperature + mu*0.5_wp*A/this%Temperature
         !call nickel_calcul(res,arg1,arg2)
         call nickelThanh(res,arg1,arg2)
         
         yp = yp + theta(abs(this%E0+mu*A)-this%E0)*res
         mu = -1
      end do
    end subroutine DI
  end function calcul
  function Eperp(this,k_perp,q_perp,typ)
    class(buble),intent(inout) ::this
    real(kind = wp),intent(in):: k_perp,q_perp
    character(*)::typ
    real ( kind = wp ) :: Eperp
    
    select case(typ)
    case('P')
       Eperp= 2*this%t_perp_ini*cos(k_perp*this%vec) + 2*this%t_perp2_ini*cos(2*k_perp*this%vec)&
            +2*this%t_perp_ini*cos( (q_perp+k_perp)*this%vec) +&
            2*this%t_perp2_ini*cos(2*(q_perp+k_perp)*this%vec)
              
    case('C')
       Eperp= 2*this%t_perp_ini*cos(k_perp*this%vec) + 2*this%t_perp2_ini*cos(2*k_perp*this%vec)&
            -2*this%t_perp_ini*cos( (q_perp-k_perp)*this%vec) -&
            2*this%t_perp2_ini*cos(2*(q_perp-k_perp)*this%vec)
    end select
  end function EPERP
end module class_bubles


