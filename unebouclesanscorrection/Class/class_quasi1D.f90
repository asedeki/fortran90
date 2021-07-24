module class_quasi1D
  use parameters
  use class_system
  use class_interaction
  use class_bubles
  use class_ResponseFunctions
  use class_FreeEnergie

  !use class_FreeEnergie
  implicit none
  !private
  type, public :: quasi1D
     integer::Nequations=0
     real ( kind = wp ) :: t_perp=0
     real ( kind = wp ) :: t_perp2=0
     real ( kind = wp)  :: g1_ini=0, g2_ini=0,g3_ini=0
     real ( kind = wp)  :: g1_perp_ini=0,g2_perp_ini=0,g3_perp_ini=0
     real ( kind = wp ) :: Temperature=0
     real ( kind = wp ) :: E0=0
     real ( kind = wp ) :: l=500
     real ( kind = wp ) :: lc=0,ic=0.0
     type(System):: sys
     type(Interaction)::interaction
     type(buble) ::b
     type(ResponseFunction)::RF
     type(ResponseFunction)::deriv_RF
     type(FreeEnergie):: FE
  
   contains
     procedure::Intdeallocate
     !procedure::TempVarflow
     procedure::flow
     procedure::TheDerivative
     procedure::pack
     procedure::unpack
     procedure::WriteOutputs
     procedure::End
     procedure::equal
     GENERIC :: ASSIGNMENT(=)=> equal
  end type quasi1D
  interface quasi1D
     module procedure new_quasi1D
  end interface quasi1D
contains
  subroutine Intdeallocate(this)
    class(quasi1D),intent(out)::this
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call this%deriv_RF%RFdeallocate()
  end subroutine Intdeallocate
  subroutine equal(this,from)
    class(quasi1D),intent(out)::this
    type(quasi1D) ,intent(in)::from
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(wp)::val
    this%Nequations=from%Nequations
    this%t_perp=from%t_perp
    this%t_perp2=from%t_perp2
    this%g1_ini=from%g1_ini;this%g2_ini=from%g2_ini;this%g3_ini=from%g3_ini
    this%g1_perp_ini=from%g1_perp_ini;this%g2_perp_ini=from%g2_perp_ini
    this%g3_perp_ini=from%g3_perp_ini
    this%Temperature=from%Temperature
    this%E0=from%E0
    this%l=from%l
    this%sys=from%sys
    this%interaction=from%interaction
    this%b=from%b
    
    this%RF=from%RF
    this%deriv_RF=from%deriv_RF
    this%FE=from%FE
  end subroutine equal
  function new_quasi1D(fileN,ch)
    type(quasi1D)::new_quasi1D
    character(*) , intent(in)   :: fileN
    character(*) ,optional::ch
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(wp)::val
    if(present(ch))then
       new_quasi1D%sys=System(fileN,ch)
    else
       new_quasi1D%sys=System(fileN)
    end if

    new_quasi1D%t_perp=new_quasi1D%sys%in%t_perp_ini
    new_quasi1D%t_perp2=new_quasi1D%sys%in%t_perp2_ini
    new_quasi1D%g1_ini=new_quasi1D%sys%in%g1_ini
    new_quasi1D%g2_ini=new_quasi1D%sys%in%g2_ini
    new_quasi1D%g3_ini=new_quasi1D%sys%in%g3_ini
    new_quasi1D%g1_perp_ini=new_quasi1D%sys%in%g1_perp_ini
    new_quasi1D%g2_perp_ini=new_quasi1D%sys%in%g2_perp_ini
    new_quasi1D%g3_perp_ini=new_quasi1D%sys%in%g3_perp_ini
    new_quasi1D%Temperature=new_quasi1D%sys%in%Temperature
    new_quasi1D%E0=new_quasi1D%sys%in%E0
    
    new_quasi1D%interaction=Interaction(new_quasi1D%sys)
    
    new_quasi1D%b=buble(new_quasi1D%sys)
    
    val=1.0
    new_quasi1D%RF=ResponseFunction(new_quasi1D%sys%in%N_patche,val)
    val=0.0
    new_quasi1D%deriv_RF=ResponseFunction(new_quasi1D%sys%in%N_patche,val)
    new_quasi1D%FE=FreeEnergie()
    new_quasi1D%Nequations=new_quasi1D%interaction%Neq+new_quasi1D%RF%Neq&
         +new_quasi1D%FE%Neq
    
  end function new_quasi1D
  logical function flow(this,l_ini,l_fin)
    class(quasi1D)::this
    real(kind=wp),intent(inout)::l_ini,l_fin
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer:: flag
    real ( kind = wp ):: abs_err
    real(kind=wp),allocatable:: y(:)
    real(kind=wp),allocatable:: yp(:)
    real ( kind = wp )::RG_error,cst=0.0
    logical::arret
    flow=.false.
    arret=.false.
    this%deriv_RF=ResponseFunction(this%sys%in%N_patche,cst)
    
    allocate(yp(this%Nequations));yp=0.0_wp
    allocate(y(this%Nequations))
    
    call this%pack(y)
    abs_err=this%sys%in%error
    RG_error=this%sys%in%error
    flag = -1
    baba:do while ( flag < 0 .or. flag==4)
       call rkf45_d(fd, this%Nequations,y,yp,l_ini,&
            l_fin,RG_error,abs_err, flag)
       if (abs(this.lc-l_ini)<=1E-6)then
          this.ic=this.ic+1.0_wp
          if (this.ic>=10)then
             arret=.true.
             this.ic=0.0_wp
             this.lc=0.0_wp
          end if
       end if
       this.lc=l_ini
       print*,"#x=",l_ini,this.lc,this.ic,arret
       if(arret) exit baba
    end do baba
    
    if (abs(flag) /= 2)then
       if(abs(flag)==6.or.arret)then
          print*,"#Critical region probably reached for Temperature=",this%Temperature
          flow=.true.
       else
          print*,"Problem in flow function -> class_quasi1D.f90"
          call error(flag)
          stop
       end if
    else
       if(arret)then
          print*,"#Critical region probably reached for Temperature=",this%Temperature
          flow=.true.
       else
          !print*,"#flag=",flag
          call this%unpack(y)
       end if
    end if
  contains
    subroutine fd( t, y1, yp1 )
      implicit none
      real ( kind = wp ) :: t
      real ( kind = wp ) :: y1(this%Nequations)
      real ( kind = wp ) :: yp1(this%Nequations)
      !print*,"In fd",t,size(y1)
      call this%unpack(y1)
      call this%TheDerivative(t)
      call this%pack(yp1,"dy")
    end subroutine fd
  end function flow
  subroutine TheDerivative(this,x)
    class(quasi1D)::this
    real(kind = wp), intent(in):: x
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = wp)::E0El
    if(abs(this%l-x)>=1E-2*this%sys%in%error)then
       if (this%t_perp2==0)then
          call this%b%values_tp20(x,this%Temperature,this%t_perp,this%t_perp2)
       else
          call this%b%values(x,this%Temperature,this%t_perp,this%t_perp2)
       end if
    else
       this%l=x
    end if
    !this%b%IC=0.0_wp
    call this%interaction%IntDerivative(this%b)
    call this%deriv_RF%RFDerivative(this%RF,this%b,this%interaction)
    E0El=exp(-x)*(1.0-exp(-x))
    call this%FE%FE_Derivative(this%interaction,this%sys%in%N_patche,E0El)
  end subroutine TheDerivative
  subroutine pack(this,y,ch)
    class(quasi1D)::this
    real ( kind = wp ):: y(:)
    character(*),optional::ch
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if(present(ch))then
       y(1:this%interaction%Neq)=this%interaction%Int_pack(ch)
       y(this%interaction%Neq+1:this%Nequations-this%FE%Neq)=this%deriv_RF%RF_pack()
       y(this%Nequations)=this%FE%FE_pack(ch)
    else
       y(1:this%interaction%Neq)=this%interaction%Int_pack()
       y(this%interaction%Neq+1:this%Nequations-this%FE%Neq)=this%RF%RF_pack()
       y(this%Nequations)=this%FE%FE_pack()
    end if
  end subroutine pack
  subroutine unpack(this,y)
    class(quasi1D)::this
    real ( kind = wp ) :: y(:)
    integer::k
    real ( kind = wp ),allocatable::yy(:)
    
    k=this%interaction%Neq
    allocate(yy(this%interaction%Neq))
    yy=y(1:this%interaction%Neq)
    call this%interaction%Int_unpack(yy) 
    deallocate(yy)
    allocate(yy(this%RF%Neq))
    yy(:)=y(this%interaction%Neq+1:this%Nequations-this%FE%Neq)
    call this%RF%getValues(yy)
    deallocate(yy)
    call this%FE%getValues(y(this%Nequations))
  end subroutine unpack
  subroutine End(this)
    class(quasi1D)::this
    call timestamp_file(this%sys%out%fg)
    call timestamp_file(this%sys%out%chiMax)
    call timestamp_file(this%sys%out%fgmoy)
    call timestamp_file(this%sys%out%chiCDW)
    call timestamp_file(this%sys%out%chiSDW)
    call timestamp_file(this%sys%out%chiSS)
    call timestamp_file(this%sys%out%chiST)
    call timestamp_file(this%sys%out%FreeE)

  end subroutine End
  subroutine WriteOutputs(this)
    class(quasi1D)::this
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call this%interaction%InteractionWrite(this%t_perp2,this%Temperature,this%sys%out)
    call this%RF%RFwrite(this%t_perp2,this%Temperature,this%sys%out)
    call this%FE%FEwrite(this%t_perp2,this%Temperature,this%sys%out)
  end subroutine WriteOutputs
end module class_quasi1D
