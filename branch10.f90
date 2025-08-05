!This is a working fortran file


!This subroutine read a input file using namelist
subroutine read_input(ny,nz,Ly,Lz,res_tolerance,A,Rm,linear_solver ,pre_cond,del_t,alpha_k)
 implicit none
 integer:: ny,nz,Rm
 real :: Ly,Lz,res_tolerance,A,del_t,alpha_k
 character :: linear_solver*10 ,pre_cond*10
 
 namelist/input_parameters/ny,nz,Ly,Lz,res_tolerance,A,linear_solver ,pre_cond,Rm,del_t,alpha_k
 open(unit=10, file='../input_file/input.txt')
 read(10, nml=input_parameters)
 close(10)
end subroutine read_input

!This subroutine generates a grid
subroutine generate_grid(ny,nz,Ly,Lz,A,gridpts_y,gridpts_z)
 implicit none
 INTEGER,intent(in) :: ny,nz
 INTEGER :: i,j,k
 REAL,intent(in) :: Ly,Lz,A
 REAL(kind=8) :: dy,dz
 REAL(kind=8),intent(out) :: gridpts_y(0:ny),gridpts_z(0:nz)
 dy=Ly/(ny)
 dz=Lz/(nz)
 

 do i=0,ny
   gridpts_y(i)= -1.d0+ dy*(i) 
   gridpts_z(i)= -1.d0+ dz*(i)  
 end do
 !This loop stretches the grid point

  gridpts_y= tanh(A*gridpts_y)/tanh(A)
  gridpts_z= tanh(A*gridpts_z)/tanh(A) 

 open(unit=20, file='../output_file/grid_pts_stretched_b10.dat')
 
    do k=0,nz
       
     do j=0,ny 
       write(20, '(1X,F10.6,2X,F10.6)'),gridpts_y(k),gridpts_z(j)
     end do 
       write(20,*), " "
     do j=0,ny 
       write(20, '(1X,F10.6,2X,F10.6)'),gridpts_y(j),gridpts_z(k)
     end do 
       write(20,*), " "
    end do
  
  
 close(20)
 print*, "Grid_generation done"
end subroutine generate_grid


module node_classes
!This is a class for node informations
type node

 real(kind=8):: y_loc , z_loc
 real(kind=8),allocatable :: coeffs(:)
 integer,allocatable :: col(:)
end type


contains
  
  subroutine fill_node_array(node_arr,grid_pty,grid_ptz,ny,nz,n,f,alpha_k)
  
    implicit none
    integer,intent(in):: ny,nz,n 
    real,intent(in) ::f,alpha_k
    real(kind=8),intent(in) :: grid_pty(0:ny),grid_ptz(0:nz)
    integer :: i,j,k,count,c
   
    
    
    type(node), allocatable, intent(inout) :: node_arr(:)
 

    count=0
    
    allocate(node_arr(n))
    print*,'size of node array',size(node_arr)
     do j=0,ny
      do k=0,nz
       
       i=k*(ny+1)+j+1
       node_arr(i)%y_loc=grid_pty(j)
       node_arr(i)%z_loc=grid_ptz(k)
       
       
       
       !print*,i,j,k
       !count=count+1                               
       if(j>0 .and. j<ny .and. k>0 .and. k<nz ) then
          allocate(node_arr(i)%coeffs(5))
          allocate(node_arr(i)%col(5))
          
          node_arr(i)%coeffs(3)= -(alpha_k**2+f+2/(grid_pty(j+1)-grid_pty(j))/(grid_pty(j)-grid_pty(j-1))+2/(grid_ptz(k+1)-grid_ptz(k))/(grid_ptz(k)-grid_ptz(k-1)))
          node_arr(i)%coeffs(5)= 2/(grid_ptz(k+1)-grid_ptz(k-1))/(grid_ptz(k+1)-grid_ptz(k))
          node_arr(i)%coeffs(1)=2/(grid_ptz(k+1)-grid_ptz(k-1))/(grid_ptz(k)-grid_ptz(k-1))
          node_arr(i)%coeffs(4)=2/(grid_pty(j+1)-grid_pty(j-1))/(grid_pty(j+1)-grid_pty(j))
          node_arr(i)%coeffs(2)=2/(grid_pty(j+1)-grid_pty(j-1))/(grid_pty(j)-grid_pty(j-1))
          
          node_arr(i)%col(3)=(k)*(ny+1)+j+1
          node_arr(i)%col(5)=(k+1)*(ny+1)+j+1
          node_arr(i)%col(1)=(k-1)*(ny+1)+j+1
          node_arr(i)%col(4)=(k)*(ny+1)+j+1+1
          node_arr(i)%col(2)=(k)*(ny+1)+j-1+1
          count=count+1 
        else if(j==0)then
             c=(grid_pty(j+2)-grid_pty(j+1))/(grid_pty(j+1)-grid_pty(j))
             allocate(node_arr(i)%coeffs(3))
             allocate(node_arr(i)%col(3))
             node_arr(i)%coeffs(1)= -(2*c+c**2)/(grid_pty(j+1)-grid_pty(j))/(c+c**2)
             node_arr(i)%coeffs(2)= (1+c)**2/(grid_pty(j+1)-grid_pty(j))/(c+c**2)
             node_arr(i)%coeffs(3)= -1.0/(grid_pty(j+1)-grid_pty(j))/(c+c**2)
             node_arr(i)%col(1)=(k)*(ny+1)+j+1
             node_arr(i)%col(2)=(k)*(ny+1)+j+1+1
             node_arr(i)%col(3)=(k)*(ny+1)+j+2+1
          
        else if(j==ny)then
            !print*,'ok1'
            c=(grid_pty(j-2)-grid_pty(j-1))/(grid_pty(j-1)-grid_pty(j))
            allocate(node_arr(i)%coeffs(3))
            allocate(node_arr(i)%col(3))
            !print*,'ok2'
            node_arr(i)%coeffs(3)= (2*c+c**2)/(grid_pty(j-1)-grid_pty(j))/(c+c**2)
            node_arr(i)%coeffs(2)= -(1+c)**2/(grid_pty(j-1)-grid_pty(j))/(c+c**2)
            node_arr(i)%coeffs(1)= 1.0/(grid_pty(j-1)-grid_pty(j))/(c+c**2)
            node_arr(i)%col(3)=(k)*(ny+1)+j+1
            node_arr(i)%col(2)=(k)*(ny+1)+j-1+1
            node_arr(i)%col(1)=(k)*(ny+1)+j-2+1
          
        else
          allocate(node_arr(i)%coeffs(1))
          allocate(node_arr(i)%col(1))
          node_arr(i)%coeffs(1)= 1.0
          node_arr(i)%col(1)=(k)*(ny+1)+j+1
        end if
       
                   
       end do
     end do
      print*,"l143"
     open(unit=30, file='../output_file/node_info_b10_coeffs.dat')
     open(unit=40, file='../output_file/node_info_b10_cols.dat')
      write(30,*),"ap       an       as        aw       ae"
 
      do j=1,n
    
  
    
        write(30,"(5F15.6)"),node_arr(j)%coeffs!,node_arr(j)%y_loc,node_arr(j)%z_loc
        write(40,*),node_arr(j)%col
      end do 
  
      close(30)
      close(40)
      print*,"Node array filled" 
    
    end subroutine fill_node_array
    
    subroutine fill_q_array(node_arr,q_array,f,alpha_k)
       type(node), allocatable, intent(inout) :: node_arr(:)
       real(kind=8),allocatable, intent(out):: q_array(:)
       real, intent(in)::f,alpha_k
       real,parameter :: PI=4*atan(1.d0)
       integer :: i
       real :: y,z,constant
       
       constant=0.01*((f+alpha_k**2)+8*(PI)**2)
       allocate(q_array(0:size(node_arr)-1))
       open(unit=30, file='../output_file/rhs_b10.dat')
       do i=1,size(node_arr)
         y=node_arr(i)%y_loc
         z=node_arr(i)%z_loc
         q_array(i)=-constant*cos(2*PI*y)*sin(2*PI*z)
         if(abs(z)>0.9999999 .or. abs(y)>0.9999999)then
          q_array(i)=0.d0
         end if
         write(30,*),y,z,q_array(i)
       end do
       close(30)
       
        
    end subroutine fill_q_array
     

end module node_classes


module matrix_vector_op
   
contains
   subroutine mat_vec_mult(mat,vec_in,vec_out)
   use node_classes
   implicit none
   type(node), allocatable, intent(inout) :: mat(:)
   real(kind=8), allocatable ,intent (in) :: vec_in(:)
   real(kind=8), allocatable ,intent (out) :: vec_out(:)
   integer :: row,i
   real(kind=8) :: sum
   allocate(vec_out(size(vec_in)))
   do row=1, size(mat)
       sum=0
       do i=1,size(mat(row)%coeffs)
         sum=sum+mat(row)%coeffs(i)*vec_in(mat(row)%col(i))
       end do
        vec_out(row)=sum
   end do
     
   end subroutine mat_vec_mult
   
   
   
   subroutine vec_vec_mult(vec_in1,vec_in2,scal_out)
   use node_classes
   implicit none
   
   real(kind=8) , allocatable ,intent (in) :: vec_in1(:),vec_in2(:)
   real(kind=8)  ,intent (out) :: scal_out
   integer :: row,i
   real(kind=8):: sum
   sum=0
  
   do row=1, size(vec_in1)
       
        sum=sum+vec_in1(row)*vec_in2(row)
        
   end do
     scal_out=sum
   end subroutine vec_vec_mult
   
   subroutine scale_vec_mult(vec_in,scal_in,vec_out)
   
   implicit none
   
   real(kind=8) , allocatable ,intent (in) :: vec_in(:)
   real(kind=8) , allocatable ,intent (out) :: vec_out(:)
   real(kind=8) ,intent (in) :: scal_in
   integer :: row,i
   
   allocate(vec_out(size(vec_in)))
   do row=1, size(vec_in)
       
        vec_out(row)=scal_in*vec_in(row)
        
   end do
     
   end subroutine scale_vec_mult
   
   subroutine vec_vec_sum(vec_in1,vec_in2,vec_out,option)
   use node_classes
   implicit none
   
   real(kind=8) , allocatable ,intent (in) :: vec_in1(:),vec_in2(:)
   real(kind=8)  ,allocatable, intent (out) :: vec_out(:)
   integer,intent(in) :: option     ! +1 for sum and -1 for difference
   integer :: row,i
   
   
   
   allocate(vec_out(size(vec_in1)))
      
      select case(option)
      
      case (1)
      do row=1, size(vec_in1)
       
         vec_out(row)=vec_in1(row)+vec_in2(row)
        
      end do
      
      case (-1)
      do row=1, size(vec_in1)
       
         vec_out(row)=vec_in1(row)-vec_in2(row)
        
      end do
      end select
   end subroutine vec_vec_sum
   
   subroutine elementwise_div(vec_in1,vec_in2,vec_out)
   use node_classes
   implicit none
   
   real(kind=8) , allocatable ,intent (in) :: vec_in1(:),vec_in2(:)
   real(kind=8)  ,allocatable, intent (out) :: vec_out(:)
  
   integer :: row,i
   
   
   
   allocate(vec_out(size(vec_in1)))
      
      
      
   
      do row=1, size(vec_in1)
         if(vec_in2(row)/=0) then
           vec_out(row)=vec_in1(row)/vec_in2(row)
         else
           vec_out(row)=vec_in1(row) 
         end if
      end do
      
      
  
   end subroutine elementwise_div
   
   subroutine get_inf_norm(vec_in,inf_norm,index)
   real(kind=8),allocatable ,intent (in) :: vec_in(:)
   real(kind=8),intent(out):: inf_norm
   real(kind=8)::max_ele
   integer:: i
   integer,intent(out):: index
   max_ele=0
   do i=1,size(vec_in)
     if(max_ele<abs(vec_in(i))) then
     max_ele=abs(vec_in(i))
     index=i
     endif
   end do
   inf_norm=max_ele
   end subroutine get_inf_norm
   
   subroutine get_2_norm(vec_in,norm_2)
   real(kind=8),allocatable ,intent (in) :: vec_in(:)
   real(kind=8),intent(out):: norm_2
   real(kind=8)::sum
   integer:: i
   sum=0
   do i=1,size(vec_in)
     sum=sum+vec_in(i)**2
   end do
   norm_2=sqrt(sum)
   end subroutine get_2_norm
    
end module matrix_vector_op

module linear_solvers
 
contains
  subroutine IL_solve(A,L,y,b,ny,nz,n)
    use node_classes
    use matrix_vector_op
    implicit none
    type(node), allocatable, intent(inout) :: A(:) 
    type(node), allocatable, intent(inout)::   L(:)
    !real(kind=8),allocatable,intent(inout) :: x(:)
    real(kind=8),allocatable,intent(in) :: b(:)
    real(kind=8),allocatable,intent(inout) :: y(:)
    integer,intent(in) :: ny,nz,n
    integer :: i,j
   
    !solve l
     do i=1,size(b)
       if(i==1) then
       y(i)=b(i)/L(i)%coeffs(1)
       else if(i>=2 .and. i<=ny+1) then
       y(i)=(b(i)-L(i)%coeffs(2)*y(i-1))/L(i)%coeffs(1)
       else if(i>=ny+2 .and. i<=n) then
       y(i)=(b(i)-L(i)%coeffs(2)*y(i-1)-L(i)%coeffs(3)*y(i-ny-1))/L(i)%coeffs(1)
       end if
     end do
   
  end subroutine IL_solve
  subroutine ILU_solve(A,L,U,x,b,ny,nz,n)
    use node_classes
    use matrix_vector_op
    implicit none
    type(node), allocatable, intent(inout) :: A(:) 
    type(node), allocatable, intent(inout)::   L(:)
    type(node), allocatable, intent(inout)::   U(:)
    real(kind=8),allocatable,intent(inout) :: x(:)
    real(kind=8),allocatable,intent(in) :: b(:)
    real(kind=8),allocatable :: y(:)
    integer,intent(in) :: ny,nz,n
    integer :: i,j
    allocate(y(size(b)))
    !solve l
     do i=1,size(b)
       if(i==1) then
       y(i)=b(i)/L(i)%coeffs(1)
       else if(i>=2 .and. i<=ny+1) then
       y(i)=(b(i)-L(i)%coeffs(2)*y(i-1))/L(i)%coeffs(1)
       else if(i>=ny+2 .and. i<=n) then
       y(i)=(b(i)-L(i)%coeffs(2)*y(i-1)-L(i)%coeffs(3)*y(i-ny-1))/L(i)%coeffs(1)
       end if
     end do
    !solve u
     do i=size(b),1,-1
       if(i==n) then
       x(i)=y(i)
       else if(i<=n-1 .and. i>=n-ny-1+1) then
       x(i)=y(i)-U(i)%coeffs(1)*x(i+1)
       else if(i<=n-ny-1 .and. i>=1) then
       x(i)=y(i)-U(i)%coeffs(1)*x(i+1)-U(i)%coeffs(2)*x(i+ny+1)
       end if
     end do
     
  end subroutine ILU_solve
  
  subroutine ILU_preconditioner(A,L,U,b,q,ny,nz,n)
    use node_classes
    use matrix_vector_op
    implicit none
    type(node), allocatable, intent(inout) :: A(:) 
    type(node), allocatable,intent(out)::   L(:)
    type(node), allocatable,intent(out)::   U(:)
    real(kind=8),allocatable,intent(inout) :: b(:)
    real(kind=8),allocatable,intent(in) :: q(:)
    integer,intent(in)::ny,nz,n
    logical, dimension(:), allocatable :: mask
    integer :: i,j,k,s,t,w,bandsize
    real::  sum
    integer*4 :: count_rate, count_max, start_count, end_count
    real*8 :: start_time, end_time
  
    
    bandsize=ny+2
    allocate(L(size(A)),U(size(A)))
    open(unit=50,file="../output_file/ILU_L10.dat")
    open(unit=60,file="../output_file/ILU_U10.dat")

    do i=1,size(A)
      
      allocate(L(i)%coeffs(3),L(i)%col(3),U(i)%coeffs(2),U(i)%col(2))
     
      !L(i)%coeffs(1)=>lp
      !L(i)%coeffs(2)=>lw
      !L(i)%coeffs(3)=>ls
      !U(i)%coeffs(1)=>ue
      !U(i)%coeffs(2)=>un
      L(i)%col(1)=i
      L(i)%col(2)=i-1
      L(i)%col(3)=i-ny-1
      U(i)%col(1)=i+1
      U(i)%col(2)=i+ny+1 
      
      if(i==1)then
        L(i)%coeffs(3)=0.d0
        L(i)%coeffs(2)=0.d0
        L(i)%coeffs(1)=A(i)%coeffs(1)
        
        mask=(A(i)%col==U(i)%col(1))
        !********************************
        U(i)%coeffs(1)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            U(i)%coeffs(1)=A(i)%coeffs(t)/L(i)%coeffs(1)
          end if
         end do
        !******************************** 
        mask=(A(i)%col==U(i)%col(2))
        U(i)%coeffs(2)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            U(i)%coeffs(2)=A(i)%coeffs(t)/L(i)%coeffs(1)
         
          end if
         end do
         !***************************************
      else if(i>=2 .and. i<=ny+1) then
        L(i)%coeffs(3)=0.d0
        
         mask=(A(i)%col==L(i)%col(2))
         L(i)%coeffs(2)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            L(i)%coeffs(2)=A(i)%coeffs(t)
          
          end if
         end do
        !*****************************************************************
        L(i)%coeffs(1)=A(i)%coeffs(1)-L(i)%coeffs(2)*U(i-1)%coeffs(1)
        !***************************************************************
         mask=(A(i)%col==U(i)%col(1))
         U(i)%coeffs(1)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            U(i)%coeffs(1)=A(i)%coeffs(t)/L(i)%coeffs(1)
          
          end if
         end do
         !************************************************************** 
        mask=(A(i)%col==U(i)%col(2))
         U(i)%coeffs(2)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            U(i)%coeffs(2)=A(i)%coeffs(t)/L(i)%coeffs(1)
          end if
         end do
         
      else if(i>=ny+2 .and. i<=n-ny-1) then
        
         mask=(A(i)%col==L(i)%col(3))
         L(i)%coeffs(3)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            L(i)%coeffs(3)=A(i)%coeffs(t)
         
          end if
         end do
        
        !********************************************
         mask=(A(i)%col==L(i)%col(2))
         L(i)%coeffs(2)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            L(i)%coeffs(2)=A(i)%coeffs(t)
          
          end if
         end do
         !*****************************************************
        L(i)%coeffs(1)=A(i)%coeffs(1)-L(i)%coeffs(2)*U(i-1)%coeffs(1)-L(i)%coeffs(3)*U(i-ny-1)%coeffs(2)
        !****************************************************************
        mask=(A(i)%col==U(i)%col(1))
        U(i)%coeffs(1)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            U(i)%coeffs(1)=A(i)%coeffs(t)/L(i)%coeffs(1)
         
          end if
         end do
        !************************************************************************  
        mask=(A(i)%col==U(i)%col(2))
         U(i)%coeffs(2)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            U(i)%coeffs(2)=A(i)%coeffs(t)/L(i)%coeffs(1)
          end if
         end do
        !**********************************************************
      else if (i>=n-ny-1+1 .and. i<=n-1) then
      
        mask=(A(i)%col==L(i)%col(3))
        L(i)%coeffs(3)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            L(i)%coeffs(3)=A(i)%coeffs(t)
          
          end if
         end do
        
        !**************************************************************
         mask=(A(i)%col==L(i)%col(2))
         L(i)%coeffs(2)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            L(i)%coeffs(2)=A(i)%coeffs(t)
          
          end if
         end do
         !*************************************************************
        L(i)%coeffs(1)=A(i)%coeffs(1)-L(i)%coeffs(2)*U(i-1)%coeffs(1)-L(i)%coeffs(3)*U(i-ny-1)%coeffs(2)
         !***************************************************************
        mask=(A(i)%col==U(i)%col(1))
        U(i)%coeffs(1)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            U(i)%coeffs(1)=A(i)%coeffs(t)/L(i)%coeffs(1)
         
          end if
         end do
         !*********************************************************** 
        
        U(i)%coeffs(2)=0.d0
        
      else if (i==n) then
        mask=(A(i)%col==L(i)%col(3))
         L(i)%coeffs(3)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            L(i)%coeffs(3)=A(i)%coeffs(t)
          
          end if
         end do
         !****************************************************
        
         mask=(A(i)%col==L(i)%col(2))
         L(i)%coeffs(2)=0.d0
         do t=1,size(A(i)%col)
          if(mask(t))then
            L(i)%coeffs(2)=A(i)%coeffs(t)
          
          end if
         end do
         !**********************************************************
        L(i)%coeffs(1)=A(i)%coeffs(1)-L(i)%coeffs(2)*U(i-1)%coeffs(1)-L(i)%coeffs(3)*U(i-ny-1)%coeffs(2)
        U(i)%coeffs(1)=0.d0
        U(i)%coeffs(2)=0.d0
      end if
      
      
      
        
      write(50,"(3F15.6)"),L(i)%coeffs
      write(60,"(2F15.6)"),U(i)%coeffs
      
    end do
    !print*,"coeff of l at 300th row",L(300)%coeffs
    !print*,"coeff of A at 300th row",A(300)%coeffs
    close(50)
    close(60)
    
    
   
    
  end subroutine ILU_preconditioner
  
  
  
  subroutine biCGSTAB(A,b,q,tol,ny,nz,PRECOND)
  use node_classes
  use matrix_vector_op
  implicit none
  type(node), allocatable, intent(inout) :: A(:)
  type(node), allocatable                :: L(:)
  type(node), allocatable                :: U(:)
  real(kind=8),allocatable,intent(in)    :: q(:)
  real(kind=8),allocatable,intent(inout) :: b(:) 
  real,intent(in)                        :: tol
  integer,intent(in)                     ::ny,nz
  character,intent(in)                   ::PRECOND*10
  real(kind=8),allocatable               :: res_ini(:)
  real(kind=8),allocatable               :: res(:)
  real(kind=8),allocatable               :: res_norm(:)
  real(kind=8),allocatable               :: poly(:)
  real(kind=8),allocatable               :: k_inv_poly(:)
  real(kind=8),allocatable               :: k1_inv_poly(:)
  real(kind=8),allocatable               :: v(:)
  real(kind=8),allocatable               :: k_inv_s(:)
  real(kind=8),allocatable               :: k1_inv_s(:)
  real(kind=8),allocatable               :: k1_inv_t(:)
  real(kind=8),allocatable               :: s(:)
  real(kind=8),allocatable               :: d_arr(:) !dummy array
  real(kind=8),allocatable               :: t(:)
  real(kind=8),allocatable               :: alpha_poly(:)
  real(kind=8),allocatable               :: omega_s(:)
  real(kind=8),allocatable               :: omega_t(:)
  real(kind=8),allocatable               :: omega_v(:)
  real(kind=8)                           :: row,alpha,omega,v_dot_res,t_dot_s,t_dot_t,beta,row_new,inf_n,inf_rhs,norm_2,norm_2_rhs
  integer                                :: n,i,j,k
  
  n=size(q)
  call get_inf_norm(q,inf_rhs,j)
  call get_2_norm(q,norm_2_rhs)
  open(unit=20,file='../output_file/ini_b_q_out_b10.dat')
   write(20,*) ,'iter        b          q'
   
     do  i=1,n
     write(20,*) ,i,b(i),q(i)
     end do
  
  close(20)
  !Allocating memories
  allocate(poly(n),s(n),t(n),k_inv_poly(n),k_inv_s(n),k1_inv_poly(n),k1_inv_s(n),k1_inv_t(n))
  
  !Initial operation
  call mat_vec_mult(A,b,d_arr)           !d_arr=Ab
  call vec_vec_sum(q,d_arr,res_ini,-1)   !res_o=q-Ab
  call vec_vec_mult(res_ini,res_ini,row) !row=res_o.res_o
  poly=res_ini
  res =res_ini
  open(unit=30,file='ini_res_b10.dat')
   write(30,*) ,'iter        res'
   do  i=1,n
     write(30,*) ,i,res(i)
    end do
  close(30)
  print*,'BiCGSTAB started'
  !Preconditioning
  if(PRECOND=="ILU")then
     call ILU_preconditioner(A,L,U,b,q,ny,nz,n)
  end if
  !Iteration start
  open(unit=10,file='../output_file/bicgstab_out_b10.dat')
  do i=1,1
   if(PRECOND=="ILU")then
   
    call ILU_solve(A,L,U,k_inv_poly,poly,ny,nz,n)   ! k_inv_poly= (K^-1)p
    call IL_solve(A,L,k1_inv_poly,poly,ny,nz,n)   ! k1_inv_poly= (K1^-1)p
    
    call mat_vec_mult(A,k_inv_poly,v)             !v=A k_inv_poly
    
   else
    call mat_vec_mult(A,poly,v) 
   end if
   call vec_vec_mult(v,res_ini,v_dot_res)  !v.res_o
   
   alpha=row/v_dot_res 
   print*,alpha
   call scale_vec_mult(v,alpha,d_arr)      ! alpha * v
   call vec_vec_sum(res,d_arr,s,-1)        ! s=r-alpha * v
   if(PRECOND=="ILU")then
    call ILU_solve(A,L,U,k_inv_s,s,ny,nz,n)   ! k_inv_s= (K^-1)s
    call mat_vec_mult(A,k_inv_s,t)          !t=A k_inv_s
    call IL_solve(A,L,k1_inv_s,s,ny,nz,n)     !K_inv s
    call IL_solve(A,L,k1_inv_t,t,ny,nz,n)     !K_inv t
    call vec_vec_mult(k1_inv_t,k1_inv_s,t_dot_s)          !(K^-1)t.(K^-1)s
    call vec_vec_mult(k1_inv_t,k1_inv_t,t_dot_t)          !(K^-1)t.(K^-1)t
    omega= t_dot_s/t_dot_t 
    call scale_vec_mult(k_inv_poly,alpha,alpha_poly)   !alpha * (K^-1)p
    call scale_vec_mult(k_inv_s,omega,omega_s)         !omega * (K^-1)s
   else
    call mat_vec_mult(A,s,t)                !t=As
    call vec_vec_mult(t,s,t_dot_s)          !t.s
    call vec_vec_mult(t,t,t_dot_t)          !t.t
    omega= t_dot_s/t_dot_t 
    call scale_vec_mult(poly,alpha,alpha_poly)   !alpha * poly
    call scale_vec_mult(s,omega,omega_s)         !omega*s
   
   end if
              
   open(unit=20,file='../output_file/kinv_poly_kinv_s.dat')
   do j=1,size(k_inv_poly)
    write(20,*),k1_inv_poly(j),k_inv_poly(j),poly(j)
   end do
   
   call vec_vec_sum(alpha_poly,omega_s,d_arr,1) !alpha * poly + omega*s=d_arr
   call vec_vec_sum(b,d_arr,b,1)                !bnew=b+d_arr
   call scale_vec_mult(t,omega,omega_t)         !omega*t
   call vec_vec_sum(s,omega_t,res,-1)           !res_new=s-omega*t 
   call vec_vec_mult(res,res_ini,row_new)       !row_new=res_new.res_o   
   beta=alpha*row_new/omega/row   
   call scale_vec_mult(v,omega,omega_v)         !omega*v
   call vec_vec_sum(poly,omega_v,d_arr,-1)      !d_arr=p -  omega*v        
   call scale_vec_mult(d_arr,beta,d_arr)        !d_arr=beta*(d_arr)
   call vec_vec_sum(res,d_arr,poly,1)           !poly=res_new+d_arr
   !call elementwise_div(res,q,res_norm)  
    
   call get_inf_norm(res,inf_n,j)
   call get_2_norm(res,norm_2)
   row=row_new                                  !row updated
   inf_n=inf_n/inf_rhs
   norm_2=norm_2/norm_2_rhs
   if(inf_n<tol) exit
   write(10,"(1X,I5,1X,4f10.6)") ,i,inf_n,A(j)%y_loc,A(j)%z_loc,norm_2
  end do
  
  close(10)
  close(20) 
  print*,"alpha=",alpha
  open(unit=40,file='../output_file/by_res_final_b10.dat')
   write(40,*) ,'iter         y        z         by       res'
    do  j=0,ny
     do k=0,nz
        i=(k)*(ny+1)+j+1
        write(40,'(1X,I5,1X,7f10.6)') ,i,A(i)%y_loc,A(i)%z_loc,b(i),res(i)
     end do
       write(40,*) " "
    end do
  close(40)
  deallocate(res,v,poly,s,t,alpha_poly,omega_s,omega_t,omega_v)
  end subroutine biCGSTAB
  
  subroutine LU_solve(A,L,U,x,b)
   use node_classes
   use matrix_vector_op
   implicit none
   type(node), allocatable, intent(inout) ::   A(:)
   type(node), allocatable, intent(in) ::   L(:)
   type(node), allocatable, intent(in) ::   U(:)
   real(kind=8),dimension(:),allocatable,intent(inout) :: x
   real(kind=8),dimension(:),allocatable :: y
   real(kind=8),dimension(:),allocatable,intent(in) ::  b
   integer :: i,j
   real:: sum
   !Forward substitution
   allocate(y(size(b)))
   do i=1,size(b) 
     
     sum=0
     do j=2,size(L(i)%coeffs)
       sum=sum+L(i)%coeffs(j)*y(i-(j-1))
     end do
     y(i)=b(i)-sum
     y(i)=y(i)/L(i)%coeffs(1)
   end do 
   
   !Backward substitution
   do i=size(b),1,-1 
    sum=0
    do j=2,size(U(i)%coeffs)
      sum=sum+U(i)%coeffs(j)*x(i+(j-1))
    end do
    x(i)=y(i)-sum
    
   end do 
   open(unit=70,file="../output_file/LU_soln_b10.dat")
     write(70,*),"y    z    by"
     do i=1,size(x)
      write(70,"(3F10.6)"),A(i)%y_loc,A(i)%z_loc,x(i)
     end do
   close(70)
  end subroutine LU_solve
  
  
  subroutine LU_decompositon(A,b,q,ny,nz,n)
    use node_classes
    use matrix_vector_op
    implicit none
    type(node), allocatable, intent(inout) :: A(:) 
    type(node), allocatable::   L(:)
    type(node), allocatable::   U(:)
    real(kind=8),allocatable,intent(inout) :: b(:)
    real(kind=8),allocatable,intent(in) :: q(:)
    integer,intent(in)::ny,nz,n
    logical, dimension(:), allocatable :: mask
    integer :: i,j,k,s,t,w,bandsize
    real::  sum
    integer*4 :: count_rate, count_max, start_count, end_count
    real*8 :: start_time, end_time
    
    bandsize=ny+2
    allocate(L(size(A)),U(size(A)))
    open(unit=50,file="../output_file/L_10.dat")
    open(unit=60,file="../output_file/U_10.dat")

    do i=1,size(A)
      if(i<=bandsize) then 
      allocate(L(i)%coeffs(i),L(i)%col(i),U(i)%coeffs(bandsize),U(i)%col(bandsize))
     
       
      else if(i>(n-bandsize)) then
      allocate(L(i)%coeffs(bandsize),L(i)%col(bandsize),U(i)%coeffs(n-i+1),U(i)%col(n-i+1)) 
      
      else
      allocate(L(i)%coeffs(bandsize),L(i)%col(bandsize),U(i)%coeffs(bandsize),U(i)%col(bandsize))
     
      end if
      U(i)%coeffs(1)=1 !Diagonal element is 1
 
      !print *,i
      !solve row elements of L
      do j=size(L(i)%coeffs),1,-1
         sum=0
         s=2 
         do k=j+1,size(L(i)%coeffs)
        
           sum=sum+L(i)%coeffs(k)*U(i-k+1)%coeffs(s)
           s=s+1
         end do
        
         mask=(A(i)%col==i-(j-1))
         L(i)%coeffs(j)=-sum
         do t=1,size(A(i)%col)
          if(mask(t))then
            L(i)%coeffs(j)=A(i)%coeffs(t)+L(i)%coeffs(j)
          end if
         end do
    
      end do

      !solve row elements of U
      do j=2,size(U(i)%coeffs)
         sum=0
         w=j+1
 
         do t=2,size(L(i)%coeffs) 
          if(w>bandsize) then
            exit
          end if
      
          sum=sum+L(i)%coeffs(t)*U(i-t+1)%coeffs(w)
          w=w+1
   
         end do
         mask=(A(i)%col==i+(j-1))
         U(i)%coeffs(j)=-sum
         do t=1,size(A(i)%col)
          if(mask(t))then
            U(i)%coeffs(j)=A(i)%coeffs(t)+U(i)%coeffs(j)
          end if
         end do
        U(i)%coeffs(j)=U(i)%coeffs(j)/L(i)%coeffs(1)
      end do
      write(50,*),L(i)%coeffs
      write(60,*),U(i)%coeffs
      
    end do
    
    close(50)
    close(60)
    
    call system_clock(start_count, count_rate, count_max)
     start_time = real(start_count)/real(count_rate)
    
    call LU_solve(A,L,U,b,q)
    
    call system_clock(end_count, count_rate, count_max)
     end_time = real(end_count)/real(count_rate)
   
    print*, "Wall clock time for substitution:", end_time - start_time, "seconds"
  end subroutine LU_decompositon

  
  subroutine solver(A,b,q,tol,ny,nz,n,solver_name,precond)
   use node_classes
   use matrix_vector_op
   implicit none
   type(node), allocatable, intent(inout) :: A(:)
   real(kind=8),allocatable,intent(in):: q(:)
   real(kind=8),allocatable,intent(inout):: b(:) 
   real,intent(in):: tol
   integer,intent(in)::ny,nz,n 
   character, intent(in):: solver_name*10
   character, intent(in):: precond*10
   if(solver_name=="BICGSTAB")then
     call biCGSTAB(A,b,q,tol,ny,nz,precond)
   else if(solver_name=="LU_METHOD") then
     call LU_decompositon(A,b,q,ny,nz,n)
   else  
     print*,"solver chosen is not available"
   end if
      
  end subroutine solver
  
end module linear_solvers





!This is the main code
PROGRAM  main

 !Calling the modules
 use node_classes
 use matrix_vector_op
 use linear_solvers
 implicit none
 
 
 integer:: ny,nz,Rm,n,j,i,k,band
 real :: Ly,Lz,res_tolerance,A,del_t,alpha_k
 real :: init_by, init_bz,f
 character :: linear_solver*10 ,pre_cond*10
 real(kind=8) ,dimension (:),allocatable:: gridpts_y,gridpts_z
 real(kind=8),dimension(:),allocatable :: q,by,act_b,error_b
 type(node), allocatable :: node_info(:),L(:),U(:)
 real,parameter :: PI=4*atan(1.d0)
 real :: y,z
 integer*4 :: count_rate, count_max, start_count, end_count
 real*8 :: start_time, end_time
 
 call read_input(ny,nz,Ly,Lz,res_tolerance,A,Rm,linear_solver ,pre_cond,del_t,alpha_k)
 print *,ny,nz,Ly,Lz,res_tolerance,A,Rm,linear_solver,del_t ,pre_cond,alpha_k
 
 
 f=3*Rm/2/del_t
 band=ny+2 !No of non-zero diagonals
 n=(ny+1)*(nz+1)
 init_by=0
 init_bz=0
 
 allocate(gridpts_y(0:ny))
 allocate(gridpts_z(0:nz))
 allocate(act_b(n))
 allocate(error_b(n))
 
 
 
 
 call generate_grid(ny,nz,Ly,Lz,A,gridpts_y,gridpts_z)
 call fill_node_array(node_info,gridpts_y,gridpts_z,ny,nz,n,f,alpha_k)
 call fill_q_array(node_info,q,f,alpha_k)
 print *,'size of node_info' , size(node_info)
 print *,'size of q' , size(q)
 deallocate(gridpts_y)!memory freed
 deallocate(gridpts_z)!memory freed
 
 open(unit=30, file='../output_file/act_by_b9.dat')
       do j=0,ny
        do k=0,nz
         i=(k)*(ny+1)+j+1
         y=node_info(i)%y_loc
         z=node_info(i)%z_loc
         act_b(i)=0.01*cos(2.0*PI*y)*sin(2.0*PI*z)
     
         write(30,*),y,z,act_b(i)
       end do
        write(30,*) " "
      end do
 close(30)
 
 allocate(by(size(node_info)))
 by=init_by

 call system_clock(start_count, count_rate, count_max)
     start_time = real(start_count)/real(count_rate)
     
 call solver(node_info,by,q,res_tolerance,ny,nz,n,linear_solver,pre_cond)
 
 call system_clock(end_count, count_rate, count_max)
   end_time = real(end_count)/real(count_rate)
   
print*, "Wall clock time:", end_time - start_time, "seconds"

 error_b=abs(act_b-by) 
 
  open(unit=40, file='../output_file/analytical_numerical_error_b10.dat')
       do j=0,ny
        do k=0,nz
         i=(k)*(ny+1)+j+1
         y=node_info(i)%y_loc
         z=node_info(i)%z_loc
         write(40,*),y,z,error_b(i)
       end do
        write(40,*) " "
      end do
 close(40)
 
  if(allocated(node_info)) deallocate(node_info)
  deallocate(by,q,act_b)
  
END PROGRAM main
