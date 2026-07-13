module linear_solvers
 
contains
  subroutine IL_solve(L_csc,y,b)
   use node_classes
   use matrix_vector_op
   implicit none
   type(node_csc), intent(in)                          ::L_csc
   real(kind=8),dimension(:),allocatable               ::denum
   real(kind=8),dimension(:),allocatable               ::temp_coeffs
   real(kind=8),dimension(:),allocatable,intent(inout) ::y
   real(kind=8),dimension(:),allocatable,intent(in)    ::b
   integer                                             ::i,j ,flop,zeros_L,zeros_U
   real                                                ::sum
   integer*4                                           ::count_rate, count_max, start_count, end_count
   real*8                                              ::start_time, end_time
   
   !Forward substitution
   !open(unit=20,file="../output_file/x_at_65971_b10_3.dat")
   !print*,"Forward substitution started"
   allocate(denum(size(b)))
   flop=0
   y=b
   denum(1)=1/L_csc%coeffs(1)
  
   call system_clock(start_count, count_rate, count_max)
   start_time = real(start_count)/real(count_rate)
   
   do i=1,size(b)-1 
    denum(i+1)=1/L_csc%coeffs(L_csc%col(i)+1) 
    
    
     y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))=y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))-L_csc%coeffs(L_csc%col(i-1)+2:L_csc%col(i))*y(i)/L_csc%coeffs(L_csc%col(i-1)+1) 
     
     
     
   end do 
   
    y=y*denum
   
   
   

   call system_clock(end_count, count_rate, count_max)
   end_time = real(end_count)/real(count_rate)
   
   
   print*, "Wall clock time for IL solve:", end_time - start_time, "seconds"
   
  !print*,"no. of operations= " ,flop
   
   deallocate(denum)
  end subroutine IL_solve
  
  subroutine ILU_solve(L_csc,U_csc,x,y,b)
   use node_classes
   use matrix_vector_op
   implicit none
   type(node_csc), intent(in)                          ::L_csc
   type(node_csc), intent(in)                          ::U_csc
   real(kind=8),dimension(:),allocatable,intent(inout) ::x
   real(kind=8),dimension(:),allocatable               ::denum
   real(kind=8),dimension(:),allocatable               ::temp_coeffs
   real(kind=8),dimension(:),allocatable,intent(inout) ::y
   real(kind=8),dimension(:),allocatable,intent(in)    ::b
   integer                                             ::i,j ,flop,zeros_L,zeros_U
   real                                                ::sum
   integer*4                                           ::count_rate, count_max, start_count, end_count
   real*8                                              ::start_time, end_time
   
   !Forward substitution
   !open(unit=20,file="../output_file/x_at_65971_b10_3.dat")
   !print*,"Forward substitution started"
   allocate(denum(size(b)))
   flop=0
   zeros_L=0
   zeros_U=0
   y=b
   denum(1)=1/L_csc%coeffs(1)
  
   call system_clock(start_count, count_rate, count_max)
   start_time = real(start_count)/real(count_rate)
   
   do i=1,size(b)-1 
    denum(i+1)=1/L_csc%coeffs(L_csc%col(i)+1) 
    
    
     y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))=y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))-L_csc%coeffs(L_csc%col(i-1)+2:L_csc%col(i))*y(i)/L_csc%coeffs(L_csc%col(i-1)+1) 
     
     
     
   end do 
   
    y=y*denum
   
   !Backward substitution
   !print*,"Backward substitution started"
   x=y
   do i=size(b),1,-1 
  
    
    
      x(U_csc%row(U_csc%col(i-1)+1:U_csc%col(i)))= x(U_csc%row(U_csc%col(i-1)+1:U_csc%col(i))) - U_csc%coeffs(U_csc%col(i-1)+1:U_csc%col(i))*x(i)
 
    
  
    
   end do 

   

   call system_clock(end_count, count_rate, count_max)
   end_time = real(end_count)/real(count_rate)
   
   
   print*, "Wall clock time for ILU solve:", end_time - start_time, "seconds"
   
  !print*,"no. of operations= " ,flop
   
   deallocate(denum)
  end subroutine ILU_solve
  subroutine extract_column(A_csc,x,col)
  !input: A in csc format, x is a array of size of A 
  !output: x contains the value of A of corresponding positions
    use node_classes 
    implicit none
    type(node_csc), intent(inout)                 :: A_csc
    real(kind=8),allocatable,intent(inout)        :: x(:)
    integer,intent(in)                            :: col
    integer                                       :: i,j,K,n
    integer                                       :: local_index,global_index
    integer                                       :: matched_row
    n=size(A_csc%col)-1
    
    allocate(x(n))
    x=0.d0
      !print*,"Extracting values of column = ",col
      do j=1,A_csc%col(col)-A_csc%col(col-1) 
    
       k=A_csc%col(col-1) + j
       x(A_csc%row(k))=A_csc%coeffs(k)

      
       
      end do
   

  end subroutine extract_column
  
  subroutine extract_LUcolumn(A_csc,L_csc,U_csc)
  !input: A in csc format, L and U in csc format with symbolic LU 
  !       factorisation done        
  !output: L and U contains the value of A of corresponding positions
    use node_classes 
    implicit none
    type(node_csc), intent(inout)                 :: A_csc
    type(node_csc), intent(inout)                 :: U_csc
    type(node_csc), intent(inout)                 :: L_csc
    integer                                       :: i,j,K,n
    integer                                       :: local_index,global_index
    integer                                       :: matched_row
    n=size(L_csc%col)-1
    
    
    !open(unit=10,file="../output_file/Extracted_L_10_3.dat")
    !open(unit=20,file="../output_file/Extracted_U_10_3.dat")
    do k=1,n
      !print*,"Extracting values for LU column = ",k
      do j=1,A_csc%col(k)-A_csc%col(k-1) 
       local_index=findloc(L_csc%row(L_csc%col(k-1)+1:L_csc%col(k)), A_csc%row(A_csc%col(k-1)+j),dim=1)  
       if(local_index>0)  then
        global_index=L_csc%col(k-1) + local_index
        L_csc%coeffs(global_index)=A_csc%coeffs(A_csc%col(k-1)+j)
       end if

       local_index=findloc(U_csc%row(U_csc%col(k-1)+1:U_csc%col(k)), A_csc%row(A_csc%col(k-1)+j),dim=1)
       if(local_index>0)  then
        global_index=U_csc%col(k-1) + local_index
        U_csc%coeffs(global_index)=A_csc%coeffs(A_csc%col(k-1)+j)
       end if
      end do
      !write(10,*),"col=",k,"values ",L_csc%coeffs(L_csc%col(k-1)+1:L_csc%col(k))
      !write(20,*),"col=",k,"values ",U_csc%coeffs(U_csc%col(k-1)+1:U_csc%col(k))
    end do
   !close(10)
   !close(20)
  end subroutine extract_LUcolumn
  
  subroutine ILU_preconditioner(A_csc,L_csc,U_csc,b,q,ny,nz,n)
    use node_classes
    use matrix_vector_op
    use array_processes
    use symbolic_factorisation_op
    implicit none
    type(node_csc), intent(inout)                 :: A_csc
    type(node_csc), intent(inout)                 :: L_csc
    type(node_csc), intent(inout)                 :: U_csc
    real(kind=8),allocatable,intent(inout)        :: b(:)
    real(kind=8),allocatable,intent(in)           :: q(:)
    real(kind=8),allocatable                      :: soln_k(:)
    real(kind=8),allocatable                      :: denom(:)
    integer,allocatable                           :: nz_row(:)
    integer,intent(in)                            :: ny,nz,n
    logical, dimension(:), allocatable            :: mask
    integer                                       :: i,j,k,s,t,w,extra_size,c
    integer                                       :: global_index, local_index
    integer*4                                     :: count_rate, count_max, start_count, end_count
    real*8                                        :: start_time, end_time
    
   
    
    
    
    
    open(unit=50,file="../output_file/IL_10_4.dat")
    open(unit=60,file="../output_file/IU_10_4.dat")
    !open(unit=70,file="../output_file/soln_k.dat")
    
    !print*,"Entering Symbolic LU"
    call symbolic_ILU(A_csc,L_csc,U_csc,10)
    allocate(L_csc%coeffs(L_csc%nz_size),U_csc%coeffs(U_csc%nz_size))
    L_csc%coeffs=0.d0
    U_csc%coeffs=0.d0
   
    
    call extract_LUcolumn(A_csc,L_csc,U_csc)
    
    do k=2,size(A_csc%col)-1
       
      ! solve Lx=b
     
      allocate(denom(1:k-1))
      call extract_column(A_csc,soln_k,k)  
      
   
      do i=1,U_csc%col(k)-U_csc%col(k-1)
       s=U_csc%row(U_csc%col(k-1)+i)
       soln_k(s)=soln_k(s)/L_csc%coeffs(L_csc%col(s-1)+1)
       do j=2,L_csc%col(s)-L_csc%col(s-1) 
        soln_k(L_csc%row(L_csc%col(s-1)+j))=soln_k(L_csc%row(L_csc%col(s-1)+j))-L_csc%coeffs(L_csc%col(s-1)+j)*soln_k(s)
     
       end do
      end do 
      
      
      
      
      
      !solve U of column k
      U_csc%coeffs(U_csc%col(k-1)+1:U_csc%col(k))=soln_k(U_csc%row(U_csc%col(k-1)+1:U_csc%col(k)))
      write(60,'(*(F15.6))'),U_csc%coeffs(U_csc%col(k-1)+1:U_csc%col(k))
      !solve L of column k
      L_csc%coeffs(L_csc%col(k-1)+1:L_csc%col(k))=soln_k(L_csc%row(L_csc%col(k-1)+1:L_csc%col(k)))  
      write(50,'(*(F15.6))'),L_csc%coeffs(L_csc%col(k-1)+1:L_csc%col(k)) 
      
   
      !print*,"Numerical LU done for col=",k
      
      if(allocated(soln_k)) deallocate(soln_k)
      if(allocated(denom)) deallocate(denom)
    end do
    
    close(50)
    close(60)
    !close(70)
    
    print*,"ILU factorisation completed"
      
      
      
      !call LU_solve(L_csc,U_csc,b,q)
      
      
    
     
      
      !write solution to a file
      !open(unit=70,file="../output_file/LU_soln_b10_5.dat")
      !write(70,*),"y          z        by"
      !do i=1,size(b)
       !write(70,"(3F10.6)"),A_csc%y_loc(i),A_csc%z_loc(i),b(i)
      !end do
      !close(70)
   
    
  end subroutine ILU_preconditioner
  !#########################################################################
  subroutine jacobi_preconditioner(A,J,ny,nz,n)
   use node_classes
   use matrix_vector_op
   implicit none
   type(node_csr),intent(inout)                       :: A
   type(node_csr),intent(out)                         :: J
   integer,intent(in)                                 :: ny,nz,n
   integer                                            :: i,k,t
   logical, dimension(:), allocatable                 :: mask
   allocate(J%row(0:n),J%col(1:n),J%coeffs(1:n))
   
   J%nz_size=0
   J%row(J%nz_size)=0
   do i=1,n
    do k=A%row(i-1)+1,A%row(i)
      if(A%col(k)==i) then
        J%nz_size=J%nz_size+1
        J%coeffs(J%nz_size)=1/A%coeffs(k) 
        J%col(J%nz_size)=i 
      end if
    end do
     J%row(i)=J%nz_size
   end do             
  end subroutine jacobi_preconditioner
  !########################################################################

  subroutine read_LU(L_csc,U_csc,col_size,read_file)
   use node_classes
   implicit none
   type(node_csc), intent(inout)                       ::L_csc,U_csc
   integer,dimension(:),allocatable                    ::col_ptr,row_indices
   real*8,dimension(:),allocatable                     ::values
   integer,intent(in)                                  ::col_size 
   integer                                             ::nz_size,error_stat
   logical,intent(out)                                 ::read_file
   logical                                             ::read_L,read_U
   
   allocate(col_ptr(0:col_size))
   allocate(L_csc%col(0:col_size))
   allocate(U_csc%col(0:col_size))
   
   namelist/Column/col_ptr
   namelist/Row_and_Values/row_indices,values
   
   read_L=.true.        
   read_U=.true.
   read_file=.false.
   
   inquire(file='../input_file/L_10_4.nml',exist=read_L)
   inquire(file='../input_file/U_10_4.nml',exist=read_U)
   if (read_L .and. read_U) then
    open(unit=10, file='../input_file/L_10_4.nml',iostat=error_stat)
     if(error_stat/=0) then
      read_L=.false.
      print*,"File cannot be opened."
      
     else 
      read(10, nml=Column)
      nz_size=col_ptr(col_size)
      allocate(row_indices(nz_size),values(nz_size))
      read(10, nml=Row_and_Values)
      L_csc%col=col_ptr
      L_csc%row=row_indices
      L_csc%coeffs=values
     end if
    close(10)
   
    if(allocated(row_indices)) deallocate(row_indices)
    if(allocated(values)) deallocate(values)
   
    open(unit=20, file='../input_file/U_10_4.nml',iostat=error_stat)
   
     if(error_stat/=0) then
      read_U=.false.
      print*,"File cannot be opened."
     else 
      read(20, nml=Column)
      nz_size=col_ptr(col_size)
      allocate(row_indices(nz_size),values(nz_size))
      read(20, nml=Row_and_Values)
      U_csc%col=col_ptr
      U_csc%row=row_indices
      U_csc%coeffs=values
     end if
    close(20)
   end if
   
   if(read_L .and. read_U) read_file=.true.
   
  end subroutine read_LU
  subroutine L_solve(L_csc,y,b)
   use node_classes
   use matrix_vector_op
   implicit none
   type(node_csc), intent(in)                          ::L_csc
   real(kind=8),dimension(:),allocatable               ::denum
   real(kind=8),dimension(:),allocatable               ::temp_coeffs
   real(kind=8),dimension(:),allocatable,intent(inout) ::y
   real(kind=8),dimension(:),allocatable,intent(in)    ::b
   integer                                             ::i,j ,flop,zeros_L,zeros_U
   real                                                ::sum
   integer*4                                           ::count_rate, count_max, start_count, end_count
   real*8                                              ::start_time, end_time
   
   !Forward substitution
   !open(unit=20,file="../output_file/x_at_65971_b10_3.dat")
   !print*,"Forward substitution started"
   allocate(denum(size(b)))
   flop=0
   y=b
   denum(1)=1/L_csc%coeffs(1)
  
   call system_clock(start_count, count_rate, count_max)
   start_time = real(start_count)/real(count_rate)
   
   do i=1,size(b)-1 
    denum(i+1)=1/L_csc%coeffs(L_csc%col(i)+1) 
    
    
     y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))=y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))-L_csc%coeffs(L_csc%col(i-1)+2:L_csc%col(i))*y(i)/L_csc%coeffs(L_csc%col(i-1)+1) 
     
     
     
   end do 
   
    y=y*denum
   
   
   

   call system_clock(end_count, count_rate, count_max)
   end_time = real(end_count)/real(count_rate)
   
   
   print*, "Wall clock time for L solve:", end_time - start_time, "seconds"
   
  !print*,"no. of operations= " ,flop
   
   deallocate(denum)
  end subroutine L_solve
  
  subroutine LU_solve(L_csc,U_csc,x,y,b)
   use node_classes
   use matrix_vector_op
   implicit none
   type(node_csc), intent(in)                          ::L_csc
   type(node_csc), intent(in)                          ::U_csc
   real(kind=8),dimension(:),allocatable,intent(inout) ::x
   real(kind=8),dimension(:),allocatable               ::denum
   real(kind=8),dimension(:),allocatable               ::temp_coeffs
   real(kind=8),dimension(:),allocatable,intent(inout) ::y
   real(kind=8),dimension(:),allocatable,intent(in)    ::b
   integer                                             ::i,j ,flop,zeros_L,zeros_U
   real                                                ::sum
   integer*4                                           ::count_rate, count_max, start_count, end_count
   real*8                                              ::start_time, end_time
   
   !Forward substitution
   !open(unit=20,file="../output_file/x_at_65971_b10_3.dat")
   !print*,"Forward substitution started"
   allocate(denum(size(b)))
   flop=0
   zeros_L=0
   zeros_U=0
   y=b
   denum(1)=1/L_csc%coeffs(1)
  
   call system_clock(start_count, count_rate, count_max)
   start_time = real(start_count)/real(count_rate)
   
   do i=1,size(b)-1 
    denum(i+1)=1/L_csc%coeffs(L_csc%col(i)+1) 
    
    
     y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))=y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))-L_csc%coeffs(L_csc%col(i-1)+2:L_csc%col(i))*y(i)/L_csc%coeffs(L_csc%col(i-1)+1) 
     
     
     
   end do 
   
    y=y*denum
   
   !Backward substitution
   !print*,"Backward substitution started"
   x=y
   do i=size(b),1,-1 
  
    
    
      x(U_csc%row(U_csc%col(i-1)+1:U_csc%col(i)))= x(U_csc%row(U_csc%col(i-1)+1:U_csc%col(i))) - U_csc%coeffs(U_csc%col(i-1)+1:U_csc%col(i))*x(i)
 
    
  
    
   end do 

   

   call system_clock(end_count, count_rate, count_max)
   end_time = real(end_count)/real(count_rate)
   
   
   print*, "Wall clock time for LU solve:", end_time - start_time, "seconds"
   
  !print*,"no. of operations= " ,flop
   
   deallocate(denum)
  end subroutine LU_solve
  
  
  

  subroutine LU_decompositon(A_csc,L_csc,U_csc)
    use node_classes
    use array_processes
    use symbolic_factorisation_op
    use matrix_vector_op
    implicit none
    type(node_csc), intent(inout)                 :: A_csc
    type(node_csc), intent(inout)                 :: L_csc
    type(node_csc), intent(inout)                 :: U_csc
    real(kind=8),allocatable                      :: soln_k(:)
    real(kind=8),allocatable                      :: denom(:)
    integer,allocatable                           :: nz_row(:)
    logical, dimension(:), allocatable            :: mask
    integer                                       :: i,j,k,s,t,w,extra_size,c
    integer                                       :: global_index, local_index
    integer*4                                     :: count_rate, count_max, start_count, end_count
    real*8                                        :: start_time, end_time
    
   
    
    
    
    
    open(unit=50,file="../input_file/L_10_4.nml")
    open(unit=60,file="../input_file/U_10_4.nml")
    !open(unit=70,file="../output_file/soln_k.dat")
    
    !print*,"Entering Symbolic LU"
    call symbolic_LU(A_csc,L_csc,U_csc)
    write(50,*) "&Column "
    write(50, '(A,"",*(I0, :, ","))') "col_ptr=", L_csc%col
    write(50,*) "/"
    write(50,*) ""
    write(50,*) "&Row_and_Values "
    write(50, '(A,"",*(I0, :, ","))') "row_indices=",L_csc%row
    write(60,*) "&Column "
    write(60, '(A,"",*(I0, :, ","))') "col_ptr=",U_csc%col
    write(60,*) "/ "
    write(60,*) ""
    write(60,*) "&Row_and_Values "
    write(60, '(A,"",*(I0, :, ","))') "row_indices=", U_csc%row
    allocate(L_csc%coeffs(L_csc%nz_size),U_csc%coeffs(U_csc%nz_size))
    L_csc%coeffs=0.d0
    U_csc%coeffs=0.d0
    
    
    call extract_LUcolumn(A_csc,L_csc,U_csc)
    
    do k=2,size(A_csc%col)-1
       
      ! solve Lx=b
     
      
      call extract_column(A_csc,soln_k,k)  
      
   
      do i=1,U_csc%col(k)-U_csc%col(k-1)
       s=U_csc%row(U_csc%col(k-1)+i)
       soln_k(s)=soln_k(s)/L_csc%coeffs(L_csc%col(s-1)+1)
       do j=2,L_csc%col(s)-L_csc%col(s-1) 
        soln_k(L_csc%row(L_csc%col(s-1)+j))=soln_k(L_csc%row(L_csc%col(s-1)+j))-L_csc%coeffs(L_csc%col(s-1)+j)*soln_k(s)
     
       end do
      end do 
      
      
      
      
      
      !solve U of column k
      U_csc%coeffs(U_csc%col(k-1)+1:U_csc%col(k))=soln_k(U_csc%row(U_csc%col(k-1)+1:U_csc%col(k)))
      !write(60,'(*(F15.6))'),U_csc%coeffs(U_csc%col(k-1)+1:U_csc%col(k))
      !solve L of column k
      L_csc%coeffs(L_csc%col(k-1)+1:L_csc%col(k))=soln_k(L_csc%row(L_csc%col(k-1)+1:L_csc%col(k)))  
      !write(50,'(*(F15.6))'),L_csc%coeffs(L_csc%col(k-1)+1:L_csc%col(k)) 
      
   
      print*,"Numerical LU done for col=",k
      
      if(allocated(soln_k)) deallocate(soln_k)
      !if(allocated(soln_k)) deallocate(nz_row)
      if(allocated(denom)) deallocate(denom)
    end do
    
    write(50, '(A,"",*(G0, :, ","))') "values=", L_csc%coeffs
    write(60, '(A,"",*(G0, :, ","))') "values=", U_csc%coeffs
    write(50,*) "/"
    write(60,*) "/"
    close(50)
    close(60)
    !close(70)
    
    print*,"LU factorisation completed"
      
      
      
   !call LU_solve(L_csc,U_csc,b,q)
      
      
    
     
      
      
  end subroutine LU_decompositon
  !########################################################################
subroutine biCGSTAB(A_csr,b,q,tol,ny,nz,PRECOND)
  use node_classes
  use matrix_vector_op
  implicit none
  type(node_csr), intent(inout)          :: A_csr
  type(node_csc)                         :: A_csc
  type(node_csc)                         :: K1
  type(node_csc)                         :: K2
  real(kind=8),allocatable,intent(in)    :: q(:)
  real(kind=8),allocatable,intent(inout) :: b(:) 
  real*8,intent(in)                      :: tol
  integer,intent(in)                     :: ny,nz
  character,intent(in)                   :: PRECOND*10
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
  real(kind=8),allocatable               :: d_arr(:),d_arr1(:) !dummy array
  real(kind=8),allocatable               :: t(:)
  real(kind=8),allocatable               :: alpha_poly(:)
  real(kind=8),allocatable               :: omega_s(:)
  real(kind=8),allocatable               :: omega_t(:)
  real(kind=8),allocatable               :: omega_v(:)
  real(kind=8)                           :: row,alpha,omega,v_dot_res,t_dot_s,t_dot_t,beta,row_new,inf_n,inf_rhs,norm_2,norm_2_rhs
  integer                                :: n,i,j,k
  integer*4                              :: count_rate, count_max, start_count, end_count
  real*8                                 :: start_time, end_time
  
  n=size(q)
  call get_inf_norm(q,inf_rhs,j)
  call get_2_norm(q,norm_2_rhs)
  open(unit=20,file='../output_file/ini_b_q_out_b30.dat')
   write(20,*) ,'iter        b          q'
   
     do  i=1,n
     write(20,*) ,i,b(i),q(i)
     end do
  
  close(20)
  !Allocating memories
  allocate(poly(n),s(n),t(n),k_inv_poly(n),k_inv_s(n),k1_inv_poly(n),k1_inv_s(n),k1_inv_t(n))
  
  !Initial operation
  call mat_vec_mult(A_csr,b,d_arr)       !d_arr=Ab
  call vec_vec_sum(q,d_arr,res_ini,-1)   !res_o=q-Ab
  call vec_vec_mult(res_ini,res_ini,row,"s") !row=res_o.res_o
  poly=res_ini
  res =res_ini
  open(unit=30,file='../output_file/ini_res_b13.dat')
   write(30,*) ,'iter        res'
   do  i=1,n
     write(30,*) ,i,res(i)
   end do
  close(30)
  print*,'BiCGSTAB started'
  !Preconditioning
  if(PRECOND=="CLU")then
     call csr_to_csc(A_csr,A_csc,n)
     call LU_decompositon(A_csc,K1,K2)
  else if(PRECOND=="ILU")then
     call csr_to_csc(A_csr,A_csc,n)
     call ILU_preconditioner(A_csc,K1,K2,b,q,ny,nz,n)
    
  else if(PRECOND=="JACOBI")then
    ! call jacobi_preconditioner(A_csr,K1,ny,nz,n)  
    
 
  end if
   call system_clock(start_count, count_rate, count_max)
     start_time = real(start_count)/real(count_rate)
     
  !Iteration start
  open(unit=10,file='../output_file/bicgstab_out_b34.dat')
  do i=1,200
  
   if(PRECOND=="CLU")then
   
  
     call LU_solve(K1,K2,k_inv_poly,k1_inv_poly,poly)    ! k_inv_poly= (K^-1)p
     call mat_vec_mult(A_csr,k_inv_poly,v)               ! v=A k_inv_poly
   else if(PRECOND=="ILU")then
   
  
     call ILU_solve(K1,K2,k_inv_poly,k1_inv_poly,poly)    ! k_inv_poly= (K^-1)p
     call mat_vec_mult(A_csr,k_inv_poly,v)                !v=A k_inv_poly  
   else if(PRECOND=="JACOBI")then 
    !call mat_vec_mult(K1,poly,k_inv_poly)               ! k_inv_poly= (K^-1)p
    !call mat_vec_mult(A_csr,k_inv_poly,v)               !v=A k_inv_poly
   else
     call mat_vec_mult(A_csr,poly,v)                     !v=Ap
   end if
   call vec_vec_mult(v,res_ini,v_dot_res,"d")  !v.res_o
   
   alpha=row/v_dot_res 
   
   call scale_vec_mult(v,alpha,d_arr)                    ! alpha * v
   call vec_vec_sum(res,d_arr,s,-1)                      ! s=r-alpha * v
   
   if(PRECOND=="CLU")then
    call LU_solve(K1,K2,k_inv_s,k1_inv_s,s)              ! k_inv_s= (K^-1)s
    call mat_vec_mult(A_csr,k_inv_s,t)                   ! t=A k_inv_s
    call L_solve(K1,k1_inv_t,t)                          !K1_inv t
    call vec_vec_mult(k1_inv_t,k1_inv_s,t_dot_s,"d")         !(K1^-1)t.(K1^-1)s
    call vec_vec_mult(k1_inv_t,k1_inv_t,t_dot_t,"s")         !(K1^-1)t.(K1^-1)t
    omega= t_dot_s/t_dot_t 
    call scale_vec_mult(k_inv_poly,alpha,alpha_poly)     !alpha * (K^-1)p
    call scale_vec_mult(k_inv_s,omega,omega_s)           !omega * (K^-1)s
   else if(PRECOND=="ILU")then
    call ILU_solve(K1,K2,k_inv_s,k1_inv_s,s)              ! k_inv_s= (K^-1)s
    call mat_vec_mult(A_csr,k_inv_s,t)                   ! t=A k_inv_s
 
    
    call IL_solve(K1,k1_inv_t,t)                          !K1_inv t
    
    call vec_vec_mult(k1_inv_t,k1_inv_s,t_dot_s,"d")         !(K1^-1)t.(K1^-1)s
    call vec_vec_mult(k1_inv_t,k1_inv_t,t_dot_t,"s")         !(K1^-1)t.(K1^-1)t
    omega= t_dot_s/t_dot_t 
    call scale_vec_mult(k_inv_poly,alpha,alpha_poly)     !alpha * (K^-1)p
    call scale_vec_mult(k_inv_s,omega,omega_s)           !omega * (K^-1)s
   else if(PRECOND=="JACOBI")then 
    !call mat_vec_mult(K1,s,k_inv_s)                  ! k_inv_s= (K^-1)s 
    !call mat_vec_mult(A_csr,k_inv_s,t)                   ! t=A k_inv_s
    !call mat_vec_mult(K1,s,k1_inv_s)     !K_inv s
    !call mat_vec_mult(K1,t,k1_inv_t)     !K_inv t
    !call vec_vec_mult(k1_inv_t,k1_inv_s,t_dot_s,"d")          !(K^-1)t.(K^-1)s
    !call vec_vec_mult(k1_inv_t,k1_inv_t,t_dot_t,"s")          !(K^-1)t.(K^-1)t
    omega= t_dot_s/t_dot_t 
    !call scale_vec_mult(k_inv_poly,alpha,alpha_poly)   !alpha * (K^-1)p
    !call scale_vec_mult(k_inv_s,omega,omega_s)         !omega * (K^-1)s
   else
    call mat_vec_mult(A_csr,s,t)                          !t=As
    call vec_vec_mult(t,s,t_dot_s,"d")                        !t.s
    call vec_vec_mult(t,t,t_dot_t,"s")                        !t.t
    omega= t_dot_s/t_dot_t 
    call scale_vec_mult(poly,alpha,alpha_poly)       !alpha * poly
    call scale_vec_mult(s,omega,omega_s)             !omega*s
   
   end if
              
   
   
   call vec_vec_sum(alpha_poly,omega_s,d_arr,1) !alpha * poly + omega*s=d_arr
   call vec_vec_sum(b,d_arr,d_arr1,1)           !d_arr1=b+d_arr
   b=d_arr1                                     !b_new=d_arr1
   call scale_vec_mult(t,omega,omega_t)         !omega*t
   call vec_vec_sum(s,omega_t,res,-1)           !res_new=s-omega*t 
   call vec_vec_mult(res,res_ini,row_new,"d")   !row_new=res_new.res_o   
   beta=alpha*row_new/(omega*row)   
   call scale_vec_mult(v,omega,omega_v)         !omega*v
   call vec_vec_sum(poly,omega_v,d_arr1,-1)      !d_arr1=p -  omega*v        
   call scale_vec_mult(d_arr1,beta,d_arr)        !d_arr=beta*(d_arr1)
   call vec_vec_sum(res,d_arr,poly,1)           !poly=res_new+d_arr
   
    
   call get_inf_norm(res,inf_n,j)
   call get_2_norm(res,norm_2)
   row=row_new                                  !row updated
   inf_n=inf_n/inf_rhs
   norm_2=norm_2/norm_2_rhs
  
   if(norm_2<tol) exit
   write(10,"(1X,I5,1X,4f10.6)") ,i,inf_n,A_csr%y_loc(j),A_csr%z_loc(j),norm_2
  
  end do
  call system_clock(end_count, count_rate, count_max)
   end_time = real(end_count)/real(count_rate)
   
  print*, "Iteration clock time:", end_time - start_time, "seconds"
  close(10)
 
  print*,"alpha=",alpha
  
  open(unit=40,file="../output_file/by_bz_res_b34.dat")
   do i=1,n
     write(40,'(1X,I6,1X,*(E14.6))'),i,A_csr%y_loc(i),A_csr%z_loc(i),b(i),res(i)
   end do
  close(40)
   
  deallocate(res,v,poly,s,t,alpha_poly,omega_s,omega_t,omega_v,d_arr)
  end subroutine biCGSTAB
    
   
  subroutine inv_mat(A_csr,A_csr_inv)
   use node_classes
   use matrix_vector_op
   implicit none
   type(node_csr),intent(inout)                ::A_csr
   type(node_csc)                              ::A_csc
   type(node_csr),intent(out)                  ::A_csr_inv
   type(node_csc)                              ::A_csc_inv
   type(node_csc)                              ::L_csc,U_csc
   integer                                     ::n,k,nz,t,chunk_size,j
   real*8,allocatable,dimension(:)             ::x,y, x_values,I
   integer,allocatable,dimension(:)            ::x_index
   real(kind=8),allocatable                    ::temp_coeffs(:)
   integer,allocatable                         ::temp_col(:)
   integer,allocatable                         ::temp_row(:)
   
   n=size(A_csr%row)-1  
   chunk_size=1000
   A_csc_inv%capacity=258 
   A_csc_inv%nz_size=0
   print*,"HERE"
   call csr_to_csc(A_csr,A_csc,n)
   print*,"before LU"
   call LU_decompositon(A_csc,L_csc,U_csc)
   print*,"LU decomposition done"
   allocate(I(n))
   allocate(A_csc_inv%col(0:n),A_csc_inv%y_loc(0:n),A_csc_inv%z_loc(0:n))
   allocate(A_csc_inv%row(A_csc_inv%capacity))
   allocate(A_csc_inv%coeffs(A_csc_inv%capacity))
   A_csc_inv%col(0)=A_csc_inv%nz_size
   I=0.d0
   do k=1,n
     I(k)=1.d0
     nz=0
     call LU_solve(L_csc,U_csc,x,y,I)
     
     print*,"k=",k
     !strip vector x of zeros and store nonzeros in x_values and indices in x_index
     do t=1,n
       if(abs(x(t))>1e-7) then
       
        !======================================================
        !increase memory if capacity reached
        if((A_csc_inv%capacity-A_csc_inv%nz_size)<258)then
          if(A_csc_inv%capacity==0)then
            allocate(A_csc_inv%coeffs(chunk_size))
            allocate(A_csc_inv%row(chunk_size))
          else
            temp_coeffs=A_csc_inv%coeffs(1:A_csc_inv%nz_size)
            temp_row   =A_csc_inv%row(1:A_csc_inv%nz_size)
            if(allocated(A_csc_inv%coeffs)) deallocate(A_csc_inv%coeffs)
            if(allocated(A_csc_inv%row))   deallocate(A_csc_inv%row)
            allocate(A_csc_inv%coeffs(A_csc_inv%capacity+chunk_size))
            allocate(A_csc_inv%row(A_csc_inv%capacity+chunk_size))
            A_csc_inv%capacity=A_csc_inv%capacity+chunk_size
            A_csc_inv%coeffs(1:A_csc_inv%nz_size)=temp_coeffs
            A_csc_inv%row(1:A_csc_inv%nz_size)   =temp_row
          
         end if
        end if
        !==================================================
        nz=nz+1  
        A_csc_inv%nz_size=A_csc_inv%nz_size + 1
        A_csc_inv%coeffs( A_csc_inv%col(k-1)+nz)=x(t)  
        A_csc_inv%row( A_csc_inv%col(k-1)+nz)=t 
       end if
     end do 
     A_csc_inv%col(k)=A_csc_inv%nz_size
     
     I(k)=0.d0
       
   end do
   call csc_to_csr(A_csc_inv,A_csr_inv,n)
   
   
   
     
     deallocate(A_csc_inv%coeffs,A_csc_inv%col,A_csc_inv%row)
     deallocate(A_csc%coeffs,A_csc%col,A_csc%row)
     
     print*,"inverse matrix created"
  end subroutine inv_mat
 
  

  
  subroutine solver(A_csr,b,q,tol,ny,nz,n,solver_name,precond)
   use node_classes
   use matrix_vector_op
   implicit none
  
   type(node_csr),intent(inout)                :: A_csr
   type(node_csc)                              :: A_csc
   type(node_csc)                              :: L_csc
   type(node_csc)                              :: U_csc
   real(kind=8),allocatable,intent(in)         :: q(:)
   real(kind=8),allocatable,intent(inout)      :: b(:)
   real(kind=8),allocatable                    :: f(:) 
   real*8,intent(in)                           :: tol
   integer,intent(in)                          :: ny,nz,n 
   character, intent(in)                       :: solver_name*10
   character, intent(in)                       :: precond*10
   if(solver_name=="BICGSTAB")then
     call biCGSTAB(A_csr,b,q,tol,ny,nz,precond)
     
   else if(solver_name=="LU_METHOD") then
      print*,"LU METHOD SELECTED"
      
      call csr_to_csc(A_csr,A_csc,n)
      call LU_decompositon(A_csc,L_csc,U_csc)
      call LU_solve(L_csc,U_csc,b,f,q)
   else  
     print*,"solver chosen is not available"
   end if
      
  end subroutine solver
  
end module linear_solvers

