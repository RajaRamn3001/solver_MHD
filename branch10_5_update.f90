!This is a working fortran file


!This subroutine read a input file using namelist
subroutine read_input(ny,nz,Ly,Lz,res_tolerance,A,Rm,linear_solver ,pre_cond,del_t,alpha_k)
 implicit none
 integer					:: ny,nz,Rm
 real 						:: Ly,Lz,res_tolerance,A,del_t,alpha_k
 character 					:: linear_solver*10 ,pre_cond*10
 
 namelist/input_parameters/ny,nz,Ly,Lz,res_tolerance,A,linear_solver ,pre_cond,Rm,del_t,alpha_k
 open(unit=10, file='../input_file/input.txt')
 read(10, nml=input_parameters)
 close(10)
end subroutine read_input

!This subroutine generates a grid
subroutine generate_grid(ny,nz,Ly,Lz,A,gridpts_y,gridpts_z)
 implicit none
 INTEGER,intent(in) 				:: ny,nz
 INTEGER 					:: i,j,k
 REAL,intent(in) 				:: Ly,Lz,A
 REAL(kind=8) 					:: dy,dz
 REAL(kind=8),intent(out) 			:: gridpts_y(0:ny),gridpts_z(0:nz)
 dy=Ly/(ny)
 dz=Lz/(nz)
 

 do i=0,ny
   gridpts_y(i)= -1.d0+ dy*(i) 
   gridpts_z(i)= -1.d0+ dz*(i)  
 end do
 !This loop stretches the grid point

  gridpts_y= tanh(A*gridpts_y)/tanh(A)
  gridpts_z= tanh(A*gridpts_z)/tanh(A) 

 open(unit=20, file='../output_file/grid_pts_stretched_b10_2.dat')
 
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

 real(kind=8)                                      :: y_loc , z_loc
 real(kind=8),allocatable                          :: coeffs(:)
 integer,allocatable                               :: col(:)
 integer             			           :: max_size ,nz_size
end type


type node_csr

 real(kind=8),allocatable                          :: y_loc (:), z_loc(:)
 real(kind=8),allocatable   			   :: coeffs(:)
 integer,allocatable        			   :: col(:)
 integer,allocatable        			   :: row(:)
 integer,allocatable        			   :: max_size(:)
 integer                    			   :: nz_size,capacity
end type

type node_csc

 real(kind=8),allocatable                          :: y_loc (:), z_loc(:)
 real(kind=8),allocatable   			   :: coeffs(:)
 integer,allocatable        			   :: col(:)
 integer,allocatable        			   :: row(:)
 integer,allocatable        			   :: max_size(:)
 integer                    			   :: nz_size,capacity
end type

contains
  subroutine csr_to_csc(csr_array,csc_array,col_num)
    type(node_csr),intent(inout)                   :: csr_array
    type(node_csc),intent(out)                     :: csc_array
    integer        ,intent(in)                     :: col_num
    integer                                        :: i,j,k
    integer                                        :: nz_size
    integer                                        :: count
    count=0
    nz_size=csr_array%row(col_num)
    allocate(csc_array%col(0:col_num),csc_array%y_loc(0:col_num),csc_array%z_loc(0:col_num))
    csc_array%y_loc(0:col_num)=csr_array%y_loc
    csc_array%z_loc(0:col_num)=csr_array%z_loc
    csc_array%nz_size=0
    csc_array%col(0)=csc_array%nz_size
    allocate(csc_array%coeffs(nz_size),csc_array%row(nz_size))
    
     do k=1,col_num
      do i=1,size(csr_array%row)-1
         do j=1,csr_array%row(i)-csr_array%row(i-1)
            if(csr_array%col(csr_array%row(i-1)+j)==k) then
               csc_array%nz_size=csc_array%nz_size+1
               count=count+1
               csc_array%coeffs(csc_array%nz_size)=csr_array%coeffs(csr_array%row(i-1)+j)
               csc_array%row(csc_array%nz_size) =i
               
            end if
         end do
       end do
        csc_array%col(k)=csc_array%nz_size
       
     end do
    
     print*,"NNZ of csc_array",count
     print*,"NNZ of csr_array",nz_size
     open(unit=50, file='../output_file/node_info_csc_b10_5_coeffs.dat')
     open(unit=60, file='../output_file/node_info_csc_b10_5_rows.dat')
      write(50,*),"ap       an       as        aw       ae"
 
      do j=1,col_num
    
  
        
        
        write(50,"(5F15.6)"),csc_array%coeffs(csc_array%col(j-1)+1:csc_array%col(j))
        write(60,*),csc_array%row(csc_array%col(j-1)+1:csc_array%col(j))
       
      end do 
  
     
      close(50)
      close(60)
  end subroutine csr_to_csc
  
  subroutine fill_node_array(node_arr,node_arr_csr,grid_pty,grid_ptz,ny,nz,n,f,alpha_k)
  
    implicit none
    integer,intent(in)                             :: ny,nz,n 
    real,intent(in)                                :: f,alpha_k
    real(kind=8),intent(in)                        :: grid_pty(0:ny),grid_ptz(0:nz)
    integer                                        :: i,j,k,count,c,chunk_size
    type(node), allocatable, intent(inout)         :: node_arr(:)
    type(node_csr),intent(inout)                   :: node_arr_csr
    real(kind=8),allocatable                       :: temp_coeffs(:)
    integer,allocatable                            :: temp_col(:)
    integer,allocatable                            :: temp_row(:)
    integer                                        :: row_ptr,col_ptr

    count=0
    chunk_size=1000
    node_arr_csr%nz_size=0
    node_arr_csr%capacity=ny+2      !ny+2=258
    
   
    row_ptr=0
    col_ptr=0
    allocate(node_arr(n))                               
    allocate(node_arr_csr%row(0:n),node_arr_csr%y_loc(0:n),node_arr_csr%z_loc(0:n))
    allocate(node_arr_csr%col(node_arr_csr%capacity))
    allocate(node_arr_csr%coeffs(node_arr_csr%capacity))
    
    
    
    node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
    node_arr_csr%y_loc(0)=0.d0
    node_arr_csr%z_loc(0)=0.d0
    print*,'size of node array',size(node_arr)
     do k=0,nz
      do j=0,ny
       
       i=k*(ny+1)+j+1
       row_ptr=row_ptr+1
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
       node_arr(i)%y_loc=grid_pty(j)
       node_arr(i)%z_loc=grid_ptz(k)
       
       node_arr_csr%y_loc(row_ptr)=grid_pty(j)
       node_arr_csr%z_loc(row_ptr)=grid_ptz(k)
       
       if((node_arr_csr%capacity-node_arr_csr%nz_size)<(ny+2))then
          if(node_arr_csr%capacity==0)then
            allocate(node_arr_csr%coeffs(chunk_size))
            allocate(node_arr_csr%col(chunk_size))
          else
            temp_coeffs=node_arr_csr%coeffs(1:node_arr_csr%nz_size)
            temp_col   =node_arr_csr%col(1:node_arr_csr%nz_size)
            if(allocated(node_arr_csr%coeffs)) deallocate(node_arr_csr%coeffs)
            if(allocated(node_arr_csr%col))   deallocate(node_arr_csr%col)
            allocate(node_arr_csr%coeffs(node_arr_csr%capacity+chunk_size))
            allocate(node_arr_csr%col(node_arr_csr%capacity+chunk_size))
            node_arr_csr%capacity=node_arr_csr%capacity+chunk_size
            node_arr_csr%coeffs(1:node_arr_csr%nz_size)=temp_coeffs
            node_arr_csr%col(1:node_arr_csr%nz_size)   =temp_col
          
        end if
       end if
       
       
       
       
       
       !count=count+1                               
       if(j>0 .and. j<ny .and. k>0 .and. k<nz ) then
          allocate(node_arr(i)%coeffs(5))
          allocate(node_arr(i)%col(5))
          
          node_arr(i)%coeffs(3)= -(alpha_k**2+f+2/(grid_pty(j+1)-grid_pty(j))/(grid_pty(j)-grid_pty(j-1))+2/(grid_ptz(k+1)-grid_ptz(k))/(grid_ptz(k)-grid_ptz(k-1)))
          node_arr(i)%coeffs(5)=2/(grid_ptz(k+1)-grid_ptz(k-1))/(grid_ptz(k+1)-grid_ptz(k))
          node_arr(i)%coeffs(1)=2/(grid_ptz(k+1)-grid_ptz(k-1))/(grid_ptz(k)-grid_ptz(k-1))
          node_arr(i)%coeffs(4)=2/(grid_pty(j+1)-grid_pty(j-1))/(grid_pty(j+1)-grid_pty(j))
          node_arr(i)%coeffs(2)=2/(grid_pty(j+1)-grid_pty(j-1))/(grid_pty(j)-grid_pty(j-1))
          
          node_arr(i)%col(3)=(k)*(ny+1)+j+1
          node_arr(i)%col(5)=(k+1)*(ny+1)+j+1
          node_arr(i)%col(1)=(k-1)*(ny+1)+j+1
          node_arr(i)%col(4)=(k)*(ny+1)+j+1+1
          node_arr(i)%col(2)=(k)*(ny+1)+j-1+1
          
          node_arr_csr%coeffs(node_arr_csr%nz_size+1)=2/(grid_ptz(k+1)-grid_ptz(k-1))/(grid_ptz(k)-grid_ptz(k-1))
          node_arr_csr%coeffs(node_arr_csr%nz_size+2)=2/(grid_pty(j+1)-grid_pty(j-1))/(grid_pty(j)-grid_pty(j-1))
          node_arr_csr%coeffs(node_arr_csr%nz_size+3)=-(alpha_k**2+f+2/(grid_pty(j+1)-grid_pty(j))/(grid_pty(j)-grid_pty(j-1))+2/(grid_ptz(k+1)-grid_ptz(k))/(grid_ptz(k)-grid_ptz(k-1)))
          node_arr_csr%coeffs(node_arr_csr%nz_size+4)=2/(grid_pty(j+1)-grid_pty(j-1))/(grid_pty(j+1)-grid_pty(j))
          node_arr_csr%coeffs(node_arr_csr%nz_size+5)=2/(grid_ptz(k+1)-grid_ptz(k-1))/(grid_ptz(k+1)-grid_ptz(k))
          
          node_arr_csr%col(node_arr_csr%nz_size+1)=(k-1)*(ny+1)+j+1
          node_arr_csr%col(node_arr_csr%nz_size+2)=(k)*(ny+1)+j-1+1
          node_arr_csr%col(node_arr_csr%nz_size+3)=(k)*(ny+1)+j+1
          node_arr_csr%col(node_arr_csr%nz_size+4)=(k)*(ny+1)+j+1+1
          node_arr_csr%col(node_arr_csr%nz_size+5)=(k+1)*(ny+1)+j+1
          
       
          
          node_arr_csr%nz_size=node_arr_csr%nz_size+5
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          
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
             
             node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-(2*c+c**2)/(grid_pty(j+1)-grid_pty(j))/(c+c**2)
             node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(1+c)**2/(grid_pty(j+1)-grid_pty(j))/(c+c**2)
             node_arr_csr%coeffs(node_arr_csr%nz_size+3)=-1.0/(grid_pty(j+1)-grid_pty(j))/(c+c**2)
             
             node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j+1
             node_arr_csr%col(node_arr_csr%nz_size+2)=(k)*(ny+1)+j+1+1
             node_arr_csr%col(node_arr_csr%nz_size+3)=(k)*(ny+1)+j+2+1
             
             node_arr_csr%nz_size=node_arr_csr%nz_size+3
             node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
             
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
            
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)=1.0/(grid_pty(j-1)-grid_pty(j))/(c+c**2)
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-(1+c)**2/(grid_pty(j-1)-grid_pty(j))/(c+c**2)
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(2*c+c**2)/(grid_pty(j-1)-grid_pty(j))/(c+c**2)
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j-2+1
            node_arr_csr%col(node_arr_csr%nz_size+2)=(k)*(ny+1)+j-1+1
            node_arr_csr%col(node_arr_csr%nz_size+3)=(k)*(ny+1)+j+1
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+3
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
            
            
        else
          allocate(node_arr(i)%coeffs(1))
          allocate(node_arr(i)%col(1))
          node_arr(i)%coeffs(1)= 1.0
          node_arr(i)%col(1)=(k)*(ny+1)+j+1
          
          node_arr_csr%coeffs(node_arr_csr%nz_size+1) =1.0
          node_arr_csr%col(node_arr_csr%nz_size+1)    =(k)*(ny+1)+j+1
          
          node_arr_csr%nz_size=node_arr_csr%nz_size+1
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          
        end if
       
                   
       end do
     end do
     
     open(unit=30, file='../output_file/node_info_b10_5_coeffs.dat')
     open(unit=40, file='../output_file/node_info_b10_5_cols.dat')
     open(unit=50, file='../output_file/node_info_csr_b10_5_coeffs.dat')
     open(unit=60, file='../output_file/node_info_csr_b10_5_cols.dat')
      write(30,*),"ap       an       as        aw       ae"
 
      do j=1,n
    
  
    
        write(30,"(5F15.6)"),node_arr(j)%coeffs!,node_arr(j)%y_loc,node_arr(j)%z_loc
        write(40,*),node_arr(j)%col
        
        
        write(50,"(5F15.6)"),node_arr_csr%coeffs(node_arr_csr%row(j-1)+1:node_arr_csr%row(j))
        write(60,*),node_arr_csr%col(node_arr_csr%row(j-1)+1:node_arr_csr%row(j))
       
      end do 
  
      close(30)
      close(40)
      close(50)
      close(60)
      print*,"Node array filled" 
    
    end subroutine fill_node_array
    
    subroutine fill_q_array(node_arr_csr,q_array,f,alpha_k)
       type(node_csr),intent(inout)                     :: node_arr_csr
       real(kind=8),allocatable, intent(out)            :: q_array(:)
       real, intent(in)                                 ::f,alpha_k
       real,parameter                                   :: PI=4*atan(1.d0)
       integer                                          :: i
       real                                             :: y,z,constant
       
       constant=0.01*((f+alpha_k**2)+8*(PI)**2)
       allocate(q_array(1:size(node_arr_csr%row)-1))
       open(unit=30, file='../output_file/rhs_b10.dat')
       do i=1,size(node_arr_csr%row)-1
         y=node_arr_csr%y_loc(i)
         z=node_arr_csr%z_loc(i)
         q_array(i)=-constant*cos(2*PI*y)*sin(2*PI*z)
         if(abs(z)>0.9999999 .or. abs(y)>0.9999999)then
          q_array(i)=0.d0
         end if
         write(30,*),y,z,q_array(i)
       end do
       close(30)
       
        
    end subroutine fill_q_array
     

end module node_classes

module array_processes
contains
  SUBROUTINE QUICKSORT(ARRAY, N)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                            :: N
    INTEGER, INTENT(INOUT)                         :: ARRAY(N)
    
    ! Stack for storing partition boundaries
    INTEGER, PARAMETER :: MAXSTACK = 100
    INTEGER :: STACK(MAXSTACK, 2)
    INTEGER :: TOP, LOW, HIGH, PIVOT
    INTEGER :: I, J, TEMP
    
    ! Initialize stack
    TOP = 1
    STACK(1, 1) = 1
    STACK(1, 2) = N
    
    ! Main sorting loop
    DO WHILE (TOP > 0)
        ! Pop from stack
        LOW = STACK(TOP, 1)
        HIGH = STACK(TOP, 2)
        TOP = TOP - 1
        
        ! Only sort if more than one element
        IF (LOW < HIGH) THEN
            ! Partition the array
            CALL PARTITION(ARRAY, LOW, HIGH, PIVOT)
            
            ! Push left partition onto stack
            IF (PIVOT - 1 > LOW) THEN
                TOP = TOP + 1
                STACK(TOP, 1) = LOW
                STACK(TOP, 2) = PIVOT - 1
            END IF
            
            ! Push right partition onto stack
            IF (PIVOT + 1 < HIGH) THEN
                TOP = TOP + 1
                STACK(TOP, 1) = PIVOT + 1
                STACK(TOP, 2) = HIGH
            END IF
        END IF
    END DO
 END SUBROUTINE QUICKSORT

 SUBROUTINE PARTITION(ARRAY, LOW, HIGH, PIVOT_INDEX)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT)                             :: ARRAY(*)
    INTEGER, INTENT(IN)                                :: LOW, HIGH
    INTEGER, INTENT(OUT)                               :: PIVOT_INDEX
    
    INTEGER :: PIVOT_VALUE, I, J, TEMP
    
    ! Choose last element as pivot
    PIVOT_VALUE = ARRAY(HIGH)
    I = LOW - 1
    
    ! Partition loop
    DO J = LOW, HIGH - 1
        IF (ARRAY(J) <= PIVOT_VALUE) THEN
            I = I + 1
            ! Swap elements
            TEMP = ARRAY(I)
            ARRAY(I) = ARRAY(J)
            ARRAY(J) = TEMP
        END IF
    END DO
    
    ! Place pivot in correct position
    TEMP = ARRAY(I + 1)
    ARRAY(I + 1) = ARRAY(HIGH)
    ARRAY(HIGH) = TEMP
    
    PIVOT_INDEX = I + 1
 END SUBROUTINE PARTITION

  subroutine append_to_array(array, new_element)
    integer, allocatable, intent(inout) :: array(:)
    integer, intent(in) :: new_element
    integer, allocatable :: temp_array(:)
    integer :: n
    
    if (.not. allocated(array)) then
        allocate(array(1))
        array(1) = new_element
    else
        n = size(array)
        allocate(temp_array(n+1))
        temp_array(1:n) = array(1:n)
        temp_array(n+1) = new_element
        call move_alloc(temp_array, array)
    endif
  end subroutine
end module array_processes



module symbolic_factorisation_op

type graph                  !col--->row 
   integer,allocatable                            :: col(:),row(:)
end type

type SCC                  !row--->col 
   integer,allocatable                            :: vertices
   integer,allocatable                            :: group
end type

contains
 
  
  subroutine dfs_reach(A_csc,L_csc,dependency_gr,filled_col,column)
   use node_classes
   use array_processes
   type(graph),intent(inout)                      :: A_csc
   type(graph),intent(in)                         :: dependency_gr
   type(node_csc),intent(in)                      :: L_csc
   integer,allocatable,intent(out)                :: filled_col(:)
   integer,    intent(in)                         :: column
   integer                                        :: i,j,k,n
   logical, allocatable                           :: visited(:)
   integer                                        :: old_size, extra_size
   integer, allocatable                           :: stack(:)
   integer                                        :: stack_top
   integer                                        :: current_vertex, neighbor, start_idx, end_idx
   integer                                        :: count
   
   n=size(dependency_gr%col) - 1
 
   allocate(visited(n))
   allocate(stack(n))
   visited(1:n) = .false.
   stack_top=0 
   count=0
   
   !print*,"dfs reach for column =",column 
   
   filled_col=A_csc%row(A_csc%col(column-1)+1:A_csc%col(column))

     stack_top=stack_top+ size(dependency_gr%row(dependency_gr%col(column-1)+1:dependency_gr%col(column)))
     stack(1:stack_top)=dependency_gr%row(dependency_gr%col(column):dependency_gr%col(column-1)+1:-1)
    
     
      do while (stack_top > 0)
        current_vertex=stack( stack_top )
        stack_top=stack_top - 1
        
        if( .not. visited(current_vertex)) then
         ! mark as visited and add to result
          visited(current_vertex)= .true. 
          
          !get neighbours and push unvisited ones to stack
          
          start_idx=L_csc%col(current_vertex-1)+1
          end_idx  =L_csc%col(current_vertex)
          
          ! Push in reverse order to maintain left-to-right traversal
          
          do i=end_idx,start_idx,-1
            neighbor= L_csc%row(i)
            if(neighbor > current_vertex .and. .not. any(filled_col== neighbor) )then
               if( .not. visited(neighbor)) then
                 if(neighbor<column)then
                   stack_top = stack_top + 1
                   stack(stack_top)=neighbor
                 end if
                 
                 call append_to_array(filled_col, neighbor) 
               end if
            end if
          end do
         
        end if
      end do
      
      call quicksort(filled_col,size(filled_col))
      deallocate(stack,visited)         
  end subroutine dfs_reach

  subroutine get_dependency(A_csc,dependency_gr)
   use array_processes
   implicit none
   type(graph),intent(inout)                      :: A_csc
   type(graph),intent(out)                        :: dependency_gr
   integer                                        :: i,j,k
   integer, allocatable                           :: temp(:)
   logical, allocatable                           :: visited(:)
   integer, allocatable                           :: stack(:)
   integer, allocatable                           :: dfs_order(:)
   integer                                        :: stack_top ,stack_size
   integer                                        :: old_size, new_size
   integer                                        :: current_vertex, neighbor, start_idx, end_idx
   integer                                        :: count
   integer, allocatable                           :: temp_row(:)
   
   allocate(dependency_gr%col(0:size(A_csc%col)-1))
   dependency_gr%col=0
      
   
   do k= 2, size(A_csc%col)-1
      
     temp=pack(A_csc%row(A_csc%col(k-1)+1:A_csc%col(k)),A_csc%row(A_csc%col(k-1)+1:A_csc%col(k))<k)
     stack_top=0 
     count=0
    
     
     if(size(temp) > 0)  then
     
      stack_size=k-temp(1)+2  !extra 2 for safety
      allocate(stack(stack_size))
      stack(:)=0
      allocate(visited(k-1))
      visited(:) = .false.
      
      stack_top=stack_top + size(temp)
      stack(1:size(temp))=temp(size(temp):1:-1)  !last element in temporaray array
      
        
      
       do while (stack_top > 0)
       
       
       
         
      
        current_vertex=stack( stack_top )
          
        stack_top= stack_top - 1
        
        
        
         if( .not. visited(current_vertex)) then
         ! mark as visited and add to result
          visited(current_vertex)= .true. 
          
          !get neighbours and push unvisited ones to stack
          
          start_idx=A_csc%col(current_vertex-1)+1
          end_idx  =A_csc%col(current_vertex)
          
         
          
          ! Push in reverse order to maintain left-to-right traversal
            do i=end_idx,start_idx,-1
            
               neighbor= A_csc%row(i)
            
                if( (neighbor < k .and. neighbor > current_vertex ).and. .not. any(temp == neighbor) )then
                  
                   if( .not. visited(neighbor)) then
                     if (stack_top < stack_size) then  ! Check stack space
                          stack_top = stack_top + 1
                          stack(stack_top) = neighbor
                         ! print*,"col=",k, "Before append: temp size=", size(temp), "adding neighbor=", neighbor

                          call append_to_array(temp, neighbor)

                          !print*, "After append: temp size=", size(temp)
                         ! if (size(temp) > 0) print*, "Last element:", temp(size(temp))

                      endif
                  
                   end if
                end if
              end do
              
            end if
           
           
          end do
          
         
         
         if( allocated(stack) )   deallocate(stack)
         if( allocated(visited) ) deallocate(visited)      
       
     end if
    
     call quicksort(temp,size(temp))
     
     old_size = size(dependency_gr%row)
     new_size = old_size + size(temp)
     allocate(temp_row(new_size))
     
     temp_row(1:old_size) = dependency_gr%row
     temp_row(old_size+1:new_size) = temp
     
     if(allocated(dependency_gr%row)) deallocate(dependency_gr%row)
     allocate(dependency_gr%row(new_size))
     
     dependency_gr%row = temp_row
     deallocate(temp_row)
     dependency_gr%col(k)=dependency_gr%col(k-1)+size(temp)
     !print*, "dependency graph Completed col=", k
    
     if( allocated(temp) )    deallocate(temp)
     
     
    
      
   end do 
   
   print*,"dependency graph done"
   
  end subroutine get_dependency

  subroutine symbolic_LU(A_csc,L,U)
   use node_classes
   use array_processes
   implicit none
   type(node_csc),intent(inout)                   :: A_csc
   type(node_csc),intent(out)                     :: L
   type(node_csc),intent(out)                     :: U
   type(graph)                                    :: A_adjacency
   type(graph)                                    :: A_dependency
   integer,allocatable                            :: temp1(:)
   integer,allocatable                            :: filled_array(:)
   integer                                        :: k,j,i,s,n,min
   integer                                        :: old_size,new_size
   integer,allocatable                            :: temp_row(:)
   integer                                        :: row_indices,col_indices
   integer                                        :: chunk_size
   
   
    n= size(A_csc%col)-1
   
    allocate(L%col(0:n))
    allocate(U%col(0:n))
    allocate(A_adjacency%col(0:n))
    
    L%col(0)=0
    U%col(0)=0
    
    A_adjacency%row=A_csc%row
    A_adjacency%col=A_csc%col
    
    call get_dependency(A_adjacency,A_dependency)
    !insert existing A[1:n,i] to L[:,i]
    do i=1,n
      
      
    
      
      call dfs_reach( A_adjacency,L, A_dependency,filled_array,i)
       
     
    
      !insert existing A[1:n,i] to L[:,i]
      temp1=pack(filled_array,filled_array>=i)
     
      old_size = size(L%row)
      new_size = old_size + size(temp1)
      allocate(temp_row(new_size))
     
      temp_row(1:old_size) = L%row
      temp_row(old_size+1:new_size) = temp1
     
      if(allocated(L%row)) deallocate(L%row)
      allocate(L%row(new_size))
     
      L%row = temp_row
      deallocate(temp_row)
      
      L%col(i)=size(temp1) + L%col(i-1) 
    
    
     
   
      if(allocated(temp1)) deallocate(temp1)
      !insert existing A[i,1:n to U[i,:]
     
       temp1=pack(filled_array,filled_array<i)
       old_size = size(U%row)
       new_size = old_size + size(temp1)
       allocate(temp_row(new_size))
      
       temp_row(1:old_size) = U%row
       temp_row(old_size+1:new_size) = temp1
     
       if(allocated(U%row)) deallocate(U%row)
       allocate(U%row(new_size))
     
       U%row = temp_row
       deallocate(temp_row)
   
       U%col(i)=size(temp1)+ U%col(i-1)

      
       if(allocated(temp1)) deallocate(temp1)
     !print*,"Symbolic LU done for col=",i
   end do
    L%nz_size=size(L%row)
    U%nz_size=size(U%row)
   
    ! deallocate(temp1,A_dependency%col,A_dependency%row,A_adjacency%col,A_adjacency%row)
  print*,"Symbolic LU done"
  print*,"NNZ in L",L%nz_size
  print*,"NNZ in U",U%nz_size
  !open(unit=10,file="../output_file/U_10_3_row.dat")
  !   do i=1,size(A_csc%col)-1
   
  !    write(10,'(*(I6))'),U%row(U%col(i-1)+1:U%col(i)) 
  ! end do
  !close(10)
  end subroutine symbolic_LU
  
  
  subroutine dfs_reach_incomplete(A_csc,L_csc,dependency_gr,filled_col,column,reach_limit)
   use node_classes
   use array_processes
   type(graph),intent(inout)                      :: A_csc
   type(graph),intent(in)                         :: dependency_gr
   type(node_csc),intent(in)                      :: L_csc
   integer,allocatable,intent(out)                :: filled_col(:)
   integer,    intent(in)                         :: column,reach_limit
   integer                                        :: i,j,k,n
   logical, allocatable                           :: visited(:)
   integer                                        :: old_size, extra_size
   integer, allocatable                           :: stack(:)
   integer                                        :: stack_top
   integer                                        :: current_vertex, neighbor, start_idx, end_idx
   integer                                        :: count
   
   n=size(dependency_gr%col) - 1
 
   allocate(visited(n))
   allocate(stack(n))
   visited(1:n) = .false.
   stack_top=0 
   count=0
   
   print*,"dfs reach for column =",column 
   
   filled_col=A_csc%row(A_csc%col(column-1)+1:A_csc%col(column))

     stack_top=stack_top+ size(dependency_gr%row(dependency_gr%col(column-1)+1:dependency_gr%col(column)))
     stack(1:stack_top)=dependency_gr%row(dependency_gr%col(column-1)+1:dependency_gr%col(column))
    
     
  outer_loop:   do while (stack_top > 0 ) 
        current_vertex=stack( stack_top )
        stack_top=stack_top - 1
        
        if(count == reach_limit) exit
        
       
        
        
        if( .not. visited(current_vertex)) then
         ! mark as visited and add to result
          visited(current_vertex)= .true. 
          
          !get neighbours and push unvisited ones to stack
          
          start_idx=L_csc%col(current_vertex-1)+1
          end_idx  =L_csc%col(current_vertex)
          
          ! Push in reverse order to maintain left-to-right traversal
          
          do i=start_idx,end_idx
            neighbor= L_csc%row(i)
            if(neighbor > current_vertex .and. .not. any(filled_col== neighbor) )then
               if( .not. visited(neighbor)) then
                 if(neighbor< column)then
                   stack_top = stack_top + 1
                   stack(stack_top)=neighbor
                 end if
                 
                 call append_to_array(filled_col, neighbor)
                 count=count + 1 
                 if(count > reach_limit) exit outer_loop
               end if
            end if
          end do
         
        end if
      end do outer_loop
      
      call quicksort(filled_col,size(filled_col))
      deallocate(stack,visited)         
  end subroutine dfs_reach_incomplete

  subroutine get_dependency_incomplete(A_csc,dependency_gr,depend_limit)
   use array_processes
   implicit none
   type(graph),intent(inout)                      :: A_csc
   type(graph),intent(out)                        :: dependency_gr
   integer                                        :: i,j,k
   integer, allocatable                           :: temp(:)
   logical, allocatable                           :: visited(:)
   integer, allocatable                           :: stack(:)
   integer, intent(in)                            :: depend_limit
   integer                                        :: stack_top ,stack_size
   integer                                        :: old_size, new_size
   integer                                        :: current_vertex, neighbor, start_idx, end_idx
   integer                                        :: count
   integer, allocatable                           :: temp_row(:)
   
   allocate(dependency_gr%col(0:size(A_csc%col)-1))
   dependency_gr%col=0
      
   
   do k= 2, size(A_csc%col)-1
      
     temp=pack(A_csc%row(A_csc%col(k-1)+1:A_csc%col(k)),A_csc%row(A_csc%col(k-1)+1:A_csc%col(k))<k)
     stack_top=0 
     count=0
    
     
     if(size(temp) > 0)  then
     
      stack_size=10  !extra 2 for safety
      allocate(stack(stack_size))
      stack(:)=0
      allocate(visited(k-1))
      visited(:) = .false.
      
      stack_top=stack_top + size(temp)
      stack(1:size(temp))=temp(size(temp):1:-1)  !last element in temporaray array
      
        
      print*,k
       do while (stack_top > 0 )
       
       
       
         
      
        current_vertex=stack( stack_top )
          
        stack_top= stack_top - 1
        
        
        
         if( .not. visited(current_vertex)) then
         ! mark as visited and add to result
          visited(current_vertex)= .true. 
          
          !get neighbours and push unvisited ones to stack
          
          start_idx=A_csc%col(current_vertex-1)+1
          end_idx  =A_csc%col(current_vertex)
          
         
          
          ! Push in reverse order to maintain left-to-right traversal
            do i=end_idx,start_idx,-1
            
               neighbor= A_csc%row(i)
            
                if( (neighbor < k .and. neighbor > current_vertex ).and. .not. any(temp == neighbor) )then
                  
                   if( .not. visited(neighbor)) then
                     if (stack_top < stack_size) then  ! Check stack space
                          stack_top = stack_top + 1
                          stack(stack_top) = neighbor
                         ! print*,"col=",k, "Before append: temp size=", size(temp), "adding neighbor=", neighbor
                          !if(k==5136) print* 
                          call append_to_array(temp, neighbor)
                          
                          !print*, "After append: temp size=", size(temp)
                         ! if (size(temp) > 0) print*, "Last element:", temp(size(temp))

                      endif
                  
                   end if
                end if
              end do
              
            end if
           
           
          end do
          
         
         
         if( allocated(stack) )   deallocate(stack)
         if( allocated(visited) ) deallocate(visited)      
       
     end if
    
     call quicksort(temp,size(temp))
     
     old_size = size(dependency_gr%row)
     new_size = old_size + size(temp)
     allocate(temp_row(new_size))
     
     temp_row(1:old_size) = dependency_gr%row
     temp_row(old_size+1:new_size) = temp
     
     if(allocated(dependency_gr%row)) deallocate(dependency_gr%row)
     allocate(dependency_gr%row(new_size))
     
     dependency_gr%row = temp_row
     deallocate(temp_row)
     dependency_gr%col(k)=dependency_gr%col(k-1)+size(temp)
     !print*, "dependency graph Completed col=", k
    
     if( allocated(temp) )    deallocate(temp)
     
     
    
      
   end do 
   
   print*,"dependency graph done"
   
  end subroutine get_dependency_incomplete

  subroutine symbolic_ILU(A_csc,L,U,limit)
   use node_classes
   use array_processes
   implicit none
   type(node_csc),intent(inout)                   :: A_csc
   type(node_csc),intent(out)                     :: L
   type(node_csc),intent(out)                     :: U
   type(graph)                                    :: A_adjacency
   type(graph)                                    :: A_dependency
   integer,intent(in)                             :: limit
   integer,allocatable                            :: temp1(:)
   integer,allocatable                            :: filled_array(:)
   integer                                        :: k,j,i,s,n,min
   integer                                        :: old_size,new_size
   integer,allocatable                            :: temp_row(:)
   integer                                        :: row_indices,col_indices
   integer                                        :: chunk_size
   
   
    n= size(A_csc%col)-1
   
    allocate(L%col(0:n))
    allocate(U%col(0:n))
    allocate(A_adjacency%col(0:n))
    
    L%col(0)=0
    U%col(0)=0
    
    A_adjacency%row=A_csc%row
    A_adjacency%col=A_csc%col
    
    call get_dependency_incomplete(A_adjacency,A_dependency,limit)
    !insert existing A[1:n,i] to L[:,i]
    do i=1,n
      
      
    
      
      call dfs_reach_incomplete( A_adjacency,L, A_dependency,filled_array,i,limit)
       
     
    
      !insert existing A[1:n,i] to L[:,i]
      temp1=pack(filled_array,filled_array>=i)
     
      old_size = size(L%row)
      new_size = old_size + size(temp1)
      allocate(temp_row(new_size))
     
      temp_row(1:old_size) = L%row
      temp_row(old_size+1:new_size) = temp1
     
      if(allocated(L%row)) deallocate(L%row)
      allocate(L%row(new_size))
     
      L%row = temp_row
      deallocate(temp_row)
      
      L%col(i)=size(temp1) + L%col(i-1) 
    
    
     
   
      if(allocated(temp1)) deallocate(temp1)
      !insert existing A[i,1:n to U[i,:]
     
       temp1=pack(filled_array,filled_array<i)
       old_size = size(U%row)
       new_size = old_size + size(temp1)
       allocate(temp_row(new_size))
      
       temp_row(1:old_size) = U%row
       temp_row(old_size+1:new_size) = temp1
     
       if(allocated(U%row)) deallocate(U%row)
       allocate(U%row(new_size))
     
       U%row = temp_row
       deallocate(temp_row)
   
       U%col(i)=size(temp1)+ U%col(i-1)

      
       if(allocated(temp1)) deallocate(temp1)
     print*,"Symbolic ILU done for col=",i
   end do
    L%nz_size=size(L%row)
    U%nz_size=size(U%row)
    
    ! deallocate(temp1,A_dependency%col,A_dependency%row,A_adjacency%col,A_adjacency%row)
  print*,"Symbolic ILU done"
  print*,"NNZ in L",L%nz_size
  print*,"NNZ in U",U%nz_size
  open(unit=10,file="../output_file/L_10_5_row.dat")
     do i=1,size(A_csc%col)-1
   
      write(10,'(*(I6))'),L%row(L%col(i-1)+1:L%col(i)) 
   end do
  close(10)
  end subroutine symbolic_ILU
end module symbolic_factorisation_op


module matrix_vector_op
 implicit none
 interface mat_vec_mult
    module procedure mat_vec_mult_cir
    module procedure mat_vec_mult_csr
 end interface   
contains
   subroutine mat_vec_mult_cir(mat,vec_in,vec_out)
   use node_classes
   implicit none
   type(node), allocatable, intent(inout)     :: mat(:)
   real(kind=8), allocatable ,intent (in)     :: vec_in(:)
   real(kind=8), allocatable ,intent (out)    :: vec_out(:)
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
     
   end subroutine mat_vec_mult_cir
   
   subroutine mat_vec_mult_csr(mat,vec_in,vec_out)
   use node_classes
   implicit none
   type(node_csr), intent(inout)              :: mat
   real(kind=8), allocatable ,intent (in)     :: vec_in(:)
   real(kind=8), allocatable ,intent (out)    :: vec_out(:)
   integer :: row_pt,i
   real(kind=8) :: sum
   allocate(vec_out(size(vec_in)))
   do row_pt=1, size(mat%row)-1
       sum=0
       do i=1,mat%row(row_pt)-mat%row(row_pt-1)
         sum=sum+mat%coeffs(mat%row(row_pt-1)+i)*vec_in(mat%col(mat%row(row_pt-1)+i))
       end do
        vec_out(row_pt)=sum
   end do
     
   end subroutine mat_vec_mult_csr
   
   
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
  
       
        vec_out=vec_in*scal_in
        
   
     
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
      
       
         vec_out=vec_in1+vec_in2
        
     
      
      case (-1)
     
       
         vec_out=vec_in1-vec_in2
        
      
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
    
   
    
    
    
    
    open(unit=50,file="../output_file/IL_10_5.dat")
    open(unit=60,file="../output_file/IU_10_5.dat")
    !open(unit=70,file="../output_file/soln_k.dat")
    
    !print*,"Entering Symbolic LU"
    call symbolic_ILU(A_csc,L_csc,U_csc,4)
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
   type(node_csr), allocatable, intent(inout)         :: A(:)
   type(node), allocatable, intent(out)               :: J(:)
   integer,intent(in)                                 :: ny,nz,n
   integer                                            :: i,k,t
   logical, dimension(:), allocatable                 :: mask
   allocate(J(n))
   do i=1,size(J)
     allocate(J(i)%col(1),J(i)%coeffs(1))
     J(i)%col(1)=i
     mask=(A(i)%col==J(i)%col(1))
        do t=1,size(A(i)%col)
         if(mask(t))then
          j(i)%coeffs(1)=1/A(i)%coeffs(t)
         end if
        end do
     
   end do             
  end subroutine jacobi_preconditioner
  !########################################################################

  
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
  
   !call system_clock(start_count, count_rate, count_max)
   !start_time = real(start_count)/real(count_rate)
   
   do i=1,size(b)-1 
    denum(i+1)=1/L_csc%coeffs(L_csc%col(i)+1) 
    
    
     y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))=y(L_csc%row(L_csc%col(i-1)+2:L_csc%col(i)))-L_csc%coeffs(L_csc%col(i-1)+2:L_csc%col(i))*y(i)/L_csc%coeffs(L_csc%col(i-1)+1) 
     
     
     
   end do 
   
    y=y*denum
   
   
   

   !call system_clock(end_count, count_rate, count_max)
   !end_time = real(end_count)/real(count_rate)
   
   
   !print*, "Wall clock time for substitution:", end_time - start_time, "seconds"
   
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
   
   
   !print*, "Wall clock time for substitution:", end_time - start_time, "seconds"
   
  !print*,"no. of operations= " ,flop
   
   deallocate(denum)
  end subroutine LU_solve
  
  

  subroutine LU_decompositon(A_csc,L_csc,U_csc,b,q,ny,nz,n)
    use node_classes
    use array_processes
    use symbolic_factorisation_op
    use matrix_vector_op
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
    
   
    
    
    
    
    !open(unit=50,file="../output_file/L_10_3.dat")
    !open(unit=60,file="../output_file/U_10_3.dat")
    !open(unit=70,file="../output_file/soln_k.dat")
    
    !print*,"Entering Symbolic LU"
    call symbolic_LU(A_csc,L_csc,U_csc)
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
      !write(60,'(*(F15.6))'),U_csc%coeffs(U_csc%col(k-1)+1:U_csc%col(k))
      !solve L of column k
      L_csc%coeffs(L_csc%col(k-1)+1:L_csc%col(k))=soln_k(L_csc%row(L_csc%col(k-1)+1:L_csc%col(k)))  
      !write(50,'(*(F15.6))'),L_csc%coeffs(L_csc%col(k-1)+1:L_csc%col(k)) 
      
   
      !print*,"Numerical LU done for col=",k
      
      if(allocated(soln_k)) deallocate(soln_k)
      !if(allocated(soln_k)) deallocate(nz_row)
      if(allocated(denom)) deallocate(denom)
    end do
    
    !close(50)
    !close(60)
    !close(70)
    
    print*,"LU factorisation completed"
      
      
      
      !call LU_solve(L_csc,U_csc,b,q)
      
      
    
     
      
      !write solution to a file
      open(unit=70,file="../output_file/LU_soln_b10_5.dat")
      write(70,*),"y          z        by"
      do i=1,size(b)
       write(70,"(3F10.6)"),A_csc%y_loc(i),A_csc%z_loc(i),b(i)
      end do
      close(70)
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
  real,intent(in)                        :: tol
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
  real(kind=8),allocatable               :: d_arr(:) !dummy array
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
  open(unit=20,file='../output_file/ini_b_q_out_b10_3.dat')
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
  call vec_vec_mult(res_ini,res_ini,row) !row=res_o.res_o
  poly=res_ini
  res =res_ini
  open(unit=30,file='../output_file/ini_res_b10_3.dat')
   write(30,*) ,'iter        res'
   do  i=1,n
     write(30,*) ,i,res(i)
    end do
  close(30)
  print*,'BiCGSTAB started'
  !Preconditioning
  if(PRECOND=="CLU")then
     call csr_to_csc(A_csr,A_csc,n)
     call LU_decompositon(A_csc,K1,K2,b,q,ny,nz,n)
  else if(PRECOND=="ILU")then
     call csr_to_csc(A_csr,A_csc,n)
     call ILU_preconditioner(A_csc,K1,K2,b,q,ny,nz,n)
  else if(PRECOND=="JACOBI")then
    ! call jacobi_preconditioner(A_csr,K1,ny,nz,n)  
    
 
  end if
   call system_clock(start_count, count_rate, count_max)
     start_time = real(start_count)/real(count_rate)
     
  !Iteration start
  open(unit=10,file='../output_file/bicgstab_out_b10_5.dat')
  do i=1,2000
   if(PRECOND=="CLU")then
   
  
     call LU_solve(K1,K2,k_inv_poly,k1_inv_poly,poly)    ! k_inv_poly= (K^-1)p
     call mat_vec_mult(A_csr,k_inv_poly,v)               !v=A k_inv_poly
   else if(PRECOND=="ILU")then
   
  
     call ILU_solve(K1,K2,k_inv_poly,k1_inv_poly,poly)    ! k_inv_poly= (K^-1)p
     call mat_vec_mult(A_csr,k_inv_poly,v)                !v=A k_inv_poly  
   else if(PRECOND=="JACOBI")then 
    !call mat_vec_mult(K1,poly,k_inv_poly)               ! k_inv_poly= (K^-1)p
    !call mat_vec_mult(A_csr,k_inv_poly,v)               !v=A k_inv_poly
   else
     call mat_vec_mult(A_csr,poly,v)                     !v=Ap
   end if
   call vec_vec_mult(v,res_ini,v_dot_res)  !v.res_o
   
   alpha=row/v_dot_res 
   !print*,alpha
   call scale_vec_mult(v,alpha,d_arr)                    ! alpha * v
   call vec_vec_sum(res,d_arr,s,-1)                      ! s=r-alpha * v
   if(PRECOND=="CLU")then
    call LU_solve(K1,K2,k_inv_s,k1_inv_s,s)              ! k_inv_s= (K^-1)s
    call mat_vec_mult(A_csr,k_inv_s,t)                   ! t=A k_inv_s
 
    
    call L_solve(K1,k1_inv_t,t)                          !K1_inv t
    
    call vec_vec_mult(k1_inv_t,k1_inv_s,t_dot_s)         !(K1^-1)t.(K1^-1)s
    call vec_vec_mult(k1_inv_t,k1_inv_t,t_dot_t)         !(K1^-1)t.(K1^-1)t
    omega= t_dot_s/t_dot_t 
    call scale_vec_mult(k_inv_poly,alpha,alpha_poly)     !alpha * (K^-1)p
    call scale_vec_mult(k_inv_s,omega,omega_s)           !omega * (K^-1)s
   else if(PRECOND=="ILU")then
    call ILU_solve(K1,K2,k_inv_s,k1_inv_s,s)              ! k_inv_s= (K^-1)s
    call mat_vec_mult(A_csr,k_inv_s,t)                   ! t=A k_inv_s
 
    
    call IL_solve(K1,k1_inv_t,t)                          !K1_inv t
    
    call vec_vec_mult(k1_inv_t,k1_inv_s,t_dot_s)         !(K1^-1)t.(K1^-1)s
    call vec_vec_mult(k1_inv_t,k1_inv_t,t_dot_t)         !(K1^-1)t.(K1^-1)t
    omega= t_dot_s/t_dot_t 
    call scale_vec_mult(k_inv_poly,alpha,alpha_poly)     !alpha * (K^-1)p
    call scale_vec_mult(k_inv_s,omega,omega_s)           !omega * (K^-1)s
   else if(PRECOND=="JACOBI")then 
    !call mat_vec_mult(K1,s,k_inv_s)                  ! k_inv_s= (K^-1)s 
    !call mat_vec_mult(A_csr,k_inv_s,t)                   ! t=A k_inv_s
    !call mat_vec_mult(K1,s,k1_inv_s)     !K_inv s
    !call mat_vec_mult(K1,t,k1_inv_t)     !K_inv t
    !call vec_vec_mult(k1_inv_t,k1_inv_s,t_dot_s)          !(K^-1)t.(K^-1)s
    !call vec_vec_mult(k1_inv_t,k1_inv_t,t_dot_t)          !(K^-1)t.(K^-1)t
    omega= t_dot_s/t_dot_t 
    !call scale_vec_mult(k_inv_poly,alpha,alpha_poly)   !alpha * (K^-1)p
    !call scale_vec_mult(k_inv_s,omega,omega_s)         !omega * (K^-1)s
   else
    call mat_vec_mult(A_csr,s,t)                          !t=As
    call vec_vec_mult(t,s,t_dot_s)                        !t.s
    call vec_vec_mult(t,t,t_dot_t)                        !t.t
    omega= t_dot_s/t_dot_t 
    call scale_vec_mult(poly,alpha,alpha_poly)       !alpha * poly
    call scale_vec_mult(s,omega,omega_s)             !omega*s
   
   end if
              
   
   
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
   write(10,"(1X,I5,1X,4f10.6)") ,i,inf_n,A_csr%y_loc(j),A_csr%z_loc(j),norm_2
  end do
  call system_clock(end_count, count_rate, count_max)
   end_time = real(end_count)/real(count_rate)
   
  print*, "Iteration clock time:", end_time - start_time, "seconds"
  close(10)
  close(20) 
  print*,"alpha=",alpha
  open(unit=40,file='../output_file/by_res_final_b10.dat')
   write(40,*) ,'iter         y        z         by       res'
    do  j=0,ny
     do k=0,nz
        i=(k)*(ny+1)+j+1
        write(40,'(1X,I5,1X,7f10.6)') ,i,A_csr%y_loc(i),A_csr%z_loc(i),b(i),res(i)
     end do
       write(40,*) " "
    end do
  close(40)
  deallocate(res,v,poly,s,t,alpha_poly,omega_s,omega_t,omega_v)
  end subroutine biCGSTAB
  
  
  

  
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
   real,intent(in)                             :: tol
   integer,intent(in)                          :: ny,nz,n 
   character, intent(in)                       :: solver_name*10
   character, intent(in)                       :: precond*10
   if(solver_name=="BICGSTAB")then
     call biCGSTAB(A_csr,b,q,tol,ny,nz,precond)
     
   else if(solver_name=="LU_METHOD") then
      print*,"LU METHOD SELECTED"
      
      call csr_to_csc(A_csr,A_csc,n)
      call LU_decompositon(A_csc,L_csc,U_csc,b,q,ny,nz,n)
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
 
 
 integer                                          :: ny,nz,Rm,n,j,i,k,band
 real                                             :: Ly,Lz,res_tolerance,A,del_t,alpha_k
 real                                             :: init_by, init_bz,f
 character                                        :: linear_solver*10 ,pre_cond*10
 real(kind=8),dimension(:),allocatable            :: gridpts_y,gridpts_z
 real(kind=8),dimension(:),allocatable            :: q,by,act_b,error_b
 type(node), allocatable                          :: node_info(:),L(:),U(:)
 type(node_csr)                                   :: node_info_csr
 type(node_csc)                                   :: node_info_csc
 real,parameter                                   :: PI=4*atan(1.d0)
 real                                             :: y,z
 integer*4                                        :: count_rate, count_max, start_count, end_count
 real*8                                           :: start_time, end_time
 
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
 call fill_node_array(node_info,node_info_csr,gridpts_y,gridpts_z,ny,nz,n,f,alpha_k)
 call fill_q_array(node_info_csr,q,f,alpha_k)
 print *,'size of node_info' , size(node_info)
 print *,'size of q' , size(q)
 deallocate(gridpts_y)!memory freed
 deallocate(gridpts_z)!memory freed
 
 open(unit=30, file='../output_file/act_by_b10_5.dat')
       do k=0,nz
        do j=0,ny
         i=(k)*(ny+1)+j+1
         y=node_info_csr%y_loc(i)
         z=node_info_csr%z_loc(i)
         act_b(i)=0.01*cos(2.0*PI*y)*sin(2.0*PI*z)
     
         write(30,"(3f10.6)"),y,z,act_b(i)
       end do
        write(30,*) " "
      end do
 close(30)
 
 allocate(by(size(node_info)))
 by=init_by

 call system_clock(start_count, count_rate, count_max)
     start_time = real(start_count)/real(count_rate)
     
 call solver(node_info_csr,by,q,res_tolerance,ny,nz,n,linear_solver,pre_cond)
 
 call system_clock(end_count, count_rate, count_max)
   end_time = real(end_count)/real(count_rate)
   
print*, "Wall clock time:", end_time - start_time, "seconds"

 error_b=abs(act_b-by) 
 
 open(unit=40, file='../output_file/analytical_numerical_error_b10_5.dat')
      do k=0,nz
        do j=0,ny
         i=(k)*(ny+1)+j+1
         y=node_info_csr%y_loc(i)
         z=node_info_csr%z_loc(i)
         write(40,"(3f10.6)"),y,z,error_b(i)
       end do
        write(40,*) " "
      end do
 close(40)
 
  if(allocated(node_info)) deallocate(node_info)
  deallocate(by,q,act_b)
  
END PROGRAM main
