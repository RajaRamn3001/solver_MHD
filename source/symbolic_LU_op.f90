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
   integer                                        :: stacksize
   integer                                        :: old_size, extra_size
   integer, allocatable                           :: stack(:)
   integer                                        :: stack_top
   integer                                        :: current_vertex, neighbor, start_idx, end_idx
   integer                                        :: count
   
   n=size(dependency_gr%col) - 1
 
   
  
   stack_top=0 
   count=0
   
   print*,"dfs reach for column =",column 
   
   filled_col=A_csc%row(A_csc%col(column-1)+1:A_csc%col(column))
   stacksize=n/2
   
   allocate(visited(n))
   allocate(stack(stacksize))
   visited = .false.
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
                   if(column==773) print*, size(stack),stack_top
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
     print*, "dependency graph Completed col=", k
    
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
     print*,"Symbolic LU done for col=",i
   end do
    L%nz_size=size(L%row)
    U%nz_size=size(U%row)
   
    ! deallocate(temp1,A_dependency%col,A_dependency%row,A_adjacency%col,A_adjacency%row)
  print*,"Symbolic LU done"
  
  open(unit=10,file="../output_file/L_10_3_row.dat")
     do i=1,size(A_csc%col)-1
   
      write(10,'(*(I6))'),L%row(U%col(i-1)+1:L%col(i)) 
   end do
  close(10)
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
