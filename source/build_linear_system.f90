module node_classes
!This is a class for node informations



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
  
  
   subroutine csc_to_csr(csc_array,csr_array,row_num)
    type(node_csr),intent(out)                     :: csr_array
    type(node_csc),intent(inout)                   :: csc_array
    integer        ,intent(in)                     :: row_num
    integer                                        :: i,j,k
    integer                                        :: nz_size
    integer                                        :: count
    count=0
    nz_size=csc_array%col(row_num)
    allocate(csr_array%row(0:row_num),csr_array%y_loc(0:row_num),csr_array%z_loc(0:row_num))
    csr_array%y_loc(0:row_num)=csc_array%y_loc
    csr_array%z_loc(0:row_num)=csc_array%z_loc
    csr_array%nz_size=0
    csr_array%row(0)=csr_array%nz_size
    allocate(csr_array%coeffs(nz_size),csr_array%col(nz_size))
   
     do k=1,row_num
      do i=1,size(csc_array%col)-1
         do j=1,csc_array%col(i)-csc_array%col(i-1)
            if(csc_array%row(csc_array%col(i-1)+j)==k) then
            
               csr_array%nz_size=csr_array%nz_size+1
               count=count+1
               csr_array%coeffs(csr_array%nz_size)=csc_array%coeffs(csc_array%col(i-1)+j)
               csr_array%col(csr_array%nz_size) =i
               
               exit
            end if
         end do
           !print*,"k=",k,"csr_array%nz_size=",csr_array%nz_size
       end do
        csr_array%row(k)=csr_array%nz_size
     
     end do
    
     print*,"NNZ of csr_array",count
     print*,"NNZ of csc_array",nz_size
     open(unit=50, file='../output_file/inv_csr_b10_5_coeffs.dat')
     open(unit=60, file='../output_file/inv_csr_b10_5_cols.dat')
      write(50,*),"ap       an       as        aw       ae"
 
      do j=1,row_num
   
  
        
        
        write(50,"(5F15.6)"),csr_array%coeffs(csr_array%row(j-1)+1:csr_array%row(j))
        write(60,*),csr_array%col(csr_array%row(j-1)+1:csr_array%row(j))
       
      end do 
  
     
      close(50)
      close(60)
  end subroutine csc_to_csr

  subroutine DH_pseudo_bound_coeff_mat(node_arr_csr,grid_pty,grid_ptz,ny,nz,n,f,alpha_k)
  
    implicit none
    integer,intent(in)                             :: ny,nz,n 
    real*8,intent(in)                              :: f,alpha_k
    real(kind=8),intent(in)                        :: grid_pty(0:ny),grid_ptz(0:nz)
    integer                                        :: i,j,k,count,c,chunk_size
    type(node_csr),intent(inout)                   :: node_arr_csr
    real(kind=8),allocatable                       :: temp_coeffs(:)
    integer,allocatable                            :: temp_col(:)
    integer,allocatable                            :: temp_row(:)
    integer                                        :: row_ptr,col_ptr,start_from

    count=0
    chunk_size=1000
    node_arr_csr%nz_size=0
    node_arr_csr%capacity=ny+2                                            !ny+2=258
    
   
    row_ptr=0
    col_ptr=0                              
    allocate(node_arr_csr%row(0:n),node_arr_csr%y_loc(0:n),node_arr_csr%z_loc(0:n))
    allocate(node_arr_csr%col(node_arr_csr%capacity))
    allocate(node_arr_csr%coeffs(node_arr_csr%capacity))
    
    
    
    node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
    node_arr_csr%y_loc(0)=0.d0
    node_arr_csr%z_loc(0)=0.d0
    print*,'size of node array',size(node_arr_csr%row)-1
    
    !Preparing the coefficients that will solve by
     do k=0,nz
      do j=0,ny
       
       i=k*(ny+1)+j+1
       row_ptr=row_ptr+1
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
    
       
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
            
            
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)=1.0/(grid_pty(j-1)-grid_pty(j))/(c+c**2)
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-(1+c)**2/(grid_pty(j-1)-grid_pty(j))/(c+c**2)
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(2*c+c**2)/(grid_pty(j-1)-grid_pty(j))/(c+c**2)
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j-2+1
            node_arr_csr%col(node_arr_csr%nz_size+2)=(k)*(ny+1)+j-1+1
            node_arr_csr%col(node_arr_csr%nz_size+3)=(k)*(ny+1)+j+1
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+3
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
            
            
        else
          
          
          node_arr_csr%coeffs(node_arr_csr%nz_size+1) =1.0
          node_arr_csr%col(node_arr_csr%nz_size+1)    =(k)*(ny+1)+j+1
          
          node_arr_csr%nz_size=node_arr_csr%nz_size+1
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          
        end if
       
                   
       end do
     end do
     
     !Preparing the coefficients that will solve bz
     !row_ptr will continue from 66049
     start_from=row_ptr
     do k=0,nz
      do j=0,ny
       
       i=k*(ny+1)+j+1
       row_ptr=row_ptr+1
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
    
       
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
          
          
         
          
          node_arr_csr%coeffs(node_arr_csr%nz_size+1)=2/(grid_ptz(k+1)-grid_ptz(k-1))/(grid_ptz(k)-grid_ptz(k-1))
          node_arr_csr%coeffs(node_arr_csr%nz_size+2)=2/(grid_pty(j+1)-grid_pty(j-1))/(grid_pty(j)-grid_pty(j-1))
          node_arr_csr%coeffs(node_arr_csr%nz_size+3)=-(alpha_k**2+f+2/(grid_pty(j+1)-grid_pty(j))/(grid_pty(j)-grid_pty(j-1))+2/(grid_ptz(k+1)-grid_ptz(k))/(grid_ptz(k)-grid_ptz(k-1)))
          node_arr_csr%coeffs(node_arr_csr%nz_size+4)=2/(grid_pty(j+1)-grid_pty(j-1))/(grid_pty(j+1)-grid_pty(j))
          node_arr_csr%coeffs(node_arr_csr%nz_size+5)=2/(grid_ptz(k+1)-grid_ptz(k-1))/(grid_ptz(k+1)-grid_ptz(k))
          
          node_arr_csr%col(node_arr_csr%nz_size+1)=start_from + (k-1)*(ny+1)+j+1
          node_arr_csr%col(node_arr_csr%nz_size+2)=start_from + (k)*(ny+1)+j-1+1
          node_arr_csr%col(node_arr_csr%nz_size+3)=start_from + (k)*(ny+1)+j+1
          node_arr_csr%col(node_arr_csr%nz_size+4)=start_from + (k)*(ny+1)+j+1+1
          node_arr_csr%col(node_arr_csr%nz_size+5)=start_from + (k+1)*(ny+1)+j+1
          
       
          
          node_arr_csr%nz_size=node_arr_csr%nz_size+5
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          
          count=count+1 
          
          !Boundary condition will change for bz 
        else if(k==0)then
             c=(grid_ptz(k+2)-grid_ptz(k+1))/(grid_ptz(k+1)-grid_ptz(k))
             
             
             node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-(2*c+c**2)/(grid_ptz(k+1)-grid_ptz(k))/(c+c**2)
             node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(1+c)**2/(grid_ptz(k+1)-grid_ptz(k))/(c+c**2)
             node_arr_csr%coeffs(node_arr_csr%nz_size+3)=-1.0/(grid_ptz(k+1)-grid_ptz(k))/(c+c**2)
             
             node_arr_csr%col(node_arr_csr%nz_size+1)=start_from + (k)*(ny+1)+j+1
             node_arr_csr%col(node_arr_csr%nz_size+2)=start_from + (k+1)*(ny+1)+j+1
             node_arr_csr%col(node_arr_csr%nz_size+3)=start_from + (k+2)*(ny+1)+j+1
             
             node_arr_csr%nz_size=node_arr_csr%nz_size+3
             node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
             
        else if(k==nz)then
            !print*,'ok1'
            c=(grid_ptz(k-2)-grid_ptz(k-1))/(grid_ptz(k-1)-grid_ptz(k))
            
            
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)=1.0/(grid_ptz(k-1)-grid_ptz(k))/(c+c**2)
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-(1+c)**2/(grid_ptz(k-1)-grid_ptz(k))/(c+c**2)
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(2*c+c**2)/(grid_ptz(k-1)-grid_ptz(k))/(c+c**2)
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=start_from + (k-2)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+2)=start_from + (k-1)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+3)=start_from + (k)*(ny+1)+j+1
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+3
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
            
            
        else
          
          
          node_arr_csr%coeffs(node_arr_csr%nz_size+1) =1.0
          node_arr_csr%col(node_arr_csr%nz_size+1)    =start_from + (k)*(ny+1)+j+1
          
          node_arr_csr%nz_size=node_arr_csr%nz_size+1
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          
        end if
       
                   
       end do
     end do
    
     open(unit=50, file='../output_file/node_info_csr_b10_4_coeffs.dat')
     open(unit=60, file='../output_file/node_info_csr_b10_4_cols.dat')
      
      do j=1,n
    
  
    
       
        
        
        write(50,"(5F15.6)"),node_arr_csr%coeffs(node_arr_csr%row(j-1)+1:node_arr_csr%row(j))
        write(60,*),node_arr_csr%col(node_arr_csr%row(j-1)+1:node_arr_csr%row(j))
       
      end do 
  
     
      close(50)
      close(60)
      print*,"Node array filled" 
    
    end subroutine DH_pseudo_bound_coeff_mat
    
    subroutine fill_DH_bem_coeff_mat(node_arr_csr,bo_int_arr_csr,grid_pty,grid_ptz,ny,nz,ns,f,alpha_k)
    
    use greens_func_subroutines
    use array_processes
    implicit none
    integer,intent(in)                             :: ny,nz,ns 
    real*8,intent(in)                              :: f,alpha_k
    real(kind=8),intent(in)                        :: grid_pty(0:ny),grid_ptz(0:nz)
    integer                                        :: i,j,k,count,chunk_size
    real*8                                         :: c,d,m,n,p,q
    real*8                                         :: h1,h2,h3
    type(node_csr),intent(in)                      :: bo_int_arr_csr
    type(node_csr),intent(inout)                   :: node_arr_csr
    real(kind=8),allocatable                       :: temp_coeffs(:)
    integer,allocatable                            :: temp_col(:)
    integer,allocatable                            :: temp_row(:)
    integer                                        :: row_ptr,col_ptr,start_from
    integer*4                                      :: yb,zb,mp,me,side_tr, side_tr1,nb,np
    integer*4                                      :: tet
    
    
    count=0
    chunk_size=1000
    node_arr_csr%nz_size=0
    node_arr_csr%capacity=ny+2                                            !ny+2=258
     
    nb=2*(ny+nz)
    np=(ny+1)*(nz+1)
    row_ptr=0
    col_ptr=0                              
    allocate(node_arr_csr%row(0:ns),node_arr_csr%y_loc(0:ns),node_arr_csr%z_loc(0:ns))
    allocate(node_arr_csr%col(node_arr_csr%capacity))
    allocate(node_arr_csr%coeffs(node_arr_csr%capacity))
    
    
    
    node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
    node_arr_csr%y_loc(0)=0.d0
    node_arr_csr%z_loc(0)=0.d0
    print*,'size of node array',size(node_arr_csr%row)-1
!-------------------------------------------------------------------------------------------    
!Preparing the coefficients that will solve b_y_hat
!----------------------------------------------------------------------------------
     do k=0,nz
      do j=0,ny
       
       i=k*(ny+1)+j+1
       row_ptr=row_ptr+1
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
    
       
       node_arr_csr%y_loc(row_ptr)=grid_pty(j)
       node_arr_csr%z_loc(row_ptr)=grid_ptz(k)
       
       
      !---------------------------------------------------------------------------------
      ! Add memory to col and coeff array if required
      !---------------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------------------
! Prepare interior element coeffs  
!-----------------------------------------------------------------------------------------                          
       if(j>0 .and. j<ny .and. k>0 .and. k<nz ) then
          
          
         
          
          node_arr_csr%coeffs(node_arr_csr%nz_size+1)=1.d0
          
          node_arr_csr%col(node_arr_csr%nz_size+1)=k*(ny+1)+j+1
          
          
          node_arr_csr%nz_size=node_arr_csr%nz_size+1
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          
        
        !----------------------------------------------------------------------
        !Prepare coeffs for normal by_hat at j=0
        !----------------------------------------------------------------------
        else if(j==0 .and. .not. (k==0 .or. k==nz))then
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
             deallocate(temp_coeffs)
             deallocate(temp_col)
            end if
           end if
             h1=grid_pty(j+1)-grid_pty(j)
             h2=grid_pty(j+2)-grid_pty(j)
             h3=grid_pty(j+3)-grid_pty(j)
            
             !f''~ q*f0 + m*f1 + n*f2 + p*f3
             m=2*h2*h3/(h1*(h2-h1)*(h3-h1))     
             n=-2*h1*h3/(h2*(h2-h1)*(h3-h2))
             p=2*h1*h2/(h3*(h3-h1)*(h3-h2))
             q=-(m+n+p) 
                         
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=1.d0
             
            node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j+1
             
             node_arr_csr%nz_size=node_arr_csr%nz_size+1
             node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
        !----------------------------------------------------------------------
        !Prepare coeffs for normal by_hat at j=ny
        !----------------------------------------------------------------------      
        else if(j==ny .and. .not. (k==0 .or. k==nz))then
            !print*,'ok1'
           !-----------------------------------------------------------------
           !    Add memory to coeff and col array if required
           !-----------------------------------------------------------------
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
             deallocate(temp_coeffs)
             deallocate(temp_col)
            end if
           end if
          !-----------------------------------------------------------------
            h1=grid_pty(j)-grid_pty(j-1)
            h2=grid_pty(j)-grid_pty(j-2)
            h3=grid_pty(j)-grid_pty(j-3)
            
             !f''~ q*f0 + m*f1 + n*f2 + p*f3
             m=2*h2*h3/(h1*(h2-h1)*(h3-h1))     
             n=-2*h1*h3/(h2*(h2-h1)*(h3-h2))
             p=2*h1*h2/(h3*(h3-h1)*(h3-h2))
             q=-(m+n+p) 
             
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)=1.d0
             
            node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j+1
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+1
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
            
        !------------------------------------------------------------------------
        !Prepare coeffs for tangential by_hat at k=0 
        !------------------------------------------------------------------------    
        else if(k==0) then
         side_tr1=4  !bottom face
         !if node is corner i.e j=0 or j=ny
         mp=ibound(j,k,ny,nz)
        
         if(corner(mp,ny,nz)) then
          !extrapolate by from neighbours
           if(j==0) then
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)= -1.d0
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(grid_pty(j) - grid_pty(j+2))*(grid_pty(j) - grid_pty(j+3))/(grid_pty(j+1) - grid_pty(j+2))/(grid_pty(j+1) - grid_pty(j+3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(grid_pty(j) - grid_pty(j+1))*(grid_pty(j) - grid_pty(j+3))/(grid_pty(j+2) - grid_pty(j+1))/(grid_pty(j+2) - grid_pty(j+3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+4)=(grid_pty(j) - grid_pty(j+1))*(grid_pty(j) - grid_pty(j+2))/(grid_pty(j+3) - grid_pty(j+1))/(grid_pty(j+3) - grid_pty(j+2))
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j+1 
            node_arr_csr%col(node_arr_csr%nz_size+2)=(k)*(ny+1)+j+1 + 1
            node_arr_csr%col(node_arr_csr%nz_size+3)=(k)*(ny+1)+j+1 + 2
            node_arr_csr%col(node_arr_csr%nz_size+4)=(k)*(ny+1)+j+1 + 3
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+4
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
           else
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)=(grid_pty(j) - grid_pty(j-2))*(grid_pty(j) - grid_pty(j-3))/(grid_pty(j-1) - grid_pty(j-2))/(grid_pty(j-1) - grid_pty(j-3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(grid_pty(j) - grid_pty(j-1))*(grid_pty(j) - grid_pty(j-3))/(grid_pty(j-2) - grid_pty(j-1))/(grid_pty(j-2) - grid_pty(j-3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(grid_pty(j) - grid_pty(j-1))*(grid_pty(j) - grid_pty(j-2))/(grid_pty(j-3) - grid_pty(j-1))/(grid_pty(j-3) - grid_pty(j-2))
            node_arr_csr%coeffs(node_arr_csr%nz_size+4)=-1.d0
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j+1 - 1
            node_arr_csr%col(node_arr_csr%nz_size+2)=(k)*(ny+1)+j+1 - 2
            node_arr_csr%col(node_arr_csr%nz_size+3)=(k)*(ny+1)+j+1 - 3  
            node_arr_csr%col(node_arr_csr%nz_size+4)=(k)*(ny+1)+j+1 
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+4
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
           end if
           
         else
          side_tr=0
          tet=0
          !---------------------------------------(-1)*b_tau---------------------------------------------
          node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-0.5d0
          node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j+1 
          node_arr_csr%nz_size=node_arr_csr%nz_size + 1
          
          !----------------------------------(+H*b_n)-------------------------------------------------
           do me=1,nb
           !-----------------------------------------------------------------
           ! Check and add memory to coeff and col array if required
           !-----------------------------------------------------------------
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
             deallocate(temp_coeffs)
             deallocate(temp_col)
            end if
           end if
          !-----------------------------------------------------------------
            if(corner(me,ny,nz)) then
            
             yb=jb(me,ny,nz)
             zb=kb(me,ny,nz)
             
             if(me==1) then  !corner 
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else if(me==nz+1) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else if(me==2*nz+ny+1) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             end if 
             node_arr_csr%nz_size=node_arr_csr%nz_size + 2
             side_tr=side_tr+1
             if(tet==1) then
              tet=0
             else
              tet=1
             end if
            else
             yb=jb(me,ny,nz)
             zb=kb(me,ny,nz)
             
             if( (me> 1 .and. me < nz+1) .or. (me> 2*nz + ny + 1 .and. me < 2*nz+2*ny +1 )) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             else
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             end if
             
             node_arr_csr%nz_size=node_arr_csr%nz_size + 1
            end if
           end do
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          !temp_coeffs=node_arr_csr%coeffs(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))
          !temp_col=node_arr_csr%col(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))
         
          !call quicksort(temp_coeffs,temp_col,size(temp_col))
         
          !node_arr_csr%coeffs(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))=temp_coeffs
          !node_arr_csr%col(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))=temp_col  
         endif
          
  
        !------------------------------------------------------------------------
        !Prepare coeffs for tangential by_hat at k=nz 
        !------------------------------------------------------------------------
        else if(k==nz) then 
         side_tr1=2   !top face
         !if node is corner i.e j=0 or j=ny
         mp=ibound(j,k,ny,nz)
        
         if(corner(mp,ny,nz)) then
           !extrapolate by from neighbours
           if(j==0) then
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)= -1.d0
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(grid_pty(j) - grid_pty(j+2))*(grid_pty(j) - grid_pty(j+3))/(grid_pty(j+1) - grid_pty(j+2))/(grid_pty(j+1) - grid_pty(j+3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(grid_pty(j) - grid_pty(j+1))*(grid_pty(j) - grid_pty(j+3))/(grid_pty(j+2) - grid_pty(j+1))/(grid_pty(j+2) - grid_pty(j+3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+4)=(grid_pty(j) - grid_pty(j+1))*(grid_pty(j) - grid_pty(j+2))/(grid_pty(j+3) - grid_pty(j+1))/(grid_pty(j+3) - grid_pty(j+2))
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j+1 
            node_arr_csr%col(node_arr_csr%nz_size+2)=(k)*(ny+1)+j+1 + 1
            node_arr_csr%col(node_arr_csr%nz_size+3)=(k)*(ny+1)+j+1 + 2
            node_arr_csr%col(node_arr_csr%nz_size+4)=(k)*(ny+1)+j+1 + 3
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+4
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
           else
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)=(grid_pty(j) - grid_pty(j-2))*(grid_pty(j) - grid_pty(j-3))/(grid_pty(j-1) - grid_pty(j-2))/(grid_pty(j-1) - grid_pty(j-3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(grid_pty(j) - grid_pty(j-1))*(grid_pty(j) - grid_pty(j-3))/(grid_pty(j-2) - grid_pty(j-1))/(grid_pty(j-2) - grid_pty(j-3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(grid_pty(j) - grid_pty(j-1))*(grid_pty(j) - grid_pty(j-2))/(grid_pty(j-3) - grid_pty(j-1))/(grid_pty(j-3) - grid_pty(j-2))
            node_arr_csr%coeffs(node_arr_csr%nz_size+4)=-1.d0
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j+1 - 1
            node_arr_csr%col(node_arr_csr%nz_size+2)=(k)*(ny+1)+j+1 - 2
            node_arr_csr%col(node_arr_csr%nz_size+3)=(k)*(ny+1)+j+1 - 3  
            node_arr_csr%col(node_arr_csr%nz_size+4)=(k)*(ny+1)+j+1 
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+4
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
           end if
         !else if any node in-between
         else
          side_tr=0
          tet=0
          !---------------------------------------(-1)*b_tau---------------------------------------------
          node_arr_csr%coeffs(node_arr_csr%nz_size+1)=0.5d0
          node_arr_csr%col(node_arr_csr%nz_size+1)=(k)*(ny+1)+j+1 
          node_arr_csr%nz_size=node_arr_csr%nz_size + 1
          
          !----------------------------------(+H*b_n)-------------------------------------------------
           do me=1,nb
           !-----------------------------------------------------------------
           ! Check and add memory to coeff and col array if required
           !-----------------------------------------------------------------
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
             deallocate(temp_coeffs)
             deallocate(temp_col)
            end if
           end if
          !-----------------------------------------------------------------
            if(corner(me,ny,nz)) then
            
             yb=jb(me,ny,nz)
             zb=kb(me,ny,nz)
             if(me==1) then  !corner 
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else if(me==nz+1) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else if(me==2*nz+ny+1) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             end if 
             
             node_arr_csr%nz_size=node_arr_csr%nz_size + 2
             side_tr=side_tr+1
             if(tet==1) then
              tet=0
             else
              tet=1
             end if
            else
             yb=jb(me,ny,nz)
             zb=kb(me,ny,nz)
             
             if( (me> 1 .and. me < nz+1) .or. (me> 2*nz + ny + 1 .and. me < 2*nz+2*ny +1 )) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             else
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             end if
             
             node_arr_csr%nz_size=node_arr_csr%nz_size + 1
            end if
           end do
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size      
          !temp_coeffs=node_arr_csr%coeffs(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))
          !temp_col=node_arr_csr%col(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))
         
          !call quicksort(temp_coeffs,temp_col,size(temp_col))
         
          !node_arr_csr%coeffs(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))=temp_coeffs
          !node_arr_csr%col(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))=temp_col 
         endif
        end if
       
                   
       end do
     end do
 !--------------------------------------------------------------------------------------    
     !Preparing the coefficients that will solve b_z_hat
 !--------------------------------------------------------------------------------------
     !row_ptr will continue from 66049
     start_from=row_ptr
     
     do k=0,nz
      do j=0,ny
       
       i=k*(ny+1)+j+1
       row_ptr=row_ptr+1
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
    
       
       node_arr_csr%y_loc(row_ptr)=grid_pty(j)
       node_arr_csr%z_loc(row_ptr)=grid_ptz(k)
       
       !-----------------------------------------------------------------
       ! Check and add memory to coeff and col array if required
       !-----------------------------------------------------------------
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
             deallocate(temp_coeffs)
             deallocate(temp_col)
            end if
           end if
        !-----------------------------------------------------------------
       
       
       
       
       
       !count=count+1
!-----------------------------------------------------------------------------------------
! Prepare interior element coeffs  
!-----------------------------------------------------------------------------------------                                
       if(j>0 .and. j<ny .and. k>0 .and. k<nz ) then
          
          
         
          
          node_arr_csr%coeffs(node_arr_csr%nz_size+1)=1.d0
         
          node_arr_csr%col(node_arr_csr%nz_size+1)=start_from + (k)*(ny+1)+j+1
          
       
          
          node_arr_csr%nz_size=node_arr_csr%nz_size+1
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          
          count=count+1 
  !----------------------------------------------------------------------
  !Prepare coeffs for normal bz_hat at k=0
  !----------------------------------------------------------------------
        else if(k==0 .and. .not. (j==0 .or. j==ny))then
             h1=grid_ptz(k+1)-grid_ptz(k)
             h2=grid_ptz(k+2)-grid_ptz(k)
             h3=grid_ptz(k+3)-grid_ptz(k)
            
             !f''~ q*f0 + m*f1 + n*f2 + p*f3
             m=2*h2*h3/(h1*(h2-h1)*(h3-h1))     
             n=-2*h1*h3/(h2*(h2-h1)*(h3-h2))
             p=2*h1*h2/(h3*(h3-h1)*(h3-h2))
             q=-(m+n+p) 
              
             node_arr_csr%coeffs(node_arr_csr%nz_size+1)=1.d0
             
             node_arr_csr%col(node_arr_csr%nz_size+1)=start_from + (k)*(ny+1)+j+1
             
             node_arr_csr%nz_size=node_arr_csr%nz_size+1
             node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
 
  !----------------------------------------------------------------------
  !Prepare coeffs for normal bz_hat at k=nz
  !----------------------------------------------------------------------            
        else if(k==nz .and. .not. (j==0 .or. j==ny))then
            !print*,'ok1'
             h1=grid_ptz(k)-grid_ptz(k-1)
             h2=grid_ptz(k)-grid_ptz(k-2)
             h3=grid_ptz(k)-grid_ptz(k-3)
            
             !f''~ q*f0 + m*f1 + n*f2 + p*f3
             m=2*h2*h3/(h1*(h2-h1)*(h3-h1))     
             n=-2*h1*h3/(h2*(h2-h1)*(h3-h2))
             p=2*h1*h2/(h3*(h3-h1)*(h3-h2))
             q=-(m+n+p)
            
             node_arr_csr%coeffs(node_arr_csr%nz_size+1)=1.d0
             
             node_arr_csr%col(node_arr_csr%nz_size+1)=start_from + (k)*(ny+1)+j+1
             
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+1
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
  !------------------------------------------------------------------------
  !Prepare coeffs for tangential bz_hat at j=0 
  !------------------------------------------------------------------------             
            
        else if (j==0) then
         side_tr1=1  !left face
         !if node is corner i.e j=0 or j=ny
         mp=ibound(j,k,ny,nz)
      
         if(corner(mp,ny,nz)) then
          !extrapolate by from neighbours
           if(k==0) then
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)= -1.d0
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(grid_ptz(k) - grid_ptz(k+2))*(grid_ptz(k) - grid_ptz(k+3))/(grid_ptz(k+1) - grid_ptz(k+2))/(grid_ptz(k+1) - grid_ptz(k+3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(grid_ptz(k) - grid_ptz(k+1))*(grid_ptz(k) - grid_ptz(k+3))/(grid_ptz(k+2) - grid_ptz(k+1))/(grid_ptz(k+2) - grid_ptz(k+3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+4)=(grid_ptz(k) - grid_ptz(k+1))*(grid_ptz(k) - grid_ptz(k+2))/(grid_ptz(k+3) - grid_ptz(k+1))/(grid_ptz(k+3) - grid_ptz(k+2))
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=start_from +(k)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+2)=start_from +(k+1)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+3)=start_from +(k+2)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+4)=start_from +(k+3)*(ny+1)+j+1 
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+4
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
           else
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)=(grid_ptz(k) - grid_ptz(k-2))*(grid_ptz(k) - grid_ptz(k-3))/(grid_ptz(k-1) - grid_ptz(k-2))/(grid_ptz(k-1) - grid_ptz(k-3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(grid_ptz(k) - grid_ptz(k-1))*(grid_ptz(k) - grid_ptz(k-3))/(grid_ptz(k-2) - grid_ptz(k-1))/(grid_ptz(k-2) - grid_ptz(k-3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(grid_ptz(k) - grid_ptz(k-1))*(grid_ptz(k) - grid_ptz(k-2))/(grid_ptz(k-3) - grid_ptz(k-1))/(grid_ptz(k-3) - grid_ptz(k-2))
            node_arr_csr%coeffs(node_arr_csr%nz_size+4)= -1.d0
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=start_from +(k-1)*(ny+1)+j+1 
            node_arr_csr%col(node_arr_csr%nz_size+2)=start_from +(k-2)*(ny+1)+j+1 
            node_arr_csr%col(node_arr_csr%nz_size+3)=start_from +(k-3)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+4)=start_from +(k)*(ny+1)+j+1    
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+4
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
           end if
           
         else
          side_tr=0
          tet=0
         !---------------------------------------(-1)*b_tau---------------------------------------------
          node_arr_csr%coeffs(node_arr_csr%nz_size+1)=0.5d0
          node_arr_csr%col(node_arr_csr%nz_size+1)=start_from + (k)*(ny+1)+j+1 
          node_arr_csr%nz_size=node_arr_csr%nz_size + 1
          
          !----------------------------------(+H*b_n)-------------------------------------------------
           do me=1,nb
            !-----------------------------------------------------------------
            ! Check and add memory to coeff and col array if required
            !-----------------------------------------------------------------
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
             deallocate(temp_coeffs)
             deallocate(temp_col)
            end if
           end if
          !-----------------------------------------------------------------
       
            if(corner(me,ny,nz)) then
            
             yb=jb(me,ny,nz)
             zb=kb(me,ny,nz)
             if(me==1) then  !corner 
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else if(me==nz+1) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else if(me==2*nz+ny+1) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             end if 
             
             node_arr_csr%nz_size=node_arr_csr%nz_size + 2
             side_tr=side_tr+1
             
             if(tet==1) then
              tet=0
             else
              tet=1
             end if
             
            else
             yb=jb(me,ny,nz)
             zb=kb(me,ny,nz)
             
             if( (me> 1 .and. me < nz+1) .or. (me> 2*nz + ny + 1 .and. me < 2*nz+2*ny +1 )) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             else
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             end if
             
             node_arr_csr%nz_size=node_arr_csr%nz_size + 1
            end if
           end do
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          !temp_coeffs=node_arr_csr%coeffs(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))
          !temp_col=node_arr_csr%col(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))
         
          !call quicksort(temp_coeffs,temp_col,size(temp_col))
         
          !node_arr_csr%coeffs(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))=temp_coeffs
          !node_arr_csr%col(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))=temp_col  
         endif
          
  
  !------------------------------------------------------------------------
  !Prepare coeffs for tangential bz_hat at j=ny 
  !------------------------------------------------------------------------ 
         else if(j==ny) then 
           side_tr1=3  !right face
         !if node is corner i.e j=0 or j=ny
         mp=ibound(j,k,ny,nz)
        
         if(corner(mp,ny,nz)) then
          !extrapolate by from neighbours
           if(k==0) then
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)= -1.d0
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(grid_ptz(k) - grid_ptz(k+2))*(grid_ptz(k) - grid_ptz(k+3))/(grid_ptz(k+1) - grid_ptz(k+2))/(grid_ptz(k+1) - grid_ptz(k+3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(grid_ptz(k) - grid_ptz(k+1))*(grid_ptz(k) - grid_ptz(k+3))/(grid_ptz(k+2) - grid_ptz(k+1))/(grid_ptz(k+2) - grid_ptz(k+3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+4)=(grid_ptz(k) - grid_ptz(k+1))*(grid_ptz(k) - grid_ptz(k+2))/(grid_ptz(k+3) - grid_ptz(k+1))/(grid_ptz(k+3) - grid_ptz(k+2))
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=start_from +(k)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+2)=start_from +(k+1)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+3)=start_from +(k+2)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+4)=start_from +(k+3)*(ny+1)+j+1 
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+4
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
           else
            node_arr_csr%coeffs(node_arr_csr%nz_size+1)=(grid_ptz(k) - grid_ptz(k-2))*(grid_ptz(k) - grid_ptz(k-3))/(grid_ptz(k-1) - grid_ptz(k-2))/(grid_ptz(k-1) - grid_ptz(k-3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+2)=(grid_ptz(k) - grid_ptz(k-1))*(grid_ptz(k) - grid_ptz(k-3))/(grid_ptz(k-2) - grid_ptz(k-1))/(grid_ptz(k-2) - grid_ptz(k-3))   
            node_arr_csr%coeffs(node_arr_csr%nz_size+3)=(grid_ptz(k) - grid_ptz(k-1))*(grid_ptz(k) - grid_ptz(k-2))/(grid_ptz(k-3) - grid_ptz(k-1))/(grid_ptz(k-3) - grid_ptz(k-2))
            node_arr_csr%coeffs(node_arr_csr%nz_size+4)= -1.d0
            
            node_arr_csr%col(node_arr_csr%nz_size+1)=start_from +(k-1)*(ny+1)+j+1 
            node_arr_csr%col(node_arr_csr%nz_size+2)=start_from +(k-2)*(ny+1)+j+1 
            node_arr_csr%col(node_arr_csr%nz_size+3)=start_from +(k-3)*(ny+1)+j+1
            node_arr_csr%col(node_arr_csr%nz_size+4)=start_from +(k)*(ny+1)+j+1    
            
            node_arr_csr%nz_size=node_arr_csr%nz_size+4
            node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
           end if
           
         else
          side_tr=0
          tet=0
          !---------------------------------------(-1)*b_tau---------------------------------------------
          node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-0.5d0
          node_arr_csr%col(node_arr_csr%nz_size+1)=start_from + (k)*(ny+1)+j+1 
          node_arr_csr%nz_size=node_arr_csr%nz_size + 1
          
          !----------------------------------(+H*b_n)-------------------------------------------------
           do me=1,nb
            !-----------------------------------------------------------------
            ! Check and add memory to coeff and col array if required
            !-----------------------------------------------------------------
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
             deallocate(temp_coeffs)
             deallocate(temp_col)
            end if
           end if
           !-----------------------------------------------------------------
       
            if(corner(me,ny,nz)) then
            
             yb=jb(me,ny,nz)
             zb=kb(me,ny,nz)
             if(me==1) then  !corner 
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else if(me==nz+1) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else if(me==2*nz+ny+1) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             else
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             
              node_arr_csr%coeffs(node_arr_csr%nz_size+2)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr+1))
              node_arr_csr%col(node_arr_csr%nz_size+2)=tet*np + (zb)*(ny+1)+yb+1
             end if 
             
             node_arr_csr%nz_size=node_arr_csr%nz_size + 2
             side_tr=side_tr+1
             if(tet==1) then
              tet=0
             else
              tet=1
             end if
            else
             yb=jb(me,ny,nz)
             zb=kb(me,ny,nz)
             
             if( (me> 1 .and. me < nz+1) .or. (me> 2*nz + ny + 1 .and. me < 2*nz+2*ny +1 )) then
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=-bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             else
              node_arr_csr%coeffs(node_arr_csr%nz_size+1)=bo_int_arr_csr%coeffs(bo_int_arr_csr%row((mp+side_tr1)-1)+(me+side_tr))
              node_arr_csr%col(node_arr_csr%nz_size+1)=(1-tet)*np + (zb)*(ny+1)+yb+1
             end if
             
             node_arr_csr%nz_size=node_arr_csr%nz_size + 1
            end if
           end do
          node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
          !temp_coeffs=node_arr_csr%coeffs(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))
          !temp_col=node_arr_csr%col(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))
         
          !call quicksort(temp_coeffs,temp_col,size(temp_col))
         
          !node_arr_csr%coeffs(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))=temp_coeffs
          !node_arr_csr%col(node_arr_csr%row(row_ptr-1)+1:node_arr_csr%row(row_ptr))=temp_col  
         endif
          
  
         
         
        end if
       
                   
       end do
     end do
    
     open(unit=50, file='../output_file/node_info_csr_b33_coeffs.dat')
     open(unit=60, file='../output_file/node_info_csr_b33_cols.dat')
      
      do j=1,ns
    
  
    
       
        
        
        write(50,"(*(5F15.6))"),node_arr_csr%coeffs(node_arr_csr%row(j-1)+1:node_arr_csr%row(j))
        write(60,"(*(1XI6))"),node_arr_csr%col(node_arr_csr%row(j-1)+1:node_arr_csr%row(j))
       
      end do 
  
     
      close(50)
      close(60)
      print*,"Node array filled" 
    
    end subroutine fill_DH_bem_coeff_mat


    subroutine fill_psi_cohfmat(node_arr_csr,ny,nz,nb,alpha_k,h_array,I1,I2)
    use greens_func_subroutines  
    implicit none
    integer,intent(in)                             :: ny,nz,nb 
    real*8,intent(in)                              :: alpha_k
    integer                                        :: i,m,j,k,count,c,chunk_size
    real*8,intent(in)                              :: I1(1:nb,1:nb),I2(1:nb,1:nb)
    real*8,intent(in)                              :: h_array(1:nb)
    type(node_csr),intent(inout)                   :: node_arr_csr
    real(kind=8),allocatable                       :: temp_coeffs(:)
    real,parameter                                 :: PI=4*atan(1.d0)
    integer,allocatable                            :: temp_col(:)
    integer,allocatable                            :: temp_row(:)
    integer                                        :: row_ptr,col_ptr
    real(kind=8)                                   :: cohf_value
    
    count=0
    chunk_size=1000
    node_arr_csr%nz_size=0
    node_arr_csr%capacity=nb                                 !nb=1024
    
   
    row_ptr=0
    col_ptr=0
                                 
    allocate(node_arr_csr%row(0:nb),node_arr_csr%y_loc(0:nb),node_arr_csr%z_loc(0:nb))
    allocate(node_arr_csr%col(node_arr_csr%capacity))
    allocate(node_arr_csr%coeffs(node_arr_csr%capacity))
    
    
   
    node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
    node_arr_csr%y_loc(0)=0.d0
    node_arr_csr%z_loc(0)=0.d0
   
     do i=1,nb
     row_ptr=row_ptr+1
     
     
      do m=1,nb
       
   
       
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
  
       
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
        if(m==i) then
         node_arr_csr%coeffs(node_arr_csr%nz_size+1)=beta_const(i,ny,nz)*(2*PI) 
         node_arr_csr%col(node_arr_csr%nz_size+1)=m
         node_arr_csr%nz_size=node_arr_csr%nz_size+1
        else
         cohf_value=-1.d0*(I2(i,m_indx(m-1,nb))- I2(i,m_indx(m,nb))+ I1(i,m_indx(m,nb))) 
         if(abs(cohf_value)>1e-6) then
           node_arr_csr%coeffs(node_arr_csr%nz_size+1)=cohf_value  
           node_arr_csr%col(node_arr_csr%nz_size+1)=m
           node_arr_csr%nz_size=node_arr_csr%nz_size+1
          end if
        end if
          
          
          
        
        
             
  
            
 
          
          
                   
       end do
       node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
       
     end do
     
  
     open(unit=50, file='../output_file/node_info_csr_b13_coeffs.dat')
     open(unit=60, file='../output_file/node_info_csr_b13_cols.dat')
      write(30,*),"ap       an       as        aw       ae"
 
      do j=1,nb
    
  
    
   
        
        
        write(50,"(*(F15.6))"),node_arr_csr%coeffs(node_arr_csr%row(j-1)+1:node_arr_csr%row(j))
        write(60,"(*(5I))"),node_arr_csr%col(node_arr_csr%row(j-1)+1:node_arr_csr%row(j))
       
      end do 
  
     
      close(50)
      close(60)
      print*,"Node array filled" 
    
    end subroutine fill_psi_cohfmat
    
    subroutine fill_psi_cohfmat_non_sparse(node_arr_csr,ny,nz,nb,alpha_k,h_array,I1,I2)
    use greens_func_subroutines  
    implicit none
    integer,intent(in)                             :: ny,nz,nb 
    real*8,intent(in)                              :: alpha_k
    integer                                        :: i,m,j,k,count,c,chunk_size
    real*8,intent(in)                              :: I1(1:nb,1:nb),I2(1:nb,1:nb)
    real*8,intent(in)                              :: h_array(1:nb)
    type(node_csr),intent(inout)                   :: node_arr_csr
    real(kind=8),allocatable                       :: temp_coeffs(:)
    real,parameter                                 :: PI=4*atan(1.d0)
    integer,allocatable                            :: temp_col(:)
    integer,allocatable                            :: temp_row(:)
    integer                                        :: row_ptr,col_ptr
    real(kind=8)                                   :: cohf_value
    
    count=0
    chunk_size=1000
    node_arr_csr%nz_size=0
    node_arr_csr%capacity=nb                                 !nb=1024
    
   
    row_ptr=0
    col_ptr=0
                                 
    allocate(node_arr_csr%row(0:nb),node_arr_csr%y_loc(0:nb),node_arr_csr%z_loc(0:nb))
    allocate(node_arr_csr%col(node_arr_csr%capacity))
    allocate(node_arr_csr%coeffs(node_arr_csr%capacity))
    
    
   
    node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
    node_arr_csr%y_loc(0)=0.d0
    node_arr_csr%z_loc(0)=0.d0
   
     do i=1,nb
     row_ptr=row_ptr+1
     
     
      do m=1,nb
       
   
       
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
  
       
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
        if(m==i) then
         node_arr_csr%coeffs(node_arr_csr%nz_size+1)=beta_const(i,ny,nz)*(2*PI) 
         node_arr_csr%col(node_arr_csr%nz_size+1)=m
         node_arr_csr%nz_size=node_arr_csr%nz_size+1
        else
         cohf_value=-1.d0*(I2(i,m_indx(m-1,nb))- I2(i,m_indx(m,nb))+ I1(i,m_indx(m,nb))) 
      
           node_arr_csr%coeffs(node_arr_csr%nz_size+1)=cohf_value  
           node_arr_csr%col(node_arr_csr%nz_size+1)=m
           node_arr_csr%nz_size=node_arr_csr%nz_size+1
      
        end if
          
          
          
        
        
             
  
            
 
          
          
                   
       end do
       node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
       
     end do
     
  
     
  
     
      close(50)
      close(60)
      print*,"Node array filled" 
    
    end subroutine fill_psi_cohfmat_non_sparse    
    
    subroutine fill_psi_cohf_E(node_arr_csr,ny,nz,nb,alpha_k,h_array,I1,I2)
    use greens_func_subroutines  
    implicit none
    integer,intent(in)                             :: ny,nz,nb 
    real*8,intent(in)                              :: alpha_k
    integer                                        :: i,m,j,k,count,c,chunk_size
    real*8,allocatable,intent(in)                  :: I1(:,:),I2(:,:)
    real*8,intent(in)                              :: h_array(1:nb)
    type(node_csr),intent(inout)                   :: node_arr_csr
    real(kind=8),allocatable                       :: temp_coeffs(:)
    real,parameter                                 :: PI=4*atan(1.d0)
    integer,allocatable                            :: temp_col(:)
    integer,allocatable                            :: temp_row(:)
    integer                                        :: row_ptr,col_ptr
    real(kind=8)                                   :: cohf_value
    integer                                        :: n_size
    integer*4                                      :: side_tr
    logical                                        :: flag                                     
    
    count=0
    chunk_size=1000
    node_arr_csr%nz_size=0
    node_arr_csr%capacity=nb                                 !nb=1024
    
    n_size=2*(ny+1)+2*(nz+1)
    row_ptr=0
    col_ptr=0
    side_tr=0
    flag=.true.                             
    allocate(node_arr_csr%row(0:n_size),node_arr_csr%y_loc(0:n_size),node_arr_csr%z_loc(0:n_size))
    allocate(node_arr_csr%col(node_arr_csr%capacity))
    allocate(node_arr_csr%coeffs(node_arr_csr%capacity))
    
    
   
    node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
    node_arr_csr%y_loc(0)=0.d0
    node_arr_csr%z_loc(0)=0.d0
   
     do i=1,nb
     
      
     
      if(corner(i,ny,nz)) then
      !-------------------------backward tangential at corner------------------------------
       row_ptr=row_ptr+1 
       do m=1,nb
       
   
       
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
  
       
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
                                     
         cohf_value=I2(i+side_tr,m_indx(m-1,nb))- I2(i+side_tr,m_indx(m,nb))+ I1(i+side_tr,m_indx(m,nb)) 
         
       
         if(m==i) then 
          
          cohf_value=I2(i+side_tr,m_indx(m-1,nb)) + I1(i+side_tr,m_indx(m,nb)) 
          
         else if(m==m_indx(i+1,nb)) then
          
          cohf_value=I2(i+side_tr,m_indx(m-1,nb))- I2(i+side_tr,m_indx(m,nb))+ I1(i+side_tr,m_indx(m,nb)) 
         
         end if
         
         
        if(abs(cohf_value)>1e-6) then
           node_arr_csr%coeffs(node_arr_csr%nz_size+1)=cohf_value  
           node_arr_csr%col(node_arr_csr%nz_size+1)=m
           node_arr_csr%nz_size=node_arr_csr%nz_size+1
          end if  
                   
       end do
       
       
       node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
       
      !---------------------------------------------------------------------------------------
      !--------------------------forward tangential at corner---------------------------------
      row_ptr=row_ptr+1
      do m=1,nb
       
   
       
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
  
       
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
                                     
         cohf_value=I2(i+side_tr+1,m_indx(m-1,nb))- I2(i+side_tr+1,m_indx(m,nb))+ I1(i+side_tr+1,m_indx(m,nb)) 
         
       
         if(m==i) then 
          
          cohf_value=I2(i+side_tr+1,m_indx(m-1,nb))- I2(i+side_tr+1,m_indx(m,nb))+ I1(i+side_tr+1,m_indx(m,nb)) 
         
         else if(m==m_indx(i-1,nb)) then
          
          cohf_value=I2(i+side_tr+1,m_indx(m-1,nb)) + I1(i+side_tr+1,m_indx(m,nb)) 
         
         end if
         
         
        if(abs(cohf_value)>1e-6) then
           node_arr_csr%coeffs(node_arr_csr%nz_size+1)=cohf_value  
           node_arr_csr%col(node_arr_csr%nz_size+1)=m
           node_arr_csr%nz_size=node_arr_csr%nz_size+1
          end if  
                   
       end do
       
       
       node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
    
      !---------------------------------------------------------------------------------------
       side_tr= side_tr + 1
      !-----------------------------------non-corner elements---------------------------------
      else
       row_ptr=row_ptr+1
       do m=1,nb
       
   
       
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
  
       
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
         cohf_value=I2(i+side_tr,m_indx(m-1,nb))- I2(i+side_tr,m_indx(m,nb))+ I1(i+side_tr,m_indx(m,nb)) 
         if(abs(cohf_value)>1e-6) then
           node_arr_csr%coeffs(node_arr_csr%nz_size+1)=cohf_value  
           node_arr_csr%col(node_arr_csr%nz_size+1)=m
           node_arr_csr%nz_size=node_arr_csr%nz_size+1
          end if
       
          
   
          
          
                   
       end do
       node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
      !--------------------------------------------------------------------------------------------
      end if
 
     end do
     
  
   
    
    end subroutine fill_psi_cohf_E
    
    subroutine fill_psi_cohf_E_non_sparse(node_arr_csr,ny,nz,nb,alpha_k,h_array,I1,I2)
    use greens_func_subroutines  
    implicit none
    integer,intent(in)                             :: ny,nz,nb 
    real*8,intent(in)                              :: alpha_k
    integer                                        :: i,m,j,k,count,c,chunk_size
    real*8,allocatable,intent(in)                  :: I1(:,:),I2(:,:)
    real*8,intent(in)                              :: h_array(1:nb)
    type(node_csr),intent(inout)                   :: node_arr_csr
    real(kind=8),allocatable                       :: temp_coeffs(:)
    real,parameter                                 :: PI=4*atan(1.d0)
    integer,allocatable                            :: temp_col(:)
    integer,allocatable                            :: temp_row(:)
    integer                                        :: row_ptr,col_ptr
    real(kind=8)                                   :: cohf_value
    integer                                        :: n_size
    
    
    
    count=0
    chunk_size=1000
    node_arr_csr%nz_size=0
    node_arr_csr%capacity=nb                                 !nb=1024
    
    n_size=2*(ny+1)+2*(nz+1)
    row_ptr=0
    col_ptr=0
                                 
    allocate(node_arr_csr%row(0:n_size),node_arr_csr%y_loc(0:n_size),node_arr_csr%z_loc(0:n_size))
    allocate(node_arr_csr%col(node_arr_csr%capacity))
    allocate(node_arr_csr%coeffs(node_arr_csr%capacity))
    
    
   
    node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
    node_arr_csr%y_loc(0)=0.d0
    node_arr_csr%z_loc(0)=0.d0
   
     do i=1,n_size
     
    
      row_ptr=row_ptr+1
      do m=1,nb
       
   
       
       
       !print*,i,node_arr_csr%nz_size ,node_arr_csr%capacity
  
       
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
         cohf_value=I2(i,m_indx(m-1,nb))- I2(i,m_indx(m,nb))+ I1(i,m_indx(m,nb)) 
        
           node_arr_csr%coeffs(node_arr_csr%nz_size+1)=cohf_value  
           node_arr_csr%col(node_arr_csr%nz_size+1)=m
           node_arr_csr%nz_size=node_arr_csr%nz_size+1
      
       
          
   
          
          
                   
       end do
       node_arr_csr%row(row_ptr)=node_arr_csr%nz_size
      
     end do
     
  
   
    
    end subroutine fill_psi_cohf_E_non_sparse
    
    subroutine fill_rhs_array(node_arr_csr,q_array,f,alpha_k) !for double helmholtz + pseudo vacuum bound
       type(node_csr),intent(inout)                     :: node_arr_csr
       real(kind=8),allocatable, intent(out)            :: q_array(:)
       real*8, intent(in)                               ::f,alpha_k
       real,parameter                                   :: PI=4*atan(1.d0)
       integer                                          :: i
       real                                             :: y,z,constant
       
       constant=0.01*((f+alpha_k**2)+8*(PI)**2)
       allocate(q_array(1:size(node_arr_csr%row)-1))
       open(unit=30, file='../output_file/rhs_y_b10.dat')
       open(unit=40, file='../output_file/rhs_z_b10.dat')
       do i=1,(size(node_arr_csr%row)-1)/2
         y=node_arr_csr%y_loc(i)
         z=node_arr_csr%z_loc(i)
         q_array(i)=-constant*cos(2*PI*y)*sin(2*PI*z)
         if(abs(z)>0.9999999 .or. abs(y)>0.9999999)then
          q_array(i)=0.d0
         end if
         write(30,*),y,z,q_array(i)
       end do
       
        do i=(size(node_arr_csr%row)-1)/2 + 1,size(node_arr_csr%row)-1
         y=node_arr_csr%y_loc(i)
         z=node_arr_csr%z_loc(i)
         q_array(i)=-constant*sin(2*PI*y)*cos(2*PI*z)
         if(abs(z)>0.9999999 .or. abs(y)>0.9999999)then
          q_array(i)=0.d0
         end if
         write(40,*),y,z,q_array(i)
       end do
       close(30)
       close(40)
       
        
    end subroutine fill_rhs_array
     
    subroutine fill_act_b_array(node_arr_csr,b_array,f,alpha_k) !test against del_K0= -b to verify the  
       type(node_csr),intent(inout)                     :: node_arr_csr
       real(kind=8),allocatable, intent(out)            :: b_array(:)
       real*8, intent(in)                               ::f,alpha_k
       real*8,parameter                                 :: PI=4*atan(1.d0)
       integer                                          :: i,t
       real*8                                           :: y,z,constant,r
       integer                                          :: n
       real*8                                           ::E0(1:6), C0(1:7),D0(1:7)
       real*8                                           ::I0, I0_dash, I0_double_dash
       real*8                                           ::G0, G0_dash, G0_double_dash
       E0(1:6) = (/3.5156229,3.0899424,1.2067492,0.2659732,0.360768e-1,0.45813e-2/)
       C0(1:7)=(/-0.57721566,0.42278420,0.23069756,0.348859e-1,0.262698e-2,0.10750e-3,0.74e-5/)
       D0(1:7)= (/1.25331414, -0.7832358e-1, 0.2189568e-1, -0.1062446e-1, 0.587872e-2, -0.251540e-2, 0.53208e-3/)
                                           !             size of by     size of bz
       n=size(node_arr_csr%row)-1         !Total size (ny+1)*(nz+1) + (ny+1)*(nz+1)
       allocate(b_array(1:size(node_arr_csr%row)-1))
       open(unit=30, file='../output_file/b_act_b34.dat')
       !----------------------by------------------------------------------
       do i=1,n/2
         y=node_arr_csr%y_loc(i) + 0.d0
         z=node_arr_csr%z_loc(i) + 0.d0
         r= sqrt(y**2+z**2)
         
         if(r*alpha_k<=2.d0) then
         b_array(i)= -1/r
         do t=1,6
            I0=  E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t)
            I0_dash= E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t-1)*(2*t)
           
            b_array(i)=b_array(i) - I0/r - I0_dash*log(alpha_k*r/2)                
         end do
              
         do t=1,7
              
            b_array(i)=b_array(i) + C0(t)*(2*t-2)*r**(2*t-3)*(alpha_k)**(2*t-2)*0.25**(t-1)
               
         end do
         
      
         if(y==-1.d0 + 0.d0 .or. y==1.d0 + 0.d0 .or. z==-1.d0 + 0.d0 .or. z==1.d0 + 0.d0) then
            b_array(i)= -b_array(i)*y/r
         else
          b_array(i)=1.d0
         end if
         
         else
          b_array(i)=0.d0
          
          do t=1,7
           G0= D0(t)*2**(t-1)/alpha**(t-0.5)*(1/r)**(t-0.5)
           G0_dash= -D0(t)*(t-0.5)*2**(t-1)*(1/alpha_k)**(t-0.5)*(1/r)**(t+0.5)
           
           K0 = K0 + exp(-alpha_k*r)*G0
           K0_dash = K0_dash -alpha_k*exp(-alpha_k*r)*G0 + exp(-alpha_k*r)*G0_dash
          end do
           b_array(i)= 1.d0
           if(y==-1.d0 + 0.d0 .or. y==1.d0 + 0.d0 .or. z==-1.d0 + 0.d0 .or. z==1.d0 + 0.d0 ) then
            b_array(i)= -K0_dash*y/r
           end if
         end if
         write(30,*),y,z,b_array(i)
       end do
   !------------------------------bz-----------------------------------    
       do i=1,n/2
         y=node_arr_csr%y_loc(n/2 + i) + 0.d0
         z=node_arr_csr%z_loc(n/2 + i) + 0.d0
         r= sqrt(y**2+z**2)
         
         if(r*alpha_k<=2.d0) then
         b_array(n/2 + i)= -1/r
         do t=1,6
            I0=  E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t)
            I0_dash= E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t-1)*(2*t)
           
            b_array(n/2 + i)=b_array(n/2 + i) - I0/r - I0_dash*log(alpha_k*r/2)                
         end do
              
         do t=1,7
              
            b_array(n/2 + i)=b_array(n/2 + i) + C0(t)*(2*t-2)*r**(2*t-3)*(alpha_k)**(2*t-2)*0.25**(t-1)
               
         end do
         
          if(y==-1.d0 + 0.d0 .or. y==1.d0 + 0.d0 .or. z==-1.d0 + 0.d0 .or. z==1.d0 + 0.d0) then
            b_array(n/2+i)= -b_array(n/2+i)*z/r
         else
          b_array(n/2+i)=1.d0
         end if
         
         else
          b_array(n/2 + i)=0.d0
          do t=1,7
           G0= D0(t)*2**(t-1)/alpha**(t-0.5)*(1/r)**(t-0.5)
           G0_dash= -D0(t)*(t-0.5)*2**(t-1)*(1/alpha_k)**(t-0.5)*(1/r)**(t+0.5)
           
           K0 = K0 + exp(-alpha_k*r)*G0
           K0_dash = K0_dash -alpha_k*exp(-alpha_k*r)*G0 + exp(-alpha_k*r)*G0_dash
          end do
           b_array(n/2 + i)= 1.d0
           if(y==-1.d0 + 0.d0 .or. y==1.d0 + 0.d0 .or. z==-1.d0 + 0.d0 .or. z==1.d0 + 0.d0) then
            b_array(n/2 + i)= -K0_dash*z/r
           end if
         end if
         write(30,*),y,z,b_array(n/2 + i)
       end do
       close(30)
       
        
    end subroutine fill_act_b_array
    
    subroutine fill_q_array(node_arr_csr,q_array,f,alpha_k)
       type(node_csr),intent(inout)                     :: node_arr_csr
       real(kind=8),allocatable, intent(out)            :: q_array(:)
       real*8, intent(in)                               :: f,alpha_k
       real*8,parameter                                 :: PI=4*atan(1.d0)
       integer                                          :: i,t
       real*8                                           :: y,z,constant,r
       integer                                          :: n
       real*8                                           :: E0(1:6), C0(1:7), D0(1:7)
       real*8                                           :: I0, I0_dash, I0_double_dash,I0_trip_dash
       real*8                                           :: P0, P0_dash, P0_double_dash,P0_trip_dash
       real*8                                           :: G0, G0_dash, G0_double_dash,G0_trip_dash
       real*8                                           :: K0, K0_dash, K0_double_dash,K0_trip_dash
       real*8                                           :: K0_tder_y, K0_tder_z, K0_dder_z_der_y, K0_dder_y_der_z
       
       ! tder_y -> triple derivative w.r.t y  
       ! tder_z -> triple derivative w.r.t z   
       ! dder_z_der_y -> double derivative w.r.t z derivative w.r.t y   
       ! dder_y_der_z -> double derivative w.r.t y derivative w.r.t z   
       
                                           
       E0(1:6) = (/3.5156229,3.0899424,1.2067492,0.2659732,0.360768e-1,0.45813e-2/)
       C0(1:7)=(/-0.57721566,0.42278420,0.23069756,0.348859e-1,0.262698e-2,0.10750e-3,0.74e-5/)
       D0(1:7)= (/1.25331414, -0.7832358e-1, 0.2189568e-1, -0.1062446e-1, 0.587872e-2, -0.251540e-2, 0.53208e-3/)
       
                                           !             size of by     size of bz
       n=size(node_arr_csr%row)-1         !Total size (ny+1)*(nz+1) + (ny+1)*(nz+1)
       allocate(q_array(1:size(node_arr_csr%row)-1))
       open(unit=30, file='../output_file/rhs_b34.dat')
       
       !*************** for b_y *****************************!
       do i=1,n/2
         y=node_arr_csr%y_loc(i) + 0.d0
         z=node_arr_csr%z_loc(i) + 0.d0
         
         
        if(z==-1.d0 + 0.d0 .or. z== 1.d0 + 0.d0) then
          q_array(i)= 0.d0
         else
         
         
         r= sqrt(y**2+z**2)
         
          if(alpha_k*r<=2.d0) then 
           !initialise with coefficient of '1+' term of I0
           !---------------------------
           K0= -log(alpha_k*r/2.0)
           K0_dash = -1/r
           K0_double_dash = (1/r)**2
           K0_trip_dash   = -2*(1/r)**3
           !---------------------------
           do t=1,6
            I0=  E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t)
            I0_dash= E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t-1)*(2*t)
            I0_double_dash= E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t-2)*(2*t)*(2*t-1)
            I0_trip_dash= E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t-3)*(2*t)*(2*t-1)*(2*t-2)
            
           
            K0 = K0 - log(alpha_k*r/2)*I0
            K0_dash = K0_dash -I0/r - log(alpha_k*r/2)*I0_dash  
            K0_double_dash = K0_double_dash + I0*(1/r)**2 - 2*I0_dash/r - log(alpha_k*r/2)*I0_double_dash
            K0_trip_dash   = K0_trip_dash  - 2*I0*(1/r)**3  + 3*I0_dash *(1/r)**2 - 3*I0_double_dash/r - log(alpha_k*r/2)*I0_trip_dash 
            
                   
           end do
           
           do t=1,7
            P0=C0(t)*r**(2*t-2)*(alpha_k)**(2*t-2)*0.25**(t-1)
            P0_dash=C0(t)*(2*t-2)*r**(2*t-3)*(alpha_k)**(2*t-2)*0.25**(t-1)  
            P0_double_dash=C0(t)*(2*t-2)*(2*t-3)*r**(2*t-4)*(alpha_k)**(2*t-2)*0.25**(t-1)
            P0_trip_dash=C0(t)*(2*t-2)*(2*t-3)*(2*t-4)*r**(2*t-5)*(alpha_k)**(2*t-2)*0.25**(t-1)
       
            K0 = K0 + P0
            K0_dash = K0_dash + P0_dash 
            K0_double_dash = K0_double_dash + P0_double_dash
            K0_trip_dash   = K0_trip_dash + P0_trip_dash
           end do
         
           K0_tder_y = K0_trip_dash*(y/r)**3 +   3.d0*K0_double_dash*(y*(1/r)**2 - y**3*(1/r)**4) + K0_dash*(-3*y*(1/r)**3 + 3.d0*y**3*(1/r)**5)               
           K0_dder_z_der_y = K0_trip_dash*(z**2*y*(1/r)**3) +   K0_double_dash*((y*(1/r)**2 - 3.d0*y*z**2*(1/r)**4) ) + K0_dash*(-y*(1/r)**3 + 3.d0*z**2*y*(1/r)**5) 
           !q_array(i)= -(-(f+alpha_k**2)*K0_dash*y/r + K0_tder_y + K0_dder_z_der_y)
            q_array(i)= 1.d0
           if(y==-1.d0 + 0.d0 .or. y==1.d0 + 0.d0 ) then
            q_array(i)= -K0_dash*y/r
           end if
          else
           
           K0 = 0.d0
           K0_dash = 0.d0  
           K0_double_dash = 0.d0
           K0_trip_dash   = 0.d0
            
           do t=1,7
            G0= D0(t)*2**(t-1)/alpha**(t-0.5)*(1/r)**(t-0.5)
            G0_dash= -D0(t)*(t-0.5)*2**(t-1)*(1/alpha_k)**(t-0.5)*(1/r)**(t+0.5)
            G0_double_dash= D0(t)*(t-0.5)*(t+0.5)*2**(t-1)*(1/alpha_k)**(t-0.5)*(1/r)**(t+1.5)
            G0_trip_dash= -D0(t)*(t-0.5)*(t+0.5)*(t+1.5)*2**(t-1)*(1/alpha_k)**(t-0.5)*(1/r)**(t+2.5)
            
           
            K0 = K0 + exp(-alpha_k*r)*G0
            K0_dash = K0_dash -alpha_k*exp(-alpha_k*r)*G0 + exp(-alpha_k*r)*G0_dash
            K0_double_dash = K0_double_dash + alpha_k**2*exp(-alpha_k*r)*G0 - 2*alpha_k*exp(-alpha_k*r)*G0_dash + exp(-alpha_k*r)*G0_double_dash
            K0_trip_dash   = K0_trip_dash - alpha_k**3*exp(-alpha_k*r)*G0 + 3*alpha_k*2*exp(-alpha_k*r)*G0_dash - 3*alpha_k*exp(-alpha_k*r)*G0_double_dash + exp(-alpha_k*r)*G0_trip_dash
            
                   
           end do
           
                    
           K0_tder_y = K0_trip_dash*(y/r)**3 +   3.d0*K0_double_dash*(y*(1/r)**2 - y**3*(1/r)**4) + K0_dash*(-3*y*(1/r)**3 + 3.d0*y**3*(1/r)**5)               
           K0_dder_z_der_y = K0_trip_dash*(z**2*y*(1/r)**3) +   K0_double_dash*((y*(1/r)**2 - 3.d0*y*z**2*(1/r)**4) ) + K0_dash*(-y*(1/r)**3 + 3.d0*z**2*y*(1/r)**5) 
           !q_array(i)= -(-(f+alpha_k**2)*K0_dash*y/r + K0_tder_y + K0_dder_z_der_y)
           q_array(i)= 1.d0
           if(y==-1.d0 + 0.d0 .or. y==1.d0 + 0.d0 ) then
            q_array(i)= -K0_dash*y/r
           end if
          end if
         end if 
         write(30,*),y,z,q_array(i)
       end do
       !*************** for b_z *****************************! 
       do i=1,n/2
         y=node_arr_csr%y_loc(n/2 + i) + 0.d0
         z=node_arr_csr%z_loc(n/2 + i) + 0.d0
         r= sqrt(y**2+z**2)
         
         if(y==-1.d0 + 0.d0 .or. y== 1.d0 + 0.d0 ) then
          q_array(n/2 + i)= 0.d0
         else
          r= sqrt(y**2+z**2)
          
          if(alpha_k*r<=2.d0) then  
           !initialise with coefficient of '1+' term of I0
           !---------------------------
           K0= -log(alpha_k*r/2)
           K0_dash = -1/r
           K0_double_dash = (1/r)**2
           K0_trip_dash   = -2*(1/r)**3
           !---------------------------
          do t=1,6
            I0=  E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t)
            I0_dash= E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t-1)*(2*t)
            I0_double_dash= E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t-2)*(2*t)*(2*t-1)
            I0_trip_dash= E0(t)*(alpha_k/3.75)**(2*t)*r**(2*t-3)*(2*t)*(2*t-1)*(2*t-2)
            
            
            K0 = K0 - log(alpha_k*r/2)*I0
            K0_dash = K0_dash -I0/r - log(alpha_k*r/2)*I0_dash  
            K0_double_dash = K0_double_dash + I0*(1/r)**2 - 2*I0_dash/r - log(alpha_k*r/2)*I0_double_dash
            K0_trip_dash   = K0_trip_dash  - 2*I0*(1/r)**3  + 3*I0_dash *(1/r)**2 - 3*I0_double_dash/r - log(alpha_k*r/2)*I0_trip_dash                       
          end do
              
          do t=1,7
              
            P0=C0(t)*r**(2*t-2)*(alpha_k)**(2*t-2)*0.25**(t-1)
            P0_dash=C0(t)*(2*t-2)*r**(2*t-3)*(alpha_k)**(2*t-2)*0.25**(t-1)  
            P0_double_dash=C0(t)*(2*t-2)*(2*t-3)*r**(2*t-4)*(alpha_k)**(2*t-2)*0.25**(t-1)
            P0_trip_dash=C0(t)*(2*t-2)*(2*t-3)*(2*t-4)*r**(2*t-5)*(alpha_k)**(2*t-2)*0.25**(t-1)
       
            K0 = K0 + P0
            K0_dash = K0_dash + P0_dash 
            K0_double_dash = K0_double_dash + P0_double_dash
            K0_trip_dash   = K0_trip_dash + P0_trip_dash  
          end do
         
          K0_tder_z = K0_trip_dash*(z/r)**3 +   3.d0*K0_double_dash*(z*(1/r)**2 - z**3*(1/r)**4) + K0_dash*(-3*z*(1/r)**3 + 3.d0*z**3*(1/r)**5)               
          K0_dder_y_der_z = K0_trip_dash*(y**2*z*(1/r)**3) +   K0_double_dash*((z*(1/r)**2 - 3.d0*z*y**2*(1/r)**4) ) + K0_dash*(-z*(1/r)**3 + 3.d0*y**2*z*(1/r)**5) 
          !q_array(n/2 + i)= -(-(f+alpha_k**2)*K0_dash*z/r + K0_tder_z + K0_dder_y_der_z)
           q_array(n/2 + i)= 1.d0
           if(z==-1.d0 + 0.d0 .or. z==1.d0 + 0.d0 ) then
            q_array(n/2 + i)= -K0_dash*z/r
           end if
          else
           K0 = 0.d0
           K0_dash = 0.d0  
           K0_double_dash = 0.d0
           K0_trip_dash   = 0.d0
            
           do t=1,7
            G0= D0(t)*2**(t-1)/alpha**(t-0.5)*(1/r)**(t-0.5)
            G0_dash= -D0(t)*(t-0.5)*2**(t-1)*(1/alpha_k)**(t-0.5)*(1/r)**(t+0.5)
            G0_double_dash= D0(t)*(t-0.5)*(t+0.5)*2**(t-1)*(1/alpha_k)**(t-0.5)*(1/r)**(t+1.5)
            G0_trip_dash= -D0(t)*(t-0.5)*(t+0.5)*(t+1.5)*2**(t-1)*(1/alpha_k)**(t-0.5)*(1/r)**(t+2.5)
            
           
            K0 = K0 + exp(-alpha_k*r)*G0
            K0_dash = K0_dash -alpha_k*exp(-alpha_k*r)*G0 + exp(-alpha_k*r)*G0_dash
            K0_double_dash = K0_double_dash + alpha_k**2*exp(-alpha_k*r)*G0 - 2*alpha_k*exp(-alpha_k*r)*G0_dash + exp(-alpha_k*r)*G0_double_dash
            K0_trip_dash   = K0_trip_dash - alpha_k**3*exp(-alpha_k*r)*G0 + 3*alpha_k*2*exp(-alpha_k*r)*G0_dash - 3*alpha_k*exp(-alpha_k*r)*G0_double_dash + exp(-alpha_k*r)*G0_trip_dash
            
                   
           end do
           
                    
           K0_tder_z = K0_trip_dash*(z/r)**3 +   3.d0*K0_double_dash*(z*(1/r)**2 - z**3*(1/r)**4) + K0_dash*(-3*z*(1/r)**3 + 3.d0*z**3*(1/r)**5)               
           K0_dder_y_der_z = K0_trip_dash*(y**2*z*(1/r)**3) +   K0_double_dash*((z*(1/r)**2 - 3.d0*z*y**2*(1/r)**4) ) + K0_dash*(-z*(1/r)**3 + 3.d0*y**2*z*(1/r)**5) 
           !q_array(n/2 + i)= -(-(f+alpha_k**2)*K0_dash*z/r + K0_tder_z + K0_dder_y_der_z)
           q_array(n/2 + i)= 1.d0
           
            if(z==-1.d0 + 0.d0 .or. z==1.d0 + 0.d0 ) then
             q_array(n/2 + i)= -K0_dash*z/r
            end if
          end if
         end if
         write(30,*),y,z,q_array(n/2 + i)
       end do
     
       close(30)
       
        
    end subroutine fill_q_array
    
    subroutine fill_bn_cohf_D(D_csr,alpha_k,ny,nz,nb,h_array,I1,I2)
    use greens_func_subroutines
    use array_processes 
    implicit none
       integer,intent(in)                               ::ny,nz,nb
       type(node_csr),intent(out)                       ::D_csr
       real(kind=8),intent(in)                          ::h_array(1:nb)
       real(kind=8),intent(in)                          ::I1(1:nb,1:nb),I2(1:nb,1:nb)   
       real*8, intent(in)                               ::alpha_k
       real,parameter                                   ::PI=4*atan(1.d0)
       integer                                          ::i,m,j,k,count,c,chunk_size
       real                                             ::y,z,constant
       real(kind=8),allocatable                         ::temp_coeffs(:)
       integer,allocatable                              ::temp_col(:)
       integer,allocatable                              ::temp_row(:)
       integer                                          ::row_ptr,col_ptr,side, n_size
       
       
       chunk_size=1000
       D_csr%nz_size=0
       D_csr%capacity=nb                                 !nb=1024
    
       n_size=2*(ny+1)+2*(nz+1)                          !n_size=1028, no. of unique b_normals
       row_ptr=0
       col_ptr=0
                                 
       allocate(D_csr%row(0:nb),D_csr%y_loc(0:nb),D_csr%z_loc(0:nb))
       allocate(D_csr%col(D_csr%capacity))
       allocate(D_csr%coeffs(D_csr%capacity))
    
    
   
       D_csr%row(row_ptr)=D_csr%nz_size
       D_csr%y_loc(0)=0.d0
       D_csr%z_loc(0)=0.d0
       
        
    
       open(unit=10, file='../output_file/rhs_b13.dat')
      
       
       do i=1,nb
         
         row_ptr=row_ptr+1 
         
        
         
         side=0
         do j=1,nb
         
          if((D_csr%capacity-D_csr%nz_size)<(ny+2))then
           if(D_csr%capacity==0)then
             allocate(D_csr%coeffs(chunk_size))
             allocate(D_csr%col(chunk_size))
            
           else
            temp_coeffs=D_csr%coeffs(1:D_csr%nz_size)
            temp_col   =D_csr%col(1:D_csr%nz_size)
            if(allocated(D_csr%coeffs)) deallocate(D_csr%coeffs)
            if(allocated(D_csr%col))   deallocate(D_csr%col)
            allocate(D_csr%coeffs(D_csr%capacity+chunk_size))
            allocate(D_csr%col(D_csr%capacity+chunk_size))
            D_csr%capacity=D_csr%capacity+chunk_size
            D_csr%coeffs(1:D_csr%nz_size)=temp_coeffs
            D_csr%col(1:D_csr%nz_size)   =temp_col
          
           end if
          end if  
       
          if( j== i )then
            
             
                       
            
            
             
             
             if(corner(j,ny,nz)) then
              
              D_csr%coeffs(D_csr%nz_size+1)=I1(i,m_indx(j-1,nb))-(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))
              D_csr%col(D_csr%nz_size+1)=m_indx(j,nb) + side
              D_csr%coeffs(D_csr%nz_size+2)=I1(i,m_indx(j,nb))-(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              D_csr%col(D_csr%nz_size+2)=m_indx(m_indx(j,nb) + side+1,n_size)
              D_csr%coeffs(D_csr%nz_size+3)=(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))+(I1(i,m_indx(j+1,nb))- I2(i,m_indx(j+1,nb)))
              D_csr%col(D_csr%nz_size+3)=m_indx(m_indx(j,nb) + side + 2,n_size)
              D_csr%coeffs(D_csr%nz_size+4)=(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I2(i,m_indx(j-2,nb))
              D_csr%col(D_csr%nz_size+4)=m_indx(m_indx(j,nb) + side - 1,n_size)
              D_csr%nz_size=D_csr%nz_size+4
              side=side+1
             !newly added condition 11/3/2025
             elseif(corner(m_indx(j+1,nb),ny,nz)) then
              D_csr%coeffs(D_csr%nz_size+1)=I1(i,m_indx(j-1,nb))-(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I1(i,m_indx(j,nb))-(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              D_csr%col(D_csr%nz_size+1)=m_indx(j,nb) + side
              D_csr%coeffs(D_csr%nz_size+2)=(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              D_csr%col(D_csr%nz_size+2)=m_indx(m_indx(j,nb) + side + 1,n_size)
              D_csr%coeffs(D_csr%nz_size+3)=(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I2(i,m_indx(j-2,nb))
              D_csr%col(D_csr%nz_size+3)=m_indx(m_indx(j,nb) + side - 1,n_size)
              D_csr%nz_size=D_csr%nz_size+3
             elseif(corner(m_indx(j-1,nb),ny,nz)) then
              D_csr%coeffs(D_csr%nz_size+1)=I1(i,m_indx(j-1,nb))-(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I1(i,m_indx(j,nb))-(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              D_csr%col(D_csr%nz_size+1)=m_indx(j,nb) + side
              D_csr%coeffs(D_csr%nz_size+2)=(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))+(I1(i,m_indx(j+1,nb))- I2(i,m_indx(j+1,nb)))
              D_csr%col(D_csr%nz_size+2)=m_indx(m_indx(j,nb) + side + 1,n_size)
              D_csr%coeffs(D_csr%nz_size+3)=(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))
              D_csr%col(D_csr%nz_size+3)=m_indx(m_indx(j,nb) + side - 1,n_size)
              D_csr%nz_size=D_csr%nz_size+3
             !newly added
             else
              D_csr%coeffs(D_csr%nz_size+1)=I1(i,m_indx(j-1,nb))-(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I1(i,m_indx(j,nb))-(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              D_csr%col(D_csr%nz_size+1)=m_indx(j,nb) + side
              D_csr%coeffs(D_csr%nz_size+2)=(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))+(I1(i,m_indx(j+1,nb))- I2(i,m_indx(j+1,nb)))
              D_csr%col(D_csr%nz_size+2)=m_indx(m_indx(j,nb) + side + 1,n_size)
              D_csr%coeffs(D_csr%nz_size+3)=(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I2(i,m_indx(j-2,nb))
              D_csr%col(D_csr%nz_size+3)=m_indx(m_indx(j,nb) + side - 1,n_size)
              D_csr%nz_size=D_csr%nz_size+3
             end if 
             
             
          elseif(j==m_indx(i-1,nb) .and. corner(j,ny,nz)) then
          
              D_csr%coeffs(D_csr%nz_size+1)=I2(i,m_indx(j-1,nb))
              D_csr%col(D_csr%nz_size+1)=j+side
              D_csr%nz_size=D_csr%nz_size+1 
              side=side+1
          elseif(j==m_indx(i-1,nb))  then
          else
       
          
              
             
              
              if(corner(j,ny,nz)) then
               
               if(m_indx(j-1,nb)==i)then
                 D_csr%coeffs(D_csr%nz_size+1)=I1(i,m_indx(j,nb))- I2(i,m_indx(j,nb))
                 D_csr%col(D_csr%nz_size+1)=j+side+1
                 D_csr%nz_size=D_csr%nz_size+1
               else
                 D_csr%coeffs(D_csr%nz_size+1)=I2(i,m_indx(j-1,nb))
                 D_csr%col(D_csr%nz_size+1)=j+side
                 D_csr%coeffs(D_csr%nz_size+2)=I1(i,m_indx(j,nb))- I2(i,m_indx(j,nb))
                 D_csr%col(D_csr%nz_size+2)=j+side+1 
                 D_csr%nz_size=D_csr%nz_size+2               !at corner we have 2 two b_normals 
               end if
               side=side+1
              else
               if(m_indx(j-1,nb)/=i)then
                D_csr%coeffs(D_csr%nz_size+1)=( I1(i,m_indx(j,nb))- I2(i,m_indx(j,nb)) + I2(i,m_indx(j-1,nb)))
                D_csr%col(D_csr%nz_size+1)=j+side
                D_csr%nz_size=D_csr%nz_size+1
               
               end if
              end if
             
          end if       
         end do
         
         D_csr%row(row_ptr)= D_csr%nz_size  
         temp_coeffs=D_csr%coeffs(D_csr%row(row_ptr-1)+1:D_csr%row(row_ptr))
         temp_col=D_csr%col(D_csr%row(row_ptr-1)+1:D_csr%row(row_ptr))
         
         call quicksort(temp_coeffs,temp_col,size(temp_col))
         
         D_csr%coeffs(D_csr%row(row_ptr-1)+1:D_csr%row(row_ptr))=temp_coeffs
         D_csr%col(D_csr%row(row_ptr-1)+1:D_csr%row(row_ptr))=temp_col
       end do
       close(10)
       
        
    end subroutine fill_bn_cohf_D
    
    
    subroutine fill_bn_cohf_F(F_csr,alpha_k,ny,nz,nb,h_array,I1,I2)
    use greens_func_subroutines
    use array_processes 
    implicit none
       integer,intent(in)                               ::ny,nz,nb
       type(node_csr),intent(out)                       ::F_csr
       real(kind=8),intent(in)                          ::h_array(1:nb)
       real(kind=8),allocatable,intent(in)              ::I1(:,:),I2(:,:)   
       real*8, intent(in)                               ::alpha_k
       real,parameter                                   ::PI=4*atan(1.d0)
       integer                                          ::i,m,j,k,count,c,chunk_size
       real                                             ::y,z,constant
       real(kind=8),allocatable                         ::temp_coeffs(:)
       integer,allocatable                              ::temp_col(:)
       integer,allocatable                              ::temp_row(:)
       integer                                          ::row_ptr,col_ptr, n_size
       integer                                          ::side,side_tr
       logical                                          ::flag
      
       chunk_size=1000
       F_csr%nz_size=0
       F_csr%capacity=nb                                 !nb=1024
    
       n_size=2*(ny+1)+2*(nz+1)                          !n_size=1028, no. of unique b_normals
       row_ptr=0
       col_ptr=0
                                 
       allocate(F_csr%row(0:n_size),F_csr%y_loc(0:n_size),F_csr%z_loc(0:n_size))
       allocate(F_csr%col(F_csr%capacity))
       allocate(F_csr%coeffs(F_csr%capacity))
       
    
   
       F_csr%row(row_ptr)=F_csr%nz_size
       F_csr%y_loc(0)=0.d0
       F_csr%z_loc(0)=0.d0
       side_tr=0
       flag=.true.  
    
       open(unit=10, file='../output_file/rhs_b13.dat')
      
       
       do i=1,n_size
       
        row_ptr=row_ptr+1 
        side=0
       
        
         do j=1,nb
         
          if((F_csr%capacity-F_csr%nz_size)<(ny+2))then
           if(F_csr%capacity==0)then
             allocate(F_csr%coeffs(chunk_size))
             allocate(F_csr%col(chunk_size))
            
           else
            temp_coeffs=F_csr%coeffs(1:F_csr%nz_size)
            temp_col   =F_csr%col(1:F_csr%nz_size)
            if(allocated(F_csr%coeffs)) deallocate(F_csr%coeffs)
            if(allocated(F_csr%col))   deallocate(F_csr%col)
            allocate(F_csr%coeffs(F_csr%capacity+chunk_size))
            allocate(F_csr%col(F_csr%capacity+chunk_size))
            F_csr%capacity=F_csr%capacity+chunk_size
            F_csr%coeffs(1:F_csr%nz_size)=temp_coeffs
            F_csr%col(1:F_csr%nz_size)   =temp_col
          
           end if
          end if  
          
          
          if( j==i-side_tr)then
            
             
           
            
             
             
           if(corner(j,ny,nz)) then
              
              F_csr%coeffs(F_csr%nz_size+1)=(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))
              F_csr%col(F_csr%nz_size+1)=m_indx(j,nb) + side
              F_csr%coeffs(F_csr%nz_size+2)=I1(i,m_indx(j,nb))-(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              F_csr%col(F_csr%nz_size+2)=m_indx(m_indx(j,nb) + side+1,n_size)
              F_csr%coeffs(F_csr%nz_size+3)=(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))+(I1(i,m_indx(j+1,nb))- I2(i,m_indx(j+1,nb)))
              F_csr%col(F_csr%nz_size+3)=m_indx(m_indx(j,nb) + side + 2,n_size)
              F_csr%coeffs(F_csr%nz_size+4)=I1(i,m_indx(j-1,nb))-(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I2(i,m_indx(j-2,nb))
              F_csr%col(F_csr%nz_size+4)=m_indx(m_indx(j,nb) + side - 1,n_size)
              F_csr%nz_size=F_csr%nz_size+4
              side=side+1
           !newly added condition 11/3/2025
            else if(corner(m_indx(j+1,nb),ny,nz)) then
              F_csr%coeffs(F_csr%nz_size+1)=(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I1(i,m_indx(j,nb))-(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              F_csr%col(F_csr%nz_size+1)=m_indx(j,nb) + side
              F_csr%coeffs(F_csr%nz_size+2)=(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              F_csr%col(F_csr%nz_size+2)=m_indx(m_indx(j,nb) + side + 1,n_size)
              F_csr%coeffs(F_csr%nz_size+3)=I1(i,m_indx(j-1,nb))-(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I2(i,m_indx(j-2,nb))
              F_csr%col(F_csr%nz_size+3)=m_indx(m_indx(j,nb) + side - 1,n_size)
              F_csr%nz_size=F_csr%nz_size+3 
             else if(corner(m_indx(j-1,nb),ny,nz)) then
              F_csr%coeffs(F_csr%nz_size+1)=(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I1(i,m_indx(j,nb))-(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              F_csr%col(F_csr%nz_size+1)=m_indx(j,nb) + side
              F_csr%coeffs(F_csr%nz_size+2)=(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))+(I1(i,m_indx(j+1,nb))- I2(i,m_indx(j+1,nb)))
              F_csr%col(F_csr%nz_size+2)=m_indx(m_indx(j,nb) + side + 1,n_size)
              F_csr%coeffs(F_csr%nz_size+3)=I1(i,m_indx(j-1,nb))-(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))
              F_csr%col(F_csr%nz_size+3)=m_indx(m_indx(j,nb) + side - 1,n_size)
              F_csr%nz_size=F_csr%nz_size+3
           !newly added
             else
              F_csr%coeffs(F_csr%nz_size+1)=(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I1(i,m_indx(j,nb))-(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))
              F_csr%col(F_csr%nz_size+1)=m_indx(j,nb) + side
              F_csr%coeffs(F_csr%nz_size+2)=(I2(i,m_indx(j,nb))/h_array(m_indx(j,nb)))+(I1(i,m_indx(j+1,nb))- I2(i,m_indx(j+1,nb)))
              F_csr%col(F_csr%nz_size+2)=m_indx(m_indx(j,nb) + side + 1,n_size)
              F_csr%coeffs(F_csr%nz_size+3)=I1(i,m_indx(j-1,nb))-(I2(i,m_indx(j-1,nb))/h_array(m_indx(j-1,nb)))+I2(i,m_indx(j-2,nb))
              F_csr%col(F_csr%nz_size+3)=m_indx(m_indx(j,nb) + side - 1,n_size)
              F_csr%nz_size=F_csr%nz_size+3
             end if 
             
             
          elseif(j==m_indx(i-side_tr-1,nb) .and. corner(j,ny,nz)) then
          
              F_csr%coeffs(F_csr%nz_size+1)=I2(i,m_indx(j-1,nb))
              F_csr%col(F_csr%nz_size+1)=j+side
              F_csr%nz_size=F_csr%nz_size+1 
              side=side+1
          elseif(j==m_indx(i-side_tr-1,nb))  then
          
          else
       
          
              
             
              
              if(corner(j,ny,nz)) then
               
               if(m_indx(j-1,nb)==i-side_tr)then
                 F_csr%coeffs(F_csr%nz_size+1)=I1(i,m_indx(j,nb))- I2(i,m_indx(j,nb))
                 F_csr%col(F_csr%nz_size+1)=j+side+1
                 F_csr%nz_size=F_csr%nz_size+1
               else
                 F_csr%coeffs(F_csr%nz_size+1)=I2(i,m_indx(j-1,nb))
                 F_csr%col(F_csr%nz_size+1)=j+side
                 F_csr%coeffs(F_csr%nz_size+2)=I1(i,m_indx(j,nb))- I2(i,m_indx(j,nb))
                 F_csr%col(F_csr%nz_size+2)=j+side+1 
                 F_csr%nz_size=F_csr%nz_size+2               !at corner we have 2 two b_normals 
                 
                 
               end if
               side=side+1
              else
               if(m_indx(j-1,nb)/=i-side_tr)then
                F_csr%coeffs(F_csr%nz_size+1)=( I1(i,m_indx(j,nb))- I2(i,m_indx(j,nb)) + I2(i,m_indx(j-1,nb)))
                F_csr%col(F_csr%nz_size+1)=j+side
                F_csr%nz_size=F_csr%nz_size+1
               
               end if
              end if
             
          end if       
         end do
        
        if(corner(i-side_tr,ny,nz) .and. flag)then
            side_tr=side_tr+1 
            flag=.false.  
        else
            flag=.true.      
        end if 
         
        
          F_csr%row(row_ptr)= F_csr%nz_size
          
         temp_coeffs=F_csr%coeffs(F_csr%row(row_ptr-1)+1:F_csr%row(row_ptr))
         temp_col=F_csr%col(F_csr%row(row_ptr-1)+1:F_csr%row(row_ptr))
         
         call quicksort(temp_coeffs,temp_col,size(temp_col))
         
         F_csr%coeffs(F_csr%row(row_ptr-1)+1:F_csr%row(row_ptr))=temp_coeffs
         F_csr%col(F_csr%row(row_ptr-1)+1:F_csr%row(row_ptr))=temp_col
        
       end do
       close(10)
       
        
    end subroutine fill_bn_cohf_F

end module node_classes
