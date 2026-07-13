module matrix_vector_op
 implicit none
 interface mat_vec_mult
    module procedure mat_vec_mult_csr
 end interface   
contains
   subroutine mat_mat_mult_csr(mat_in1,mat_in2,mat_out)
   use node_classes
   implicit none
   type(node_csr), intent(inout)              :: mat_in1, mat_in2
   type(node_csr), intent(out)                :: mat_out

   integer                                    :: k,j,i,l,end_indx, r_indx
   integer,parameter                          :: buffer_size=5
   integer,parameter                          :: chunk_size=1000
   real*8,allocatable                         :: temp_coeffs(:)
   integer*4,allocatable                      :: temp_col(:)
   mat_out%capacity=size(mat_in2%row)
   mat_out%nz_size=0
  
   allocate(mat_out%row(0:size(mat_in1%row)-1))
   allocate(mat_out%coeffs( mat_out%capacity),mat_out%col( mat_out%capacity))
   mat_out%row(0)=mat_out%nz_size
   mat_out%coeffs=0.d0
   mat_out%col=0 
   
   do k=1, size(mat_in1%row)-1
      
     
      j=mat_in1%row(k-1)+1
      end_indx=mat_in1%row(k)
      
      
      do i=1,size(mat_in2%row)-1
        if((mat_in1%col(j)==i) .and. j<=end_indx)then
         
         
         r_indx=mat_out%row(k-1)+1
         do l=1,mat_in2%row(i)-mat_in2%row(i-1)
         
         !===============================================================
         !Allocate new memory size, if required
         !===============================================================
          if((mat_out%capacity-mat_out%nz_size)<buffer_size)then
           if(mat_out%capacity==0)then
             allocate(mat_out%coeffs(chunk_size))
             allocate(mat_out%col(chunk_size))
           else
             temp_coeffs=mat_out%coeffs(1:mat_out%nz_size)
             temp_col   =mat_out%col(1:mat_out%nz_size)
             if(allocated(mat_out%coeffs)) deallocate(mat_out%coeffs)
             if(allocated(mat_out%col))   deallocate(mat_out%col)
             allocate(mat_out%coeffs(mat_out%capacity+chunk_size))
             allocate(mat_out%col(mat_out%capacity+chunk_size))
             mat_out%coeffs=0.d0 
             mat_out%col=0            
             mat_out%capacity=mat_out%capacity+chunk_size
             mat_out%coeffs(1:mat_out%nz_size)=temp_coeffs
             mat_out%col(1:mat_out%nz_size)   =temp_col
          
           end if
          end if
         !===================================================================
          
      
          
          
          if(mat_out%col(r_indx) == mat_in2%col(mat_in2%row(i-1)+l)) then
            mat_out%coeffs(r_indx)=mat_out%coeffs(r_indx) + mat_in1%coeffs(j)*mat_in2%coeffs(mat_in2%row(i-1)+l)
            r_indx=r_indx+1
            !if found already existing column
           
          else if( .not. any(mat_out%col(r_indx:) == mat_in2%col(mat_in2%row(i-1)+l))) then
            mat_out%nz_size=mat_out%nz_size+1
            mat_out%coeffs(mat_out%nz_size)=mat_out%coeffs(mat_out%nz_size) + mat_in1%coeffs(j)*mat_in2%coeffs(mat_in2%row(i-1)+l)
            mat_out%col(mat_out%nz_size)=mat_in2%col(mat_in2%row(i-1)+l)
            !if new column is found, increase nz_size by 1
          else 
            r_indx=r_indx+1
            !Not increase r_indx if new column is found
          end if
           
         end do
         
         j=j+1
        end if 
         mat_out%row(k)=mat_out%nz_size
         
      end do
      
      
   print*,"mat_out%row",mat_out%row(k)   
   end do
    print*,"okay"   
   end subroutine mat_mat_mult_csr

   
   subroutine mat_vec_mult_csr(mat,vec_in,vec_out)
   use node_classes
   implicit none
   type(node_csr), intent(inout)              :: mat
   real(kind=8), allocatable ,intent (in)     :: vec_in(:)
   real(kind=8), allocatable ,intent (out)    :: vec_out(:)
   integer :: row_pt,i
   real(kind=8) :: sum
   allocate(vec_out(size(mat%row)-1))
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

