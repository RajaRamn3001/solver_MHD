
PROGRAM  main

 !Calling the modules
 use node_classes
 use matrix_vector_op
 use linear_solvers
 use greens_func_subroutines
 implicit none
 
 
 integer                                          :: ny,nz,Rm,n,nb,j,i,k,side_tr,n_size
 real*8                                           :: Ly,Lz,res_tolerance,A,del_t,alpha_k
 real*8                                           :: init_by, init_bz,f
 character                                        :: linear_solver*10 ,pre_cond*10
 real(kind=8),dimension(:),allocatable            :: gridpts_y,gridpts_z
 real(kind=8),dimension(:),allocatable            :: q,b,act_b,error_b,p,p_new,error_p
 real(kind=8),dimension(:),allocatable            :: psi 
 real(kind=8),dimension(:),allocatable            :: psi_act 
 real(kind=8),dimension(:),allocatable            :: rad_bn,Ia1,Ia2 
 real(kind=8),dimension(:),allocatable            :: b_normal,b_tau,E_psi,F_bn,E_Cinv_D_bn,b_tau_act
 real(kind=8),dimension(:),allocatable            :: b_radf, b_radb
 real(kind=8),dimension(:,:),allocatable          :: rad,rad1,rad2,rad3,rad4,rad5,rad6
 real(kind=8),dimension(:,:),allocatable          :: rad1h,rad2h,rad3h,rad4h,rad5h,rad6h
 real(kind=8),dimension(:,:),allocatable          :: rad1h_2,rad2h_2,rad3h_2,rad4h_2,rad5h_2,rad6h_2
 real(kind=8),dimension(:,:),allocatable          :: I1,I2,I_nor1,I_nor2,I_tan1,I_tan2,I_nt1,I_nt2,K02
 real(kind=8),dimension(:,:),allocatable          :: dK0_nor,K0
 real(kind=8),dimension(:),allocatable            :: h
 type(node_csr)                                   :: node_info_csr
 type(node_csr)                                   :: psi_cohf_mat_C ,psi_cohf_mat_inv_C, bn_cohf,psi_cohf_mat_E,psi_cohf_mat_C_ns,psi_cohf_mat_E_ns
 type(node_csr)                                   :: C_inv_D,bn_cohf_F,E_C_inv_D,H_mat
 type(node_csc)                                   :: node_info_csc
 real,parameter                                   :: PI=4*atan(1.d0)
 real                                             :: y,z
 integer*4                                        :: count_rate, count_max, start_count, end_count
 real*8                                           :: start_time, end_time
 real*8,allocatable                               :: F_row2(:)
 real*8,parameter                                 :: Sx=0.0d0,Sy=0.0d0     !shift of K0 function origin in x,y
 call read_input(ny,nz,Ly,Lz,res_tolerance,A,Rm,linear_solver ,pre_cond,del_t,alpha_k)
 print *,ny,nz,Ly,Lz,res_tolerance,A,Rm,linear_solver,del_t ,pre_cond,alpha_k
 
 
 f=3*Rm/2/del_t
 n=(ny+1)*(nz+1)*2
 nb=2*(ny+nz)
 n_size=2*(ny+1)+2*(nz+1)
 init_by=0
 init_bz=0
 
 
 allocate(gridpts_y(0:ny))
 allocate(gridpts_z(0:nz))
 allocate(act_b(n))
 allocate(b(n))
 allocate(error_b(n))
 allocate(psi(nb))
 allocate(Ia1(nb),Ia2(nb))
 allocate(rad_bn(nb))
 allocate(h(nb))
 allocate(dK0_nor(1:nb,1:nb),K0(1:nb,1:nb))
 allocate(rad(1:nb,1:nb),rad1(1:nb,1:nb),rad2(1:nb,1:nb),rad3(1:nb,1:nb),rad4(1:nb,1:nb),rad5(1:nb,1:nb),rad6(1:nb,1:nb))
 allocate(rad1h(1:nb+8,1:nb),rad2h(1:nb+8,1:nb),rad3h(1:nb+8,1:nb),rad4h(1:nb+8,1:nb),rad5h(1:nb+8,1:nb),rad6h(1:nb+8,1:nb))
 allocate(rad1h_2(1:nb+16,1:nb),rad2h_2(1:nb+16,1:nb),rad3h_2(1:nb+16,1:nb),rad4h_2(1:nb+16,1:nb),rad5h_2(1:nb+16,1:nb),rad6h_2(1:nb+16,1:nb))
 allocate(I1(1:nb,1:nb),I2(1:nb,1:nb))
 allocate(I_nor1(1:nb,1:nb),I_nor2(1:nb,1:nb))
 
 !allocate(I_tan1(1:n_size,1:nb),I_tan2(1:n_size,1:nb))
 !allocate(I_nt1(1:n_size,1:nb),I_nt2(1:n_size,1:nb))
 
 
 call generate_grid(ny,nz,Ly,Lz,A,gridpts_y,gridpts_z)
 !call fill_q_array(node_info_csr,q,f,alpha_k)
 
 !Calling Greens Function operation subroutines
 call create_rad_array(ny,nz,nb,gridpts_y,gridpts_z,Sx,Sy,rad_bn)
 call create_rad_mat(ny,nz,nb,gridpts_y,gridpts_z,rad)
 call create_h_array(ny,nz,nb,gridpts_y,gridpts_z,h)
 call create_6rad_mat(ny,nz,nb,gridpts_y,gridpts_z,rad1,rad2,rad3,rad4,rad5,rad6)
 call create_6rad_mat_extra(ny,nz,nb,h,gridpts_y,gridpts_z,rad1h,rad2h,rad3h,rad4h,rad5h,rad6h)
 call create_6rad_mat_extra2(ny,nz,nb,h,gridpts_y,gridpts_z,rad1h_2,rad2h_2,rad3h_2,rad4h_2,rad5h_2,rad6h_2)
 call integrals3(nb,h,alpha_k,Ia1,Ia2)
 call integral_K0(ny,nz,nb,alpha_k,gridpts_y,gridpts_z,h,rad1,rad2,rad3,rad4,rad5,rad6,I1,I2)
 call integral_dK0nor(ny,nz,nb,gridpts_y,gridpts_z,alpha_k,h,rad1,rad2,rad3,rad4,rad5,rad6,I_nor1,I_nor2)
 call integral_dK0dtau(ny,nz,nb,alpha_k,gridpts_y,gridpts_z,h,rad1,rad2,rad3,rad4,rad5,rad6,rad1h_2,rad2h_2,rad3h_2,rad4h_2,rad5h_2,rad6h_2,I_tan1,I_tan2)
 call integral_dK0nor_tau(ny,nz,nb,gridpts_y,gridpts_z,alpha_k,h,rad1,rad2,rad3,rad4,rad5,rad6,rad1h_2,rad2h_2,rad3h_2,rad4h_2,rad5h_2,rad6h_2,I_nt1,I_nt2)
 call dk0rad1d(ny,nz,nb,alpha_k,rad_bn,b_radf, b_radb)
 call dk0nor1d(ny,nz,nb,rad_bn,gridpts_y,gridpts_z,Sx,Sy,b_radf, b_radb,b_normal)
 
 call fill_psi_cohfmat(psi_cohf_mat_C,ny,nz,nb,alpha_k,h,I_nor1,I_nor2)
 call fill_psi_cohfmat_non_sparse(psi_cohf_mat_C_ns,ny,nz,nb,alpha_k,h,I_nor1,I_nor2)
 call fill_psi_cohf_E(psi_cohf_mat_E,ny,nz,nb,alpha_k,h,I_nt1,I_nt2)
 call fill_psi_cohf_E_non_sparse(psi_cohf_mat_E_ns,ny,nz,nb,alpha_k,h,I_nt1,I_nt2)
 call fill_bn_cohf_D(bn_cohf,alpha_k,ny,nz,nb,h,I1,I2)
 call fill_bn_cohf_F(bn_cohf_F,alpha_k,ny,nz,nb,h,I_tan1,I_tan2)

 !call modify_corner_coeffs(psi_cohf_mat_E,bn_cohf_F,ny,nz)
 
 !*********************************************
 print *,'size of node_info' , psi_cohf_mat_C%nz_size
 print *,'size of p' ,         size(b_normal)
 !deallocate(gridpts_y)!memory freed
 !deallocate(gridpts_z)!memory freed
 
 
 
 b_normal=-b_normal
 psi=0.65

 
   
call inv_mat(psi_cohf_mat_C,psi_cohf_mat_inv_C)
call mat_mat_mult_csr(psi_cohf_mat_inv_C,bn_cohf,C_inv_D)
print*,"C_inv_D created"
call mat_vec_mult(C_inv_D,b_normal,psi)
print*,"psi computed"
call act_func(nb,alpha_k,rad_bn,psi_act)
call mat_vec_mult(psi_cohf_mat_E,psi,E_psi)
print*,"E_psi computed"
call mat_vec_mult(bn_cohf_F,b_normal,F_bn)
print*,"F_bn computed"
call mat_mat_mult_csr(psi_cohf_mat_E,C_inv_D,E_C_inv_D)
print*,"E_C_inv_D computed"
call mat_vec_mult(E_C_inv_D,b_normal,E_Cinv_D_bn)
print*,"E_Cinv_D_bn computed"
!call vec_vec_sum(E_Cinv_D_bn,F_bn,b_tau,1)
call mat_mat_sum(E_C_inv_D,bn_cohf_F,H_mat,1)
print*,"H_mat computed"
call mat_vec_mult(H_mat,b_normal,b_tau)
print*,"b_tau computed"
call fill_DH_bem_coeff_mat(node_info_csr,H_mat,gridpts_y,gridpts_z,ny,nz,n,f,alpha_k)

b_tau=-2*b_tau

print*,"b_tau computed"
print*,"mat_vec multiplication done"  

!call act_func(nb,alpha_k,rad_bn,psi_act)
!call dk0tau1d(ny,nz,nb,rad_bn,gridpts_y,gridpts_z,b_radf, b_radb,b_tau_act)
call btau_func(nb,ny,nz,alpha_k,rad_bn, psi_act,b_tau_act,h)
!call mat_vec_mult(node_info_csr,psi,p_new)
b_tau_act=-b_tau_act
!--------------------------------------------------- 
call fill_act_b_array(node_info_csr,act_b,f,alpha_k)
call fill_q_array(node_info_csr,q,f,alpha_k)
!--------------------------------------------------
!------------------------------------------------------
!initialise
b=120.d0
!-----------------------------------------------------

call system_clock(start_count, count_rate, count_max)
     start_time = real(start_count)/real(count_rate)
     
 call solver(node_info_csr,b,q,res_tolerance,ny,nz,n,linear_solver,pre_cond)
 
 call system_clock(end_count, count_rate, count_max)
 end_time = real(end_count)/real(count_rate)
 print*, "Wall clock time:", end_time - start_time, "seconds"
 !----------------------------------------------------------------
 open(unit=50, file='../output_file/F_csr_b29_coeffs.dat')

  
 
      do j=1,1028
    
  
    
   
        
        
       ! write(50,"(*(F15.6))"),bn_cohf_F%coeffs(bn_cohf_F%row(j-1)+1:bn_cohf_F%row(j))
       
       
      end do 
  
     
  close(50)
  close(60)                     
!F_row2=bn_cohf_F%coeffs(bn_cohf_F%row(3-1)+1:bn_cohf_F%row(3))
!call vec_vec_mult(F_row2,b_normal,answ,"d")
!print*,"answ",answ
 
             
 open(unit=10,file='../output_file/psi_value_b35.dat')
 do i=1,nb
   write(10,'(1X,I5,2X,*(E14.6))'), i,psi(i),psi_act(i),b_radf(i)
   
 end do

 close(10)
side_tr=0 
open(unit=60,file='../output_file/b_tau_b35.dat')
 do i=1,nb
   if(corner(i,ny,nz)) then
    
    write(60,'(1X,I5,2X,*(E14.6))'), i+side_tr,b_tau(i+side_tr),b_normal(i+side_tr),b_tau_act(i+side_tr),E_psi(i+side_tr),F_bn(i+side_tr),psi(i)
    write(60,'(1X,I5,2X,*(E14.6))'), i+side_tr+1,b_tau(i+side_tr+1),b_normal(i+side_tr+1),b_tau_act(i+side_tr+1),E_psi(i+side_tr+1),F_bn(i+side_tr+1),psi(i)
    side_tr=side_tr+1
   else
    write(60,'(1X,I5,2X,*(E14.6))'), i+side_tr,b_tau(i+side_tr),b_normal(i+side_tr),b_tau_act(i+side_tr),E_psi(i+side_tr),F_bn(i+side_tr),psi(i)
   end if
 end do

open(unit=20,file='../output_file/I1_nt_b35.dat')
open(unit=30,file='../output_file/I2_nt_b35.dat')
open(unit=40,file='../output_file/I1_tan_b35.dat')
open(unit=50,file='../output_file/I2_tan_b35.dat')

open(unit=22,file='../output_file/I1nor_b35.dat')
open(unit=32,file='../output_file/I2nor_b35.dat')
open(unit=42,file='../output_file/I1_b35.dat')
open(unit=52,file='../output_file/I2_b35.dat')

open(unit=21,file='../output_file/H_mat_col_b35.dat')
open(unit=31,file='../output_file/H_mat_coeff_b35.dat')
open(unit=41,file='../output_file/H_mat_coeff_b35.dat')
open(unit=51,file='../output_file/error_rhs_b34.dat')

open(unit=23,file='../output_file/abs_error_b_b35.dat')
!open(unit=31,file='../output_file/E_corner_b27.dat')

 do i=1,nb+4

   write(20,'(*(E14.6))'),I_nt1(i,:)
   write(30,'(*(E14.6))'),I_nt2(i,:)

   write(40,'(*(E14.6))'),I_tan1(i,:)
   write(50,'(*(E14.6))'),I_tan2(i,:)
   
 end do
 
 
 do i=1,nb

   write(22,'(*(E14.6))'),I_nor1(i,:)
   write(32,'(*(E14.6))'),I_nor2(i,:)

   write(42,'(*(E14.6))'),I1(i,:)
   write(52,'(*(E14.6))'),I2(i,:)
   
 end do 
  do i=1,nb+4
   
      write(21,'(*(1XI5))'),H_mat%col(H_mat%row(i-1)+1:H_mat%row(i))
      write(31,'(*(E14.6))'),H_mat%coeffs(H_mat%row(i-1)+1:H_mat%row(i))           
  end do 
  !do i=1,1024
   
      !write(31,'(*(E14.6))'),bn_cohf%coeffs(bn_cohf%row(i-1)+1:bn_cohf%row(i))
                 
  !end do
  
  
!-----------------------------------------------------
call mat_vec_mult(node_info_csr,act_b,p)
error_p=abs(p-q)
error_b=abs(b-act_b)
do j=1,size(error_p)
 write(51,"(1X,I6,1X,E14.6)"),j,error_p(j)
 write(23,"(1X,*(E14.6))"),node_info_csr%y_loc(j),node_info_csr%z_loc(j),error_b(j)
end do
!-----------------------------------------------------
 close(20)  
 close(30)
 close(21)
 close(51)
 close(40)  
 close(50)
 close(60)
 close(22)  
 close(32)
 close(42)
 close(52) 
 close(23)
  deallocate(psi)
  
END PROGRAM main
