
PROGRAM  main

 !Calling the modules
 use node_classes
 use matrix_vector_op
 use linear_solvers
 use input_sub
 use grid_generate_sub
 use greens_func_subroutines
 
 implicit none
 
 
 integer                                          :: ny,nz,Rm,n,nb,j,i,k,band
 real*8                                           :: Ly,Lz,res_tolerance,A,del_t,alpha_k
 real*8                                           :: init_by, init_bz,f
 character                                        :: linear_solver*10 ,pre_cond*10
 real(kind=8),dimension(:),allocatable            :: gridpts_y,gridpts_z
 real(kind=8),dimension(:),allocatable            :: q,b,act_b,error_b,p,p_new,error_p
 real(kind=8),dimension(:),allocatable            :: psi,x,y_vec 
 real(kind=8),dimension(:),allocatable            :: psi_act 
 real(kind=8),dimension(:),allocatable            :: rad_bn,Ia1,Ia2 
 real(kind=8),dimension(:),allocatable            :: b_normal_f , b_normal_b,b_tau_f , b_tau_b
 real(kind=8),dimension(:),allocatable            :: b_radf, b_radb
 real(kind=8),dimension(:,:),allocatable          :: rad,rad1,rad2,rad3,rad4,rad5,rad6
 real(kind=8),dimension(:,:),allocatable          :: I1,I2,I_nor1,I_nor2
 real(kind=8),dimension(:,:),allocatable          :: dK0_nor,K0
 real(kind=8),dimension(:),allocatable            :: h
 type(node_csr)                                   :: node_info_csr,A_mat,B_mat,C_mat,A_mat_inv
 type(node_csc)                                   :: node_info_csc
 real,parameter                                   :: PI=4*atan(1.d0)
 real                                             :: y,z
 integer*4                                        :: count_rate, count_max, start_count, end_count
 real*8                                           :: start_time, end_time
 
 call read_input(ny,nz,Ly,Lz,res_tolerance,A,Rm,linear_solver ,pre_cond,del_t,alpha_k)
 print *,ny,nz,Ly,Lz,res_tolerance,A,Rm,linear_solver,del_t ,pre_cond,alpha_k
 
 
 f=3*Rm/2/del_t
 band=ny+2 !No of non-zero diagonals
 n=(ny+1)*(nz+1)*2
 nb=2*(ny+nz)
 init_by=0
 init_bz=0
 
 allocate(gridpts_y(0:ny))
 allocate(gridpts_z(0:nz))
 allocate(act_b(n))
 allocate(error_b(n))
 allocate(b(n))
 allocate(psi(nb))
 allocate(Ia1(nb),Ia2(nb))
 allocate(rad_bn(nb))
 allocate(h(nb))
 allocate(dK0_nor(1:nb,1:nb),K0(1:nb,1:nb))
 allocate(rad(1:nb,1:nb),rad1(1:nb,1:nb),rad2(1:nb,1:nb),rad3(1:nb,1:nb),rad4(1:nb,1:nb),rad5(1:nb,1:nb),rad6(1:nb,1:nb))
 allocate(I1(1:nb,1:nb),I2(1:nb,1:nb))
 allocate(I_nor1(1:nb,1:nb),I_nor2(1:nb,1:nb))
 print*, "memory allocation complete"
 b(1:n/2)=init_by
 b(n/2 +1 : )=init_bz
 
 !call generate_grid(ny,nz,Ly,Lz,A,gridpts_y,gridpts_z)
 !call fill_node_array(node_info_csr,gridpts_y,gridpts_z,ny,nz,n,f,alpha_k)
 !call fill_q_array(node_info_csr,q,f,alpha_k)
 
!Calling Greens Funcion operation subroutines
 !call create_rad_array(ny,nz,nb,gridpts_y,gridpts_z,rad_bn)
 !call create_rad_mat(ny,nz,nb,gridpts_y,gridpts_z,rad)
 !call create_h_array(ny,nz,nb,gridpts_y,gridpts_z,h)
 !call create_6rad_mat(ny,nz,nb,gridpts_y,gridpts_z,rad1,rad2,rad3,rad4,rad5,rad6)
 !call integrals2(nb,h,alpha_k,Ia1,Ia2)
 !call integral_K0(ny,nz,nb,alpha_k,gridpts_y,gridpts_z,h,rad1,rad2,rad3,rad4,rad5,rad6,I1,I2)
 !call integral_dK0nor(ny,nz,nb,gridpts_y,gridpts_z,alpha_k,h,rad1,rad2,rad3,rad4,rad5,rad6,I_nor1,I_nor2)
 !call dk0rad1d(ny,nz,nb,alpha_k,rad_bn,b_radf, b_radb)
 !call dk0nor1d(ny,nz,nb,rad_bn,gridpts_y,gridpts_z,b_radf, b_radb,b_normal_f,b_normal_b)
 
 !call fill_psi_coeffmat(node_info_csr,ny,nz,nb,alpha_k,h,I_nor1,I_nor2)
 !call fill_P_array(p,alpha_k,ny,nz,nb,b_normal_f,b_normal_b,h,I1,I2)
 
 
 !*********************************************
 deallocate(gridpts_y)!memory freed
 deallocate(gridpts_z)!memory freed
 
 
 allocate(A_mat%row(0:3),A_mat%y_loc(0:3),A_mat%z_loc(0:3))
 allocate(B_mat%row(0:3),B_mat%y_loc(0:3),B_mat%z_loc(0:3))
 A_mat%coeffs=(/7,4,8,5,7,3,7,8,5/)
 A_mat%col=(/1,2,3,1,2,3,1,2,3/)
 A_mat%row=(/0,3,6,9/)
 
 B_mat%coeffs=(/4,4,4,3,2,1,2,4,4,2,2,2,4,4,1/)
 B_mat%col=(/1,2,3,4,5,1,2,3,4,5,1,2,3,4,5/)
 B_mat%row=(/0,5,10,15/)
 
 y_vec=(/1,2,3,4,5/)
 
 print*,"yyy"
 call inv_mat(A_mat,A_mat_inv)
 
 call mat_mat_mult_csr(A_mat_inv,B_mat,C_mat)
 call mat_vec_mult(C_mat,y_vec,x)
 !print*,size(C_mat%coeffs)
 !print*,size(C_mat%col)
 print*,C_mat%row
 print*,x
END PROGRAM main
