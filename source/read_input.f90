!This subroutine read a input file using namelist
module input_sub
contains
 subroutine read_input(ny,nz,Ly,Lz,res_tolerance,A,Rm,linear_solver ,pre_cond,del_t,alpha_k)
  implicit none
  integer					:: ny,nz,Rm
  real*8 				        :: Ly,Lz,res_tolerance,A,del_t,alpha_k
  character 					:: linear_solver*10 ,pre_cond*10
 
  namelist/input_parameters/ny,nz,Ly,Lz,res_tolerance,A,linear_solver ,pre_cond,Rm,del_t,alpha_k
  open(unit=10, file='../input_file/input.txt')
  read(10, nml=input_parameters)
  close(10)
 end subroutine read_input

end module input_sub
