module grid_generate_sub

contains
 subroutine generate_grid(ny,nz,Ly,Lz,A,gridpts_y,gridpts_z)
 implicit none
 INTEGER,intent(in) 				:: ny,nz
 INTEGER 					:: i,j,k
 REAL*8,intent(in) 				:: Ly,Lz,A
 REAL(kind=8) 					:: dy,dz
 REAL(kind=8),intent(out) 			:: gridpts_y(0:ny),gridpts_z(0:nz)
 dy=Ly/(ny)
 dz=Lz/(nz)
 

 do i=0,ny
   gridpts_y(i)= -1.d0+ dy*(i)
 end do
 
 do i=0,nz 
   gridpts_z(i)= -1.d0+ dz*(i)  
 end do
 !This loop stretches the grid point

  gridpts_y= tanh(A*gridpts_y)/tanh(A)
  gridpts_z= tanh(A*gridpts_z)/tanh(A) 

 open(unit=20, file='../output_file/grid_pts_stretched_b30.dat')
 
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
end module grid_generate_sub
