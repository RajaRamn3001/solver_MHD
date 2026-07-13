module greens_func_subroutines

 contains

! ***************************************************
!    To determine if a boundary point is a corner
! ***************************************************

!DEC$ ATTRIBUTES FORCEINLINE :: corner
    integer*4 function corner(indx,ny,nz)

      implicit none
      integer*4,intent(in)::indx,ny,nz

      if((indx==1) .OR. (indx==nz+1) .OR. (indx==nz+ny+1) .OR. (indx==2*nz+ny+1)) then
        corner = 1
      else
        corner = 0
      endif
  
    return
    end function corner

! *************************************
!    To compute j index given indx
! *************************************

!DEC$ ATTRIBUTES FORCEINLINE :: jb
    integer*4 function jb(indx,ny,nz)

      implicit none
      integer*4,intent(in)::indx,ny,nz

      if(indx<nz+1) then
         jb=0 
      elseif(indx<(nz+ny+1)) then
         jb=indx-nz-1          
      elseif(indx<(2*nz+ny+1)) then
         jb=ny          
      else
         jb=2*ny+2*nz-indx+1     
      endif
  
    return
    end function jb


! *************************************
!    To compute k index given i
! *************************************

!DEC$ ATTRIBUTES FORCEINLINE :: kb
    integer*4 function kb(indx,ny,nz)

      implicit none
      integer*4,intent(in)::indx,ny,nz

      if(indx<nz+1) then
         kb=indx-1
      elseif(indx<(nz+ny+1)) then 
         kb=nz          
      elseif(indx<(2*nz+ny+1)) then
	 kb=2*nz+ny-indx+1          
      else
	 kb=0         
      endif

    return
    end function kb
    
! *************************************
!    To compute cyclic index given m
! *************************************
integer*4 function m_indx(indx,nb)

  implicit none
  integer*4,intent(in)::nb
  integer*4,intent(in)::indx

  if(indx <= 0) then
     m_indx = nb + indx
  elseif(indx > nb) then 
     m_indx = indx - nb            
  else  ! This handles 0 < indx <= nb
     m_indx = indx
  endif


    return
  end function m_indx 
    
    ! *************************************
!    To compute  beta
! *************************************

    real*8 function beta_const(indx,ny,nz)

      implicit none
      integer*4,intent(in)::ny,nz
      integer*4,intent(in)::indx

      if(corner(indx,ny,nz)==1) then
         beta_const=3.d0/4.d0
      elseif(corner(indx,ny,nz)==0) then 
          beta_const=0.5           
      endif

    return
    end function beta_const  
    
! ***************************************************
! To compute the radial coordinate for each boundary node in an array
!****************************************************
subroutine create_rad_array(ny,nz,nb,gridpts_y,gridpts_z,sx,sy,rad)
  integer , intent(in)                          ::ny,nz,nb
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  , intent(in)                          ::sx,sy
  real*8  , intent(out)                         ::rad(1:nb)

  integer*4                                     ::i
  real*8                                        ::y_disp,z_disp
 
 
 do i=1,nb
    y_disp=gridpts_y(jb(i,ny,nz)) + sx
    z_disp=gridpts_z(kb(i,ny,nz)) + sy
    rad(i)=sqrt(y_disp**2 + z_disp**2) 
 end do
  
  
  print*,'Radial coordinate array generated'
end subroutine create_rad_array
! ***************************************************
! To compute the radial coordinate matrix
!****************************************************
subroutine create_rad_mat(ny,nz,nb,gridpts_y,gridpts_z,rad)
  integer , intent(in)                          ::ny,nz,nb
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  , intent(out)                         ::rad(1:nb,1:nb)

  integer*4                                     ::i,m 
  real*8                                        ::y_disp,z_disp
 
  do i=1,nb
    do m=1,nb
      y_disp=gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m,ny,nz))
      z_disp=gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m,ny,nz))
      rad(i,m)=sqrt(y_disp**2 + z_disp**2)
    end do
  end do 
  
  print*,'Radial coordinate 2d array generated'
end subroutine create_rad_mat


! ***************************************************
! To compute the array of lengths of each elements 
!****************************************************
subroutine create_h_array(ny,nz,nb,gridpts_y,gridpts_z,h)
  integer , intent(in)                          ::ny,nz,nb
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  , intent(out)                         ::h(1:nb)

  integer*4                                     ::i,m 
  real*8                                        ::y_disp,z_disp
  
  do i=1,nb
     y_disp=gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz))
     z_disp=gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz))
     h(i)=sqrt(y_disp**2 + z_disp**2)
  end do 
  print*,'h array created'
end subroutine create_h_array

! **********************************************************************************
!  To get the corresponding boundary index from 2D coordinates of boundary points
! **********************************************************************************
  integer*4 function ibound(j,k,ny,nz)

      implicit none
      !input variables
      integer*4 :: j,k,ny,nz

      !local variables
     
      
     

         if(j==0) then 
           ibound = k+1
         elseif(k==nz) then
           ibound = nz+1+j
         elseif(j==ny) then
           ibound = nz+ny+1 + (nz-k)
         else
           ibound = 2*nz+ny+1+ (ny-j)
         endif

      return
      end function ibound

! **************************
!     Bessel function I0 
! **************************
!DEC$ ATTRIBUTES FORCEINLINE :: bessi0
    real*8 function bessi0(x)

      implicit none
      real*8,intent(in)::x

      ! local variables
      real*8 :: ax,y

	  ax=abs(x)
          if(ax<3.75) then
            y=x/3.75
            y=y**2
            bessi0 = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))))
          else
            y=3.75/ax   
	    bessi0 = (exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2 +  & 
		     y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1 + y*0.392377e-2))))))))
          endif

    return
    end function bessi0


! *********************************************************
!     Evaluation of radial derrivative of I0 w.r.t to x
! *********************************************************
!DEC$ ATTRIBUTES FORCEINLINE :: dbessi0
    real*8 function dbessi0(x,alpha)

      implicit none
      real*8,intent(in)::x,alpha

      ! local variables
      real*8 :: ax,const1(1:6),const2(1:9)
      integer*4::k
      
      const1(1:6) = (/3.5156229,3.0899424,1.2067492,0.2659732,0.360768e-1,0.45813e-2/)

      const2(1:9)=(/0.39894228,0.1328592e-1,0.225319e-2,-0.157565e-2,0.916281e-2,-0.2057706e-1,0.2635537e-1,-0.1647633e-1,0.392377e-2/)

	  ax=abs(x*alpha)

          dbessi0=0.0

          if(ax<3.75) then

            do k=1,6
              dbessi0 = dbessi0 + (const1(k)*2.0*k*(alpha**(2*k))*(x**(2*k-1))/(3.75**(2*k)))           
            enddo 

          else 

            do k=1,9 
              dbessi0 = dbessi0 + exp(ax)*const2(k)*(3.75**(k-1))* (alpha*sqrt(x) + (0.5-k)/sqrt(x))/((alpha**(k-0.5))*(x**k)) 
            enddo 

          endif

    return
    end function dbessi0



! ******************************************************************************
!     Evaluation of radial derrivative of Mcdonald function(K0) for a 2D array
! ******************************************************************************
      subroutine dk0rad(nb,alpha,rad,output1)
      
      implicit none
      integer*4,intent(in)                        :: nb
      real*8,intent(in)                           :: alpha
      real*8,intent(in)                           :: rad(1:nb,1:nb)
      real*8,allocatable,intent(out)              :: output1(:,:)

      ! local variables
      real*8, pointer :: sum1(:,:),sum2(:,:),sum3(:,:)
      real*8 :: const1(1:7),const2(1:7)
      integer*4 :: i,m,k

      allocate(sum1(1:nb,1:nb),sum2(1:nb,1:nb),sum3(1:nb,1:nb))
      allocate(output1(1:nb,1:nb))
      ! Constants for series expansion of K0 function
      const1(1:7)=(/-0.57721566,0.42278420,0.23069756,0.348859e-1,0.262698e-2,0.10750e-3,0.74e-5/)
      const2(1:7)=(/1.25331414,-0.7832358e-1,0.2189568e-1,-0.1062446e-1,0.587872e-2,-0.251540e-2,0.53208e-3/)

      if(alpha==0) then

        do i=1,nb
          do m=1,nb

           if(i/=m) then
             output1(i,m) =  -1.0/rad(i,m)
           endif

          enddo
        enddo
   
      else

        do i=1,nb
          do m=1,nb

           output1(i,m) = 0.0d0

           if(i/=m) then

              if((alpha*rad(i,m))<=2.0d0) then

                 sum1(i,m)=-1.0*bessi0(alpha*rad(i,m))/rad(i,m) + (log(2.0d0)-log(alpha*rad(i,m)))*dbessi0(rad(i,m),alpha)             
                 sum2(i,m)=0.0
       
                 do k=1,6
                   sum2(i,m) = sum2(i,m) + (const1(k+1)*2.0*k*(alpha**(2*k))*(rad(i,m)**(2*k-1)))/(4.0**k)         
                 enddo

                 output1(i,m) = sum1(i,m) + sum2(i,m) 
 
              else
             
                 sum3(i,m)=0.0

                 do k=1,7
                   sum3(i,m) = sum3(i,m) + exp(-1.0d0*alpha*rad(i,m)) * const2(k)*(2.0**(k-1))*(alpha**(0.5-k))*(rad(i,m)**(0.5-k)) * (-1.0*alpha + (0.5-k)/rad(i,m))
                 enddo 

	         output1(i,m) = sum3(i,m) 

              endif

           endif

          enddo
        enddo

      endif

      deallocate(sum1,sum2,sum3)
     
      return
      end subroutine dk0rad

! ******************************************************************************
!     Evaluation of radial derrivative of Mcdonald function(K0) for a 2D array
! ******************************************************************************
!dk0radh is used to find delG/delr for the sub-pole radial coordinates
      subroutine dk0radh(nb,alpha,rad,output1)
      
      implicit none
      integer*4,intent(in)                        :: nb
      real*8,intent(in)                           :: alpha
      real*8,allocatable,intent(in)               :: rad(:,:)
      real*8,allocatable,intent(out)              :: output1(:,:)

      ! local variables
      real*8, pointer :: sum1(:,:),sum2(:,:),sum3(:,:)
      real*8 :: const1(1:7),const2(1:7)
      integer*4 :: i,m,k,n_size
      
      n_size=size(rad,1)
      allocate(sum1(1:n_size,1:nb),sum2(1:n_size,1:nb),sum3(1:n_size,1:nb))
      if(allocated(output1)) deallocate(output1)
      allocate(output1(1:n_size,1:nb))
      ! Constants for series expansion of K0 function
      const1(1:7)=(/-0.57721566,0.42278420,0.23069756,0.348859e-1,0.262698e-2,0.10750e-3,0.74e-5/)
      const2(1:7)=(/1.25331414,-0.7832358e-1,0.2189568e-1,-0.1062446e-1,0.587872e-2,-0.251540e-2,0.53208e-3/)

      if(alpha==0) then

        do i=1,n_size
          do m=1,nb

           if(i/=m) then
             output1(i,m) =  -1.0/rad(i,m)
           endif

          enddo
        enddo
   
      else

        do i=1,n_size
          
          do m=1,nb

           output1(i,m) = 0.d0

           if(rad(i,m)/=0.d0) then

              if((alpha*rad(i,m))<=2.0d0) then

                 sum1(i,m)=-1.0*bessi0(alpha*rad(i,m))/rad(i,m) + (log(2.0d0)-log(alpha*rad(i,m)))*dbessi0(rad(i,m),alpha)             
                 sum2(i,m)=0.0
       
                 do k=1,6
                   sum2(i,m) = sum2(i,m) + (const1(k+1)*2.0*k*(alpha**(2*k))*(rad(i,m)**(2*k-1)))/(4.0**k)         
                 enddo

                 output1(i,m) = sum1(i,m) + sum2(i,m) 
 
              else
             
                 sum3(i,m)=0.0

                 do k=1,7
                   sum3(i,m) = sum3(i,m) + exp(-1.0d0*alpha*rad(i,m)) * const2(k)*(2.0**(k-1))*(alpha**(0.5-k))*(rad(i,m)**(0.5-k)) * (-1.0*alpha + (0.5-k)/rad(i,m))
                 enddo 

	         output1(i,m) = sum3(i,m) 

              endif

           endif

          enddo
        enddo

      endif

      deallocate(sum1,sum2,sum3)
     
      return
      end subroutine dk0radh

! ******************************************************************************
!     Evaluation of radial derrivative of Mcdonald function(K0(r)) at boundary nodes 
!     given radial coordinate array r.
! ******************************************************************************
      subroutine dk0rad1d(ny,nz,nb,alpha,rad,output_f,output_b)
      
      implicit none
      integer*4,intent(in)                        :: nb,ny,nz
      real*8,intent(in)                           :: alpha
      real*8,intent(in)                           :: rad(1:nb)
      real*8,allocatable,intent(out)              :: output_f(:),output_b(:)

      ! local variables
      real*8, pointer :: sum1(:),sum2(:),sum3(:)
      real*8 :: const1(1:7),const2(1:7)
      integer*4 :: i,m,k

      allocate(sum1(1:nb),sum2(1:nb),sum3(1:nb))
      allocate(output_f(1:nb),output_b(1:nb))
      ! Constants for series expansion of K0 function
      const1(1:7)=(/-0.57721566,0.42278420,0.23069756,0.348859e-1,0.262698e-2,0.10750e-3,0.74e-5/)
      const2(1:7)=(/1.25331414,-0.7832358e-1,0.2189568e-1,-0.1062446e-1,0.587872e-2,-0.251540e-2,0.53208e-3/)

      if(alpha==0) then

        do i=1,nb
          
          if(i<nz+1)            then
            output_f(i)=-1.0/rad(i)
            output_b(i)=-1.0/rad(i)
          elseif(i<(nz+ny+1))   then
            output_f(i)=-1.0/rad(i)
            output_b(i)=-1.0/rad(i)       
          elseif(i<(2*nz+ny+1)) then
            output_f(i)=-1.0/rad(i) 
            output_b(i)=-1.0/rad(i)       
          elseif (i<=nb)         then
            output_f(i)=-1.0/rad(i)
            output_b(i)=-1.0/rad(i)      
          endif
          
         
           

         
        enddo
   
      else
      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
        do i=1,nb
          

           output_f(i) = 0.0d0
           output_b(i) = 0.0d0
           

              if((alpha*rad(i))<=2.0d0) then

                 sum1(i)=-1.0*bessi0(alpha*rad(i))/rad(i) + (log(2.0d0)-log(alpha*rad(i)))*dbessi0(rad(i),alpha)             
                 sum2(i)=0.0
       
                 do k=1,6
                   sum2(i) = sum2(i) + (const1(k+1)*2.0*k*(alpha**(2*k))*(rad(i)**(2*k-1)))/(4.0**k)         
                 enddo

                  
                  
                 
                  output_f(i)=sum1(i) + sum2(i)
                  output_b(i)=sum1(i) + sum2(i)
               
          
                
              else
             
                 

                 do k=1,7
                   sum3(i) = sum3(i) + exp(-1.0d0*alpha*rad(i)) * const2(k)*(2.0**(k-1))*(alpha**(0.5-k))*(rad(i)**(0.5-k)) * (-1.0*alpha + (0.5-k)/rad(i))
                 enddo 

	   
              
                  output_f(i)=sum3(i)
                  output_b(i)=sum3(i)
                 
              endif

        

          
        enddo
             

      endif
   
      deallocate(sum1,sum2,sum3)
     
      return
      end subroutine dk0rad1d



! **************************************************************************************************************
!     Evaluation of normal derrivative of Mcdonald function(K0) for a 2D array given the radial derrivative
! **************************************************************************************************************
      subroutine dk0nor(ny,nz,nb,rad,grid_pty,grid_ptz,dk0dr,output)
      
      implicit none
      integer*4,intent(in)                           ::nz,ny,nb
      real*8,intent(in)                              ::grid_pty(0:ny),grid_ptz(0:nz)
      real*8,intent(in)                              ::rad(1:nb,1:nb),dk0dr(1:nb,1:nb)
      real*8,intent(out)                             ::output(1:nb,1:nb)
      real*8,allocatable                             ::temp(:,:)
      
      ! local variables
      integer*4 :: i,m,j1,j2,k1,k2     
      
      temp=dk0dr
      do i=1,nb

        j1 = jb(i,ny,nz)
        k1 = kb(i,ny,nz)

        do m=1,nb

          j2 = jb(m,ny,nz)
          k2 = kb(m,ny,nz)

          
          
          if((m>=nz+1 .AND. m<nz+ny+1) .OR. (m>=2*nz+ny+1)) then
            output(i,m) = temp(i,m)*abs(grid_ptz(k1)-grid_ptz(k2))/rad(i,m)
          else
            output(i,m) = temp(i,m)*abs(grid_pty(j1)-grid_pty(j2))/rad(i,m)
          endif

        enddo
      enddo

      return
       
      end subroutine dk0nor

! **************************************************************************************************************
!     Evaluation of normal derrivative of Mcdonald function(K0) for a 2D array given the radial derrivative
! **************************************************************************************************************
!normal derivative from sub-poles with 2 subpoles in corner elements
      subroutine dk0norh(ny,nz,nb,rad,grid_pty,grid_ptz,dk0dr,output)
      
      implicit none
      integer*4,intent(in)                           ::nz,ny,nb
      real*8,intent(in)                              ::grid_pty(0:ny),grid_ptz(0:nz)
      real*8,allocatable,intent(in)                  ::rad(:,:),dk0dr(:,:)
      real*8,allocatable,intent(out)                 ::output(:,:)
      real*8,allocatable                             ::temp(:,:)
      
      ! local variables
      integer*4 :: i,m,j1,j2,k1,k2,n_size 
      real*8    ::y1,y12,y13,z1,z12,z13    
      integer   :: side_tr
      
      n_size=size(rad,1)
  
     
      side_tr=0
      
      if (allocated(output)) deallocate(output)
      allocate(output(1:n_size,1:nb))
       
      do i=1,nb 
        if(corner(i,ny,nz))then
         y1 = grid_pty(jb(i,ny,nz))-0.25*(grid_pty(jb(i,ny,nz))-grid_pty(jb(m_indx(i-1,nb),ny,nz)))
         z1 = grid_ptz(kb(i,ny,nz))-0.25*(grid_ptz(kb(i,ny,nz))-grid_ptz(kb(m_indx(i-1,nb),ny,nz)))
         
         y12 = grid_pty(jb(i,ny,nz))+ 0.25*(grid_pty(jb(m_indx(i+1,nb),ny,nz))-grid_pty(jb(i,ny,nz)))
         z12 = grid_ptz(kb(i,ny,nz))+ 0.25*(grid_ptz(kb(m_indx(i+1,nb),ny,nz))-grid_ptz(kb(i,ny,nz)))
         
         y13 = grid_pty(jb(i,ny,nz))+ 0.5*(grid_pty(jb(m_indx(i+1,nb),ny,nz))-grid_pty(jb(i,ny,nz)))
         z13 = grid_ptz(kb(i,ny,nz))+ 0.5*(grid_ptz(kb(m_indx(i+1,nb),ny,nz))-grid_ptz(kb(i,ny,nz)))
         do m=1,nb
         
          j2 = jb(m,ny,nz)
          k2 = kb(m,ny,nz)

          
          
          if((m>=nz+1 .AND. m<nz+ny+1) .OR. (m>=2*nz+ny+1)) then
            output(i+side_tr,m) = dk0dr(i+side_tr,m)*abs(z1-grid_ptz(k2))/rad(i+side_tr,m)
            output(i+side_tr+1,m) = dk0dr(i+side_tr+1,m)*abs(z12-grid_ptz(k2))/rad(i+side_tr+1,m)
            output(i+side_tr+2,m) = dk0dr(i+side_tr+2,m)*abs(z13-grid_ptz(k2))/rad(i+side_tr+2,m)
          else
            output(i+side_tr,m) = dk0dr(i+side_tr,m)*abs(y1-grid_pty(j2))/rad(i+side_tr,m)
            output(i+side_tr+1,m) = dk0dr(i+side_tr+1,m)*abs(y12-grid_pty(j2))/rad(i+side_tr+1,m)
            output(i+side_tr+2,m) = dk0dr(i+side_tr+2,m)*abs(y13-grid_pty(j2))/rad(i+side_tr+2,m)
          endif

         enddo
          
         side_tr=side_tr+2
        else
         y1 = grid_pty(jb(i,ny,nz))+ 0.5*(grid_pty(jb(m_indx(i+1,nb),ny,nz))-grid_pty(jb(i,ny,nz)))
         z1 = grid_ptz(kb(i,ny,nz))+ 0.5*(grid_ptz(kb(m_indx(i+1,nb),ny,nz))-grid_ptz(kb(i,ny,nz)))
        do m=1,nb
         
          j2 = jb(m,ny,nz)
          k2 = kb(m,ny,nz)

          
          
          if((m>=nz+1 .AND. m<nz+ny+1) .OR. (m>=2*nz+ny+1)) then
            output(i+side_tr,m) = dk0dr(i+side_tr,m)*abs(z1-grid_ptz(k2))/rad(i+side_tr,m)
          else
            output(i+side_tr,m) = dk0dr(i+side_tr,m)*abs(y1-grid_pty(j2))/rad(i+side_tr,m)
          endif

        enddo
       end if
      enddo
      
      
      return
       
      end subroutine dk0norh

 !normal derivative from sub-poles with 3 subpoles in corner elements
   subroutine dk0norh2(ny,nz,nb,rad,grid_pty,grid_ptz,dk0dr,output)
      
      implicit none
      integer*4,intent(in)                           ::nz,ny,nb
      real*8,intent(in)                              ::grid_pty(0:ny),grid_ptz(0:nz)
      real*8,allocatable,intent(in)                  ::rad(:,:),dk0dr(:,:)
      real*8,allocatable,intent(out)                 ::output(:,:)
      real*8,allocatable                             ::temp(:,:)
      
      ! local variables
      integer*4                                      :: i,m,j1,j2,k1,k2,n_size 
      real*8                                         ::y1,y12,y13,y14,y15,z1,z12,z13,z14,z15    
      integer                                        :: side_tr
      real*8,parameter                               ::c1=0.25,c2=0.001,c3=0.5
    
      n_size=size(rad,1)
  
     
      side_tr=0
      
      if (allocated(output)) deallocate(output)
      allocate(output(1:n_size,1:nb))
       
      do i=1,nb 
        if(corner(i,ny,nz))then
         y1 = grid_pty(jb(i,ny,nz))-c1*(grid_pty(jb(i,ny,nz))-grid_pty(jb(m_indx(i-1,nb),ny,nz)))
         z1 = grid_ptz(kb(i,ny,nz))-c1*(grid_ptz(kb(i,ny,nz))-grid_ptz(kb(m_indx(i-1,nb),ny,nz)))
         
         y12 = grid_pty(jb(i,ny,nz))-c2*(grid_pty(jb(i,ny,nz))-grid_pty(jb(m_indx(i-1,nb),ny,nz)))
         z12 = grid_ptz(kb(i,ny,nz))-c2*(grid_ptz(kb(i,ny,nz))-grid_ptz(kb(m_indx(i-1,nb),ny,nz)))
         
         y13 = grid_pty(jb(i,ny,nz))+c2*(grid_pty(jb(i+1,ny,nz))-grid_pty(jb(m_indx(i,nb),ny,nz)))
         z13 = grid_ptz(kb(i,ny,nz))+c2*(grid_ptz(kb(i+1,ny,nz))-grid_ptz(kb(m_indx(i,nb),ny,nz)))
         
         y14 = grid_pty(jb(i,ny,nz))+c1*(grid_pty(jb(m_indx(i+1,nb),ny,nz))-grid_pty(jb(i,ny,nz)))
         z14 = grid_ptz(kb(i,ny,nz))+c1*(grid_ptz(kb(m_indx(i+1,nb),ny,nz))-grid_ptz(kb(i,ny,nz)))
         
         y15 = grid_pty(jb(i,ny,nz))+c3*(grid_pty(jb(m_indx(i+1,nb),ny,nz))-grid_pty(jb(i,ny,nz)))
         z15 = grid_ptz(kb(i,ny,nz))+c3*(grid_ptz(kb(m_indx(i+1,nb),ny,nz))-grid_ptz(kb(i,ny,nz)))
         do m=1,nb
         
          j2 = jb(m,ny,nz)
          k2 = kb(m,ny,nz)

          
          
          if((m>=nz+1 .AND. m<nz+ny+1) .OR. (m>=2*nz+ny+1)) then
            output(i+side_tr,m) = dk0dr(i+side_tr,m)*abs(z1-grid_ptz(k2))/rad(i+side_tr,m)
            output(i+side_tr+1,m) = dk0dr(i+side_tr+1,m)*abs(z12-grid_ptz(k2))/rad(i+side_tr+1,m)
            output(i+side_tr+2,m) = dk0dr(i+side_tr+2,m)*abs(z13-grid_ptz(k2))/rad(i+side_tr+2,m)
            output(i+side_tr+3,m) = dk0dr(i+side_tr+3,m)*abs(z14-grid_ptz(k2))/rad(i+side_tr+3,m)
            output(i+side_tr+4,m) = dk0dr(i+side_tr+4,m)*abs(z15-grid_ptz(k2))/rad(i+side_tr+4,m)
          else
            output(i+side_tr,m) = dk0dr(i+side_tr,m)*abs(y1-grid_pty(j2))/rad(i+side_tr,m)
            output(i+side_tr+1,m) = dk0dr(i+side_tr+1,m)*abs(y12-grid_pty(j2))/rad(i+side_tr+1,m)
            output(i+side_tr+2,m) = dk0dr(i+side_tr+2,m)*abs(y13-grid_pty(j2))/rad(i+side_tr+2,m)
            output(i+side_tr+3,m) = dk0dr(i+side_tr+3,m)*abs(y14-grid_pty(j2))/rad(i+side_tr+3,m)
            output(i+side_tr+4,m) = dk0dr(i+side_tr+4,m)*abs(y15-grid_pty(j2))/rad(i+side_tr+4,m)
          endif

         enddo
          
         side_tr=side_tr+4
        else
         y1 = grid_pty(jb(i,ny,nz))+ c3*(grid_pty(jb(m_indx(i+1,nb),ny,nz))-grid_pty(jb(i,ny,nz)))
         z1 = grid_ptz(kb(i,ny,nz))+ c3*(grid_ptz(kb(m_indx(i+1,nb),ny,nz))-grid_ptz(kb(i,ny,nz)))
        do m=1,nb
         
          j2 = jb(m,ny,nz)
          k2 = kb(m,ny,nz)

          
          
          if((m>=nz+1 .AND. m<nz+ny+1) .OR. (m>=2*nz+ny+1)) then
            output(i+side_tr,m) = dk0dr(i+side_tr,m)*abs(z1-grid_ptz(k2))/rad(i+side_tr,m)
          else
            output(i+side_tr,m) = dk0dr(i+side_tr,m)*abs(y1-grid_pty(j2))/rad(i+side_tr,m)
          endif

        enddo
       end if
      enddo
      
      
      return
       
      end subroutine dk0norh2
! **************************************************************************************************************
!     Evaluation of normal derrivative of Mcdonald function(K0) on boundary nodes 
!     given the radial derrivative
! **************************************************************************************************************
!Evaluating normal mag. field at boundary with sourcr shifting enabled
      subroutine dk0nor1d(ny,nz,nb,rad,grid_pty,grid_ptz,sx,sy,dk0drf1d,dk0drb1d,output)
      
      implicit none
      integer*4,intent(in)                           ::nz,ny,nb
      real*8,intent(in)                              ::grid_pty(0:ny),grid_ptz(0:nz)
      real*8,intent(in)                              ::rad(1:nb),dk0drf1d(1:nb),dk0drb1d(1:nb)
      real*8,intent(in)                              ::sx,sy
      real*8,allocatable,intent(out)                 ::output(:)
      integer                                        ::n_size, side_tr,tet
      ! local variables
      integer*4 :: i,m,j1,j2,k1,k2     
      real*8    :: cos_theta,rad_unshifted,dirn_coeff
      
      n_size=(ny+1)*2 + (nz+1)*2
      side_tr=0
      tet=0
      allocate(output(1:n_size))
       do i=1,nb

        j1 = jb(i,ny,nz)
        k1 = kb(i,ny,nz)

        rad_unshifted=sqrt(grid_pty(j1)**2 + grid_ptz(k1)**2)
        
        cos_theta= (grid_pty(j1)*(grid_pty(j1) + sx) + grid_ptz(k1)*(grid_ptz(k1) + sy) ) &
                       /(rad_unshifted * rad(i))
      
         !cos_theta=1.d0 
          if(i>1 .and. i<nz+1)              then
            dirn_coeff=(grid_pty(j1)*(grid_pty(j1) + sx))/(abs(grid_pty(j1)*abs(grid_pty(j1) + sx)))
            
            print*,i,dirn_coeff
            !output(i+side_tr)=dk0drf1d(i)*cos_theta*abs(grid_pty(j1))/rad_unshifted
            output(i+side_tr)=dk0drf1d(i)*dirn_coeff*abs(grid_pty(j1)+sx)/rad(i)
            
          elseif(i>nz+1 .and. i<(nz+ny+1))   then
            dirn_coeff=(grid_ptz(k1)*(grid_ptz(k1) + sy))/(abs(grid_ptz(k1))*abs(grid_ptz(k1) + sy))
            !output(i+side_tr)=dk0drf1d(i)*cos_theta*abs(grid_ptz(k1))/rad_unshifted
            output(i+side_tr)=dk0drf1d(i)*dirn_coeff*abs(grid_ptz(k1)+sy)/rad(i)
            print*,i,dirn_coeff       
          elseif(i>(nz+ny+1) .and. i<(2*nz+ny+1))  then
             
            dirn_coeff=(grid_pty(j1)*(grid_pty(j1) + sx))/(abs(grid_pty(j1))*abs(grid_pty(j1) + sx))
            !output(i+side_tr)=dk0drf1d(i)*cos_theta*abs(grid_pty(j1))/rad_unshifted 
            output(i+side_tr)=dk0drf1d(i)*dirn_coeff*abs(grid_pty(j1)+sx)/rad(i)
            print*,i,dirn_coeff       
          elseif (i>(2*nz+ny+1) .and. i<=nb)     then
          
            dirn_coeff=(grid_ptz(k1)*(grid_ptz(k1) + sy))/(abs(grid_ptz(k1))*abs(grid_ptz(k1) + sy))
            !output(i+side_tr)=dk0drf1d(i)*cos_theta*abs(grid_ptz(k1))/rad_unshifted
            output(i+side_tr)=dk0drf1d(i)*dirn_coeff*abs(grid_ptz(k1)+sy)/rad(i)
            print*,i,dirn_coeff
          else if(corner(i,ny,nz)) then
            
                       
            !output(i+side_tr)=dk0drb1d(i)*cos_theta*((1 - tet)*abs(grid_ptz(kb(i,ny,nz))) + tet*abs(grid_pty(jb(i,ny,nz))))/rad_unshifted 
            dirn_coeff=tet*(grid_pty(j1)*(grid_pty(j1) + sx))/(abs(grid_pty(j1))*abs(grid_pty(j1) + sx)) &
                       +(1-tet)*(grid_ptz(k1)*(grid_ptz(k1) + sy))/(abs(grid_ptz(k1))*abs(grid_ptz(k1) + sy))
            print*,i,dirn_coeff           
            output(i+side_tr)=dk0drb1d(i)*((1 - tet)*dirn_coeff*abs(grid_ptz(kb(i,ny,nz))+sy) + tet*dirn_coeff*abs(grid_pty(jb(i,ny,nz))+sx))/rad(i)
            
            !output(i+side_tr +1)=dk0drf1d(i)*cos_theta*(tet*abs(grid_ptz(kb(i,ny,nz))) + (1 - tet)*abs(grid_pty(jb(i,ny,nz))))/rad_unshifted
            dirn_coeff=(1-tet)*(grid_pty(j1)*(grid_pty(j1) + sx))/(abs(grid_pty(j1))*abs(grid_pty(j1) + sx)) &
                       +(tet)*(grid_ptz(k1)*(grid_ptz(k1) + sy))/(abs(grid_ptz(k1))*abs(grid_ptz(k1) + sy))
            print*,i,dirn_coeff           
            output(i+side_tr +1)=dk0drf1d(i)*(tet*dirn_coeff*abs(grid_ptz(kb(i,ny,nz))+sy) + (1 - tet)*dirn_coeff*abs(grid_pty(jb(i,ny,nz))+sx))/rad(i)
            
            side_tr=side_tr+1
            
            if(tet==0) then
             tet=1
            else
             tet=0
            end if
            
          endif

        
       enddo
         
          
   
      return
       
      end subroutine dk0nor1d
!    
! ****************************************************
!    Mcdonald function K0 evaluation for a 2D array
! ****************************************************
      subroutine bessk0(n,alpha,rad,output)
      
      implicit none
      integer*4,intent(in)                ::n
      real*8,intent(in)                   ::alpha
      real*8,intent(in)                   ::rad(1:n,1:n)
      real*8,allocatable,intent(out)      ::output(:,:)

      ! local variables
      real*8::r,z
      integer*4::i,m

      allocate(output(1:n,1:n))

      if(alpha==0.d0) then

        do i=1,n
          do m=1,n
            output(i,m) = -1.0*log(rad(i,m))
          enddo
        enddo
 
      else

        do i=1,n
          do m=1,n

	    r=alpha*rad(i,m)

            if(r<=2.0d0) then
              z=(r*r)/4.0
              output(i,m)=(-log(r/2.0)*bessi0(r))+(-0.57721566 + z*(0.42278420 + z*(0.23069756 + z*(0.348859e-1 + z*(0.262698e-2 + z*(0.10750e-3 + z*0.74e-5))))))
            else
              z=2.0/r   
	      output(i,m)=(exp(-r)/sqrt(r))*(1.25331414 + z*(-0.7832358e-1 + z*(0.2189568e-1 + z*(-0.1062446e-1 + z*(0.587872e-2 + z*(-0.251540e-2 + z*0.53208e-3))))))
            endif

          enddo
        enddo

      endif

      return
      end subroutine bessk0

! ****************************************************
!    Mcdonald function K0(r) evaluation for node r as radial coordinate of boundary nodes
! ****************************************************
      subroutine act_func(n,alpha,rad,output)
      
      implicit none
      integer*4,intent(in)                ::n
      real*8,intent(in)                   ::rad(1:n),alpha
      real*8,allocatable,intent(out)      ::output(:)

      ! local variables
      real*8::r,z
      integer*4::i

      allocate(output(1:n))

     
 
    

        do i=1,n
         

	    r=alpha*rad(i)

            if(r<=2.0d0) then
              z=(r*r)/4.0
              output(i)=(-log(r/2.0)*bessi0(r))+(-0.57721566 + z*(0.42278420 + z*(0.23069756 + z*(0.348859e-1 + z*(0.262698e-2 + z*(0.10750e-3 + z*0.74e-5))))))
            else
              z=2.0/r   
	      output(i)=(exp(-r)/sqrt(r))*(1.25331414 + z*(-0.7832358e-1 + z*(0.2189568e-1 + z*(-0.1062446e-1 + z*(0.587872e-2 + z*(-0.251540e-2 + z*0.53208e-3))))))
            endif

         
        enddo

     

      return
      end subroutine act_func



! ****************************************************
!    Mcdonald function K0(r) evaluation for node r as radial coordinate of boundary nodes
! ****************************************************
     subroutine btau_func(n,ny,nz,alpha,rad,psi,output,h)
      
      implicit none
      integer*4,intent(in)                ::n,ny,nz
      real*8,intent(in)                   ::rad(1:n),h(1:n)
      real*8,intent(in)                   ::alpha
      real*8,allocatable,intent(out)      ::output(:)
      real*8,allocatable,intent(in)       ::psi(:)
      ! local variables
      real*8                              ::r,z,e
      integer*4                           ::i
      integer                             ::n_size,side
      
      n_size=2*(ny+1)+2*(nz+1)
      allocate(output(1:n_size))
      side=0
     
 
    

        do i=1,n
          if(corner(i,ny,nz)) then
           e=h(m_indx(i-2,n))/h(m_indx(i-1,n))
           output(i+side)=-1.d0*(-(2*e+e**2)/(h(m_indx(i-1,n)))/(e+e**2)*psi(i) + (1+e)**2/(h(m_indx(i-1,n)))/(e+e**2)*psi(m_indx(i-1,n))-1.0/(h(m_indx(i-1,n)))/(e+e**2)*psi(m_indx(i-2,n)))
           e=h(m_indx(i+1,n))/h(m_indx(i,n))
           output(i+side+1)=1.d0*(-(2*e+e**2)/(h(m_indx(i,n)))/(e+e**2)*psi(i) + (1+e)**2/(h(m_indx(i,n)))/(e+e**2)*psi(m_indx(i+1,n))-1.0/(h(m_indx(i,n)))/(e+e**2)*psi(m_indx(i+2,n)))

        
           side=side+1
           
          else
           output(i+side)=(psi(m_indx(i+1,n))-psi(i))/h(i)
	   
	  end if
        enddo

     

      return
      end subroutine btau_func




!**************************************************************************
!Subroutine for fetching Gauss legendre x and w given number of points
!************************************************************************** 
subroutine legendre_gauss(n, x, w)

integer, intent(in) :: n
real(kind=8),allocatable, intent(out) :: x(:)
real(kind=8),allocatable, intent(out) :: w(:)



allocate(x(1:n),w(1:n))








select case(n)
case(2)
 x(1) = -0.57735026918962576451d0; x(2) = -x(1)
 w(1) = 1.0d0; w(2) = 1.0d0

return
case(3)
 x(1) = -0.77459666924148337704d0; x(2) = 0.0d0; x(3) = -x(1)
 w(1) = 0.55555555555555555556d0; w(2) = 0.88888888888888888889d0; w(3) = w(1)

return
case(4)
 x(1) = -0.86113631159405257522d0; x(2) = -0.33998104358485626480d0
 x(3) = -x(2); x(4) = -x(1)
 w(1) = 0.34785484513745385737d0; w(2) = 0.65214515486254614263d0
 w(3) = w(2); w(4) = w(1)


return
case(5)
 x(1) = -0.90617984593866399280d0; x(2) = -0.53846931010568309103d0
 x(3) = 0.0d0; x(4) = -x(2); x(5) = -x(1)
 w(1) = 0.23692688505618908751d0; w(2) = 0.47862867049936646804d0
 w(3) = 0.56888888888888888889d0; w(4) = w(2); w(5) = w(1)


return
case(6)
 x(1) = -0.93246951420315202781d0; x(2) = -0.66120938646626451366d0
 x(3) = -0.23861918608319690863d0; x(4) = -x(3); x(5) = -x(2); x(6) = -x(1)
 w(1) = 0.17132449237917034504d0; w(2) = 0.36076157304813860757d0
 w(3) = 0.46791393457269104739d0; w(4) = w(3); w(5) = w(2); w(6) = w(1)


return
case(7)
 x(1) = -0.94910791234275852453d0; x(2) = -0.74153118559939443986d0
 x(3) = -0.40584515137739716691d0; x(4) = 0.0d0
 x(5) = -x(3); x(6) = -x(2); x(7) = -x(1)
 w(1) = 0.12948496616886969327d0; w(2) = 0.27970539148927666790d0
 w(3) = 0.38183005050511894495d0; w(4) = 0.41795918367346938776d0
 w(5) = w(3); w(6) = w(2); w(7) = w(1)

return
case(8)
 x(1) = -0.96028985649753623168d0; x(2) = -0.79666647741362673959d0
 x(3) = -0.52553240991632898581d0; x(4) = -0.18343464249564980494d0
 x(5) = -x(4); x(6) = -x(3); x(7) = -x(2); x(8) = -x(1)
 w(1) = 0.10122853629037625915d0; w(2) = 0.22238103445337447054d0
 w(3) = 0.31370664587788728734d0; w(4) = 0.36268378337836198297d0
 w(5) = w(4); w(6) = w(3); w(7) = w(2); w(8) = w(1)

return
case default
print*,"data for given no. of points is not available"
end select



end subroutine legendre_gauss

! *************************************************************************
! To compute the radial coordinate matrix for gauss legendre point in each element
!*************************************************************************
subroutine create_6rad_mat(ny,nz,nb,gridpts_y,gridpts_z,rad1,rad2,rad3,rad4,rad5,rad6)
  integer , intent(in)                          ::ny,nz,nb
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  , intent(out)                         ::rad1(1:nb,1:nb),rad2(1:nb,1:nb),rad3(1:nb,1:nb),rad4(1:nb,1:nb),rad5(1:nb,1:nb),rad6(1:nb,1:nb)

  integer*4                                     ::i,m 
  real*8                                        ::y_disp1,z_disp1
  real*8                                        ::y_disp2,z_disp2
  real*8                                        ::y_disp3,z_disp3
  real*8                                        ::y_disp4,z_disp4
  real*8                                        ::y_disp5,z_disp5
  real*8                                        ::y_disp6,z_disp6
  real*8,allocatable                            ::x(:),t(:),w(:)         ! x =[-1,1] and t=[r(m),r(m+1)]
  
  call legendre_gauss(6, x, w)
  allocate(t(1:size(x)))
  do i=1,nb
    do m=1,nb
      if((m<nz+1) .or. (m>=ny+nz+1 .and. m<2*nz+ny+1)) then
       t(1:)=(gridpts_z(kb(m+1,ny,nz))-gridpts_z(kb(m,ny,nz)))/2*x(1:) + (gridpts_z(kb(m+1,ny,nz))+gridpts_z(kb(m,ny,nz)))/2
       y_disp=gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))-t(6)
       rad1(i,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i,m)=sqrt(y_disp**2 + z_disp6**2)
      else if (( m>=nz+1 .and. m< ny+nz+1) .or. (m>=2*nz+ny+1  .and. m<=nb)) then
       t(1:)=(gridpts_y(jb(m+1,ny,nz))-gridpts_y(jb(m,ny,nz)))/2*x(1:) + (gridpts_y(jb(m+1,ny,nz))+gridpts_y(jb(m,ny,nz)))/2
       y_disp1=gridpts_y(jb(i,ny,nz))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m,ny,nz))
       
       rad1(i,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i,m)=sqrt(y_disp6**2 + z_disp**2)
      end if
      
      
    end do
  end do 
  
  print*,"Radial coordinate 2d array for 4 interior gauss points created"
end subroutine create_6rad_mat



! *************************************************************************
! To compute the radial coordinate matrix for gauss legendre point in each element
!*************************************************************************
!Creating radial coordinate matrix for subpoles with 2 subpoles in the corner element
subroutine create_6rad_mat_extra(ny,nz,nb,h,gridpts_y,gridpts_z,rad1,rad2,rad3,rad4,rad5,rad6)
  integer , intent(in)                          ::ny,nz,nb
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  , intent(out)                         ::rad1(1:nb+8,1:nb),rad2(1:nb+8,1:nb),rad3(1:nb+8,1:nb),rad4(1:nb+8,1:nb),rad5(1:nb+8,1:nb),rad6(1:nb+8,1:nb)
  real*8  , intent(in)                          ::h(1:nb)
  integer*4                                     ::i,m 
  real*8                                        ::y_disp1,z_disp1
  real*8                                        ::y_disp2,z_disp2
  real*8                                        ::y_disp3,z_disp3
  real*8                                        ::y_disp4,z_disp4
  real*8                                        ::y_disp5,z_disp5
  real*8                                        ::y_disp6,z_disp6
  real*8,allocatable                            ::x(:),t(:),w(:)         ! x =[-1,1] and t=[r(m),r(m+1)]
  integer                                       ::side_tr
  call legendre_gauss(6, x, w)
  allocate(t(1:size(x)))
  side_tr=0
  do i=1,nb
    if(corner(m_indx(i,nb),ny,nz)) then
      do m=1,nb
      if((m<nz+1) .or. (m>=ny+nz+1 .and. m<2*nz+ny+1)) then
       t(1:)=(gridpts_z(kb(m+1,ny,nz))-gridpts_z(kb(m,ny,nz)))/2*x(1:) + (gridpts_z(kb(m+1,ny,nz))+gridpts_z(kb(m,ny,nz)))/2
       
       y_disp=gridpts_y(jb(i,ny,nz))- 0.25*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))- 0.25*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))- 0.25*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))- 0.25*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))- 0.25*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))- 0.25*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))- 0.25*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(6)
       rad1(i+side_tr,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr,m)=sqrt(y_disp**2 + z_disp6**2)
       
       y_disp=gridpts_y(jb(i,ny,nz))+ 0.25*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))+ 0.25*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))+ 0.25*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))+ 0.25*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))+ 0.25*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))+ 0.25*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))+ 0.25*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(6)
       rad1(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp6**2)
       
       y_disp=gridpts_y(jb(i,ny,nz))+ 0.5*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))+ 0.5*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))+ 0.5*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))+ 0.5*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))+ 0.5*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))+ 0.5*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))+ 0.5*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(6)
       rad1(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp6**2)
       
      else if (( m>=nz+1 .and. m< ny+nz+1) .or. (m>=2*nz+ny+1  .and. m<=nb)) then
       t(1:)=(gridpts_y(jb(m+1,ny,nz))-gridpts_y(jb(m,ny,nz)))/2*x(1:) + (gridpts_y(jb(m+1,ny,nz))+gridpts_y(jb(m,ny,nz)))/2

       y_disp1=gridpts_y(jb(i,ny,nz))-0.25*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))-0.25*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))-0.25*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))-0.25*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))-0.25*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))-0.25*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))-0.25*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-gridpts_z(kb(m,ny,nz))
       rad1(i+side_tr,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr,m)=sqrt(y_disp6**2 + z_disp**2)
       
       y_disp1=gridpts_y(jb(i,ny,nz))+0.25*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))+0.25*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))+0.25*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))+0.25*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))+0.25*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))+0.25*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))+0.25*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-gridpts_z(kb(m,ny,nz))
       rad1(i+side_tr+1,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr+1,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr+1,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr+1,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr+1,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr+1,m)=sqrt(y_disp6**2 + z_disp**2)
       
       y_disp1=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))+0.5*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-gridpts_z(kb(m,ny,nz))
       rad1(i+side_tr+2,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr+2,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr+2,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr+2,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr+2,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr+2,m)=sqrt(y_disp6**2 + z_disp**2)
      end if
      
      
    end do
     side_tr=side_tr+2
    else
    do m=1,nb
      if((m<nz+1) .or. (m>=ny+nz+1 .and. m<2*nz+ny+1)) then
       t(1:)=(gridpts_z(kb(m+1,ny,nz))-gridpts_z(kb(m,ny,nz)))/2*x(1:) + (gridpts_z(kb(m+1,ny,nz))+gridpts_z(kb(m,ny,nz)))/2
       y_disp=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))+0.5*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))+0.5*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))+0.5*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))+0.5*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))+0.5*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))+0.5*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(6)
       rad1(i+side_tr,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr,m)=sqrt(y_disp**2 + z_disp6**2)
      else if (( m>=nz+1 .and. m< ny+nz+1) .or. (m>=2*nz+ny+1  .and. m<=nb)) then
       t(1:)=(gridpts_y(jb(m+1,ny,nz))-gridpts_y(jb(m,ny,nz)))/2*x(1:) + (gridpts_y(jb(m+1,ny,nz))+gridpts_y(jb(m,ny,nz)))/2
       y_disp1=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))+0.5*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))+0.5*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-gridpts_z(kb(m,ny,nz))
       
       rad1(i+side_tr,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr,m)=sqrt(y_disp6**2 + z_disp**2)
      end if
      
      
    end do
    end if
  end do 
  
  print*,"Radial coordinate 2d array for 4 interior gauss points created"
end subroutine create_6rad_mat_extra

!Creating radial coordinate matrix for subpoles with 3 subpoles in the corner element
subroutine create_6rad_mat_extra2(ny,nz,nb,h,gridpts_y,gridpts_z,rad1,rad2,rad3,rad4,rad5,rad6)
  integer , intent(in)                          ::ny,nz,nb
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  , intent(out)                         ::rad1(1:nb+16,1:nb),rad2(1:nb+16,1:nb),rad3(1:nb+16,1:nb),rad4(1:nb+16,1:nb),rad5(1:nb+16,1:nb),rad6(1:nb+16,1:nb)
  real*8  , intent(in)                          ::h(1:nb)
  integer*4                                     ::i,m 
  real*8                                        ::y_disp1,z_disp1
  real*8                                        ::y_disp2,z_disp2
  real*8                                        ::y_disp3,z_disp3
  real*8                                        ::y_disp4,z_disp4
  real*8                                        ::y_disp5,z_disp5
  real*8                                        ::y_disp6,z_disp6
  real*8,allocatable                            ::x(:),t(:),w(:)         ! x =[-1,1] and t=[r(m),r(m+1)]
  integer                                       ::side_tr
  real*8,parameter                              ::c1=0.25,c2=0.001,c3=0.5
 
  call legendre_gauss(6, x, w)
  allocate(t(1:size(x)))
  side_tr=0
  do i=1,nb
    if(corner(m_indx(i,nb),ny,nz)) then
      do m=1,nb
      if((m<nz+1) .or. (m>=ny+nz+1 .and. m<2*nz+ny+1)) then
       t(1:)=(gridpts_z(kb(m+1,ny,nz))-gridpts_z(kb(m,ny,nz)))/2*x(1:) + (gridpts_z(kb(m+1,ny,nz))+gridpts_z(kb(m,ny,nz)))/2
       
       y_disp=gridpts_y(jb(i,ny,nz))- c1*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))- c1*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))- c1*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))- c1*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))- c1*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))- c1*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))- c1*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(6)
       rad1(i+side_tr,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr,m)=sqrt(y_disp**2 + z_disp6**2)
       
       
       y_disp=gridpts_y(jb(i,ny,nz))- c2*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))- c2*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))- c2*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))- c2*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))- c2*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))- c2*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))- c2*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-t(6)
       rad1(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr+1,m)=sqrt(y_disp**2 + z_disp6**2)
       
       y_disp=gridpts_y(jb(i,ny,nz))+ c2*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))+ c2*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))+ c2*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))+ c2*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))+ c2*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))+ c2*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))+ c2*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(6)
       rad1(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr+2,m)=sqrt(y_disp**2 + z_disp6**2)
       
       y_disp=gridpts_y(jb(i,ny,nz))+ c1*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))+ c1*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))+ c1*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))+ c1*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))+ c1*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))+ c1*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))+ c1*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(6)
       rad1(i+side_tr+3,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr+3,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr+3,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr+3,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr+3,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr+3,m)=sqrt(y_disp**2 + z_disp6**2)
       
       y_disp=gridpts_y(jb(i,ny,nz))+ c3*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))+ c3*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))+ c3*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))+ c3*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))+ c3*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))+ c3*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))+ c3*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-t(6)
       rad1(i+side_tr+4,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr+4,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr+4,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr+4,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr+4,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr+4,m)=sqrt(y_disp**2 + z_disp6**2)
       
      else if (( m>=nz+1 .and. m< ny+nz+1) .or. (m>=2*nz+ny+1  .and. m<=nb)) then
       t(1:)=(gridpts_y(jb(m+1,ny,nz))-gridpts_y(jb(m,ny,nz)))/2*x(1:) + (gridpts_y(jb(m+1,ny,nz))+gridpts_y(jb(m,ny,nz)))/2

       y_disp1=gridpts_y(jb(i,ny,nz))-c1*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))-c1*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))-c1*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))-c1*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))-c1*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))-c1*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))-c1*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-gridpts_z(kb(m,ny,nz))
       rad1(i+side_tr,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr,m)=sqrt(y_disp6**2 + z_disp**2)
       
       y_disp1=gridpts_y(jb(i,ny,nz))-c2*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))-c2*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))-c2*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))-c2*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))-c2*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))-c2*(gridpts_y(jb(i,ny,nz))-gridpts_y(jb(m_indx(i-1,nb),ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))-c2*(gridpts_z(kb(i,ny,nz))-gridpts_z(kb(m_indx(i-1,nb),ny,nz)))-gridpts_z(kb(m,ny,nz))
       rad1(i+side_tr+1,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr+1,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr+1,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr+1,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr+1,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr+1,m)=sqrt(y_disp6**2 + z_disp**2)
       
       y_disp1=gridpts_y(jb(i,ny,nz))+c2*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))+c2*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))+c2*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))+c2*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))+c2*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))+c2*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))+c2*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-gridpts_z(kb(m,ny,nz))
       rad1(i+side_tr+2,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr+2,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr+2,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr+2,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr+2,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr+2,m)=sqrt(y_disp6**2 + z_disp**2)
       
       y_disp1=gridpts_y(jb(i,ny,nz))+c1*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))+c1*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))+c1*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))+c1*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))+c1*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))+c1*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))+c1*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-gridpts_z(kb(m,ny,nz))
       rad1(i+side_tr+3,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr+3,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr+3,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr+3,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr+3,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr+3,m)=sqrt(y_disp6**2 + z_disp**2)
       
       y_disp1=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(m_indx(i+1,nb),ny,nz))-gridpts_y(jb(i,ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))+c3*(gridpts_z(kb(m_indx(i+1,nb),ny,nz))-gridpts_z(kb(i,ny,nz)))-gridpts_z(kb(m,ny,nz))
       rad1(i+side_tr+4,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr+4,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr+4,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr+4,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr+4,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr+4,m)=sqrt(y_disp6**2 + z_disp**2)
      end if
      
      
    end do
     side_tr=side_tr+4
    else
    do m=1,nb
      if((m<nz+1) .or. (m>=ny+nz+1 .and. m<2*nz+ny+1)) then
       t(1:)=(gridpts_z(kb(m+1,ny,nz))-gridpts_z(kb(m,ny,nz)))/2*x(1:) + (gridpts_z(kb(m+1,ny,nz))+gridpts_z(kb(m,ny,nz)))/2
       y_disp=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-gridpts_y(jb(m,ny,nz))
       z_disp1=gridpts_z(kb(i,ny,nz))+c3*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(1)
       z_disp2=gridpts_z(kb(i,ny,nz))+c3*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(2)
       z_disp3=gridpts_z(kb(i,ny,nz))+c3*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(3)
       z_disp4=gridpts_z(kb(i,ny,nz))+c3*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(4)
       z_disp5=gridpts_z(kb(i,ny,nz))+c3*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(5)
       z_disp6=gridpts_z(kb(i,ny,nz))+c3*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-t(6)
       rad1(i+side_tr,m)=sqrt(y_disp**2 + z_disp1**2)
       rad2(i+side_tr,m)=sqrt(y_disp**2 + z_disp2**2)
       rad3(i+side_tr,m)=sqrt(y_disp**2 + z_disp3**2)
       rad4(i+side_tr,m)=sqrt(y_disp**2 + z_disp4**2)
       rad5(i+side_tr,m)=sqrt(y_disp**2 + z_disp5**2)
       rad6(i+side_tr,m)=sqrt(y_disp**2 + z_disp6**2)
      else if (( m>=nz+1 .and. m< ny+nz+1) .or. (m>=2*nz+ny+1  .and. m<=nb)) then
       t(1:)=abs(gridpts_y(jb(m+1,ny,nz))-gridpts_y(jb(m,ny,nz)))/2*x(1:) + (gridpts_y(jb(m+1,ny,nz))+gridpts_y(jb(m,ny,nz)))/2
       y_disp1=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(1)
       y_disp2=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(2)
       y_disp3=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(3)
       y_disp4=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(4)
       y_disp5=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(5)
       y_disp6=gridpts_y(jb(i,ny,nz))+c3*(gridpts_y(jb(i+1,ny,nz))-gridpts_y(jb(i,ny,nz)))-t(6)
       z_disp=gridpts_z(kb(i,ny,nz))+c3*(gridpts_z(kb(i+1,ny,nz))-gridpts_z(kb(i,ny,nz)))-gridpts_z(kb(m,ny,nz))
       
       rad1(i+side_tr,m)=sqrt(y_disp1**2 + z_disp**2)
       rad2(i+side_tr,m)=sqrt(y_disp2**2 + z_disp**2)
       rad3(i+side_tr,m)=sqrt(y_disp3**2 + z_disp**2)
       rad4(i+side_tr,m)=sqrt(y_disp4**2 + z_disp**2)
       rad5(i+side_tr,m)=sqrt(y_disp5**2 + z_disp**2)
       rad6(i+side_tr,m)=sqrt(y_disp6**2 + z_disp**2)
      end if
      
      
    end do
    end if
  end do 
  
  print*,"Radial coordinate 2d array for 4 interior gauss points created"
end subroutine create_6rad_mat_extra2
! ************************************************************************************************************************************************************************************************
!    To compute i1=integral of k0dl and i2=integral of l.k0dl that is used in the computation of bterm2, accounting for the function switch that can occur within the element for high wavenumbers
! ********************************************************************************************************************************************************************************************************
      subroutine integrals2(n,h,alpha,i1,i2) !output in array form
      
      implicit none
      integer*4,intent(in) :: n
      real*8,intent(in) :: alpha,h(1:n)
      real*8,intent(out) :: i1(1:n),i2(1:n)
      ! local variables
      real,parameter :: PI=4*atan(1.d0)
      integer*4 :: i,k
      real*8 :: hlft
      real*8 :: r,z
      real*8 :: k0lft,k0hi
      real*8 :: const0(1:6),const1(1:7)

      const0(1:6) = (/3.5156229,3.0899424,1.2067492,0.2659732,0.360768e-1,0.45813e-2/)
      const1(1:7)=(/-0.57721566,0.42278420,0.23069756,0.348859e-1,0.262698e-2,0.10750e-3,0.74e-5/)

      if(alpha==0) then

        do i=1,n
          i1(i) = h(i)*(1.0d0-log(h(i))) 
          i2(i) = 0.25*(h(i)**2)*(1.0d0-2.0*log(h(i)))
        enddo

      else

        do i=1,n

          if(alpha*h(i) > 2.0d0) then
            hlft = (2.0/alpha)
          else
            hlft = h(i)
          endif

          i1(i) = hlft*(1.0d0-log(0.5*alpha*hlft)) 
          i2(i) = 0.25*(hlft**2)*(1.0d0-2.0*log(0.5*alpha*hlft))
  
          do k=1,6
            i1(i) = i1(i) + const0(k) * ((alpha/3.75)**(2*k)) * (hlft**(2*k+1)) * (1.0d0 - (2*k+1)*log(0.5*alpha*hlft)) / ((2*k+1)**2)
	    i2(i) = i2(i) + const0(k) * ((alpha/3.75)**(2*k)) * (hlft**(2*k+2)) * (1.0d0 - (2*k+2)*log(0.5*alpha*hlft)) / ((2*k+2)**2)
          enddo

          do k=1,7
            i1(i) = i1(i) + const1(k) * (0.5*alpha)**(2*k-2) * (hlft**(2*k-1)) / (2*k-1)
            i2(i) = i2(i) + const1(k) * (0.5*alpha)**(2*k-2) * (hlft**(2*k)) / (2*k)
          enddo

          ! Remaining part of the integral is added only if the function changed form within the element h(i)
          if(alpha*h(i) > 2.0d0) then

	     r = alpha*hlft
             z = 2.0/r   
	     k0lft = (exp(-r)/sqrt(r))*(1.25331414 + z*(-0.7832358e-1 + z*(0.2189568e-1 + z*(-0.1062446e-1 + z*(0.587872e-2 + z*(-0.251540e-2 + z*0.53208e-3))))))

	     r = alpha*h(i)
             z = 2.0/r   
	     k0hi = (exp(-r)/sqrt(r))*(1.25331414 + z*(-0.7832358e-1 + z*(0.2189568e-1 + z*(-0.1062446e-1 + z*(0.587872e-2 + z*(-0.251540e-2 + z*0.53208e-3))))))

             i1(i) = i1(i) + 0.5*(h(i)-hlft)*(k0lft + k0hi)
             i2(i) = i2(i) + 0.5*(h(i)-hlft)*(hlft*k0lft + h(i)*k0hi)

          endif

        
        enddo

      endif

      end subroutine integrals2
 
     subroutine integrals2h(alpha,i1,i2,u_lim) !output in single real value form
      
      implicit none
      real*8,intent(in)                          :: alpha,u_lim
      real*8,intent(out)                         :: i1,i2
      ! local variables
      real,parameter                             :: PI=4*atan(1.d0)
      integer*4                                  :: i,k
      real*8                                     :: hlft
      real*8                                     :: r,z
      real*8                                     :: k0lft,k0hi
      real*8                                     :: const0(1:6),const1(1:7)

      const0(1:6) = (/3.5156229,3.0899424,1.2067492,0.2659732,0.360768e-1,0.45813e-2/)
      const1(1:7)=(/-0.57721566,0.42278420,0.23069756,0.348859e-1,0.262698e-2,0.10750e-3,0.74e-5/)

      if(alpha==0) then

    
          i1 = u_lim*(1.0d0-log(u_lim)) 
          i2 = 0.25*(u_lim**2)*(1.0d0-2.0*log(u_lim))
   

      else

      

          if(alpha*u_lim > 2.0d0) then
            hlft = (2.0/alpha)
          else
            hlft = u_lim
          endif

          i1 = hlft*(1.0d0-log(0.5*alpha*hlft)) 
          i2 = 0.25*(hlft**2)*(1.0d0-2.0*log(0.5*alpha*hlft))
  
          do k=1,6
            i1 = i1 + const0(k) * ((alpha/3.75)**(2*k)) * (hlft**(2*k+1)) * (1.0d0 - (2*k+1)*log(0.5*alpha*hlft)) / ((2*k+1)**2)
	    i2 = i2 + const0(k) * ((alpha/3.75)**(2*k)) * (hlft**(2*k+2)) * (1.0d0 - (2*k+2)*log(0.5*alpha*hlft)) / ((2*k+2)**2)
          enddo

          do k=1,7
            i1 = i1 + const1(k) * (0.5*alpha)**(2*k-2) * (hlft**(2*k-1)) / (2*k-1)
            i2 = i2 + const1(k) * (0.5*alpha)**(2*k-2) * (hlft**(2*k)) / (2*k)
          enddo

          ! Remaining part of the integral is added only if the function changed form within the element h(i)
          if(alpha*u_lim > 2.0d0) then

	     r = alpha*hlft
             z = 2.0/r   
	     k0lft = (exp(-r)/sqrt(r))*(1.25331414 + z*(-0.7832358e-1 + z*(0.2189568e-1 + z*(-0.1062446e-1 + z*(0.587872e-2 + z*(-0.251540e-2 + z*0.53208e-3))))))

	     r = alpha*u_lim
             z = 2.0/r   
	     k0hi = (exp(-r)/sqrt(r))*(1.25331414 + z*(-0.7832358e-1 + z*(0.2189568e-1 + z*(-0.1062446e-1 + z*(0.587872e-2 + z*(-0.251540e-2 + z*0.53208e-3))))))

             i1 = i1 + 0.5*(u_lim-hlft)*(k0lft + k0hi)
             i2 = i2 + 0.5*(u_lim-hlft)*(hlft*k0lft + u_lim*k0hi)

          endif

        
    

      endif

      end subroutine integrals2h

 !recursive subroutine to calculate analytical integration of 
 !functions in the form x^n/sqrt(x^2-a^2) and ln(x)x^n/sqrt(x^2-a^2)
  recursive real*8 function analytical_int(xf,xi,const,alpha,n,funct_type)

      implicit none
      real*8,intent(in)                          ::xf,xi,alpha,const
      integer*4,intent(in)                       ::n
      integer*4,intent(in)                       ::funct_type
      real*8                                     ::tf,ti
      select case(funct_type)
      
      case (0) !x^n/sqrt(x^2-a^2)
         
       if(n==0)then 
         analytical_int= log(xf+sqrt(xf**2 - const**2)) - log(xi+sqrt(xi**2 - const**2))
       else if(n==1)then  
         analytical_int= sqrt(xf**2 - const**2)-sqrt(xi**2 - const**2) 
       else
         analytical_int= xf**(n-1)*sqrt(xf**2 - const**2)/n - xi**(n-1)*sqrt(xi**2 - const**2)/n + (n-1)/n*const**2*analytical_int(xf,xi,const,alpha,n-2,0)
       end if 
      
      case (1)  !ln(x)x^n/sqrt(x^2-a^2)
       
       if(n==0)then 
         tf=acosh(xf/const)
         ti=acosh(xi/const)
         analytical_int=(tf-ti)*log(const)+(tf**3-ti**3)/6-(tf**5-ti**5)/60   
       else if(n==1)then  
         analytical_int=sqrt(xf**2 - const**2)*(log(xf)-1)-sqrt(xi**2 - const**2)*(log(xi)-1)+const*atan(sqrt(xf**2 - const**2)/const)-const*atan(sqrt(xi**2 - const**2)/const)
       else
         analytical_int=xf**(n-1)*sqrt(xf**2 - const**2)/n*(log(xf)-1/n)-xi**(n-1)*sqrt(xi**2 - const**2)/n*(log(xi)-1/n) + (const/n)**2*analytical_int(xf,xi,const,alpha,n-2,0) + const**2*(n-1)/n*analytical_int(xf,xi,const,alpha,n-2,1)
       end if 
         
     end select

    return
    end function analytical_int       
! ************************************************************************************************************************************************************************************************
!    To compute i1=integral of (dk0/dtau)dl and i2=integral of l(dk0/dtau)dl that is used in the computation of bterm2, accounting for the function switch that can occur within the element for high wavenumbers
! ********************************************************************************************************************************************************************************************************
      subroutine integrals4(alpha,xi,xf,const,i1,i2)
      !analytical integration of delG/deln from subpoles in an element1 over another element2
      !such that both element1 and element2 are corner adjacent and perpendicular to ech other.
      implicit none
      real*8,intent(in)                               :: alpha,xi,xf,const
      real*8,intent(out)                              :: i1,i2
      ! local variables
      real,parameter                                  :: PI=4*atan(1.d0)
      integer*4                                       :: i,k,j
      real*8                                          :: hlft
      real*8                                          :: r,z
      real*8                                          :: k0lft,k0hi
      real*8                                          :: E0(1:6),C0(1:7),D0(1:7)

      E0(1:6) = (/3.5156229,3.0899424,1.2067492,0.2659732,0.360768e-1,0.45813e-2/)
      C0(1:7)=(/-0.57721566,0.42278420,0.23069756,0.348859e-1,0.262698e-2,0.10750e-3,0.74e-5/)
      D0(1:7)=(/1.25331414,-0.7832358e-1,0.2189568e-1,-0.1062446e-1,0.587872e-2,-0.251540e-2,0.53208e-3/)

      if(alpha==0) then

       
          !????????????????????????????
    

      else


        
       

          i1 = -(acos(const/xf)-acos(const/xi))/const
          !i2 = -(log((xf*alpha/3.75)**2)-log((xi*alpha/3.75)**2))/2
          i2 = -(log(alpha*xf/2)-log(alpha*xi/2))
          
          !integration of I0(alpha_x)/x dx
          do k=1,6
            i1 = i1 - E0(k)*(alpha/3.75)**(2*k)*analytical_int(xf,xi,const,alpha,2*k-1,0)
	    !i2 = i2 - (1/2)*E0(k) *(((xf*alpha/3.75)**(2*k))/k-((xi*alpha/3.75)**(2*k))/k)
	    !i2 = i2 - E0(k) *(((xf*alpha/3.75)**(2*k))-((xi*alpha/3.75)**(2*k))/k) 
          enddo
          
          
          
          
          do k=1,6
            i1 = i1 - E0(k)*(2*k)*(alpha/3.75)**(2*k)*log(alpha/2)*analytical_int(xf,xi,const,alpha,2*k-1,0)-E0(k)*(2*k)*(alpha/3.75)**(2*k)*analytical_int(xf,xi,const,alpha,2*k-1,1)
	    !i2 = i2 - E0(k)*(((xf*alpha/3.75)**(2*k))*(log(alpha*xf/2)-1/(2*k))-((xi*alpha/3.75)**(2*k))*(log(alpha*xi/2)-1/(2*k)))
	    i2 = i2 - E0(k)*(((xf*alpha/3.75)**(2*k))*(log(alpha*xf/2))-((xi*alpha/3.75)**(2*k))*(log(alpha*xi/2)))
          enddo
          
          !integration of sum of Cn 
          do k=2,7
            i1 = i1 + C0(k) * (k-1)*(alpha**(2*k-2))/4**(k-2) * analytical_int(xf,xi,const,alpha,2*k-3,0)/2
            i2 = i2 + C0(k) * (((alpha)**2)/4)**(k-1)*(xf**(2*k-2)-xi**(2*k-2))
          enddo

          ! Remaining part of the integral is added only if the function changed form within the element h(i)
          
          
            i1=i1*const
            i2=(i2-i1*sqrt(xi**2-const**2))*const
            
        

      endif

      end subroutine integrals4
      
   subroutine integrals5(alpha,xi,xf,const,i1,i2)
   !analytical integration of G from subpoles in an element1 over another element2
   !such that both element1 and element2 are corner adjacent and perpendicular to ech other.   
      implicit none
      real*8,intent(in)                               :: alpha,xi,xf,const
      real*8,intent(out)                              :: i1,i2
      ! local variables
      real,parameter                                  :: PI=4*atan(1.d0)
      integer*4                                       :: i,k,j
      real*8                                          :: hlft
      real*8                                          :: r,z
      real*8                                          :: k0lft,k0hi
      real*8                                          :: E0(1:6),C0(1:7),D0(1:7)

      E0(1:6) = (/3.5156229,3.0899424,1.2067492,0.2659732,0.360768e-1,0.45813e-2/)
      C0(1:7)=(/-0.57721566,0.42278420,0.23069756,0.348859e-1,0.262698e-2,0.10750e-3,0.74e-5/)
      D0(1:7)=(/1.25331414,-0.7832358e-1,0.2189568e-1,-0.1062446e-1,0.587872e-2,-0.251540e-2,0.53208e-3/)

      if(alpha==0) then

       
          !????????????????????????????
    
    
      else


        
       

          i1 = -log(alpha/2)*analytical_int(xf,xi,const,alpha,1,0)-analytical_int(xf,xi,const,alpha,1,1)
          i2 =  -(0.5*xf**2*(log(alpha*xf/2)-0.5)-0.5*xi**2*(log(alpha*xi/2)-0.5))
          
          !integration of I0(alpha_x)/x dx
          do k=1,6
            i1 = i1 - E0(k)*(alpha/3.75)**(2*k)*(log(alpha/2)*analytical_int(xf,xi,const,alpha,2*k+1,0)+analytical_int(xf,xi,const,alpha,2*k+1,1))
	    i2 = i2 - E0(k) *((alpha/3.75)**(2*k))*(xf**(2*(k+1))-xi**(2*(k+1)))/(2*(k+1)) 
          enddo
          
          
          
          
          
          
          !integration of sum of Cn 
          do k=1,7
            i1 = i1 + C0(k) * (alpha**(2*k-2))/4**(k-1) * analytical_int(xf,xi,const,alpha,2*k-1,0)
            i2 = i2 + C0(k) * (((alpha)**2)/4)**(k-1)*(xf**(2*k)-xi**(2*k))/(2*k)
          enddo

          ! Remaining part of the integral is added only if the function changed form within the element h(i)
            
            
          
            i2=(i2-i1*sqrt(xi**2-const**2))
            
            if(const==0.d0) then
             i1=0.d0
             i2=0.d0
            end if

      endif

 end subroutine integrals5      
! *************************************************************************
! To compute the Integral I1=K0dl and I2=K0dl in each element
!*************************************************************************
subroutine integral_K0(ny,nz,nb,alpha,gridpts_y,gridpts_z,h_array,rad1,rad2,rad3,rad4,rad5,rad6,I1,I2)
  implicit none
  integer , intent(in)                          ::ny,nz,nb
  real*8  , intent(in)                          ::rad1(1:nb,1:nb),rad2(1:nb,1:nb),rad3(1:nb,1:nb),rad4(1:nb,1:nb),rad5(1:nb,1:nb),rad6(1:nb,1:nb)
  real*8  , intent(in)                          ::h_array(1:nb)
  real*8  , intent(in)                          ::alpha
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  ,allocatable                          ::K01(:,:),K02(:,:),K03(:,:),K04(:,:),K05(:,:),K06(:,:)
  real*8  , intent(out)                         ::I1(1:nb,1:nb),I2(1:nb,1:nb)
  real*8                                        ::Ia_1(1:nb),Ia_2(1:nb)
  integer*4                                     ::i,m 
  real,parameter                                ::PI=4*atan(1.d0)
  real*8                                        ::y_disp1,z_disp1
  real*8                                        ::y_disp2,z_disp2
  real*8                                        ::y_disp3,z_disp3
  real*8                                        ::y_disp4,z_disp4
  real*8                                        ::y_disp5,z_disp5
  real*8                                        ::y_disp6,z_disp6
  real*8,allocatable                            ::x(:),t(:),w(:)         ! x =[-1,1] and t=[r(m),r(m+1)]
  
  call integrals2(nb,h_array,alpha,Ia_1,Ia_2)
  call bessk0(nb,alpha,rad1,K01)
  call bessk0(nb,alpha,rad2,K02)
  call bessk0(nb,alpha,rad3,K03)
  call bessk0(nb,alpha,rad4,K04)
  call bessk0(nb,alpha,rad5,K05)
  call bessk0(nb,alpha,rad6,K06)
  call legendre_gauss(6, x, w)
  allocate(t(1:size(x)))
  
      
  do i=1,nb
    do m=1,nb
      if((m<nz+1) .or. (m>=ny+nz+1 .and. m<2*nz+ny+1)) then
       if(m==m_indx(i,nb) .or. m==m_indx(i-1,nb)) then
        I1(i,m)=Ia_1(m)
        I2(i,m)=Ia_2(m)
       else
         t(1:)=(gridpts_z(kb(m+1,ny,nz))-gridpts_z(kb(m,ny,nz)))*0.5*(x(1:) + 1.0d0)
         
         I1(i,m)=0.5*h_array(m)*(w(1)*K01(i,m)+w(2)*K02(i,m)+w(3)*K03(i,m)+w(4)*K04(i,m)+w(5)*K05(i,m)+w(6)*K06(i,m))
    
         I2(i,m)=0.5*(w(1)*K01(i,m)*t(1)+w(2)*K02(i,m)*t(2)+w(3)*K03(i,m)*t(3)+w(4)*K04(i,m)*t(4)+w(5)*K05(i,m)*t(5)+w(6)*K06(i,m)*t(6))
        
        end if
      else if (( m>=nz+1 .and. m< ny+nz+1) .or. (m>=2*nz+ny+1  .and. m<=nb)) then
       if(m==m_indx(i,nb) .or. m==m_indx(i-1,nb)) then
        I1(i,m)=Ia_1(m)
        I2(i,m)=Ia_2(m)
       else
        t(1:)=(gridpts_y(jb(m+1,ny,nz))-gridpts_y(jb(m,ny,nz)))*0.5*(x(1:) + 1.0d0)
      
        I1(i,m)=h_array(m)*0.5*(w(1)*K01(i,m)+w(2)*K02(i,m)+w(3)*K03(i,m)+w(4)*K04(i,m)+w(5)*K05(i,m)+w(6)*K06(i,m))
        
       
        I2(i,m)=0.5*(w(1)*K01(i,m)*t(1)+w(2)*K02(i,m)*t(2)+w(3)*K03(i,m)*t(3)+w(4)*K04(i,m)*t(4)+w(5)*K05(i,m)*t(5)+w(6)*K06(i,m)*t(6))
     
       end if
      end if
    end do
  end do
  
  print*,'Integral of K0dl and K0ldl'
end subroutine integral_K0     

!*************************************************************************
! To compute the Integral I1=del/del_tau(int K0dl) and I2=del/del_tau(int K0ldl) in each element
!*************************************************************************
subroutine integral_dK0dtau(ny,nz,nb,alpha,gridpts_y,gridpts_z,h_array,rad1,rad2,rad3,rad4,rad5,rad6,rad1h,rad2h,rad3h,rad4h,rad5h,rad6h,I1,I2)
  implicit none
  integer , intent(in)                          ::ny,nz,nb
  real*8  , allocatable, intent(in)             ::rad1(:,:),rad2(:,:),rad3(:,:),rad4(:,:),rad5(:,:),rad6(:,:)
  real*8  , allocatable,intent(in)              ::rad1h(:,:),rad2h(:,:),rad3h(:,:),rad4h(:,:),rad5h(:,:),rad6h(:,:)
  real*8  , intent(in)                          ::h_array(1:nb)
  real*8  , intent(in)                          ::alpha
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  , allocatable                         ::K01(:,:),K02(:,:),K03(:,:),K04(:,:),K05(:,:),K06(:,:)
  real*8  , allocatable                         ::K01h(:,:),K02h(:,:),K03h(:,:),K04h(:,:),K05h(:,:),K06h(:,:)
  real*8  , allocatable                         ::dK01_dt(:,:),dK02_dt(:,:),dK03_dt(:,:),dK04_dt(:,:),dK05_dt(:,:),dK06_dt(:,:)
  real*8  , allocatable, intent(out)            ::I1(:,:),I2(:,:)
  real*8                                        ::Ia_1,Ia_2,Ia_1hlft,Ia_2hlft,Ia_1hrt,Ia_2hrt
  integer*4                                     ::i,m 
  real,parameter                                ::PI=4*atan(1.d0)
  real*8                                        ::y_disp1,z_disp1
  real*8                                        ::y_disp2,z_disp2
  real*8                                        ::y_disp3,z_disp3
  real*8                                        ::y_disp4,z_disp4
  real*8                                        ::y_disp5,z_disp5
  real*8                                        ::y_disp6,z_disp6
  real*8,allocatable                            ::x(:),t(:),w(:)         ! x =[-1,1] and t=[r(m),r(m+1)]
  real*8                                        ::e,a,b,c,h_arr                !dummy_variables 
  integer                                       ::n_size                 !no. of b_tau
  integer                                       ::side_tr,side_tr2
  real*8                                        ::I1_k0h,I2_k0h,I1_k0,I2_k0,I1_k0h_prev,I2_k0h_prev
  real*8,parameter                              ::cn=0.00001,cm=0.22,cf=0.44
  n_size=2*(ny+1)+2*(nz+1)
  side_tr=0
  side_tr2=0
 
  call bessk0(nb,alpha,rad1,K01)
  call bessk0(nb,alpha,rad2,K02)
  call bessk0(nb,alpha,rad3,K03)
  call bessk0(nb,alpha,rad4,K04)
  call bessk0(nb,alpha,rad5,K05)
  call bessk0(nb,alpha,rad6,K06)
  
  call bessk0(nb,alpha,rad1h,K01h)
  call bessk0(nb,alpha,rad2h,K02h)
  call bessk0(nb,alpha,rad3h,K03h)
  call bessk0(nb,alpha,rad4h,K04h)
  call bessk0(nb,alpha,rad5h,K05h)
  call bessk0(nb,alpha,rad6h,K06h)
  
 
  
  call legendre_gauss(6, x, w)
  allocate(t(1:size(x)))
  allocate(dK01_dt(n_size,nb),dK02_dt(n_size,nb),dK03_dt(n_size,nb),dK04_dt(n_size,nb),dK05_dt(n_size,nb),dK06_dt(n_size,nb))
  allocate(I1(1:n_size,1:nb),I2(1:n_size,1:nb))    
    
do i=1,nb
   if(corner(i,ny,nz))then
    do m=1,nb
     e=1.d0
     
       if(m==m_indx(i,nb))then
        call integrals2h(alpha,Ia_1,Ia_2,h_array(m))
        I1_k0=Ia_1
        I2_k0=Ia_2
        !call integrals2h(alpha,Ia_1hlft,Ia_2hlft,h_array(m)+0.25*h_array(m_indx(m-1,nb)))
        !call integrals2h(alpha,Ia_1hrt,Ia_2hrt,0.25*h_array(m_indx(m-1,nb)))
        !I1_k0h=Ia_1hlft-Ia_1hrt
        !I2_k0h=Ia_2hlft-Ia_2hrt
        !call integrals2h(alpha,Ia_1hlft,Ia_2hlft,h_array(m)+0.5*h_array(m_indx(m-1,nb)))
        !call integrals2h(alpha,Ia_1hrt,Ia_2hrt,0.5*h_array(m_indx(m-1,nb)))
        !I1_k0h_prev=Ia_1hlft-Ia_1hrt
        !I2_k0h_prev=Ia_2hlft-Ia_2hrt
        
        call integrals5(alpha,cm*h_array(m_indx(m-1,nb)),sqrt(h_array(m)**2+(cm*h_array(m_indx(m-1,nb)))**2),cm*h_array(m_indx(m-1,nb)),Ia_1,Ia_2)
        I1_k0h=Ia_1
        I2_k0h=Ia_2
        call integrals5(alpha,cf*h_array(m_indx(m-1,nb)),sqrt(h_array(m)**2+(cf*h_array(m_indx(m-1,nb)))**2),cf*h_array(m_indx(m-1,nb)),Ia_1,Ia_2)
        I1_k0h_prev=Ia_1
        I2_k0h_prev=Ia_2
        
     
        
        I1(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h_prev/(2*PI))
        I2(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h_prev/(2*PI))
        !I2(i+side_tr,m)=0
      
    
        call integrals2h(alpha,Ia_1,Ia_2,h_array(m))
        I1_k0=Ia_1
        I2_k0=Ia_2
        
        call integrals2h(alpha,Ia_1hlft,Ia_2hlft,(1-cm)*h_array(m_indx(m,nb)))
        call integrals2h(alpha,Ia_1hrt,Ia_2hrt,cm*h_array(m_indx(m,nb)))
        
        I1_k0h=Ia_1hlft+Ia_1hrt
        I2_k0h=Ia_2hlft+(cm*h_array(m_indx(m,nb))*Ia_1hrt-Ia_2hrt)
        
        call integrals2h(alpha,Ia_1hlft,Ia_2hlft,(1-cf)*h_array(m_indx(m,nb)))
        call integrals2h(alpha,Ia_1hrt,Ia_2hrt,cf*h_array(m_indx(m,nb)))
        
        I1_k0h_prev=Ia_1hlft+Ia_1hrt
        I2_k0h_prev=Ia_2hlft+(cf*h_array(m_indx(m,nb))*Ia_1hrt-Ia_2hrt)
        
        I1(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h_prev/(2*PI))
        I2(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h_prev/(2*PI))
      
       else if(m==m_indx(i-1,nb))then
        
        call integrals2h(alpha,Ia_1,Ia_2,h_array(m))
        I1_k0=Ia_1
        I2_k0=(h_array(m)*Ia_1-Ia_2)
        call integrals2h(alpha,Ia_1hlft,Ia_2hlft,cm*h_array(m_indx(m,nb)))
        call integrals2h(alpha,Ia_1hrt,Ia_2hrt,(1-cm)*h_array(m_indx(m,nb)))
        I1_k0h=Ia_1hlft+Ia_1hrt
        I2_k0h=Ia_2hlft+((1-cm)*h_array(m_indx(m,nb))*Ia_1hrt-Ia_2hrt)
        call integrals2h(alpha,Ia_1hlft,Ia_2hlft,cf*h_array(m_indx(m,nb)))
        call integrals2h(alpha,Ia_1hrt,Ia_2hrt,(1-cf)*h_array(m_indx(m,nb)))
        I1_k0h_prev=Ia_1hlft+Ia_1hrt
        I2_k0h_prev=Ia_2hlft+((1-cf)*h_array(m_indx(m,nb))*Ia_1hrt-Ia_2hrt)
        
        I1(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h_prev/(2*PI))
        !I1(i+side_tr,m)=1.d0*alpha/(2*PI)
        I2(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h_prev/(2*PI))
        !I2(i+side_tr,m)=1.d0*alpha/(2*PI)
        print*,"-----------------"
        print*,I1_k0,I1_k0h,I1_k0h_prev
        print*,"-----------------"
        call integrals2h(alpha,Ia_1,Ia_2,h_array(m_indx(m,nb)))
        I1_k0=Ia_1
        I2_k0=(h_array(m_indx(m,nb))*Ia_1-Ia_2)
        !call integrals2h(alpha,Ia_1hlft,Ia_2hlft,0.25*h_array(m_indx(m+1,nb))+h_array(m_indx(m,nb)))
        !call integrals2h(alpha,Ia_1hrt,Ia_2hrt,0.25*h_array(m_indx(m+1,nb)))
        !I1_k0h=Ia_1hlft-Ia_1hrt
        !I2_k0h=(Ia_1hlft*(0.25*h_array(m_indx(m+1,nb))+h_array(m_indx(m,nb)))-Ia_1hrt*(0.25*h_array(m_indx(m+1,nb))))-(Ia_2hlft-Ia_2hrt)
        !call integrals2h(alpha,Ia_1hlft,Ia_2hlft,0.5*h_array(m_indx(m+1,nb))+h_array(m_indx(m,nb)))
        !call integrals2h(alpha,Ia_1hrt,Ia_2hrt,0.5*h_array(m_indx(m+1,nb)))
        !I1_k0h_prev=Ia_1hlft-Ia_1hrt
        !I2_k0h_prev=(Ia_1hlft*(0.5*h_array(m_indx(m+1,nb))+h_array(m_indx(m,nb)))-Ia_1hrt*(0.5*h_array(m_indx(m+1,nb))))-(Ia_2hlft-Ia_2hrt)
        
        call integrals5(alpha,cm*h_array(m_indx(m+1,nb)),sqrt(h_array(m)**2+(cm*h_array(m_indx(m+1,nb)))**2),cm*h_array(m_indx(m+1,nb)),Ia_1,Ia_2)
        I1_k0h=Ia_1
        I2_k0h=(h_array(m_indx(m,nb))*Ia_1-Ia_2)
        call integrals5(alpha,cf*h_array(m_indx(m+1,nb)),sqrt(h_array(m)**2+(cf*h_array(m_indx(m+1,nb)))**2),cf*h_array(m_indx(m+1,nb)),Ia_1,Ia_2)
        I1_k0h_prev=Ia_1
        I2_k0h_prev=(h_array(m_indx(m,nb))*Ia_1-Ia_2)
        I1(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h_prev/(2*PI))
        I2(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h_prev/(2*PI))
      
       else
        t(1:)=h_array(m)*0.5*(x(1:) + 1.0d0)
        
        I1_k0h=h_array(m)*0.5*(w(1)*K01h(i+side_tr2,m)+w(2)*K02h(i+side_tr2,m)+w(3)*K03h(i+side_tr2,m)+w(4)*K04h(i+side_tr2,m)+w(5)*K05h(i+side_tr2,m)+w(6)*K06h(i+side_tr2,m))
        I2_k0h=0.5*(w(1)*K01h(i+side_tr2,m)*t(1)+w(2)*K02h(i+side_tr2,m)*t(2)+w(3)*K03h(i+side_tr2,m)*t(3)+w(4)*K04h(i+side_tr2,m)*t(4)+w(5)*K05h(i+side_tr2,m)*t(5)+w(6)*K06h(i+side_tr2,m)*t(6))
        I1_k0=h_array(m)*0.5*(w(1)*K01(m_indx(i,nb),m)+w(2)*K02(m_indx(i,nb),m)+w(3)*K03(m_indx(i,nb),m)+w(4)*K04(m_indx(i,nb),m)+w(5)*K05(m_indx(i,nb),m)+w(6)*K06(m_indx(i,nb),m))
        I2_k0=0.5*(w(1)*K01(m_indx(i,nb),m)*t(1)+w(2)*K02(m_indx(i,nb),m)*t(2)+w(3)*K03(m_indx(i,nb),m)*t(3)+w(4)*K04(m_indx(i,nb),m)*t(4)+w(5)*K05(m_indx(i,nb),m)*t(5)+w(6)*K06(i+side_tr,m)*t(6))
        I1_k0h_prev=h_array(m)*0.5*(w(1)*K01h(m_indx(i+side_tr2-1,nb+8),m)+w(2)*K02h(m_indx(i+side_tr2-1,nb+8),m)+w(3)*K03h(m_indx(i+side_tr2-1,nb+8),m)+w(4)*K04h(m_indx(i+side_tr2-1,nb+8),m)+w(5)*K05h(m_indx(i+side_tr2-1,nb+8),m)+w(6)*K06h(m_indx(i+side_tr2-1,nb+8),m))
        I2_k0h_prev=0.5*(w(1)*K01h(m_indx(i+side_tr2-1,nb+8),m)*t(1)+w(2)*K02h(m_indx(i+side_tr2-1,nb+8),m)*t(2)+w(3)*K03h(m_indx(i+side_tr2-1,nb+8),m)*t(3)+w(4)*K04h(m_indx(i+side_tr2-1,nb+8),m)*t(4)+w(5)*K05h(m_indx(i+side_tr2-1,nb+8),m)*t(5)+w(6)*K06h(m_indx(i+side_tr2-1,nb+8),m)*t(6))
      
    
        
     
        I1(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0/(2*PI) + (1+e)**2/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h/(2*PI) -1.0/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h_prev/(2*PI))
        I2(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0/(2*PI) + (1+e)**2/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h/(2*PI) -1.0/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h_prev/(2*PI))
        
        
        I1_k0h=h_array(m)*0.5*(w(1)*K01h(i+side_tr2+1,m)+w(2)*K02h(i+side_tr2+1,m)+w(3)*K03h(i+side_tr2+1,m)+w(4)*K04h(i+side_tr2+1,m)+w(5)*K05h(i+side_tr2+1,m)+w(6)*K06h(i+side_tr2+1,m))
        I2_k0h=0.5*(w(1)*K01h(i+side_tr2+1,m)*t(1)+w(2)*K02h(i+side_tr2+1,m)*t(2)+w(3)*K03h(i+side_tr2+1,m)*t(3)+w(4)*K04h(i+side_tr2+1,m)*t(4)+w(5)*K05h(i+side_tr2+1,m)*t(5)+w(6)*K06h(i+side_tr2+1,m)*t(6))
        I1_k0=h_array(m)*0.5*(w(1)*K01(m_indx(i,nb),m)+w(2)*K02(m_indx(i,nb),m)+w(3)*K03(m_indx(i,nb),m)+w(4)*K04(m_indx(i,nb),m)+w(5)*K05(m_indx(i,nb),m)+w(6)*K06(m_indx(i,nb),m))
        I2_k0=0.5*(w(1)*K01(m_indx(i,nb),m)*t(1)+w(2)*K02(m_indx(i,nb),m)*t(2)+w(3)*K03(m_indx(i,nb),m)*t(3)+w(4)*K04(m_indx(i,nb),m)*t(4)+w(5)*K05(m_indx(i,nb),m)*t(5)+w(6)*K06(i+side_tr,m)*t(6))
        I1_k0h_prev=h_array(m)*0.5*(w(1)*K01h(m_indx(i+side_tr2+2,nb+8),m)+w(2)*K02h(m_indx(i+side_tr2+2,nb+8),m)+w(3)*K03h(m_indx(i+side_tr2+2,nb+8),m)+w(4)*K04h(m_indx(i+side_tr2+2,nb+8),m)+w(5)*K05h(m_indx(i+side_tr2+2,nb+8),m)+w(6)*K06h(m_indx(i+side_tr2+2,nb+8),m))
        I2_k0h_prev=0.5*(w(1)*K01h(m_indx(i+side_tr2+2,nb+8),m)*t(1)+w(2)*K02h(m_indx(i+side_tr2+2,nb+8),m)*t(2)+w(3)*K03h(m_indx(i+side_tr2+2,nb+8),m)*t(3)+w(4)*K04h(m_indx(i+side_tr2+2,nb+8),m)*t(4)+w(5)*K05h(m_indx(i+side_tr2+2,nb+8),m)*t(5)+w(6)*K06h(m_indx(i+side_tr2+2,nb+8),m)*t(6))
      
    
        I1(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0/(2*PI) + (1+e)**2/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h/(2*PI) -1.0/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I1_k0h_prev/(2*PI))
        I2(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0/(2*PI) + (1+e)**2/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h/(2*PI) -1.0/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I2_k0h_prev/(2*PI))
      
       end if
       
    
        
        
      
    end do
     !if(i==1) print*,dK02_dt(i+side_tr,:) 
     !I1(i+side_tr,:)=I1(i+side_tr,:)/(beta_const(i,ny,nz)*2*PI)
     !I1(i+side_tr + 1,:)=I1(i+side_tr+1,:)/(beta_const(i,ny,nz)*2*PI)
     !I2(i+side_tr ,:)=I2(i+side_tr,:)/(beta_const(i,ny,nz)*2*PI)
     !I2(i+side_tr + 1,:)=I2(i+side_tr + 1,:)/(beta_const(i,ny,nz)*2*PI)
     
     side_tr=side_tr+1
     side_tr2=side_tr2+2
   else
    a=(0.5*h_array(i-1))**2/(0.5*h_array(i-1)*0.5*h_array(i)*(0.5*h_array(i-1)+0.5*h_array(i)))
    c=-(0.5*h_array(i))**2/(0.5*h_array(i-1)*0.5*h_array(i)*(0.5*h_array(i-1)+0.5*h_array(i)))
    b=((0.5*h_array(i))**2 - (0.5*h_array(i-1))**2)/(0.5*h_array(i-1)*0.5*h_array(i)*(0.5*h_array(i-1)+0.5*h_array(i)))
    do m=1,nb
      
      t(1:)=h_array(m)*0.5*(x(1:) + 1.0d0)
      
      I1_k0h=h_array(m)*0.5*(w(1)*K01h(i+side_tr2,m)+w(2)*K02h(i+side_tr2,m)+w(3)*K03h(i+side_tr2,m)+w(4)*K04h(i+side_tr2,m)+w(5)*K05h(i+side_tr2,m)+w(6)*K06h(i+side_tr2,m))
      I2_k0h=0.5*(w(1)*K01h(i+side_tr2,m)*t(1)+w(2)*K02h(i+side_tr2,m)*t(2)+w(3)*K03h(i+side_tr2,m)*t(3)+w(4)*K04h(i+side_tr2,m)*t(4)+w(5)*K05h(i+side_tr2,m)*t(5)+w(6)*K06h(i+side_tr2,m)*t(6))
      I1_k0=h_array(m)*0.5*(w(1)*K01(m_indx(i,nb),m)+w(2)*K02(m_indx(i,nb),m)+w(3)*K03(m_indx(i,nb),m)+w(4)*K04(m_indx(i,nb),m)+w(5)*K05(m_indx(i,nb),m)+w(6)*K06(m_indx(i,nb),m))
      I2_k0=0.5*(w(1)*K01(m_indx(i,nb),m)*t(1)+w(2)*K02(m_indx(i,nb),m)*t(2)+w(3)*K03(m_indx(i,nb),m)*t(3)+w(4)*K04(m_indx(i,nb),m)*t(4)+w(5)*K05(m_indx(i,nb),m)*t(5)+w(6)*K06(i,m)*t(6))
      I1_k0h_prev=h_array(m)*0.5*(w(1)*K01h(m_indx(i+side_tr2-1,nb+8),m)+w(2)*K02h(m_indx(i+side_tr2-1,nb+8),m)+w(3)*K03h(m_indx(i+side_tr2-1,nb+8),m)+w(4)*K04h(m_indx(i+side_tr2-1,nb+8),m)+w(5)*K05h(m_indx(i+side_tr2-1,nb+8),m)+w(6)*K06h(m_indx(i+side_tr2-1,nb+8),m))
      I2_k0h_prev=0.5*(w(1)*K01h(m_indx(i+side_tr2-1,nb+8),m)*t(1)+w(2)*K02h(m_indx(i+side_tr2-1,nb+8),m)*t(2)+w(3)*K03h(m_indx(i+side_tr2-1,nb+8),m)*t(3)+w(4)*K04h(m_indx(i+side_tr2-1,nb+8),m)*t(4)+w(5)*K05h(m_indx(i+side_tr2-1,nb+8),m)*t(5)+w(6)*K06h(m_indx(i+side_tr2-1,nb+8),m)*t(6))
      

      
       if(m==m_indx(i,nb))then
        call integrals2h(alpha,Ia_1,Ia_2,h_array(m))
        I1_k0=Ia_1
        I2_k0=Ia_2
        call integrals2h(alpha,Ia_1hlft,Ia_2hlft,0.5*h_array(m_indx(m,nb)))
        call integrals2h(alpha,Ia_1hrt,Ia_2hrt,0.5*h_array(m_indx(m,nb)))
        I1_k0h=Ia_1hlft+Ia_1hrt
        I2_k0h=Ia_2hlft+(0.5*h_array(m_indx(m,nb))*Ia_1hrt-Ia_2hrt)
        call integrals2h(alpha,Ia_1hlft,Ia_2hlft,h_array(m)+0.5*h_array(m_indx(m-1,nb)))
        call integrals2h(alpha,Ia_1hrt,Ia_2hrt,0.5*h_array(m_indx(m-1,nb)))
        I1_k0h_prev=Ia_1hlft-Ia_1hrt
        I2_k0h_prev=Ia_2hlft-Ia_2hrt
        
        I1(i+side_tr,m)=a*I1_k0h/(2*PI) + b*I1_k0/(2*PI) + c*I1_k0h_prev/(2*PI)
        I2(i+side_tr,m)=a*I2_k0h/(2*PI) + b*I2_k0/(2*PI) + c*I2_k0h_prev/(2*PI)
       elseif(m==m_indx(i-1,nb))then
        call integrals2h(alpha,Ia_1,Ia_2,h_array(m))
        I1_k0=Ia_1
        I2_k0=(h_array(m)*Ia_1-Ia_2)
        call integrals2h(alpha,Ia_1hlft,Ia_2hlft,0.5*h_array(m_indx(m+1,nb))+h_array(m))
        call integrals2h(alpha,Ia_1hrt,Ia_2hrt,0.5*h_array(m_indx(m+1,nb)))
        I1_k0h=Ia_1hlft-Ia_1hrt
        I2_k0h=(Ia_1hlft*(0.5*h_array(m_indx(m+1,nb))+h_array(m_indx(m,nb)))-Ia_1hrt*(0.5*h_array(m_indx(m+1,nb))))-(Ia_2hlft-Ia_2hrt)
        call integrals2h(alpha,Ia_1hlft,Ia_2hlft,0.5*h_array(m_indx(m,nb)))
        call integrals2h(alpha,Ia_1hrt,Ia_2hrt,0.5*h_array(m_indx(m,nb)))
        I1_k0h_prev=Ia_1hlft+Ia_1hrt
        I2_k0h_prev=Ia_2hlft+(0.5*h_array(m_indx(m,nb))*Ia_1hrt-Ia_2hrt)
        
        I1(i+side_tr,m)=a*I1_k0h/(2*PI) + b*I1_k0/(2*PI) + c*I1_k0h_prev/(2*PI)
        I2(i+side_tr,m)=a*I2_k0h/(2*PI) + b*I2_k0/(2*PI) + c*I2_k0h_prev/(2*PI)
       else
       
     
        I1(i+side_tr,m)=a*I1_k0h/(2*PI) + b*I1_k0/(2*PI) + c*I1_k0h_prev/(2*PI)
        I2(i+side_tr,m)=a*I2_k0h/(2*PI) + b*I2_k0/(2*PI) + c*I2_k0h_prev/(2*PI)
       end if
       
      
      
      
    end do
     
   end if
  end do
  
  print*,'Integral of dK0dt and dK0dtl'
  
  
end subroutine integral_dK0dtau     

!*************************************************************************
! To compute the Integral I1=dK0nor_dl and I2=dK0nor_dl in each element
!*************************************************************************
subroutine integral_dK0nor(ny,nz,nb,gridpts_y,gridpts_z,alpha,h_array,rad1,rad2,rad3,rad4,rad5,rad6,I1,I2)
  integer , intent(in)                          ::ny,nz,nb
  real*8  , intent(in)                          ::alpha
  real*8  , intent(in)                          ::h_array(1:nb)
  real*8  ,allocatable, intent(in)              ::rad1(:,:),rad2(:,:),rad3(:,:),rad4(:,:),rad5(:,:),rad6(:,:)
  real*8  , intent(out)                         ::I1(1:nb,1:nb),I2(1:nb,1:nb)
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  ,allocatable                          ::dK01r(:,:),dK02r(:,:),dK03r(:,:),dK04r(:,:),dK05r(:,:),dK06r(:,:)
  real*8  ,allocatable                          ::dK01(:,:),dK02(:,:),dK03(:,:),dK04(:,:),dK05(:,:),dK06(:,:)
  integer*4                                     ::i,m 
  real,parameter                                ::PI=4*atan(1.d0)
  real*8                                        ::y_disp1,z_disp1
  real*8                                        ::y_disp2,z_disp2
  real*8                                        ::y_disp3,z_disp3
  real*8                                        ::y_disp4,z_disp4
  real*8                                        ::y_disp5,z_disp5
  real*8                                        ::y_disp6,z_disp6
  real*8,allocatable                            ::x(:),t(:),w(:)         ! x =[-1,1] and t=[r(m),r(m+1)]
  

  call dk0rad(nb,alpha,rad1,dK01r)
  call dk0nor(ny,nz,nb,rad1,gridpts_y,gridpts_z,dK01r,dK01)
  call dk0rad(nb,alpha,rad2,dK02r)
  call dk0nor(ny,nz,nb,rad2,gridpts_y,gridpts_z,dK02r,dK02)
  call dk0rad(nb,alpha,rad3,dK03r)
  call dk0nor(ny,nz,nb,rad3,gridpts_y,gridpts_z,dK03r,dK03)
  call dk0rad(nb,alpha,rad4,dK04r)
  call dk0nor(ny,nz,nb,rad4,gridpts_y,gridpts_z,dK04r,dK04)
  call dk0rad(nb,alpha,rad5,dK05r)
  call dk0nor(ny,nz,nb,rad5,gridpts_y,gridpts_z,dK05r,dK05)
  call dk0rad(nb,alpha,rad6,dK06r)
  call dk0nor(ny,nz,nb,rad6,gridpts_y,gridpts_z,dK06r,dK06)
  
  
  call legendre_gauss(6, x, w)
  allocate(t(1:size(x)))
  
      
  do i=1,nb
    do m=1,nb
      if((m<nz+1) .or. (m>=ny+nz+1 .and. m<2*nz+ny+1)) then
       t(1:)=h_array(m)*0.5*(x(1:) + 1.0d0)
     
       I1(i,m)=h_array(m)*0.5*(w(1)*dK01(i,m)+w(2)*dK02(i,m)+w(3)*dK03(i,m)+w(4)*dK04(i,m)+w(5)*dK05(i,m)+w(6)*dK06(i,m))
     
       I2(i,m)=0.5*(w(1)*dK01(i,m)*t(1)+w(2)*dK02(i,m)*t(2)+w(3)*dK03(i,m)*t(3)+w(4)*dK04(i,m)*t(4)+w(5)*dK05(i,m)*t(5)+w(6)*dK06(i,m)*t(6))
      
      else if (( m>=nz+1 .and. m< ny+nz+1) .or. (m>=2*nz+ny+1  .and. m<=nb)) then
       t(1:)=h_array(m)*0.5*(x(1:) + 1.0d0)
   
       
       I1(i,m)=h_array(m)*0.5*(w(1)*dK01(i,m)+w(2)*dK02(i,m)+w(3)*dK03(i,m)+w(4)*dK04(i,m)+w(5)*dK05(i,m)+w(6)*dK06(i,m))
     
       
       I2(i,m)=0.5*(w(1)*dK01(i,m)*t(1)+w(2)*dK02(i,m)*t(2)+w(3)*dK03(i,m)*t(3)+w(4)*dK04(i,m)*t(4)+w(5)*dK05(i,m)*t(5)+w(6)*dK06(i,m)*t(6))
    
      end if
    end do
  end do
  print*,'Integral of dK0nordl and dK0norldl'
end subroutine integral_dK0nor


!*************************************************************************
! To compute the Integral I1=dK0nor_tau and I2=dK0nor_tau_dl in each element
!*************************************************************************
subroutine integral_dK0nor_tau(ny,nz,nb,gridpts_y,gridpts_z,alpha,h_array,rad1,rad2,rad3,rad4,rad5,rad6,rad1h,rad2h,rad3h,rad4h,rad5h,rad6h,I1,I2)
  integer , intent(in)                          ::ny,nz,nb
  real*8  , intent(in)                          ::alpha
  real*8  , intent(in)                          ::h_array(1:nb)
  real*8  ,allocatable, intent(in)              ::rad1(:,:),rad2(:,:),rad3(:,:),rad4(:,:),rad5(:,:),rad6(:,:)
  real*8  , allocatable,intent(in)              ::rad1h(:,:),rad2h(:,:),rad3h(:,:),rad4h(:,:),rad5h(:,:),rad6h(:,:)
  real*8  ,allocatable, intent(out)             ::I1(:,:),I2(:,:)
  real*8  , intent(in)                          ::gridpts_y(0:ny),gridpts_z(0:nz)
  real*8  ,allocatable                          ::dK01r(:,:),dK02r(:,:),dK03r(:,:),dK04r(:,:),dK05r(:,:),dK06r(:,:)
  real*8  ,allocatable                          ::dK01(:,:),dK02(:,:),dK03(:,:),dK04(:,:),dK05(:,:),dK06(:,:)
  real*8  ,allocatable                          ::dK01h(:,:),dK02h(:,:),dK03h(:,:),dK04h(:,:),dK05h(:,:),dK06h(:,:)
  real*8  ,allocatable                          ::dK01_nt(:,:),dK02_nt(:,:),dK03_nt(:,:),dK04_nt(:,:),dK05_nt(:,:),dK06_nt(:,:)
  integer*4                                     ::i,m 
  real,parameter                                ::PI=4*atan(1.d0)
  real*8,allocatable                            ::x(:),t(:),w(:)         ! x =[-1,1] and t=[r(m),r(m+1)]
  real*8                                        ::e,a,b,c                !dummy_variables  
  integer                                       ::n_size
  integer                                       ::side_tr, side_tr2
  real*8                                        ::I1_dk0h,I2_dk0h,I1_dk0,I2_dk0,I1_dk0h_prev,I2_dk0h_prev   
  real*8                                        ::xi,xf,int1,int2,h_arr
  real*8,parameter                              ::cf=0.44,cm=0.22,cn=0.00001
   
  call dk0rad(nb,alpha,rad1,dK01r)
  call dk0nor(ny,nz,nb,rad1,gridpts_y,gridpts_z,dK01r,dK01)
  call dk0rad(nb,alpha,rad2,dK02r)
  call dk0nor(ny,nz,nb,rad2,gridpts_y,gridpts_z,dK02r,dK02)
  call dk0rad(nb,alpha,rad3,dK03r)
  call dk0nor(ny,nz,nb,rad3,gridpts_y,gridpts_z,dK03r,dK03)
  call dk0rad(nb,alpha,rad4,dK04r)
  call dk0nor(ny,nz,nb,rad4,gridpts_y,gridpts_z,dK04r,dK04)
  call dk0rad(nb,alpha,rad5,dK05r)
  call dk0nor(ny,nz,nb,rad5,gridpts_y,gridpts_z,dK05r,dK05)
  call dk0rad(nb,alpha,rad6,dK06r)
  call dk0nor(ny,nz,nb,rad6,gridpts_y,gridpts_z,dK06r,dK06)
  
 
  call dk0radh(nb,alpha,rad1h,dK01r)
  call dk0norh2(ny,nz,nb,rad1h,gridpts_y,gridpts_z,dK01r,dK01h)
  call dk0radh(nb,alpha,rad2h,dK02r)
  call dk0norh2(ny,nz,nb,rad2h,gridpts_y,gridpts_z,dK02r,dK02h)
  call dk0radh(nb,alpha,rad3h,dK03r)
  call dk0norh2(ny,nz,nb,rad3h,gridpts_y,gridpts_z,dK03r,dK03h)
  call dk0radh(nb,alpha,rad4h,dK04r)
  call dk0norh2(ny,nz,nb,rad4h,gridpts_y,gridpts_z,dK04r,dK04h)
  call dk0radh(nb,alpha,rad5h,dK05r)
  call dk0norh2(ny,nz,nb,rad5h,gridpts_y,gridpts_z,dK05r,dK05h)
  call dk0radh(nb,alpha,rad6h,dK06r)
  call dk0norh2(ny,nz,nb,rad6h,gridpts_y,gridpts_z,dK06r,dK06h)
  
  call legendre_gauss(6, x, w)
  allocate(t(1:size(x)))
  n_size=2*(ny+1)+2*(nz+1)
  side_tr=0
  side_tr2=0
  
  
  allocate(I1(1:n_size,1:nb),I2(1:n_size,1:nb))
      
  do i=1,nb
   if(corner(i,ny,nz))then
    do m=1,nb
     e=1.d0
     
     
     t(1:)=h_array(m)*0.5*(x(1:) + 1.0d0)
     
     I1_dk0h=h_array(m)*0.5*(w(1)*dK01h(i+side_tr2,m)+w(2)*dK02h(i+side_tr2,m)+w(3)*dK03h(i+side_tr2,m)+w(4)*dK04h(i+side_tr2,m)+w(5)*dK05h(i+side_tr2,m)+w(6)*dK06h(i+side_tr2,m))
     I2_dk0h=0.5*(w(1)*dK01h(i+side_tr2,m)*t(1)+w(2)*dK02h(i+side_tr2,m)*t(2)+w(3)*dK03h(i+side_tr2,m)*t(3)+w(4)*dK04h(i+side_tr2,m)*t(4)+w(5)*dK05h(i+side_tr2,m)*t(5)+w(6)*dK06h(i+side_tr2,m)*t(6))
     I1_dk0=h_array(m)*0.5*(w(1)*dK01h(m_indx(i+side_tr2+1,nb+16),m)+w(2)*dK02h(m_indx(i+side_tr2+1,nb+16),m)+w(3)*dK03h(m_indx(i+side_tr2+1,nb+16),m)+w(4)*dK04h(m_indx(i+side_tr2+1,nb+16),m)+w(5)*dK05h(m_indx(i+side_tr2+1,nb+16),m)+w(6)*dK06h(m_indx(i+side_tr2+1,nb+16),m))
     I2_dk0=0.5*(w(1)*dK01h(m_indx(i+side_tr2+1,nb+16),m)*t(1)+w(2)*dK02h(m_indx(i+side_tr2+1,nb+16),m)*t(2)+w(3)*dK03h(m_indx(i+side_tr2+1,nb+16),m)*t(3)+w(4)*dK04h(m_indx(i+side_tr2+1,nb+16),m)*t(4)+w(5)*dK05h(m_indx(i+side_tr2+1,nb+16),m)*t(5)+w(6)*dK06h(m_indx(i+side_tr2+1,nb+16),m)*t(6))
     I1_dk0h_prev=h_array(m)*0.5*(w(1)*dK01h(m_indx(i+side_tr2-1,nb+16),m)+w(2)*dK02h(m_indx(i+side_tr2-1,nb+16),m)+w(3)*dK03h(m_indx(i+side_tr2-1,nb+16),m)+w(4)*dK04h(m_indx(i+side_tr2-1,nb+16),m)+w(5)*dK05h(m_indx(i+side_tr2-1,nb+16),m)+w(6)*dK06h(m_indx(i+side_tr2-1,nb+16),m))
     I2_dk0h_prev=0.5*(w(1)*dK01h(m_indx(i+side_tr2-1,nb+16),m)*t(1)+w(2)*dK02h(m_indx(i+side_tr2-1,nb+16),m)*t(2)+w(3)*dK03h(m_indx(i+side_tr2-1,nb+16),m)*t(3)+w(4)*dK04h(m_indx(i+side_tr2-1,nb+16),m)*t(4)+w(5)*dK05h(m_indx(i+side_tr2-1,nb+16),m)*t(5)+w(6)*dK06h(m_indx(i+side_tr2-1,nb+16),m)*t(6))
 
     !I1_dk0h=h_array(m)*0.5*(w(1)*dK01h(i+side_tr2,m)+w(2)*dK02h(i+side_tr2,m)+w(3)*dK03h(i+side_tr2,m)+w(4)*dK04h(i+side_tr2,m)+w(5)*dK05h(i+side_tr2,m)+w(6)*dK06h(i+side_tr2,m))
     !I2_dk0h=0.5*(w(1)*dK01h(i+side_tr2,m)*t(1)+w(2)*dK02h(i+side_tr2,m)*t(2)+w(3)*dK03h(i+side_tr2,m)*t(3)+w(4)*dK04h(i+side_tr2,m)*t(4)+w(5)*dK05h(i+side_tr2,m)*t(5)+w(6)*dK06h(i+side_tr2,m)*t(6))
     !I1_dk0=h_array(m)*0.5*(w(1)*dK01(m_indx(i,nb),m)+w(2)*dK02(m_indx(i,nb),m)+w(3)*dK03(m_indx(i,nb),m)+w(4)*dK04(m_indx(i,nb),m)+w(5)*dK05(m_indx(i,nb),m)+w(6)*dK06(m_indx(i,nb),m))
     !I2_dk0=0.5*(w(1)*dK01(m_indx(i,nb),m)*t(1)+w(2)*dK02(m_indx(i,nb),m)*t(2)+w(3)*dK03(m_indx(i,nb),m)*t(3)+w(4)*dK04(m_indx(i,nb),m)*t(4)+w(5)*dK05(m_indx(i,nb),m)*t(5)+w(6)*dK06(m_indx(i,nb),m)*t(6))
     !I1_dk0h_prev=h_array(m)*0.5*(w(1)*dK01(m_indx(i-1,nb),m)+w(2)*dK02(m_indx(i-1,nb),m)+w(3)*dK03(m_indx(i-1,nb),m)+w(4)*dK04(m_indx(i-1,nb),m)+w(5)*dK05(m_indx(i-1,nb),m)+w(6)*dK06(m_indx(i-1,nb),m))
     !I2_dk0h_prev=0.5*(w(1)*dK01(m_indx(i-1,nb),m)*t(1)+w(2)*dK02(m_indx(i-1,nb),m)*t(2)+w(3)*dK03(m_indx(i-1,nb),m)*t(3)+w(4)*dK04(m_indx(i-1,nb),m)*t(4)+w(5)*dK05(m_indx(i-1,nb),m)*t(5)+w(6)*dK06(m_indx(i-1,nb),m)*t(6))
 
 
     I1(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(0.25*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0/(2*PI) + (1+e)**2/(0.25*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0h/(2*PI) -1.0/(0.25*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0h_prev/(2*PI))
     I2(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(0.25*h_array(m_indx(i-1,nb)))/(e+e**2)*I2_dk0/(2*PI) + (1+e)**2/(0.25*h_array(m_indx(i-1,nb)))/(e+e**2)*I2_dk0h/(2*PI) -1.0/(0.25*h_array(m_indx(i-1,nb)))/(e+e**2)*I2_dk0h_prev/(2*PI))
     
     !I1(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(0.5*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0/(2*PI) + (1+e)**2/(0.5*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0h/(2*PI) -1.0/(0.5*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0h_prev/(2*PI))
     !I2(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(0.5*h_array(m_indx(i-1,nb)))/(e+e**2)*I2_dk0/(2*PI) + (1+e)**2/(0.5*h_array(m_indx(i-1,nb)))/(e+e**2)*I2_dk0h/(2*PI) -1.0/(0.5*h_array(m_indx(i-1,nb)))/(e+e**2)*I2_dk0h_prev/(2*PI))
 
     if(m==1) print*,"m=1",I1_dk0h,I1_dk0,I1_dk0h_prev,I1(i+side_tr,m)    
     e=1.0d0
     
        
         
     !I1_dk0h=h_array(m)*0.5*(w(1)*dK01h(i+side_tr2+3,m)+w(2)*dK02h(i+side_tr2+3,m)+w(3)*dK03h(i+side_tr2+3,m)+w(4)*dK04h(i+side_tr2+3,m)+w(5)*dK05h(i+side_tr2+3,m)+w(6)*dK06h(i+side_tr2+3,m))
     !I2_dk0h=0.5*(w(1)*dK01h(i+side_tr2+3,m)*t(1)+w(2)*dK02h(i+side_tr2+3,m)*t(2)+w(3)*dK03h(i+side_tr2+3,m)*t(3)+w(4)*dK04h(i+side_tr2+3,m)*t(4)+w(5)*dK05h(i+side_tr2+3,m)*t(5)+w(6)*dK06h(i+side_tr2+3,m)*t(6))
     !I1_dk0=h_array(m)*0.5*(w(1)*dK01(m_indx(i,nb),m)+w(2)*dK02(m_indx(i,nb),m)+w(3)*dK03(m_indx(i,nb),m)+w(4)*dK04(m_indx(i,nb),m)+w(5)*dK05(m_indx(i,nb),m)+w(6)*dK06(m_indx(i,nb),m))
     !I2_dk0=0.5*(w(1)*dK01(m_indx(i,nb),m)*t(1)+w(2)*dK02(m_indx(i,nb),m)*t(2)+w(3)*dK03(m_indx(i,nb),m)*t(3)+w(4)*dK04(m_indx(i,nb),m)*t(4)+w(5)*dK05(m_indx(i,nb),m)*t(5)+w(6)*dK06(m_indx(i,nb),m)*t(6))
     !I1_dk0h_prev=h_array(m)*0.5*(w(1)*dK01(m_indx(i+1,nb),m)+w(2)*dK02(m_indx(i+1,nb),m)+w(3)*dK03(m_indx(i+1,nb),m)+w(4)*dK04(m_indx(i+1,nb),m)+w(5)*dK05(m_indx(i+1,nb),m)+w(6)*dK06(m_indx(i+1,nb),m))
     !I2_dk0h_prev=0.5*(w(1)*dK01(m_indx(i+1,nb),m)*t(1)+w(2)*dK02(m_indx(i+1,nb),m)*t(2)+w(3)*dK03(m_indx(i+1,nb),m)*t(3)+w(4)*dK04(m_indx(i+1,nb),m)*t(4)+w(5)*dK05(m_indx(i+1,nb),m)*t(5)+w(6)*dK06(m_indx(i+1,nb),m)*t(6))
     
     I1_dk0h=h_array(m)*0.5*(w(1)*dK01h(i+side_tr2+3,m)+w(2)*dK02h(i+side_tr2+3,m)+w(3)*dK03h(i+side_tr2+3,m)+w(4)*dK04h(i+side_tr2+3,m)+w(5)*dK05h(i+side_tr2+3,m)+w(6)*dK06h(i+side_tr2+3,m))
     I2_dk0h=0.5*(w(1)*dK01h(i+side_tr2+3,m)*t(1)+w(2)*dK02h(i+side_tr2+3,m)*t(2)+w(3)*dK03h(i+side_tr2+3,m)*t(3)+w(4)*dK04h(i+side_tr2+3,m)*t(4)+w(5)*dK05h(i+side_tr2+3,m)*t(5)+w(6)*dK06h(i+side_tr2+3,m)*t(6))
     I1_dk0=h_array(m)*0.5*(w(1)*dK01h(m_indx(i+side_tr2+2,nb+16),m)+w(2)*dK02h(m_indx(i+side_tr2+2,nb+16),m)+w(3)*dK03h(m_indx(i+side_tr2+2,nb+16),m)+w(4)*dK04h(m_indx(i+side_tr2+2,nb+16),m)+w(5)*dK05h(m_indx(i+side_tr2+2,nb+16),m)+w(6)*dK06h(m_indx(i+side_tr2+2,nb+16),m))
     I2_dk0=0.5*(w(1)*dK01h(m_indx(i+side_tr2+2,nb+16),m)*t(1)+w(2)*dK02h(m_indx(i+side_tr2+2,nb+16),m)*t(2)+w(3)*dK03h(m_indx(i+side_tr2+2,nb+16),m)*t(3)+w(4)*dK04h(m_indx(i+side_tr2+2,nb+16),m)*t(4)+w(5)*dK05h(m_indx(i+side_tr2+2,nb+16),m)*t(5)+w(6)*dK06h(m_indx(i+side_tr2+2,nb+16),m)*t(6))
     I1_dk0h_prev=h_array(m)*0.5*(w(1)*dK01h(m_indx(i+side_tr2+4,nb+16),m)+w(2)*dK02h(m_indx(i+side_tr2+4,nb+16),m)+w(3)*dK03h(m_indx(i+side_tr2+4,nb+16),m)+w(4)*dK04h(m_indx(i+side_tr2+4,nb+16),m)+w(5)*dK05h(m_indx(i+side_tr2+4,nb+16),m)+w(6)*dK06h(m_indx(i+side_tr2+4,nb+16),m))
     I2_dk0h_prev=0.5*(w(1)*dK01h(m_indx(i+side_tr2+4,nb+16),m)*t(1)+w(2)*dK02h(m_indx(i+side_tr2+4,nb+16),m)*t(2)+w(3)*dK03h(m_indx(i+side_tr2+4,nb+16),m)*t(3)+w(4)*dK04h(m_indx(i+side_tr2+4,nb+16),m)*t(4)+w(5)*dK05h(m_indx(i+side_tr2+4,nb+16),m)*t(5)+w(6)*dK06h(m_indx(i+side_tr2+4,nb+16),m)*t(6))
     
     
     I1(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I1_dk0/(2*PI) + (1+e)**2/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I1_dk0h/(2*PI) -1.0/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I1_dk0h_prev/(2*PI))
     I2(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I2_dk0/(2*PI) + (1+e)**2/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I2_dk0h/(2*PI) -1.0/(0.25*h_array(m_indx(i,nb)))/(e+e**2)*I2_dk0h_prev/(2*PI))
     
     !I1(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(0.5*h_array(m_indx(i,nb)))/(e+e**2)*I1_dk0/(2*PI) + (1+e)**2/(0.5*h_array(m_indx(i,nb)))/(e+e**2)*I1_dk0h/(2*PI) -1.0/(0.5*h_array(m_indx(i,nb)))/(e+e**2)*I1_dk0h_prev/(2*PI))
     !I2(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(0.5*h_array(m_indx(i,nb)))/(e+e**2)*I2_dk0/(2*PI) + (1+e)**2/(0.5*h_array(m_indx(i,nb)))/(e+e**2)*I2_dk0h/(2*PI) -1.0/(0.5*h_array(m_indx(i,nb)))/(e+e**2)*I2_dk0h_prev/(2*PI))
     if(m==1023) print*,"I2_dk0h=",i,I2_dk0h 
     if(m==1023) print*,"I2_dk0h_prev=",i,I2_dk0h_prev  
     if(m==1023) print*,"I2_dk0=",i,I2_dk0 
     if(m==1023) print*,"I2(i+side_tr+1,m)=",i,I2(i+side_tr+1,m) 
  
    
    
     
       
      if(m==m_indx(i,nb))then
       xi=sqrt((cm*h_array(m_indx(i-1,nb)))**2 + 0*h_array(m_indx(i,nb))**2)
       xf=sqrt((cm*h_array(m_indx(i-1,nb)))**2 + h_array(m_indx(i,nb))**2)
       call integrals4(alpha,xi,xf,cm*h_array(m_indx(i-1,nb)),int1,int2)
       I1_dk0h=int1
       I2_dk0h=int2
       
       xi=sqrt((cf*h_array(m_indx(i-1,nb)))**2 + 0*h_array(m_indx(i,nb))**2)
       xf=sqrt((cf*h_array(m_indx(i-1,nb)))**2 + h_array(m_indx(i,nb))**2)
       call integrals4(alpha,xi,xf,cf*h_array(m_indx(i-1,nb)),int1,int2)
       I1_dk0h_prev=int1
       I2_dk0h_prev=int2
       
       
       xi=sqrt((cn*h_array(m_indx(i-1,nb)))**2 + 0*h_array(m_indx(i,nb))**2)
       xf=sqrt((cn*h_array(m_indx(i-1,nb)))**2 + h_array(m_indx(i,nb))**2)
       call integrals4(alpha,xi,xf,cn*h_array(m_indx(i-1,nb)),int1,int2)
       I1_dk0=int1
       I2_dk0=int2
     
       e=(cf-cm)/(cm-cn)
       
       
       I1(i+side_tr,m)=-1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0h/(2*PI) -1.0/(cm*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0h_prev/(2*PI))
       !I1(i+side_tr,m)=-1.d0/h_array(m_indx(i,nb))/(2*PI)
       I2(i+side_tr,m)=-1.d0/h_array(m)*(-(2*e+e**2)/(cm*h_array(m_indx(i-1,nb)))/(e+e**2)*I2_dk0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i-1,nb)))/(e+e**2)*I2_dk0h/(2*PI) -1.0/(cm*h_array(m_indx(i-1,nb)))/(e+e**2)*I2_dk0h_prev/(2*PI))
       !I2(i+side_tr,m)=-1.d0/h_array(m)/(2*PI)
       if(m==1) print*,"I1_dk0h=",i,I1_dk0h 
       if(m==1) print*,"I1_dk0h_prev=",i,I1_dk0h_prev  
       if(m==1) print*,"I1_dk0=",i,I1_dk0 
       if(m==1) print*,"I1(i+side_tr,m)=",i,I1(i+side_tr,m) 
       
       if(m==1) print*,"I1(i+side_tr,m)=",i,I1(i+side_tr,m),-(2*e+e**2)/(cm*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0/(2*PI)&
                                            ,(1+e)**2/(cm*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0h/(2*PI)&
                                            ,-1.0/(cm*h_array(m_indx(i-1,nb)))/(e+e**2)*I1_dk0h_prev/(2*PI) 
      else if(m==m_indx(i-1,nb) )then
       !I2(i+side_tr,m)=-1/2/PI*(PI/4+1/1.4142)/h_array(m)
       
       xi=sqrt((cm*h_array(m_indx(i,nb)))**2 + 0*h_array(m_indx(i-1,nb))**2)
       xf=sqrt((cm*h_array(m_indx(i,nb)))**2 + h_array(m_indx(i-1,nb))**2)
       call integrals4(alpha,xi,xf,cm*h_array(m_indx(i,nb)),int1,int2)
       I1_dk0h=int1
       I2_dk0h=h_array(m_indx(i-1,nb))*int1-int2
       
       xi=sqrt((cf*h_array(m_indx(i,nb)))**2 + 0*h_array(m_indx(i-1,nb))**2)
       xf=sqrt((cf*h_array(m_indx(i,nb)))**2 + h_array(m_indx(i-1,nb))**2)
       call integrals4(alpha,xi,xf,cf*h_array(m_indx(i,nb)),int1,int2)
       I1_dk0h_prev=int1
       I2_dk0h_prev=h_array(m_indx(i-1,nb))*int1-int2 
       
       xi=sqrt((cn*h_array(m_indx(i,nb)))**2 + 0*h_array(m_indx(i-1,nb))**2)
       xf=sqrt((cn*h_array(m_indx(i,nb)))**2 + h_array(m_indx(i-1,nb))**2)
       call integrals4(alpha,xi,xf,cn*h_array(m_indx(i,nb)),int1,int2)
       I1_dk0=int1
       I2_dk0=h_array(m_indx(i-1,nb))*int1-int2
       
       e=(cf-cm)/(cm-cn)
       I1(i+side_tr+1,m)=1.d0*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_dk0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_dk0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I1_dk0h_prev/(2*PI))
       
                                  
       
       !I1(i+side_tr+1,m)=0.d0
       I2(i+side_tr+1,m)=1.d0/h_array(m)*(-(2*e+e**2)/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_dk0/(2*PI) + (1+e)**2/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_dk0h/(2*PI) -1.0/(cm*h_array(m_indx(i,nb)))/(e+e**2)*I2_dk0h_prev/(2*PI))
       !I2(i+side_tr+1,m)=0.d0
       if(m==1024) print*,"I2_dk0h=",i,I1_dk0h 
       if(m==1024) print*,"I1_dk0h_prev=",i,I1_dk0h_prev  
       if(m==1024) print*,"I1_dk0=",i,I1_dk0 
       if(m==1024) print*,"I1(i+side_tr+1,m)=",i,I1(i+side_tr+1,m) 
      end if
       
       
 
    end do
     !I1(i+side_tr,:)=I1(i+side_tr,:)/(beta_const(i,ny,nz)*2*PI)
     !I1(i+side_tr + 1,:)=I1(i+side_tr+1,:)/(beta_const(i,ny,nz)*2*PI)
     !I2(i+side_tr ,:)=I2(i+side_tr,:)/(beta_const(i,ny,nz)*2*PI)
     !I2(i+side_tr + 1,:)=I2(i+side_tr + 1,:)/(beta_const(i,ny,nz)*2*PI)
     
     side_tr=side_tr+1
     side_tr2=side_tr2+4
     
 
     
   else 
   
    a=(0.5*h_array(i-1))**2/(0.5*h_array(i-1)*0.5*h_array(i)*(0.5*h_array(i-1)+0.5*h_array(i)))
    c=-(0.5*h_array(i))**2/(0.5*h_array(i-1)*0.5*h_array(i)*(0.5*h_array(i-1)+0.5*h_array(i)))
    b=((0.5*h_array(i))**2 - (0.5*h_array(i-1))**2)/(0.5*h_array(i-1)*0.5*h_array(i)*(0.5*h_array(i-1)+0.5*h_array(i))) 
    
    do m=1,nb
      if(m<ny+nz+1)then
       h_arr=h_array(m)
      else
       h_arr=h_array(m)
      end if
      t(1:)=h_arr*0.5*(x(1:) + 1.0d0)
      I1_dk0h=h_arr*0.5*(w(1)*dK01h(i+side_tr2,m)+w(2)*dK02h(i+side_tr2,m)+w(3)*dK03h(i+side_tr2,m)+w(4)*dK04h(i+side_tr2,m)+w(5)*dK05h(i+side_tr2,m)+w(6)*dK06h(i+side_tr2,m))
      I2_dk0h=0.5*(w(1)*dK01h(i+side_tr2,m)*t(1)+w(2)*dK02h(i+side_tr2,m)*t(2)+w(3)*dK03h(i+side_tr2,m)*t(3)+w(4)*dK04h(i+side_tr2,m)*t(4)+w(5)*dK05h(i+side_tr2,m)*t(5)+w(6)*dK06h(i+side_tr2,m)*t(6))
      I1_dk0=h_arr*0.5*(w(1)*dK01(m_indx(i,nb),m)+w(2)*dK02(m_indx(i,nb),m)+w(3)*dK03(m_indx(i,nb),m)+w(4)*dK04(m_indx(i,nb),m)+w(5)*dK05(m_indx(i,nb),m)+w(6)*dK06(m_indx(i,nb),m))
      I2_dk0=0.5*(w(1)*dK01(m_indx(i,nb),m)*t(1)+w(2)*dK02(m_indx(i,nb),m)*t(2)+w(3)*dK03(m_indx(i,nb),m)*t(3)+w(4)*dK04(m_indx(i,nb),m)*t(4)+w(5)*dK05(m_indx(i,nb),m)*t(5)+w(6)*dK06(m_indx(i,nb),m)*t(6))
      I1_dk0h_prev=h_arr*0.5*(w(1)*dK01h(i+side_tr2-1,m)+w(2)*dK02h(i+side_tr2-1,m)+w(3)*dK03h(i+side_tr2-1,m)+w(4)*dK04h(i+side_tr2-1,m)+w(5)*dK05h(i+side_tr2-1,m)+w(6)*dK06h(i+side_tr2-1,m))
      I2_dk0h_prev=0.5*(w(1)*dK01h(i+side_tr2-1,m)*t(1)+w(2)*dK02h(i+side_tr2-1,m)*t(2)+w(3)*dK03h(i+side_tr2-1,m)*t(3)+w(4)*dK04h(i+side_tr2-1,m)*t(4)+w(5)*dK05h(i+side_tr2-1,m)*t(5)+w(6)*dK06h(i+side_tr2-1,m)*t(6))
      
      I1(i+side_tr,m)=a*I1_dk0h/(2*PI) + b*I1_dk0/(2*PI) + c*I1_dk0h_prev/(2*PI)
      I2(i+side_tr,m)=a*I2_dk0h/(2*PI) + b*I2_dk0/(2*PI) + c*I2_dk0h_prev/(2*PI)     
        
     
       
       
    end do
      !I1(i+side_tr, : )=I1(i+side_tr,:)/(beta_const(i,ny,nz)*2*PI)
      !I2(i+side_tr ,:)=I2(i+side_tr,:)/(beta_const(i,ny,nz)*2*PI)
   end if
  end do
  
  print*,'Integral of dK0nor_taudl and dK0nor_tauldl'
end subroutine integral_dK0nor_tau



    
end module greens_func_subroutines
