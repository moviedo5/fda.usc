      Subroutine rp_stat(proj_X_ord,residuals,n_proj,n,rp_stat_proj)
      implicit none
	  
	  ! Arguments
      Integer n_proj, n
      Integer proj_X_ord(n,n_proj)
      Double Precision residuals(n_proj,n), rp_stat_proj(n_proj,2)

	  ! Local variables
      Integer i
      Double Precision y(n), cvm, ks
  	  
	  ! Interface for function cumsum
      interface
	  
        Function cumsum(x)
        Double Precision, intent(in) :: x(:)
        Double Precision :: cumsum(size(x))
        end function cumsum	  
	  
      end interface
	  
      do i=1,n_proj
	  
        ! Compute the empirical process	  
        y=cumsum(residuals(i,proj_X_ord(:,i)))
        
        ! Statistics (CVM and KS, rows)
        cvm=sum(y*y)/(n*n)
        ks=maxval(abs(y))/sqrt(real(n))

		! Store in rp_stat_proj
        rp_stat_proj(i,1)=cvm
        rp_stat_proj(i,2)=ks
		
      end do
	
      Return
      End
	  
	  ! Cumulative sum
      Function cumsum(x)
      implicit none
	 
    ! Arguments
      Integer j, sx
      Double Precision, intent(in) :: x(:)
      Double Precision :: cumsum(size(x))
	  
	  ! Size of x
      sx=size(x)	  
	  
	  ! Cumulative sum
      cumsum(1)=x(1)
      if(sx>1) then
        do j=2,sx
           cumsum(j)=cumsum(j-1)+x(j)
        end do
      end if
	  
      End Function cumsum

	  	    
