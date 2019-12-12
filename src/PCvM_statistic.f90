      
      !----------------------------------------------------------------------!
      !    Subroutine to compute the Projected Cramér-von Mises statistic    !
      !----------------------------------------------------------------------!
    
      !----------------------------------------------------------------------!
      !               File created by Eduardo García-Portugués               !
      !----------------------------------------------------------------------!

      Subroutine pcvm_statistic(n,Adot_vec,residuals,statistic)
      
      ! Arguments
      Double Precision statistic, Adot_vec((n*n-n+2)/2), residuals(n)
      Integer n
       
      ! Local variables
      Double Precision sums
      Integer i,j
    
      ! Sum for the symmetric part
      sums=0
      do i=2,n
        do j=1,i-1
          sums=sums+residuals(i)*Adot_vec(1+((i-1)*(i-2)/2)+j)*residuals(j)
        end do
      end do
      
      ! Statistic is computed as the sum of the diagonal and the symmetric part
      statistic=Adot_vec(1)*dot_product(residuals,residuals)+2*sums
    
      Return
      End