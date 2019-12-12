      
      !----------------------------------------------------------------------!
      !      Subroutine to compute the symmetric part of the matrix Adot     !
      !----------------------------------------------------------------------!
    
      !----------------------------------------------------------------------!
      !               File created by Eduardo García-Portugués               !
      !----------------------------------------------------------------------!
    
    
      Subroutine adot(n,inprod,Adot_vec)
      
      ! Arguments
      Double Precision inprod((n*n+n)/2), Adot_vec((n*n-n+2)/2)
      Integer n
      
      ! Local variables
      Double Precision num, den, quo, aux1, aux2, sumr
      Integer i, j, r, ij, ir, ii, jr, jj, rj, rr, auxi, auxj, auxr
      Real, parameter :: Pi = 3.1415926536
     
      ! Calculus of Adot
      
      ! First element of Adot_vec is common the diagonal element
      Adot_vec(1)=Pi*(n+1)
      
      ! Rest of elements are the lower triangle matrix of Adot
      do i=2,n
        do j=1,i-1
    
          ! Sum on the r index
          sumr=0
          do r=1,n
            
            ! From the definition of Aijr0
            if((i==r) .or. (j==r)) then
        
              ! Sum variable
              sumr=sumr+Pi
      
            else
        
              ! Auxiliar variables for the indexes
              auxi=i*(i-1)/2
              auxj=j*(j-1)/2
              auxr=r*(r-1)/2
        
              ! Indexes
              ij=auxi+j
              ii=auxi+i
              jj=auxj+j
              rr=auxr+r
              
              if(i>r) then !We are in the lower triangle of the matrix
                ir=auxi+r
              else !Upper triangle
                ir=auxr+i
              end if
              
              if(r>j) then !We are in the lower triangle of the matrix
                rj=auxr+j
              else !Upper triangle
                rj=auxj+r
              end if
              
              ! Symmetry
              jr=rj
    
              ! Calculus of the quotient
              num=inprod(ij)-inprod(ir)-inprod(rj)+inprod(rr)
              aux1=sqrt(inprod(ii)-2*inprod(ir)+inprod(rr))
              aux2=sqrt(inprod(jj)-2*inprod(jr)+inprod(rr))
              den=aux1*aux2
              quo=num/den
              
              ! Avoid numerical problems on acos
              if(quo< -1) then
                quo=-1
              else if(quo>1) then
                quo=1
              end if
        
            ! Sum
            sumr=sumr+abs(Pi-acos(quo))
        
            end if
  
          end do
      
          ! Enter the element
          Adot_vec(1+((i-1)*(i-2)/2)+j)=sumr
      
        end do
      end do
           
      Return
      End