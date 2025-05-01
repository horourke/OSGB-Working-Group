        subroutine code3(n,A,C)
        implicit double precision (a-h,p-z)
        integer n
        
        double precision A(n,n)
        
        C=0.0
        do 100 i=1,n-3 
         do 200 j=i+1,n
    
            do 300 k=i+1,n-1
              if (k .NE. j) then
                do 400 l=k+1,n                   
           if ( l .NE. j) then
                 C=C+A(i,j)*A(k,l)                 
                 end if
                                              
400    continue 
            end if               
300    continue
200    continue        
100    continue        
        
       end