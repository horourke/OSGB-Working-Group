        subroutine code2(n,A,C)
        implicit double precision (a-h,p-z)
        integer n
        
        double precision A(n,n)
        
        C=0.0
        do 100 i=1,n 
         do 200 j=1,(n-1)
          if (j .NE. i) then        
            do 300 k=j+1,n   
               if (k .NE. i) then
                 C=C+A(i,j)*A(i,k)
               end if
                
300    continue
            end if
200    continue        
100    continue        
        
       end