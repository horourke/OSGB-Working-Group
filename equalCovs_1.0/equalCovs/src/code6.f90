        subroutine code6(n,m,A,C)
        implicit double precision (a-h,p-z)
        integer n,m
        
        double precision A(n,m)
        
        C=0.0
        do 100 i=1,n-1 
         do 200 j=1,m
    
            do 300 k=i+1,n
                do 400 l=1,m   
               if ( l .NE. j) then
                 C=C+A(i,j)*A(k,l)
               end if                                
400    continue              
300    continue
200    continue        
100    continue        
        
       end