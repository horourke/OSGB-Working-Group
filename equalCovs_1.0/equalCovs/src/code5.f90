        subroutine code5(n,m,A,C)
        implicit double precision (a-h,p-z)
        integer n,m
        
        double precision A(n,m)
        
        C=0.0
        do 100 i=1,n 
         do 200 j=1,m-1     
            do 300 k=j+1,m   
                 C=C+A(i,j)*A(i,k)               
300    continue
200    continue        
100    continue        
        
       end