        subroutine code1(n,a,c)
        implicit double precision (a-h,p-z)
        integer n
        
        double precision a(n,n)
        
        c=0.0
        do 100 i=1,n-1 
         do 200 j=i+1,n
     
            c=c+a(i,j)**2.d0   
          
200    continue        
100    continue        
       return 
       end