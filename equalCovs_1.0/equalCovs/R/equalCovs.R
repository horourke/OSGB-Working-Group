# R code for testing the equality of two covariance matrices
# sam1:first sample; sam2:second sample
# size1 and size 2 are sample sizes
# sam1 and sam2 must be array with structure size1 x p or size2 x p
# p is the dimension of data

equalCovs<-function(sam1,sam2,size1,size2){

#########################################################
# obtain test statistic given in eqn (2.1) of the paper #
#########################################################

A_mat<- sam1 %*% t(sam1)
out<-0
storage.mode(A_mat)<-"double"
storage.mode(out)<-"double"
nr<-as.integer(size1)

# find A1

result1<-.Fortran(code1,nr,A_mat,out=out)
A1<-2/(size1*(size1-1))*result1[[3]]
              
# find A2

result2<-.Fortran(code2,nr,A_mat, out=out)
A2<-4/(size1*(size1-1)*(size1-2))*result2[[3]]           

# find A3

result3<-.Fortran(code3,nr,A_mat,out=out)
A3<-8/(size1*(size1-1)*(size1-2)*(size1-3))*result3[[3]]

# obtain the statistic given by eqn (2.1) in our paper
A_n1<-A1-A2+A3

# consider the sample 2

B_mat<- sam2 %*% t(sam2)
out<-0
storage.mode(B_mat)<-"double"
storage.mode(out)<-"double"
nrr<-as.integer(size2)

# find B1

result4<-.Fortran(code1,nrr,B_mat,out=out)
B1<-2/(size2*(size2-1))*result4[[3]]
              
# find B2

result5<-.Fortran(code2,nrr,B_mat,out=out)
B2<-4/(size2*(size2-1)*(size2-2))*result5[[3]]          

# find B3

result6<-.Fortran(code3,nrr,B_mat,out=out)
B3<-8/(size2*(size2-1)*(size2-2)*(size2-3))*result6[[3]]

B_n2<-B1-B2+B3

#############################################################
# obtain the test statistic given in eqn (2.2) of the paper #
#############################################################

C_mat1<- sam1 %*% t(sam2)
C_mat2<- sam2 %*% t(sam1)
out<-0
storage.mode(C_mat1)<-"double"
storage.mode(C_mat2)<-"double"
storage.mode(out)<-"double"
nrrr<-as.integer(size1)
nl<-as.integer(size2)

# find C1

result7<-.Fortran(code4,nrrr,nl,C_mat1,out=out)
C1<--2/(size1*size2)*result7[[4]]
                            
# find C2

result8<-.Fortran(code5,nl,nrrr,C_mat2,out=out)
C2<-4/(size1*size2*(size1-1))*result8[[4]]            

# find C3

result9<-.Fortran(code5,nrrr,nl,C_mat1,out=out)
C3<-4/(size1*size2*(size2-1))*result9[[4]]            

# find C4

result10<-.Fortran(code6,nrrr,nl,C_mat1,out=out)
C4<--4/(size1*size2*(size1-1)*(size2-1))*result10[[4]] 

# test statistic given by eqn (2.2) of the paper
C_n<-C1+C2+C3+C4

# the estimator 
T_n<-A_n1+B_n2+C_n

# the standard deviation
Sd_prime<-2*(1/size1+1/size2)*((size1/(size1+size2))*A_n1+(size2/(size1+size2))*B_n2)

test_stat<-T_n/Sd_prime
pvalue<-1-pnorm(test_stat)
test<-c(test_stat,pvalue)	
return(test)	
}