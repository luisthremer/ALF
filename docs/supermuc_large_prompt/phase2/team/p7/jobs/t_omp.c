#include <stdio.h>
#include <omp.h>
int main(){ int n=1<<20; double *a=(double*)omp_target_alloc(n*8,0); 
 for(int r=0;r<5;r++){
 #pragma omp target teams distribute parallel for is_device_ptr(a)
 for(int i=0;i<n;i++) a[i]=i*0.5;}
 double s=0;
 #pragma omp target teams distribute parallel for reduction(+:s) is_device_ptr(a)
 for(int i=0;i<n;i++) s+=a[i];
 printf("s=%g\n",s); return 0;}
