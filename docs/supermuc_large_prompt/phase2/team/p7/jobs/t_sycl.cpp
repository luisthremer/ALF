#include <sycl/sycl.hpp>
#include <cstdio>
int main(){ sycl::queue q{sycl::gpu_selector_v, sycl::property::queue::in_order{}}; int n=1<<20; double*a=sycl::malloc_device<double>(n,q);
 for(int r=0;r<5;r++) q.parallel_for(n,[=](sycl::id<1> i){a[i]=i[0]*0.5;}); q.wait(); printf("ok\n"); return 0;}
