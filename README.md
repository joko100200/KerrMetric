This code is an geodesic integrator on the kerr metric. Simply it takes in initial parameters for an orbiting particle, such as its coordinate position and velocity. It also takes in angular momentum J and mass M of the blackhole.   

The code is in geometric units where c=G=1

Code Organization: <br>
In the propagator folder there are two files main.f95 and paramater.inp. main.f95 holds all the code for running the simulation while paramater.inp hold all input variables the simulation needs to run.
For an explination on variables look at variable description comment in main. <br>

Code Outputs:<br>
Code will output two files rthetaphi.out and xyz.out. rthetaphi.out are the outputted in Boyer-Lindquist coordinates (r,theta,phi). These are the coordinates used in the integration. The other output file, xyz.out, are the sudo cartisan coordinates which are only there for plotting reasons. The conversion happens in main.f95 and follow these equations. <br>
![equation](https://latex.codecogs.com/gif.image?%5Csmall%20%5Cdpi%7B120%7D%5Cbg%7Bwhite%7D%5Cbegin%7Bbmatrix%7Dx%5C%5Cy%5C%5Cz%5C%5C%5Cend%7Bbmatrix%7D=%5Cbegin%7Bbmatrix%7D%5Csqrt%7Br%5E2&plus;a%5E2%7Dcos(%5Cphi)sin(%5Ctheta)%5C%5C%5Csqrt%7Br%5E2&plus;a%5E2%7Dsin(%5Cphi)sin(%5Ctheta)%5C%5Crcos(%5Ctheta)%5C%5C%5Cend%7Bbmatrix%7D)

Refrences:<br>
[Wikipedia for Kerr Metric](https://en.wikipedia.org/wiki/Kerr_metric) Used for cordinate transformation from Boyer-Lindquist to "Cartisian" and metric refrence. <br>
[Müller, T., & Grave, F. (2009). Catalogue of spacetimes. arXiv preprint arXiv:0904.4184.](https://arxiv.org/abs/0904.4184) (page 52). Used for their derivation of the christoffel symbols for the kerr metric. 
