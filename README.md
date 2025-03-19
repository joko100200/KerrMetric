## Important Notes
This code is an geodesic integrator on the kerr metric. Simply it takes in initial parameters for an orbiting particle, such as its coordinate position and velocity. It also takes in angular momentum J and mass M of the blackhole.<br>

The main chunk of this code is held in the simplification of the geodesic equation by removing all connections which are equal to 0. Using the non-zero and unique connections, which are evaluated in the connections subroutine, to continue to propagate the geodesic equation. <br>

The Geodesic equation that is being propagated in this code if from the perspective of a frame at r->\infinity away. Because of this the standard geodesic equation needs to be rewritten from proper time of the particle to coordinate time of the simulation. The equation using simple chain rule is thus: <br>
![equation](https://latex.codecogs.com/gif.image?%5Csmall%20%5Cdpi%7B120%7D%5Cbg%7Bwhite%7D%5Cfrac%7Bd%5E2x%5Ei%7D%7Bdt%5E2%7D=(-%5CGamma%5Ei_%7B%5Calpha%5Cbeta%7D&plus;%5CGamma%5E0_%7B%5Calpha%5Cbeta%7D%5Cfrac%7Bdx%5Ei%7D%7Bdt%7D)%5Cfrac%7Bdx%5E%5Calpha%7D%7Bdt%7D%5Cfrac%7Bdx%5E%5Cbeta%7D%7Bdt%7D)
<br>
Where the latin index i only goes from 1,2,3
The code is in geometric units where c=G=1

## Code Organization
In the propagator folder there are two files main.f95 and paramater.inp. main.f95 holds all the code for running the simulation while paramater.inp hold all input variables the simulation needs to run.
For an explination on variables look at variable description comment in main. <br>

Code Outputs:<br>
Code will output two files rthetaphi.out and xyz.out. rthetaphi.out are the outputted in Boyer-Lindquist coordinates (r,theta,phi). These are the coordinates used in the integration. The other output file, xyz.out, are the sudo cartisan coordinates which are only there for plotting reasons. The conversion happens in main.f95 and follow these equations. <br>
![equation](https://latex.codecogs.com/gif.image?%5Csmall%20%5Cdpi%7B120%7D%5Cbg%7Bwhite%7D%5Cbegin%7Bbmatrix%7Dx%5C%5Cy%5C%5Cz%5C%5C%5Cend%7Bbmatrix%7D=%5Cbegin%7Bbmatrix%7D%5Csqrt%7Br%5E2&plus;a%5E2%7Dcos(%5Cphi)sin(%5Ctheta)%5C%5C%5Csqrt%7Br%5E2&plus;a%5E2%7Dsin(%5Cphi)sin(%5Ctheta)%5C%5Crcos(%5Ctheta)%5C%5C%5Cend%7Bbmatrix%7D)

## Refrences:
[Wikipedia for Kerr Metric](https://en.wikipedia.org/wiki/Kerr_metric) Used for cordinate transformation from Boyer-Lindquist to "Cartisian" and metric refrence. <br>
[Müller, T., & Grave, F. (2009). Catalogue of spacetimes. arXiv preprint arXiv:0904.4184.](https://arxiv.org/abs/0904.4184) (page 52). Used for their derivation of the christoffel symbols for the kerr metric. 
