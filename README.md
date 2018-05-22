# KdV_equation
Numerical solution of Korteweg-de Vries equation

This programm is an illustration to Zaborsky and Kruskal article.
( https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.15.240#fulltext )
It computes solutions of Korteweg-de Vries equation with a period of 2, so the x-axis range is set to [0;2].
The template is 
          * 
          | 
*----*----*----*----*
          |
          *
The initial condition at t = 0 is cos (pi * x). On the next timestep, the data required for the computational template 
are taken from an equation described in the article.

Zaborsky & Kruskal pointed out that if the time is small (t < 1/pi), the evolution of the system can be described by an implicit 
equation u = cos (pi* (x - u * t)). Since the timestep is typically tiny for that template, one can obtain the explicit 
equation for u on early steps using Taylor series: 
  u = cos(pi*x)*cos(pi*u*t) - sin(pi*x)*sin(pi*u*t)
  u = cos(pi*x) - sin(pi*x)*pi*u*t
  u * (1 + sin(pi*x)*t*pi) = cos(pi*x)
  u = cos(pi*x) / (1 + sin(pi*x)*t*pi)  
Cases with x = 0, 1 and x = N - 2, N - 1 do not require a special template since we are looking for periodical solutions.

Please, be careful while using the programm: it tends to consume quite a lot of memory. It is a good idea not to allocate the memory
for all mesh points at once, but to compute them block by block, dumping each block into a .txt file and using constant memory for 
any required number of points. Maybe this feature will be added later.
