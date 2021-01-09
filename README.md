## CONJUGATE GRADIENT

This code implements CONJUGATE GRADIENT to minimize the potential energy of a system of 
soft-spheres. The code is divided into two parts where the conj_grad.cpp is the conjugate 
gradient part which contains no code specific to the system that is being optimized and 
soft_sphere.cpp	which contains the system specific code. 

Conjugate gradient algorithm is described in detail in 

[Vetterling, William T., et al. Numerical Recipes Example Book (C++): The Art of Scientific Computing. Cambridge University Press, 2002.][1]


The conjugate gradient algorithm implemented in this code is the Polak-Ribiere variant and the line minimization implemented is not 
described in the book.


Author : Varghese Babu ( varghese@jncasr.ac.in )


[1]: http://numerical.recipes/
