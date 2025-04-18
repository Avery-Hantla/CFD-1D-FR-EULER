# CFD 1D Flux Recontruction Euler Solver

This is a 1D Flux Reconstruction Euler code capable of up to the 4th order of accuracy. For solutions with discontinuities, the squeeze limiter is implemented. Both Roe and Rusanov fluxes are implemented.  

## Inputs

Xbounds - Domain Size  
num_points -  Number of solutions point with equal spacing  
sigma - CFL Number  
gamma - Gamma  
P - Polynomial Order  
N -  Number of time steps    
output_step - How many iterations between plots  
islimteron - Is the squeeze limiter on?  
isplot - Plotting during sim?  true/false  
flow - Flow problem 1/2  
plotall - A Debugging Feature to plot on a small grid

There are two flow problems implemented into the code.

1) Subsonic flow with an exit Mach number of 0.4. Assuming that the density at the exit is 1, and pressure at the exit is 1/gamma and the speed of sound at the exit is 1, and the exit velocity is 0.4.

2) Transonic flow with a shock wave downstream of the throat. Assuming that the inlet Mach number is 0.2006533, inlet density 1, and pressure 1/gamma, and exit pressure 0.6071752.
