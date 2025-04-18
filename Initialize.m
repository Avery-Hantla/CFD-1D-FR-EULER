%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compressive flow through a nozzle   
%    flow == 1: subsonic
%    flow == 2: transonic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = Initialize(flow)

% If transonic, the shock location
global x_shock

% Store the p0 and A* before and after the shock 
global p0 rho0 astar p02 astar2

% Set up some constants such as the ratio of specific heats
gam=1.4;
gamm1 = gam-1.;
gamp1 = gam+1.;

if flow == 2
  % transonic flow
  % Inlet condition and exit pressure
  mach_i = 0.2006533;
  rho_i = 1.;
  p_i = 1./gam;
  
  % Exit pressure
  p_e = 0.6071752;
  
  % Compute the total pressure, density at the inlet
  temp = 1+0.5*gamm1*mach_i^2;
  p0=p_i*temp^(gam/gamm1);
  rho0=rho_i*temp^(1./gamm1);

  % Need to find the location of the shock. Guess first
  x_shock = 0.5;
  
  % Help to converge the shock location faster
  x_shock_small = 0.;
  x_shock_large = 1.;
  
  % Guess an exit press
  pe = 0.;

  astar = area(0);  % Throat area is astar

  while abs(p_e - pe) > 0.00001
    % The area at the shock location   
    A = area(x_shock);

    aas2 = (A/astar)^2;

    % Find the pre-shock Mach no.
    machs = fzero(@(mach)MachArea(mach, aas2), 2.);

    % Find the pre-shock pressure
    temp = 1+0.5*gamm1*machs^2;
    p1 = p0/temp^(gam/gamm1);
     
    % after shock properties
    p2 = p1*(1.+2*gam/gamp1*(machs^2-1.));
    mach2 = sqrt((1+0.5*gamm1*machs^2)/(gam*machs^2-0.5*gamm1));   
    aas2 = 1./mach2^2*(2/gamp1*(1+0.5*gamm1*mach2^2))^(gamp1/gamm1);
    astar2 = sqrt(A^2/aas2);
  
    % After shock total pressure
    temp = 1+0.5*gamm1*mach2^2;
    p02 = p2*temp^(gam/gamm1);
    
    % Find the exit Mach
    aas2 = (area(4)/astar2)^2;
    mache = fzero(@(mach)MachArea(mach, aas2), 0.05);

    % Find the exit pressure
    temp = 1+0.5*gamm1*mache^2;
    pe = p02/temp^(gam/gamm1);
    
    % Revise shock location
    if pe > p_e
        x_shock_small = x_shock;
        x_shock = 0.5*(x_shock + x_shock_large);
    else
        x_shock_large = x_shock;
        x_shock = 0.5*(x_shock + x_shock_small);
    end
  end
  display(x_shock);
elseif flow == 1  
  % subsonic
  % exit Mach, pressure, density and area
  area_e = area(5.0);
  mach_e = 0.4;
  rho_e=1.;
  p_e = 1./gam;

  % Total pressure and density
  temp = 1+0.5*gamm1*mach_e^2;
  p0=p_e*temp^(gam/gamm1);
  rho0=rho_e*temp^(1./gamm1);

  % A*
  temp=1./mach_e^2*(2/gamp1*(1+0.5*gamm1*mach_e^2))^(gamp1/gamm1);
  astar = area_e/sqrt(temp);
end