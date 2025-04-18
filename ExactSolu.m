%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute the exact solution at location x               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = ExactSolu(x, flow)

% Some global variables determined in Initialize(flow)
global x_shock
global p0 rho0 astar p02 astar2

gam=1.4;
gamm1 = gam-1.;
gamp1 = gam+1.;

% The conserved variables to be returned
q=zeros(3,1);

if flow == 2
  % transonic sonic

  T0 = p0/rho0;

  if x <= x_shock   
    A = area(x);
    aas2 = (A/astar)^2;
    
    if aas2 < 1.01
      if x >= 0
        machs = fzero(@(mach)MachArea(mach, aas2), 1.1);
      else
        machs = fzero(@(mach)MachArea(mach, aas2), 0.99);
      end
    else
      % Find the local Mach number
      if x > 0
          machs = fzero(@(mach)MachArea(mach, aas2), 2);
      else
          machs = fzero(@(mach)MachArea(mach, aas2), 0.1);
      end
    end
    
    temp = 1+0.5*gamm1*machs^2;
    p1 = p0/temp^(gam/gamm1);
    rho1 = rho0/temp^(1.0/gamm1);
    c = sqrt(gam*p1/rho1);
    u1 = abs(machs)*c;
  else   
    % Find the local Mach
    aas2 = (area(x)/astar2)^2;    
    mache = fzero(@(mach)MachArea(mach, aas2), 0.05);
    
    temp = 1+0.5*gamm1*mache^2;
    p1 = p02/temp^(gam/gamm1);
    T1 = T0/temp;
    rho1 = p1/T1;
    c = sqrt(gam*T1);
    u1 = abs(mache)*c;    
  end

  q(1)=rho1;
  q(2)=rho1*u1;
  q(3)=p1/gamm1+0.5*rho1*u1^2;

elseif flow == 1
  % subsonic
  % find the Mach at x
  temp = (area(x)/astar)^2;
  mach = fzero(@(mach)MachArea(mach,temp), 0.05);

  % Now the density, pressure
  temp = 1+0.5*gamm1*mach^2;
  rho=rho0/temp^(1./gamm1);
  p = p0/temp^(gam/gamm1);
  c = sqrt(gam*p/rho);
  u=abs(mach)*c;
  
  q(1)=rho;
  q(2)=rho*u;
  q(3)=p/gamm1+0.5*rho*u^2;
end