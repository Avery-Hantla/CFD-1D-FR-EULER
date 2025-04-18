%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function releating (A/A*)^2 and local Mach number               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function aas = MachArea(mach, aas2)
gam=1.4;
gamp1=gam+1;
gamm1=gam-1;
aas = 1./mach^2*(2/gamp1*(1+0.5*gamm1*mach^2))^(gamp1/gamm1)-aas2;