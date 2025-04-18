%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Function to Calculate Flow Variables from Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho,u,E,P,c] = flowvariables(Q,gamma)
    rho = Q(:,:,1);
    u = Q(:,:,2)./rho;
    E = Q(:,:,3);
    P = (E-0.5.*rho.*(u.^2)).*(gamma-1);
    c = sqrt(gamma.*(P./rho));
end