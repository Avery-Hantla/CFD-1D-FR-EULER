%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Matlab Code 
%           Exact Solution for 1D Euler equations 
%               (Z.J. Wang, October 2021)
%
% Notes:
%   1.  flow == 1: subsonic
%       flow == 2: transonic with a shock wave after the throat
%   2. Domain size [xmin,xmax];
%   3. Uniform mesh is used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [x,q] = ExactNozzle(Xbounds,flow,np,X_sub)
    % Initialize the flow    
    Initialize(flow);

    x = X_sub(:);
    q=zeros(3,length(x));

    for i=1:length(x)
       q(:,i) = ExactSolu(x(i), flow);
    end
end
% function [x,q] = ExactNozzle(Xbounds,flow,np,X_sub) 
%     % Domain size
%     xmin=Xbounds(1);
%     xmax=Xbounds(2);
% 
%     % Initialize the flow
%     Initialize(flow);
% 
%     % no of points
%     dx = (xmax - xmin)/(np-1);
%     x = xmin:dx:xmax;
% 
%     % Solution array for (rho, rho*u, E_t)
%     q=zeros(3,np);
%     %
%     % exact solution
%     %
%     for i=1:np
%        q(:,i) = ExactSolu(x(i), flow);
%     end
% end