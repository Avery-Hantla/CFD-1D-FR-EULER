%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Function to compute the flux distontinous flux, 
%                      flux jump, and correction flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dFd_dEta,fi_m,fi_p,F] = flux(Qbar,L,Lprime,gamma,F_iphalf,F_imhalf,k,X,X_sub,plotall)
    % Solve for flux at solution points
    [rho,u,E,P,~] = flowvariables(Qbar,gamma);

    F = cat(3,rho.*u,(rho.*u.^2)+P,u.*(E+P));

    % Find the discontinous flux derivative
    dFd_dEta = zeros(size(Qbar));
    for idx = 1:k+1
        temp = 0;
        for jdx = 1:k+1
            temp = temp+F(jdx,:,:).*Lprime(idx,jdx);
        end
        dFd_dEta(idx,:,:) = reshape(temp,[],3);
        % dFd_dEta(idx,:,:) = reshape(F(1,:,:).*Lprime(idx,1)+F(2,:,:)*Lprime(idx,2)+F(3,:,:)*Lprime(idx,3),[],3);
    end

    m_temp = 0;
    p_temp = 0;
    for idx = 1:k+1
        m_temp = m_temp+F(idx,:,:).*L(1,idx);
        p_temp = p_temp+F(idx,:,:).*L(end,idx);
    end
    fiD_m = reshape(m_temp,[],3);
    fiD_p = reshape(p_temp,[],3);

    % Compute flux jump 
    % fiD_m = reshape(F(1,:,:).*L(1,1)+F(2,:,:)*L(1,2)+F(3,:,:)*L(1,3),[],3);
    fi_m = F_imhalf - fiD_m;
    
    % fiD_p = reshape(F(1,:,:).*L(end,1)+F(2,:,:)*L(end,2)+F(3,:,:)*L(end,3),[],3);
    fi_p = F_iphalf - fiD_p;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotall == true
    scatter(X(1:end-1)+(X(2)-X(1))*0.025,fiD_m(:,1),'red')
    scatter(X(2:end)-(X(2)-X(1))*0.025,fiD_p(:,1),'red')

    syms x l eta xi x1 x2 x3
    dX = X(2)-X(1);
    eta = ((2*(x-xi))/dX)-1;

    x0(1)=-sqrt(3/5);
    x0(2)=0;
    x0(3)=-x0(1);

    % the Lagrange basis
    for i=1:3
        l(i)=1;
        for j=1:3
            if i ~= j
                l(i) = l(i)*(eta-x0(j))/(x0(i)-x0(j));
            end
        end
    end

    L = x1*l(1)+x2*l(2)+x3*l(3);

    for idx = 1:length(X)-1 % number of cells
        scatter(X_sub(:,idx),F(:,idx,1),'filled','red')
        xi = X(idx); x1 = F(1,idx,1); x2 = F(2,idx,1); x3 = F(3,idx,1); L_sub = subs(L);
        fplot(L_sub,[X(idx),X(idx+1)],'red')

        scatter(X_sub(:,idx),dFd_dEta(:,idx,1),'red')
        xi = X(idx); x1 = dFd_dEta(1,idx,1); x2 = dFd_dEta(2,idx,1); x3 = dFd_dEta(3,idx,1); L_sub = subs(L);
        fplot(L_sub,[X(idx),X(idx+1)],'red','linestyle','-.')
    end
    end
end