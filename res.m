%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Function to compute the residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dQdt, Qbar] = res(Qbar_o,L,Lprime,QBC,gamma,A,dAdX,dX,alphaL,alphaR,k,X,X_sub,islimteron, W, isBCfixed, plotall)
    % Reconstruction for face values and solve for riemann flux
    [QL_iphalf, QR_iphalf, Qbar] = face_reconstruction(Qbar_o,L,QBC,X,X_sub,islimteron,W, isBCfixed,plotall); % VALIDATED
    % [F_iphalf,F_imhalf] = riemann(gamma,QL_iphalf,QR_iphalf,X,X_sub,plotall); % VALIDATED
    [F_iphalf,F_imhalf] = roes_reimann(gamma,QL_iphalf,QR_iphalf,X,X_sub,plotall);
    
    % Find discountinous flux 
    [dFd_dEta,fi_m,fi_p,F_temp] = flux(Qbar,L,Lprime,gamma,F_iphalf,F_imhalf,k,X,X_sub,plotall); % VALIDATED
    
    % Compute G Matrix 
    [rho,u,E,P,~] = flowvariables(Qbar,gamma);
    
    Gbar = cat(3,-rho.*u.*(1./A).*dAdX,-rho.*u.^2.*(1./A).*dAdX,-((u.*(E+P))./A).*dAdX);
    
    fiCL = zeros(size(Qbar));
    fiCR = zeros(size(Qbar));
    for idx = 1:3 
        fiCL(:,:,idx) = alphaL'.*fi_m(:,idx)';
        fiCR(:,:,idx) = alphaR'.*fi_p(:,idx)';
    end
    
    % figure(3); plot(X_sub,Gbar(:,:,3),X_sub,Qbar(:,:,3));
    % temp = fiCL+fiCR;
    % figure(4); plot(X_sub,temp(:,:,3));

    % Compute Residual 
    dQdt = Gbar - (2./dX).*(dFd_dEta+fiCL+fiCR);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotall == true
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

    % Define Legendre polynomial
    syms n
    fun = (eta^2-1)^n;

    P_m1 = subs(1/(2^n*factorial(n))*diff(fun, x, 3-1), n, 3-1);
    P_ = subs(1/(2^n*factorial(n))*diff(fun, x, 3), n, 3);

    rad = -1/2*(P_-P_m1);
    radR = subs(rad,x,-(x-(xi-dX)));
    for idx = 1:3 
        fiCL1(:,idx) = rad.*fi_m(:,idx)';
        fiCR1(:,idx) = rad.*fi_p(:,idx)';
    end

    % funx = (x^2-1)^n;
    % 
    % P_m1x = subs(1/(2^n*factorial(n))*diff(funx, x, 3-1), n, 3-1);
    % P_x = subs(1/(2^n*factorial(n))*diff(funx, x, 3), n, 3);
    % 
    % radx = -1/2*(P_x-P_m1x);
    % L1 = double(subs(radx, x, x0));
    % R1 = -flip(L1);
    % 
    % fiCL1 = zeros(size(Qbar));
    % fiCR1 = zeros(size(Qbar));
    % for idx = 1:3 
    %     fiCL1(:,:,idx) = L1'.*fi_m(:,idx)';
    %     fiCR1(:,:,idx) = R1'.*fi_p(:,idx)';
    % end

    for idx = 1:length(X)-1 % number of cells
        % scatter(X_sub(:,idx),fiCL1(:,idx,1),'magenta')
        xi = X(idx); fiCL1_sub = subs(fiCL1(idx,1));
        fplot(fiCL1_sub,[X(idx),X(idx+1)],'magenta','LineStyle',':')

        % scatter(X_sub(:,idx),fiCR1(:,idx,1),'magenta')
        xi = X(idx); fiCR1_sub = subs(fiCR1(idx,1));
        fplot(fiCR1_sub,[X(idx),X(idx+1)],'magenta','LineStyle','-.')

        scatter(X_sub(:,idx),F_temp(:,idx,1),'red')
        xi = X(idx); x1 = F_temp(1,idx,1); x2 = F_temp(2,idx,1); x3 = F_temp(3,idx,1); F_sub = subs(L);
        fplot(F_sub,[X(idx),X(idx+1)],'red','linestyle','-.')

        Flux_total = F_sub + fiCL1 + fiCR1;
        % scatter(X_sub(:,idx),Flux_total(:,idx,1),'filled','magenta')
        xi = X(idx); Flux_total = subs(Flux_total(idx,1));
        fplot(Flux_total,[X(idx),X(idx+1)],'magenta')
    end
    end
end