%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Function that reconstructs the left and right face values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ql, Qr, Qbar] = face_reconstruction(Qbar,L,QBC,X,X_sub,islimteron,W, isBCfixed,plotall)

    Ll = L(1,:);
    Lr = L(end,:);
   
    L_temp = 0;
    R_temp = 0;
    for idx = 1:length(Ll)
        L_temp = L_temp+Qbar(idx,:,:)*Lr(idx);
        R_temp = R_temp+Qbar(idx,:,:)*Ll(idx);
    end
    Ql = reshape(L_temp,[],3);
    Qr = reshape(R_temp,[],3);

    % Ql = reshape(Qbar(1,:,:)*Lr(1)+Qbar(2,:,:)*Lr(2)+Qbar(3,:,:)*Lr(3),[],3);
    % Qr = reshape(Qbar(1,:,:)*Ll(1)+Qbar(2,:,:)*Ll(2)+Qbar(3,:,:)*Ll(3),[],3);
        
    if islimteron == true
        [Ql, Qr, Qbar] = squeeze_lim(Qbar,W,Ql,Qr);
        % [Ql, Qr, Qbar] = minmod(Qbar,W,Ql,Qr,X(2)-X(1));
    end
    
    % Apply boundary conditions
    if isBCfixed == true 
        Ql = cat(1,QBC(:,1)',Ql); % Ql(i, Q)
        Qr = cat(1,Qr,QBC(:,2)');
    else
        
    end

    % figure(3); plot(X_sub,Qbar(:,:,3))
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotall == true
    % Plot solution points and reconstruction
    scatter(X_sub,Qbar(:,:,1),'filled','blue')
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
        xi = X(idx); x1 = Qbar(1,idx,1); x2 = Qbar(2,idx,1); x3 = Qbar(3,idx,1); L_sub = subs(L);
        fplot(L_sub,[X(idx),X(idx+1)],'blue')
    end
    end
end