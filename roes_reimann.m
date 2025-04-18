%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Function that calulcates the Roe flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F_iphalf,F_imhalf] = roes_reimann(gamma,QL_iphalf,QR_iphalf,X,X_sub,plotall)
    F_temp = zeros(length(X),3);
    for idx = 1:length(X)
        QL = QL_iphalf(idx,:);
        QR = QR_iphalf(idx,:);

        % Find value for left side
        [rhoL,uL,EL,PL,~] = flowvariables2D(QL,gamma);
        HL = (EL+PL);
        FL = [rhoL.*uL,(rhoL.*uL.^2)+PL,uL.*(EL+PL)];
    
        % Find values for right side
        [rhoR,uR,ER,PR,~] = flowvariables2D(QR,gamma);
        HR = (ER+PR);
        FR = [rhoR.*uR,(rhoR.*uR.^2)+PR,uR.*(ER+PR)];

        % Find roes average u and c
        rho_bar = sqrt(rhoL.*rhoR);
        u_bar = (sqrt(rhoL).*uL+sqrt(rhoR).*uR)/(sqrt(rhoL)+sqrt(rhoR));
        H_bar = (sqrt(rhoL).*HL+sqrt(rhoR).*HR)/(sqrt(rhoL)+sqrt(rhoR));
        c_bar = sqrt((gamma-1).*(H_bar-(u_bar.^2)./2));
        
        P_bar = (c_bar.^2)*rho_bar/gamma;
        GAMMA = abs([u_bar, 0, 0;
                 0, u_bar+c_bar, 0;
                 0, 0, u_bar-c_bar]);

        R = [2/u_bar^2, (2*rho_bar*(35*P_bar + 5*rho_bar*u_bar^2 - 2*35^(1/2)*u_bar*(P_bar*rho_bar)^(1/2)))/(245*P_bar^2 + 42*P_bar*rho_bar*u_bar^2 + 5*rho_bar^2*u_bar^4),  (2*rho_bar*(35*P_bar + 5*rho_bar*u_bar^2 + 2*35^(1/2)*u_bar*(P_bar*rho_bar)^(1/2)))/(245*P_bar^2 + 42*P_bar*rho_bar*u_bar^2 + 5*rho_bar^2*u_bar^4);
            2/u_bar,  (10*rho_bar^2*u_bar^3 + 14*7^(1/2)*P_bar*(5*P_bar*rho_bar)^(1/2) + 42*P_bar*rho_bar*u_bar - 2*7^(1/2)*rho_bar*u_bar^2*(5*P_bar*rho_bar)^(1/2))/(245*P_bar^2 + 42*P_bar*rho_bar*u_bar^2 + 5*rho_bar^2*u_bar^4),   (10*rho_bar^2*u_bar^3 - 14*7^(1/2)*P_bar*(5*P_bar*rho_bar)^(1/2) + 42*P_bar*rho_bar*u_bar + 2*7^(1/2)*rho_bar*u_bar^2*(5*P_bar*rho_bar)^(1/2))/(245*P_bar^2 + 42*P_bar*rho_bar*u_bar^2 + 5*rho_bar^2*u_bar^4); 
            1, 1, 1];
        A = R*GAMMA*(R^(-1));
        
        % Solve for F i plus half
        F_temp(idx,:) = (((FL+FR)./2)' - (A./2)*(QR-QL)')';
    end
    F_iphalf = F_temp(2:end,:);
    F_imhalf = F_temp(1:end-1,:);
end