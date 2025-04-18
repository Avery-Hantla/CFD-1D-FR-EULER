%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               MATLAB 1D Euler Code With Flux Reconstruction
%                              Avery Hantla
%                               March 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% Inputs 
Xbounds = [-4, 4];                  % Domain Size
num_points = [21,41,81,161];        % Mesh Size
sigma = 0.2;                        % CFL Number
gamma = 1.4;                        % Gamma
P = 2;                              % Polynomial Order
N = 35000;                          % Number of time steps  
output_step = 500;                  % How many iterations between plots
islimteron = false;                 % Is limiter on?
isplot = false;                     % Plotting during sim?
isBCfixed = true;                   % Not implemented
flow = 1;                           % Flow problem 1/2
plotall = false;                    % Debugging Feature

for udx = 1:length(num_points)
    clear Qbar % Clear variables between mesh simulations

    % Initilize the domain 
    dX = (Xbounds(2) - Xbounds(1))/(num_points(udx)-1);
    X = Xbounds(1):dX:Xbounds(2);
    % dX = ones(1,num_points(udx)-1)*dX;
    
    % Subdivide the cells
    [X_sub, L, Lprime, alphaL, alphaR, W] = coeff(P,X,dX);
    
    % Find analytical solution
    [x_analytical, q_analytical] = ExactNozzle(Xbounds,flow,num_points(udx),X_sub);
    
    % Find the area and its derivative 
    [A,dAdX] = find_area(X,X_sub,P);
    
    % Plot subdivided domain
    % figure
    % xline(X); hold on; scatter(X_sub,0.5*ones(1,length(X_sub))); plot(X_sub(:),A(:)); plot(X_sub(:),dAdX(:)); 
    
    % Guess Initial Conditions
    rho1 = q_analytical(1,end)*ones(P+1,length(X)-1);
    u1 = (q_analytical(2,end)/rho1(1))*ones(P+1,length(X)-1);
    P1 = (q_analytical(3,end)-0.5*rho1(1)*u1(1)^2)*(gamma-1)*ones(P+1,length(X)-1);
    E1 = ((P1./(gamma-1))+ 0.5.*rho1.*u1.^2);
    
    Qbar(:,:,:) = cat(3,rho1,rho1.*u1,(P1./(gamma-1))+0.5*rho1.*u1.^2); % Qbar(j,i,Q)
    
    % Specify the boundary conditions
    QBC = [q_analytical(:,1),q_analytical(:,end)];
    
    max_res = zeros(N,3);
    
    for n = 1:N
        if plotall == true; figure(1); clf; xline(X); hold on; end %scatter(X_sub,0.5*ones(1,length(X_sub)),'red');
        
        [~,u,~,~,c] = flowvariables(Qbar,gamma);
        % Find current iteration time step
        dt = ((X(2)-X(1)).*sigma)./(abs(u(:))+c(:));
        dt = min(dt);
    
        [dQdt, Qbar] = res(Qbar,L,Lprime,QBC,gamma,A,dAdX,dX,alphaL,alphaR,P,X,X_sub,islimteron, W, isBCfixed, plotall);
        Q_star = Qbar+dt.*dQdt;

        [dQdt, ~] = res(Q_star,L,Lprime,QBC,gamma,A,dAdX,dX,alphaL,alphaR,P,X,X_sub,islimteron, W, isBCfixed, false);
        Q_star = 3/4.*Qbar+1/4.*Q_star+1/4.*dt.*dQdt;

        [dQdt, ~] = res(Q_star,L,Lprime,QBC,gamma,A,dAdX,dX,alphaL,alphaR,P,X,X_sub,islimteron, W, isBCfixed, false);
        Qbar = 1/3.*Qbar+2/3.*Q_star+2/3*dt.*dQdt;
        
       
        max_res(n,:) = reshape(max(abs(dQdt),[],1:2),1,3);
        if mod(n,output_step) == 0
            fprintf('n: %d, res: %d, %d, %d \n',n, (max_res(n,1)), (max_res(n,2)), (max_res(n,3)))
            if isplot == true
                if isplot == true
                    figure(2)
                    subplot(3,1,1); plot(X_sub(:),q_analytical(1,:),'-.'); hold on;
                    subplot(3,1,2); plot(X_sub(:),q_analytical(2,:),'-.'); hold on;
                    subplot(3,1,3); plot(X_sub(:),q_analytical(3,:),'-.'); hold on;
                    plotQ(X_sub(:),Qbar,'-','off'); hold off;
                    drawnow 

                    % figure(3) 
                    % [~,u,~,~,c] = flowvariables(Qbar,gamma);
                    % Mach = u./c;
                    % plot(X_sub(:),Mach(:));
                end
            end
        end
    end

    Q_save{:,:,:,udx} = Qbar;
    X_save{:,:,udx} = X_sub;
    res_save{:,udx} = max_res;
    EL1(udx) = mean(reshape(Qbar(:,:,1),1,[])-q_analytical(1,:))/(length(X_sub(:)));
    EL2(udx) = sqrt((sum((reshape(Qbar(:,:,1),1,[])-q_analytical(1,:)).^2))/(length(X_sub(:))));
end
%% Post Process
close all;
% Find order of error
for zdx = 1:length(num_points)-1
    P = (log(EL2(zdx)/EL2(zdx+1)))/log(2);
    fprintf('Order of error between the mesh with %d points and %d points is: %d \n',num_points(zdx),num_points(zdx+1),P)
end

% Plot the conserved variables
subplot(3,1,1); plot(X_save{:,:,end}(:),q_analytical(1,:),'-.'); hold on;
subplot(3,1,2); plot(X_save{:,:,end}(:),q_analytical(2,:),'-.'); hold on;
subplot(3,1,3); plot(X_save{:,:,end}(:),q_analytical(3,:),'-.'); hold on;
for zdx = 1:length(num_points)
    plotQ(X_save{zdx}(:),Q_save{:,:,:,zdx},'-',false)
end

% Add legend
leg{1} = 'Analytical Solution';
for jdx = 1:length(num_points)
    leg{jdx+1} = sprintf('%d Cells',num_points(jdx)-1);
end
legend(leg,'Location','southeast')

[rho,u,~,P,~] = flowvariables2D(q_analytical',gamma);
figure(length(num_points)+3); plot(X_save{:,:,end}(:),rho,'-.'); hold on;
figure(length(num_points)+4); plot(X_save{:,:,end}(:),u,'-.'); hold on;
figure(length(num_points)+5); plot(X_save{:,:,end}(:),P,'-.'); hold on;

for idx = 1:length(num_points)
    figure(idx+1); xlabel('Iterations, n'); ylabel('Residual'); %saveas(gcf,'res.png')
    for jdx = 1:3
        figure(idx+1); semilogy(1:length(res_save{idx}(2:end,jdx)),res_save{idx}(2:end,jdx)); hold on
    end
    legend({'Continuty','Momentum','Energy'}); 
    title(sprintf('%d Cells',num_points(idx)-1))

    % Plot Residuals
    figure(length(num_points)+2); semilogy(1:length(res_save{idx}(2:end,1)),res_save{idx}(2:end,1)); hold on

    % for jdx = 1:3 % Plot Conservared variables on seperate plots
    [rho,u,~,P,~] = flowvariables(Q_save{:,:,:,idx},gamma);
    figure(length(num_points)+3); plot(X_save{:,:,idx}(:),rho(:)); hold on
    figure(length(num_points)+4); plot(X_save{:,:,idx}(:),u(:)); hold on
    figure(length(num_points)+5); plot(X_save{:,:,idx}(:),P(:)); hold on
    % end
    % Plot Mach number
    % [~,u,~,~,c] = flowvariables(Q_save{:,:,:,idx},gamma);
    % Mach = u./c;
    % figure(length(num_points)+2); plot(X_save{idx}(:),Mach(:)); hold on
end
% figure(length(num_points)+2); legend(leg{2:end},'Location','southeast'); xlabel('x'); ylabel('Mach')
figure(length(num_points)+2); legend(leg{2:end},'Location','northeast'); xlabel('Iterations'); ylabel('Density Residual')

figure(length(num_points)+3); legend(leg{1:end},'Location','southeast'); xlabel('x'); ylabel('Density')
figure(length(num_points)+4); legend(leg{1:end},'Location','northeast'); xlabel('x'); ylabel('Velocity')
figure(length(num_points)+5); legend(leg{1:end},'Location','southeast'); xlabel('x'); ylabel('Pressure')