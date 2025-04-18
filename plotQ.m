%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Function to plot Q matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotQ(X,Q,linetype,ishold)
    % Plot the analytical solution
    subplot(3,1,1);
    plot(X,reshape(Q(:,:,1),1,[]),linetype); if ishold == true; hold on; end; if ishold == 'off'; hold off; end
    xlabel('x');
    ylabel('rho');
    ylim([0.4,1.3]);
    
    subplot(3,1,2);
    plot(X,reshape(Q(:,:,2),1,[]),linetype); if ishold == true; hold on; end; if ishold == 'off'; hold off; end
    xlabel('x');
    ylabel('rho*u');
    ylim([-0.,0.8]);
    
    subplot(3,1,3);
    plot(X,reshape(Q(:,:,3),1,[]),linetype); if ishold == true; hold on; end; if ishold == 'off'; hold off; end
    xlabel('x');
    ylabel('E_t');
    ylim([1,2.2]);
end