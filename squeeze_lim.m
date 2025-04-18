%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Function to calculate squeeze limiter ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ql, Qr, Qbar] = squeeze_lim(Qbar,W,Ql,Qr)
    % Find Qbar in each cell
    temp = 0;
    for idx = 1:length(W)
        temp = temp+Qbar(idx,:,:)*W(idx);
    end
    Qibar = reshape(temp,[],3)./2;
    
    Qibar_p2 = Qibar(3:end,:);
    Qibar_i = Qibar(2:end-1,:);
    Qibar_m2 = Qibar(1:end-2,:);
    
    % Find max and min Qbar 
    Qi_max = max(max(Qibar_p2,Qibar_i),Qibar_m2);
    Qi_min = min(min(Qibar_p2,Qibar_i),Qibar_m2);
    
    % Find phi left
    phi_m = ones(size(Qibar_i));
    
    Ql_i = Qr(2:end-1,:);
    max_m_idx = Ql_i>Qi_max;
    min_m_idx = Ql_i<Qi_min;
    
    phi_m(max_m_idx) = (Qi_max(max_m_idx)-Qibar_i(max_m_idx))./(Ql_i(max_m_idx)-Qibar_i(max_m_idx));
    phi_m(min_m_idx) = (Qi_min(min_m_idx)-Qibar_i(min_m_idx))./(Ql_i(min_m_idx)-Qibar_i(min_m_idx));
    
    % Find phi right
    phi_p = ones(size(Qibar_i));
    
    Qr_i = Ql(2:end-1,:);
    max_p_idx = Qr_i>Qi_max;
    min_p_idx = Qr_i<Qi_min;
    
    phi_p(max_p_idx) = (Qi_max(max_p_idx)-Qibar_i(max_p_idx))./(Qr_i(max_p_idx)-Qibar_i(max_p_idx));
    phi_p(min_p_idx) = (Qi_min(min_p_idx)-Qibar_i(min_p_idx))./(Qr_i(min_p_idx)-Qibar_i(min_p_idx));
    
    % Find phi 
    phi = min(phi_m,phi_p);
    phi = cat(1,ones(1,3),phi,ones(1,3));

    for idx = 1:length(phi(:,1))
        % min_temp = min(phi(idx,:));
        % phi(idx,:) = [min_temp,min_temp,min_temp];
        if phi(idx,1) == 1 || phi(idx,3)  == 1
            phi(idx,:) = [1,1,1];
        elseif phi(idx,1) == 0 || phi(idx,3)  == 0
            phi(idx,:) = [0,0,0];
        end
       % if idx >= 45 && idx <= 50 
       %      phi(idx,:) = [0,0,0];
       % end
    end
    
    % phi = zeros(size(Qibar)); Qbar_test = Qbar; Ql_test = Ql; Qr_test = Qr;
    % Find the new reconstructed Qbar
    for idx = 1:length(W)
        Qbar(idx,:,:) = Qibar + phi.*(reshape(Qbar(idx,:,:),[],3)-Qibar);
    end
    Ql = Qibar + phi.*(Ql-Qibar);
    Qr = Qibar + phi.*(Qr-Qibar);

end