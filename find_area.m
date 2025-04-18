%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Function to calculate the area and its derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,dAdX] = find_area(X,X_sub,P)
    syms x 
    % aa1 = 0.536572-0.198086*exp(-log(2.)*(x/0.6)^2);
    % aa1d = diff(aa1);
    % % aa1d2 = diff(aa1d);
    % 
    % aa2 = 1.0 - 0.661514*exp(-log(2.)*(x/0.6)^2);
    % aa2d = diff(aa2);
    % aa2d2 = diff(aa2d);

    aa = 1+(0.34/(1+exp(-2.1*(x-0.4670713487163696))))-(0.8/(1+exp(-3*(x+0.9))));
    aad = diff(aa);
    % aad2 = diff(add);

    A = zeros(P+1,length(X)-1);
    dAdX = zeros(P+1,length(X)-1);
    % d2AdX2 = zeros(P+1,length(X)-1);
    for idx = 1:P+1
        for jdx = 1:length(X)-1
            x = X_sub(idx,jdx);
            A(idx,jdx) = subs(aa);
            dAdX(idx,jdx) = subs(aad);
            % d2AdX2(idx,jdx) = subs(aad2);
            % if X_sub(idx,jdx) > 0 
            %     x = X_sub(idx,jdx);
            %     A(idx,jdx) = subs(aa1);
            %     dAdX(idx,jdx) = subs(aa1d);
            %     % d2AdX2(idx,jdx) = subs(aa1d2);
            % else
            %     x = X_sub(idx,jdx);
            %     A(idx,jdx) = subs(aa2);
            %     dAdX(idx,jdx) = subs(aa2d);
            %     % d2AdX2(idx,jdx) = subs(aa2d2);
            % end
        end
    end
end