%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Derive the coefficients for a FR-DG scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_sub, L, dmat, alphaL, alphaR, W] = coeff(P,X,dX)
    % Find the Gauss Nodes
    syms x l dl;
    if P == 0
        x0(1)=0;            W(1) = 2;
        X_sub = (((x0+1)/2).*dX+X(1:end-1));
    elseif P==1
        x0(1)=-1/sqrt(3);   W(1) = 1;
        x0(2)=-x0(1);       W(2) = W(1);
        X_sub = [((((x0(1))+1)/2).*dX+X(1:end-1));((((x0(2))+1)/2).*dX+X(1:end-1))];
    elseif P==2
        x0(1)=-sqrt(3/5);   W(1) = 5/9;
        x0(2)=0;            W(2) = 8/9;
        x0(3)=-x0(1);       W(3) = W(1);
        X_sub = [((((x0(1))+1)/2).*dX+X(1:end-1));((1/2).*dX+X(1:end-1));((((x0(3))+1)/2).*dX+X(1:end-1))];
    elseif P == 3
        x0(1)=-sqrt(3/7+2/7*sqrt(6/5));     W(1) = (18-sqrt(30))/36;
        x0(2)=-sqrt(3/7-2/7*sqrt(6/5));     W(2) = (18+sqrt(30))/36;
        x0(3)=-x0(2);                       W(3) = W(2);
        x0(4)=-x0(1);                       W(4) = W(1);
        X_sub = [(((x0(1)+1)/2).*dX+X(1:end-1));((((x0(2))+1)/2).*dX+X(1:end-1));((((x0(3))+1)/2).*dX+X(1:end-1));((((x0(4))+1)/2).*dX+X(1:end-1))];
    elseif P==4
        x0(1)=-1/3*sqrt(5+2*sqrt(10/7));
        x0(2)=-1/3*sqrt(5-2*sqrt(10/7));
        x0(3)=0;
        x0(4)=-x0(2);
        x0(5)=-x0(1);
        % X_sub = [(((x0(1)+1)/2).*dX+X(1:end-1));((((x0(2))+1)/2).*dX+X(1:end-1));((((x0(3))+1)/2).*dX+X(1:end-1));((((x0(4))+1)/2).*dX+X(1:end-1))];
    else
        disp('Polynomial Order Not Implemeted');
        return;
    end
    
    % the Lagrange basis
    for i=1:P+1
        l(i)=1;
        for j=1:P+1
            if i ~= j
                l(i) = l(i)*(x-x0(j))/(x0(i)-x0(j));
            end
        end
    end
    
    % differential matrix
    dl = diff(l);
    dmat = zeros(P+1,P+1);
    for i=1:P+1
        dmat(i,:)=subs(dl,x,x0(i));
    end
    
    % reconstruction coefficients for the values at -1 and 1
    L = zeros(2,P+1);
    for i=1:P+1
        L(1,i)=subs(l(i),x,-1);
        L(2,i)=subs(l(i),x,1);
    end
    
    % Define Legendre polynomial
    syms n
    fun = (x^2-1)^n;
    
    % L_Pm1 = subs(1/(2^n*factorial(n))*diff(fun, x, P-1), n, P-1);
    L_P = subs(1/(2^n*factorial(n))*diff(fun, x, P), n, P);
    L_Pp1 = subs(1/(2^n*factorial(n))*diff(fun, x, P+1), n, P+1);

    % P1 = subs(1/(2^n*factorial(n))*diff(fun, x, 1), n, 1);
    % P2 = subs(1/(2^n*factorial(n))*diff(fun, x, 2), n, 2);
    % P3 = subs(1/(2^n*factorial(n))*diff(fun, x, 3), n, 3);
    % P4 = subs(1/(2^n*factorial(n))*diff(fun, x, 4), n, 4);
    
    rad = (((-1)^(P+1))/2)*(L_Pp1-L_P); % DG Scheme  
    % rad = (((-1)^P)/2)*(L_P-(((P+1)*L_Pm1 + P*L_Pp1)/(2*P+1))); % Huynh Scheme 

    drad = diff(rad);
    alphaL = double(subs(drad, x, x0));
    alphaR = -flip(alphaL);
end