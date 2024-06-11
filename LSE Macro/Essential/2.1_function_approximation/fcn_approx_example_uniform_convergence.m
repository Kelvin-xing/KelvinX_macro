clear
clc

for i = 1:4    % i loops over the number of nodes
        if i==1
            nnodes = 5;
        end
        if i==2
            nnodes = 11;
        end
        if i==3
            nnodes = 15;
        end
        if i==4
            nnodes = 21;
        end

% generates the Chebyshev nodes at which functions will be fitted        
chebnodes = chebnode(nnodes);  

% calculate function values and matrix with explanatory variables
fvalues = 1./(1+25*chebnodes.^2);
X = cheb_basis(chebnodes,nnodes-1);   % cheb_basis is an external function calculating Cheb. basis functions

% do projection to calculate coefficients of approximating polynomial
%beta = inv(X'*X)*X'*fvalues;
beta = (X'*X)\(X'*fvalues);

% repeat the procedure for equidistant nodes
xnodes = (-1:2/(nnodes-1):1)';
fvalues = 1./(1+25*xnodes.^2);
X = pol_basis(xnodes,nnodes-1);
beta2 = (X'*X)\(X'*fvalues);

% do comparison

T = 501;   % use a much finer grid to evaluate the two approximations
temp = zeros(T,3);
xvalues = (-1: 2/(T-1):1)';
temp(:,1) = 1./(1+25*xvalues.^2);
X = cheb_basis(xvalues,nnodes-1);
temp(:,2) = X*beta;
X = pol_basis(xvalues,nnodes-1);
temp(:,3) = X*beta2;

figure(i)
plot(xvalues,temp(:,1),xvalues,temp(:,2),xvalues,temp(:,3))

    pause
    
end