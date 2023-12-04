function v = polyinterp(x,y,u) 
% The first two input arguments, x and y, are vectors of the same length that define 
% the interpolating points. The third input argument, u, is a vector of points where
% the function is to be evaluated. 
n = length(x); 
v = zeros(size(u)); 
for k = 1:n 
    w = ones(size(u)); 
    for j = [1:k-1 k+1:n] 
        w = (u-x(j))./(x(k)-x(j)).*w; 
    end 
    v = v + w*y(k); 
end