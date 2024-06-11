function y = polyn(x,coef,order)

% y = repmat(x,1,order+1).^repmat(0:order,length(x),1)*coef;

if order == 1
    y = coef(1) + coef(2)*x;
elseif order == 2
    y = coef(1) + coef(2)*x + coef(3)*x.^2;
elseif order == 3
    y = coef(1) + coef(2)*x + coef(3)*x.^2 + coef(4)*x.^3;
elseif order == 4
    y = coef(1) + coef(2)*x + coef(3)*x.^2 + coef(4)*x.^3 ...
        + coef(5)*x.^4;
elseif order == 5
    y = coef(1) + coef(2)*x + coef(3)*x.^2 + coef(4)*x.^3 ...
        + coef(5)*x.^4 + coef(6)*x.^5;
elseif order == 6
    y = coef(1) + coef(2)*x + coef(3)*x.^2 + coef(4)*x.^3 ...
        + coef(5)*x.^4 + coef(6)*x.^5 + coef(7)*x.^6;  
elseif order == 7
   y = coef(1) + coef(2)*x + coef(3)*x.^2 + coef(4)*x.^3 ...
        + coef(5)*x.^4 + coef(6)*x.^5 + coef(7)*x.^6 + ...
        coef(8)*x.^7;
elseif order == 8
    y = coef(1) + coef(2)*x + coef(3)*x.^2 + coef(4)*x.^3 ...
        + coef(5)*x.^4 + coef(6)*x.^5 + coef(7)*x.^6 + ...
        coef(8)*x.^7 + coef(9)*x.^8;
end

