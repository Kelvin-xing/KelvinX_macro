clear
clc

icheck    = 0; % if = 1 plot function, o.w. plot function subtracted by linear
iregular  = 0; % if = 1 use regular polynomials, o.w. use chebyshev polynomials

for i = 1:4
        if i==1
            nnodes = 2;
        end
        if i==2
            nnodes = 5;
        end
        if i==3
            nnodes = 10;
        end
        if i==4
            nnodes = 25;
        end
        fprintf(' number of nodes is %4.0f \n',nnodes)

chebnodes = chebnode(nnodes);

% evaluate function values

fvalues = zeros(nnodes,1);

for j = 1:nnodes
    if chebnodes(j) < 0
    fvalues(j,1)=0.90*chebnodes(j,1);
    else
    fvalues(j,1)=0.92*chebnodes(j,1);
    end
end

% construct explanatory values

X = ones(nnodes,nnodes);
X(:,2) = chebnodes;

if iregular == 1
    for j = 3:nnodes
        X(:,j) = chebnodes.^(j-1);
    end
else
    for j = 3:nnodes
        X(:,j) = 2*chebnodes.*X(:,j-1)-X(:,j-2);
    end
end

% do projection

beta = (X'*X)\(X'*fvalues);


% do comparison

T = 201;
X = ones(T,nnodes);
for t = 1:T
    X(t,2)= -1+(t-1)*2/(T-1);
end
if nnodes > 2
    if iregular == 1
        for j = 3:nnodes
            X(:,j) = X(:,2).^(j-1);
        end
    else
        for j = 3:nnodes
            X(:,j) = 2*X(:,2).*X(:,j-1)-X(:,j-2);
        end        
    end
end

temp = zeros(T,2);
temp2 = zeros(T,2);

x = -1;
for t = 1:T
    if x < 0
    temp(t,1)=0.90*x;
    else
    temp(t,1)=0.92*x;
    end
    x=x + 2/(T-1);
end

temp(:,2) = X*beta;

temp2(:,1) = temp(:,1)-0.9*X(:,2);
temp2(:,2) = temp(:,2)-0.9*X(:,2);


if icheck ==1
    plot(X(:,2),temp)
    pause
else
    plot(X(:,2),temp2)
    pause
%    plot(temp(:,1)-temp(:,2))
%    pause
end


% do comparison outside range

fprintf('extrapolation \n')

T = 301;
upper = 3;
X = ones(T,nnodes);
for t = 1:T
    X(t,2)= -1+(t-1)*upper/(T-1);
end
if nnodes > 2
    if iregular == 1
        for j = 3:nnodes
            X(:,j) = X(:,2).^(j-1);
        end
    else
        for j = 3:nnodes
            X(:,j) = 2*X(:,2).*X(:,j-1)-X(:,j-2);
        end        
    end
end

temp = zeros(T,2);
temp2 = zeros(T,2);

x = -1;
for t = 1:T
    if x < 0
    temp(t,1)=0.90*x;
    else
    temp(t,1)=0.92*x;
    end
    x=x + upper/(T-1);
end

temp(:,2) = X*beta;

temp2(:,1) = temp(:,1)-0.9*X(:,2);
temp2(:,2) = temp(:,2)-0.9*X(:,2);


if icheck ==1
    plot(X(:,2),temp)
    pause
else
    plot(X(:,2),temp2)
    pause
%    plot(temp(:,1)-temp(:,2))
%    pause
end

%[X(:,2) temp(:,1) temp(:,2)]


end


