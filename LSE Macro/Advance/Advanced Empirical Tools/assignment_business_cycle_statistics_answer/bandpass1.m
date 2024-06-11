function output = bandpass1(input,omega1,omega2,T_lost)

% input:  vector (or matrix) with observation in columns    
% omega1: lower bound for frequency of series included
% omega2: upper bound for frequency of series included
% T_lost: number of observations lost on each side
%
% band-pass filter uses 2*T_lost + 1 observations, but because of symmetry
% we only have to calcuate T_lost + 1 coefficients

T      = size(input,1);
T_used = 2*T_lost+1;

BB   = zeros(T_lost+1,1);

BB(1,1)  = (omega2-omega1)/pi; % coefficient on current observation

for j = 1:T_lost
    BB(j+1,1) = (sin(omega2*j)-sin(omega1*j))/(pi*j);
end

sum_bb  = 2*sum(BB(2:end,1))  +BB(1,1);

BB      = BB-sum_bb/T_used; % adjust band-pass coefficients to ensure they add up to zero

output      = 0*input;      % ensures that output has the same size as the input series

for t = T_lost+1:T-T_lost
    output(t,:)   = BB(1,1) * input(t,:);
    for j = 1:T_lost
        output(t,:)  = output(t,:)  + BB(j+1,1) *(input(t-j,:)+input(t+j,:));
    end
end


end

