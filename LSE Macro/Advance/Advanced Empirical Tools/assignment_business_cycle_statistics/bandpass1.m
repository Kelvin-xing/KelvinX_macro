function output = bandpass1(input,omega1,omega2,T_lost)

% input:  vector (or matrix) with observation in columns    
% omega1: lower bound for frequency of series included
% omega2: upper bound for frequency of series included
% T_lost: number of observations lost on each side
%
% band-pass filter uses 2*T_lost + 1 observations, but because of symmetry
% we only have to calcuate T_lost + 1 coefficients of the band-pass filter

T      = size(input,1);
T_used = 2*T_lost+1;

%% first construct the coefficients of the band-pass filter (elements of BB)

BB   = zeros(T_lost+1,1);

BB(1,1)  = XXX; % coefficient on current observation

for j = 1:T_lost
    BB(j+1,1) = XXX;
end

% correct the elements of BB such that they add up to zero (make sure you
% know why we do this)

sum_bb  = XXX;

BB      = XXX; % 


%% filter the data

output      = 0*input;      % ensures that output has the same size as the input series

for t = T_lost+1:T-T_lost
    output(t,:)   = BB(1,1) * input(t,:);
    for j = 1:T_lost
        output(t,:)  = XXX;
    end
end


end

