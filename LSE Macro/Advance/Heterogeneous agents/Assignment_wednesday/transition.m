function M = transition(apg,apb)

global a_grid tau mu T alpha delta beta sigma phi N

% This programs provides a transition matrix M associated with the policy
% functions apg and apb.

a_prime = [apg,max(apb,phi)];

% Look-up index.

I = zeros(N,2);

for j=1:2
    for i=1:N
        I(i,j) = find(a_prime(i,j)>=a_grid,1,'last');
    end
end

I = min(I,N-1);

% Linear interpolation of probabilities

for j=1:2
    rho(:,j) = (a_prime(:,j)-a_grid(I(:,j)+1))./(a_grid(I(:,j))-a_grid(I(:,j)+1));
end

% Transition matrix

M1 = sparse(zeros(N,N));

for i=1:N
    M1(i,I(i,1):I(i,1)+1) = [rho(i,1),(1-rho(i,1))];
end
M1 = sparse(M1);

M2 = sparse(zeros(N,N));

for i=1:N
    M2(i,I(i,2):I(i,2)+1) = [rho(i,2),(1-rho(i,2))];
end
M2 = sparse(M2);

% Transition matrix

M = [T(1,1)*M1, T(1,2)*M1;T(2,1)*M2, T(2,2)*M2];

end

