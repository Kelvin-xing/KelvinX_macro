function E = expec(V, P, I, i)
    % Calculate the expectation 
    % given value function V at designated iteration I and state i
    % and transition matrix P from i to j
    E = 0;
    for j = 1:n   
        E = E + V(I, j) * P(i, j);
    end
    E = E - V(I,i) * P(i,i);
end