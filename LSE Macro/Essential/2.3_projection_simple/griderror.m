function ssr = griderror(coef)
global beta alpha depr gamma rho sigma q_nodes q_weights k_grid z_grid

k_number =  size(k_grid,1);
z_number =  size(z_grid,1);
q_number =  size(q_nodes,1);

ssr      =  0;

for i_k = 1:k_number
    for i_z = 1:z_number
        k       = k_grid(i_k);
        z       = z_grid(i_z);
        c       = consfun(k,z,coef);
        k_prime = z*k^alpha+(1-depr)*k-c;
        
        expec  = 0;
            for i_q = 1:q_number
                e_prime = sqrt(2)*sigma*q_nodes(i_q);
                z_prime = exp(rho*log(z)+e_prime);
                c_prime = consfun(k_prime,z_prime,coef);

                expec = expec + q_weights(i_q) * ...
                beta*c_prime^(-gamma)*(alpha*z_prime*k_prime^(alpha-1)+(1-depr));            
            end
            expec = expec/sqrt(pi);  %(constant in normal density & Jacobian combined)
        
        %add error to ssr
        ssr = ssr +  (expec-c^(-gamma))^2;
    end
end


disp(ssr)






