% Coefficients defining linear system for obtaining steady states of LM variables 

global A_SS C_SS h_SS pie_SS R_SS; 
  
global lmult1_SS lmult2_SS lmult3_SS; 
  
lmss_mat = zeros(4,3); 
lmss_vec = zeros(4,1); 
lmss_vec(1) = 1/C_SS; 
lmss_mat(1,1) = (C_SS*beta)/(C_SS^2*nbeta*pie_SS) - beta/(C_SS*pie_SS); 
lmss_mat(1,2) = (chi*epsil*h_SS*((phi*(pie_SS - 1)^2)/2 + 1))/A_SS; 
lmss_mat(1,3) = (phi*(pie_SS - 1)^2)/2 + 1; 
lmss_vec(2) = -chi*h_SS; 
lmss_mat(2,1) = 0; 
lmss_mat(2,2) = (C_SS*chi*epsil*((phi*(pie_SS - 1)^2)/2 + 1))/A_SS; 
lmss_mat(2,3) = -A_SS; 
lmss_vec(3) = 0; 
lmss_mat(3,1) = (C_SS*beta)/(C_SS*nbeta*pie_SS^2); 
lmss_mat(3,2) = (beta*phi*(pie_SS - 1) + beta*phi*pie_SS)/nbeta - phi*pie_SS - phi*(pie_SS - 1) - (phi*(2*pie_SS - 2)*((epsil - 1)*(nu + 1) - (C_SS*chi*epsil*h_SS)/A_SS))/2; 
lmss_mat(3,3) = (C_SS*phi*(2*pie_SS - 2))/2; 
lmss_vec(4) = 0; 
lmss_mat(4,1) = -1/(R_SS + 1)^2; 
lmss_mat(4,2) = 0; 
lmss_mat(4,3) = 0; 

lmss_vals = -lmss_mat \ lmss_vec; 
errcheck = max(abs(lmss_mat*lmss_vals + lmss_vec)); 
if errcheck > 1e-08, 
  disp('Warning: steady states of lagrange multipliers cannot be accurately determined'); 
  disp(['         errcheck = ',num2str(errcheck,12)]); 
end; 

  lmult1_SS = lmss_vals(1); 
  lmult2_SS = lmss_vals(2); 
  lmult3_SS = lmss_vals(3); 
