%
% stst.m
%
% program to update initial values with loop
%


init_k = 79.6248;
init_c =  6.47881;
init_h =  2.80879;

N_steps = 41;

for i = 1:N_steps
    
    nu = 0.1 + (i-1)*(4-0.1)/(N_steps-1);

    save parameterfile nu init_c init_k init_h    
    dynare ststmodel noclearall

    %load dynarerocks
    %init_k = decision(1,1);
    %init_c = decision(1,3);
    %init_h = decision(1,4);
   
    init_k = oo_.steady_state(2);
    init_c = oo_.steady_state(1);
    init_h = oo_.steady_state(3);
    
    [nu init_k init_c init_h]
    pause
end
