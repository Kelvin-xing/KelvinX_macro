function [sigma,omegabar,Gamma,Gam_muG,ixx] = get_omega_cond_Fomegabar(risk_spread,Fomegabar,mu)

%given Fomega, mu and risk spread, then optimal omegabar is returned;
%sigma is also returned which is associated with omegabar;
%an indicator is also returned to indicate possible various errors;
ixx=0;
%omega can not be zero, so here it choose a very small number as a start
omegaArr=.0000001:.1:.99;

position_of_sign_change=0;
omega=omegaArr(1);
[ff(1),ix] = find_foc_difference(omega,risk_spread,Fomegabar,mu);
for ii = 2:length(omegaArr)
    omega=omegaArr(ii);
    [ff(ii),ix] = find_foc_difference(omega,risk_spread,Fomegabar,mu);
    
    %find two continuous values which has opposite sign, this means that there is omega in the interval that solve the foc 
    if ff(ii)*ff(ii-1) < 0 && ix == 0
        position_of_sign_change=position_of_sign_change+1;
        II(position_of_sign_change)=ii;
    end    
end

if position_of_sign_change > 1 
    ixx=2;
    disp('(getomega) multiple solutions to Foc of utility maximization problem')
    sigma=[];
    omegabar=[];
    Gamma=[];
    Gam_muG=[];
    return
end
if position_of_sign_change < 1 
    ixx=1;
    disp('(getomega) no solution to  Foc of utility maximization problem')    
    sigma=[];
    omegabar=[];
    Gamma=[];
    Gam_muG=[];
    return
end
% if reach here, ij should equal to one,
% that is to say that we reduce the interval of [0,1] to a small one
omega1=omegaArr(II(1)-1);
omega2=omegaArr(II(1));
                                     
opt=optimset('diagnostics','on','TolX',1e-16,'TolFun',1e-16);
[omegabar,fval,exitflag,~] =   fzero(@find_foc_difference,[omega1,omega2],opt,risk_spread,Fomegabar,mu);
if abs(fval ) > .1e-9 || abs(imag(omega))>.1e-10 || exitflag <= 0
    error('fatal (getomega) failed to find omega')
end
%after finding the omegabar, recompute Gamma and Gam_muG and sigma that
%satisfis logncdf(omega,-sigma^2/2,sigma)==Fomegabar
[fff,ix,Gamma,Gam_muG,sigma] = find_foc_difference(omegabar,risk_spread,Fomegabar,mu);
if ix ~= 0
    ixx=3;
    error('fatal(getom(getoega) optimal omega and calibration at odd with FOC ')
end