function [omega,G,F,Gamma,Gam_muG,Fprime,er] = get_omega_cond_on_sigma(sig,mu,sp1)
%when given sigma, mu and spread, find an omega that solves the efficiency
%condition, i.e., the foc of the entrepreneurial utility maximization
%problem.
oomega=.001:.001:1.4;
er=0;
%help reduce the above large interval of omega to a small one, in the
%samller interval, we can use @fzero function to find a root.
[fu,ix,II]=reduce_to_small_interval(sig,oomega,mu,sp1);

if ix > 1
    error('fatal (get_omega_cond_on_sigma) multiple solutions')
end

%if fail to find a root in the given interval above, try to find it in a
%larger interval:from oomega=".001:.001:1.4;" to "1.4:.001:3;"
if ix == 0
    %if it looks like there might be a zero at a higher level of omega, try looking there
    if (fu(end)-fu(end-1)) < 0
        oomega=1.4:.001:3;
        [fu,ix,II]=reduce_to_small_interval(sig,oomega,mu,sp1);
    else
        omega=0;
        er=1;
        G=[];F=[];Gamma=[];Gam_muG=[];Fprime=[];
        return
    end
    %if still not found, fails completely.
    if ix == 0
        omega=0;
        er=1;
        G=[];F=[];Gamma=[];Gam_muG=[];Fprime=[];
        return
    end
end

%else, ix ==1, there is an unqie root.
x=[oomega(II(1)-1) oomega(II(1))];
[xx,fval,exitflag] = fzero('find_foc_difference',x,optimset,sig,mu,sp1);
if exitflag ~= 1 | abs(fval) > .1e-10 | abs(imag(xx)) > .1e-10
    error(' fatal (get_omega_cond_on_sigma) failed to find a zero')
end
%after obtaining omega, recomputing the basic elements of the computation
omega=xx;
[ff,G,F,Gamma,Gam_muG,Fprime]=find_foc_difference(omega,sig,mu,sp1);
