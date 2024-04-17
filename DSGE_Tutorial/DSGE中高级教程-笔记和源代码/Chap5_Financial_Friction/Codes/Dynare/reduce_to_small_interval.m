function [fu,ix,II]=reduce_to_small_interval(sig,oomega,mu,sp1)

%this function help to reduce omega interval to a small one and identified
%which is the small interval contains a root for the efficiency condition
%indicated by ix. This info will be returned to decide the interval. 
II=[];
%find_foc_difference.m serves  only as calculation of the basic elements of
%Gomega,Fomega,Gamma, etc 
[ff,G,F,Gamma,Gam_muG,Fprime]=find_foc_difference(oomega(1),sig,mu,sp1);

fu(1)=ff;
ix=0;
for ii = 2:length(oomega)   
    [ff,G,F,Gamma,Gam_muG,Fprime]=find_foc_difference(oomega(ii),sig,mu,sp1);
    fu(ii)=ff;
    if fu(ii)*fu(ii-1) < 0
        if ii > 2 % >=
            if (fu(ii)-fu(ii-1))*(fu(ii-1)-fu(ii-2))>0
                ix=ix+1;
                II(ix)=ii;
            end
        end
    end
end
