%%
% this file should be appended to newpeaproi0.m (when you have that part
% running)
% 
%%
% part 6: (is in separate file newpeaproi0B.m)
% uses the solution calculated above and calculates the solution for bond prices
% it calculates them B_iter times and then calculates some summary
% statistics. 
% finally it plots those summary statistics. On the x-axis you have the
% number of the sample; in blue you have the outcome for the statistic for
% the different maturities for the non-quad PEA and in black the
% corresponding numbers for the quad PEA

%% 6 Bond Prices
% =====================================================



B_iter   = 25; %numer of samples considered
B_mat    =  3; %maximum maturity considered

psi_bond = zeros(size(psi,1),B_iter,B_mat,2);

for iter = 1:B_iter
    
    epsi      = randn(T,1)*sigma;
    lnts(1)   = 0;
    ks(1)   = k_ss;
    for ti  = 2:T
        lnts(ti)   = rho*lnts(ti-1) + epsi(ti);
        ts(ti)     = exp(lnts(ti));
        ts_S(ti)   = lnts(ti);
        %construct scaled state variables
        ks_S(ti-1) = log(ks(ti-1));
        xx         = makepoly([po_k po_t],[ks_S(ti-1) ts_S(ti)]);
    	cs(ti)     = exp(-xx*psi/nu);
        ks(ti)     = ts(ti)*ks(ti-1)^alpha+(1-delta)*ks(ti-1)-cs(ti);
    end;

    X = makepoly([po_k po_t],[ks_S(T1-1:T-B_mat-1) ts_S(T1:T-B_mat)]);

    for j = 1:B_mat
        %XXX complete the following expression
        Y = 
        psi_b  = fminsearch(@(coef) rssfcn(coef,Y,X),psi,options);    
        disp([iter 1 j])
        psi_bond(:,iter,j,1)=psi_b;
        maturity = j;
        %XXX in RHSrealization_bond complete expressions
        Y = numi(@RHSrealization_bond,N_herm);                    
        psi_b = (X'*X)\X'*log(Y(1:end-B_mat+1));            
        psi_bond(:,iter,j,2)=psi_b;
        disp([iter 2 j])
    end
    
end    

% now do some simulations usings the same state variables but the different
% solutions for bond prices

yield = zeros(size(X,1),B_iter,B_mat,2);
CC1 = zeros(B_iter,B_mat);
CC2 = zeros(B_iter,B_mat);


for iter = 1:B_iter
    for j = 1:B_mat
        for jj = 1:2
            yield(:,iter,j,jj) = 400*(exp(-X*psi_bond(:,iter,j,jj)/j)-1);
        end
    end
end

% calculate means and standard deviations for one-period bonds

MM1 = mean(yield(:,:,1,1))';
SS1 =  std(yield(:,:,1,1))';
MM2 = mean(yield(:,:,1,2))';
SS2 =  std(yield(:,:,1,2))';

% calculate means and standard deviations for bonds with maturity > 1

for j = 2:B_mat    
    MM1 = [MM1 mean(yield(:,:,j,1))'];
    SS1 = [SS1  std(yield(:,:,j,1))'];
    MM2 = [MM2 mean(yield(:,:,j,2))'];
    SS2 = [SS2  std(yield(:,:,j,2))'];
end

%calculate autocorrelation coefficients

for iter = 1:B_iter
    for j = 1:B_mat
        xxx1 = cov(yield(2:end,iter,j,1),yield(1:end-1,iter,j,1));
        CC1(iter,j)= xxx1(1,2)/(sqrt(xxx1(1,1))*sqrt(xxx1(2,2)));
        xxx2 = cov(yield(2:end,iter,j,2),yield(1:end-1,iter,j,2));
        CC2(iter,j)= xxx2(1,2)/(sqrt(xxx2(1,1))*sqrt(xxx2(2,2)));        
    end
end

figure(3)
plot(MM1,'color','b')
hold
plot(MM2,'color','k')
title('means')

figure(4)
plot(SS1,'color','b')
hold
plot(SS2,'color','k')
title('standard deviations')

figure(5)
plot(CC1,'color','b')
hold
plot(CC2,'color','k')
title('autocorrelation')

