%%% Sign restrictions, a summary
% Michele Piffer, m.b.piffer@gmail.com
% 8 / 7 / 2014

clc
clear
close all


%% SETTING UP THE DGP

T = 50;      
k = 2; % variables
p_true = 4; % lags

% parameters
A1 = [0.6362, -0.0012; 0.0190, 0.5782]; 
A2 = [-0.0168, -0.0285; 0.5211, -0.3041]; 
A3 = [0.0273, -0.0028; 0.1568, 0.2229]; 
A4 = [0.1517, -0.0198; -0.7600, -0.3168]; 
A_true = [A1, A2, A3, A4];

% check that the VAR is stationary
companion = [A_true; kron(eye(p_true-1,p_true-1),eye(k)),zeros(k*(p_true-1),k)];

[V, C] = eig(companion);
eigenvalues = max(C); 
dominant_root = max(max(abs(C)));
assert( abs(dominant_root) <1)

% shocks
reset(RandStream.getDefaultStream);
S_true = mvnrnd(zeros(k,1),eye(k),T)'; % assume structural shocks have unitary variance

P_true = [0.025, 0.009; 0.009, -0.387]; % maps structural shocks into reduced form shocks. Shocks have asymmetric effect on impact

R_true = P_true*S_true; % reduced form shocks

% generate the data. Initial values are zero, which is the expected value of both variables
Y = zeros(k,T); 
Y(:,1) =  R_true(:,1);
for j=2:T
    step1 = [zeros(k,p_true-1), Y];
    step2 = step1(:,j-1:j-1+p_true-1);
    step3 = fliplr(step2);
    step4 = reshape(step3,k*p_true,1);
    Y(:,j) =  A_true*step4 + R_true(:,j);
end

figure(1)
subplot(1,2,1)
plot(Y(1,:)), title('Variable 1')
axis([1 T -Inf Inf])
subplot(1,2,2)
plot(Y(2,:)), title('Variable 2')
axis([1 T -Inf Inf])

%% TRUE IRFs (use true P matrix and true VAR parameters)

T_irf = 20;
shock = 1; % magnitude of structural shock given as impulse (one standard deviation)

IRF_true = zeros(k*k,T_irf); % contains impulse responses of all k variables to impulse to variable 1, then to variable 2

shocked_variab = 1;
IRF_true(1:2,:)  = zFunction_IRFs(A_true, P_true, T_irf, shocked_variab, shock);
shocked_variab = 2;
IRF_true(3:4,:)  = zFunction_IRFs(A_true, P_true, T_irf, shocked_variab, shock);

figure(1)
for ii = 1:k*k
    subplot(k,k,ii)
    plot([1:1:T_irf], IRF_true(ii,:), '-- b', 'Linewidth', 2), hold on 
    axis([1 T_irf -Inf Inf]),
    set(gca,'box','off'), 
    set(gcf, 'PaperPositionMode', 'auto'); 
    xlabel('time', 'FontAngle','italic'),

    if ii == 1 
        title('Response of Variable 1')
        ylabel('Impulse to variable 1')
        legend('True')
    elseif ii == 2
        title('Response of Variable 2')
    elseif ii == 3
        ylabel('Impulse to Variable 2')
    end
end  


%% ESTIMATION

constantD = 0;
p = p_true;
[Ahat,Consthat,Residuals,T_used,Sigmahat] = zFunction_RVAR(Y,p,constantD);           

%% HALF TRUE IRFs (use true P matrix and estimated VAR parameters)

T_irf = 20;
shock = 1; % structural shock given

IRF_halftrue = zeros(k*k,T_irf); 

shocked_variab = 1;
IRF_halftrue(1:2,:)  = zFunction_IRFs(Ahat, P_true, T_irf, shocked_variab, shock);
shocked_variab = 2;
IRF_halftrue(3:4,:)  = zFunction_IRFs(Ahat, P_true, T_irf, shocked_variab, shock);

figure(2)
for ii = 1:k*k
    subplot(k,k,ii)
    plot([1:1:T_irf], IRF_true(ii,:), '-- b', 'Linewidth', 2), hold on 
    plot([1:1:T_irf], IRF_halftrue(ii,:), 'b', 'Linewidth', 2), hold on 
    axis([1 T_irf -Inf Inf]),
    set(gca,'box','off'), 
    set(gcf, 'PaperPositionMode', 'auto'); 
    xlabel('time', 'FontAngle','italic'),

    if ii == 1 
        title('Response of Variable 1')
        ylabel('Impulse to variable 1')
        legend('True', 'Half true')
    elseif ii == 2
        title('Response of Variable 2')
    elseif ii == 3
        ylabel('Impulse to Variable 2')
    end
end  

% notes: if T is big, get the true IRF pretty well


%% CHOLESKY DECOMPOSITION

P_chol = chol(Sigmahat)';

IRF_chol = zeros(k*k,T_irf); 

shocked_variab = 1;
IRF_chol(1:2,:)  = zFunction_IRFs(Ahat, P_chol, T_irf, shocked_variab, shock);
shocked_variab = 2;
IRF_chol(3:4,:)  = zFunction_IRFs(Ahat, P_chol, T_irf, shocked_variab, shock);

figure(3)
for ii = 1:k*k
    subplot(k,k,ii)
    plot([1:1:T_irf], IRF_true(ii,:), '-- b', 'Linewidth', 2), hold on 
    plot([1:1:T_irf], IRF_chol(ii,:), 'r', 'Linewidth', 2), hold on 
    axis([1 T_irf -Inf Inf]),
    set(gca,'box','off'), 
    set(gcf, 'PaperPositionMode', 'auto'); 
    xlabel('time', 'FontAngle','italic'),

    if ii == 1 
        title('Response of Variable 1')
        ylabel('Impulse to variable 1')
        legend('True', 'Cholesky')
    elseif ii == 2
        title('Response of Variable 2')
    elseif ii == 3
        ylabel('Impulse to Variable 2')
    end
end  

% notes: even with very big T performs poorly, since imposes that variable
% 1 does not respond to structural shock to 2

%% SIGN RESTRICTIONS

% generate orthogonal matrices
M = 200;

% method = 1 % uses rotations 
% method = 2 % uses reflections
% method = 3 % uses QR decompositions
% Q_all = zFunction_OrthogonalMat(method, M);

% alternatively, use all methods above:
Q_all = cat(3, zFunction_OrthogonalMat(1, M), zFunction_OrthogonalMat(2, M), zFunction_OrthogonalMat(3, M)); 
M = 3*M;


% choose starting decomposition
P_start = P_chol;
% can also use eigenvalue-eigenvector decomposition

% compute set of IRFs
IRF_signrest_all = zeros(k*k,T_irf,M); 

for i = 1:M
    
    Q = Q_all(:,:,i);
    assert(max(max(Q'*Q - eye(2))) < 0.001)
    
    shocked_variab = 1;
    IRF_signrest_all(1:2,:,i)  = zFunction_IRFs(Ahat, Q*P_start, T_irf, shocked_variab, shock);
    shocked_variab = 2;
    IRF_signrest_all(3:4,:,i)  = zFunction_IRFs(Ahat, Q*P_start, T_irf, shocked_variab, shock);
    
end

% plot all the unrestricted IRFs
figure(4)
for ii = 1:k*k
    subplot(k,k,ii)
    for i = 1:M
        plot([1:1:T_irf], IRF_signrest_all(ii,:,i), '-- k', 'Linewidth', 1), hold on 
    end
    a = plot([1:1:T_irf], IRF_true(ii,:), '-- b', 'Linewidth', 2), hold on 
    axis([1 T_irf -Inf Inf]),
    set(gca,'box','off'), 
    set(gcf, 'PaperPositionMode', 'auto'); 
    xlabel('time', 'FontAngle','italic'),

    if ii == 1 
        title('Response of Variable 1')
        ylabel('Impulse to variable 1')
        legend(a, 'True')
    elseif ii == 2
        title('Response of Variable 2')
    elseif ii == 3
        ylabel('Impulse to Variable 2')
    end
end  

IRF_signrest_restricted = IRF_signrest_all;


% impose restriction that after a positive shock to variable 1, variable 1 increases in the first period
whichvariable = 1;                       % variables on which restrictions are imposed
shockwhere = 1;                          % shock whose IRF is restricted
when = 1;                                % time period in which restrictions are imposed
condition = IRF_signrest_restricted((shockwhere-1)*2+whichvariable,when,:) >= 0; % inequality controls for the condition imposed
IRF_signrest_restricted = IRF_signrest_restricted(:,:,condition);
M_kept = size(IRF_signrest_restricted,3) % number of identified models that meet the restrictions

% impose restriction that after a positive shock to variable 2, variable 2 decreases in the first period
whichvariable = 2;                       
shockwhere = 2;                          
when = 1;                                
condition = IRF_signrest_restricted((shockwhere-1)*2+whichvariable,when,:) <= 0; 
IRF_signrest_restricted = IRF_signrest_restricted(:,:,condition);
M_kept = size(IRF_signrest_restricted,3) 

% impose restriction that after a positive shock to variable 2, variable 2 decreases in the second period
whichvariable = 2;                       
shockwhere = 2;                          
when = 2;                                
condition = IRF_signrest_restricted((shockwhere-1)*2+whichvariable,when,:) <= 0;
IRF_signrest_restricted = IRF_signrest_restricted(:,:,condition);
M_kept = size(IRF_signrest_restricted,3) 

% impose restriction that after a positive shock to variable 2, variable 2 increases in the seventh period
whichvariable = 2;                       
shockwhere = 2;                          
when = 7;                                
condition = IRF_signrest_restricted((shockwhere-1)*2+whichvariable,when,:) >= 0;
IRF_signrest_restricted = IRF_signrest_restricted(:,:,condition);
M_kept = size(IRF_signrest_restricted,3) 

% impose restriction that after a positive shock to variable 2, variable 1 increases in the first period
whichvariable = 1;                       
shockwhere = 2;                          
when = 1;                                
condition = IRF_signrest_restricted((shockwhere-1)*2+whichvariable,when,:) >= 0;
IRF_signrest_restricted = IRF_signrest_restricted(:,:,condition);
M_kept = size(IRF_signrest_restricted,3) 

% plot only the restricted IRFs
figure(5)
for ii = 1:k*k
    subplot(k,k,ii)
    for i = 1:M_kept
        b = plot([1:1:T_irf], IRF_signrest_restricted(ii,:,i), '-- k', 'Linewidth', 1), hold on 
    end
    a = plot([1:1:T_irf], IRF_true(ii,:), '-- b', 'Linewidth', 2), hold on 
    axis([1 T_irf -Inf Inf]),
    set(gca,'box','off'), 
    set(gcf, 'PaperPositionMode', 'auto'); 
    xlabel('time', 'FontAngle','italic'),

    if ii == 1 
        title('Response of Variable 1')
        ylabel('Impulse to variable 1')
        legend([a b], 'True', 'Sign restricted')
    elseif ii == 2
        title('Response of Variable 2')
    elseif ii == 3
        ylabel('Impulse to Variable 2')
    end
end  

% compute median and median target
interval = 99;
[IRF_signrest_median, IRF_signrest_median_target, IRF_signrest_low, IRF_signrest_high, mt] = zFunction_MedianTarget(IRF_signrest_restricted,interval);

grey = [0.8,0.8,0.8]; 


figure(6)
for ii = 1:k*k
    subplot(k,k,ii)
    fill([1:1:T_irf,fliplr(1:1:T_irf)],...
        [IRF_signrest_low(ii,:),fliplr(IRF_signrest_high(ii,:))], grey,'EdgeColor','none'), hold on
    a = plot([1:1:T_irf], IRF_true(ii,:), '-- b', 'Linewidth', 2), hold on 
    b = plot([1:1:T_irf], IRF_signrest_median(ii,:), 'y', 'Linewidth', 2), hold on 
    c = plot([1:1:T_irf], IRF_signrest_median_target(ii,:), 'k', 'Linewidth', 2), hold on 
    axis([1 T_irf -Inf Inf]),
    set(gca,'box','off'), 
    set(gcf, 'PaperPositionMode', 'auto'); 
    xlabel('time', 'FontAngle','italic'),

    if ii == 1 
        title('Response of Variable 1')
        ylabel('Impulse to variable 1')
        legend([a b c], 'True', 'Median', 'Median target')
    elseif ii == 2
        title('Response of Variable 2')
    elseif ii == 3
        ylabel('Impulse to Variable 2')
    end
end  

% the median and the upper/lower bounds are not indicative of one single
% identified model, since they jump from one model to another. The median
% target is defined as the single identification that is closest to the
% median. 


