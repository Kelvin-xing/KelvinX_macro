%%**************************************************************************
%  LSE Macroeconomics Summer Program
%  Instructor: Wouter J. Den Haan
%
%  use of this program in any fee-based program requires
%  explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
%**************************************************************************
%  
%  Assignment to investigate correlation of inventories and sales
%  using different filters and to calculate standard errors using
%  heteroskedastic and autocorrelation consistent (HAC) estimators
%--------------------------------------------------------------------------

%% Section 1: read data and plot HP-filtered residuals

clear all
close all

% inventory_sales_data.mat contains the following two series 
% (i)  total final sales    called "sales"
% (ii) non-farm inventories called "inven"
% there are 261 observations from (1947Q1 to 2012Q1)

% 1967Q1 = 81 (in series starting in 47Q1)
T_67Q1 = 81;

load inventory_sales_data

% daten3.mat contains Matlab dates from 1947Q1 to 2012Q1 which are used to
% make graphs with NBER recessions

load daten3

% create hp filtered series and plot them

sales_hp      = log(sales)      - hpfilter2(log(sales),1600);
inven_hp      = log(inven)      - hpfilter2(log(inven),1600);

%  these graphs plot HP-filtered series (programs used are in the
%  subdirectory shadenber, which should be included in the path)

figure(1)
subplot('position',[0.075 0.6 0.85 0.33])
plot(daten3(81:end),sales_hp(81:end),'LineWidth',2)
title('Cyclical Sales Component')
datetick('x','YY')
set(gcf,'Renderer','zbuffer');
shadenber
axis([daten3(81) daten3(end) -0.08 0.08])
set(gca,'Ytick',(-0.08:0.02:0.08))
set(gca,'YtickLabel',{'-0.08','-0.06','-0.04','-0.02','0','0.02','0.04','0.06','0.08'})
box off
subplot('position',[0.075 0.2 0.85 0.33])
plot(daten3(81:end),inven_hp(81:end),'LineStyle','--','LineWidth',2)
datetick('x','YY')
set(gcf,'Renderer','zbuffer');
shadenber
axis([daten3(82) daten3(end-1) -0.08 0.08])
set(gca,'Ytick',(-0.08:0.02:0.08))
set(gca,'YtickLabel',{'-0.08','-0.06','-0.04','-0.02','0','0.02','0.04','0.06','0.08'})
title('Cyclical Inventory Component')
box off
maximize   % "blows" up the figure to full screen

%% Section 2: decompose series with band-pass filter and compare with HP filter
%  starting date should be 1967Q1 (81st observation)
%  end date should be as recent as possible (but some obs are lost due to
%  using band-pass filter)

% T_lost =XXX;
% omega1 =XXX;  % period of 32 quarters
% omega2 =XXX;
% 
% sales_bb_32Q_8 = bandpass1(log(sales(XXX:XXX)),omega1,omega2,T_lost);
% inven_bb_32Q_8 = bandpass1(log(inven(XXX:XXX)),omega1,omega2,T_lost);
% 
% figure(2)
% subplot('position',[0.075 0.6 0.85 0.33])
% plot(daten3(T_67Q1:end-T_lost),sales_hp(T_67Q1:end-T_lost),'--','LineWidth',2)
% hold all
% plot(daten3(T_67Q1:end-T_lost),sales_bb_32Q_8(T_lost+1:end-T_lost),'LineWidth',2)
% datetick('x','YY')
% set(gcf,'Renderer','zbuffer');
% shadenber
% axis([daten3(82) daten3(end-1) -inf +inf])
% title('Cyclical Sales - Band Pass & HP')
% box off
% subplot('position',[0.075 0.2 0.85 0.33])
% plot(daten3(T_67Q1:end-T_lost),inven_hp(T_67Q1:end-T_lost),'--','LineWidth',2)
% hold all
% plot(daten3(T_67Q1:end-T_lost),inven_bb_32Q_8(T_lost+1:end-T_lost),'LineWidth',2)
% datetick('x','YY')
% set(gcf,'Renderer','zbuffer');
% shadenber
% axis([daten3(82) daten3(end-1) -inf +inf])
% title('Cyclical Inventories - Band Pass & HP')
% box off
% maximize   % "blows" up the figure to full screen

%% Section 3: decompose series with band-pass filter focusing on high frequencies
%  starting date should be 1967Q1 (81st observation)
%  end date should be as recent as possible (but some obs are lost due to
%  using band-pass filter)

% T_lost = XXX;
% omega1 = XXX;  % period of 4 quarters
% omega2 = XXX;
% 
% sales_bb_4Q_8 = bandpass1(log(sales(XXX:XXX)),omega1,omega2,T_lost);
% inven_bb_4Q_8 = bandpass1(log(inven(XXX:XXX)),omega1,omega2,T_lost);
% 
% figure(3)
% plot(daten3(T_67Q1:end-T_lost),sales_bb_4Q_8(T_lost+1:end-T_lost),'--','LineWidth',2)
% hold all
% plot(daten3(T_67Q1:end-T_lost),inven_bb_4Q_8(T_lost+1:end-T_lost),'LineWidth',2)
% datetick('x','YY')
% set(gcf,'Renderer','zbuffer');
% shadenber
% axis([daten3(82) daten3(end-1) -inf +inf])
% title('Sales & Inventories - High Frequency Component')
% box off
% maximize   % "blows" up the figure to full screen
% 

% %% Section 4: calculate standard errors for the correlation coefficient
% %  you can impose that the mean equals zero
% 
% I_data = 1;
% % I_data = 1: use high-frequency data
% % I_data = 2: use business-cycle frequency data
% % I_data = 3: use some completely unrelated data (which are more persistent)
% 
% 
% if I_data == 1;
%     data     = [inven_bb_4Q_8(T_lost+1:end-T_lost) sales_bb_4Q_8(T_lost+1:end-T_lost)];
% end
% 
% if I_data == 2;
%     data     = [inven_bb_32Q_8(T_lost+1:end-T_lost) sales_bb_32Q_8(T_lost+1:end-T_lost)];
% end
%     
% if I_data == 3;
%     TT = 100000;
%     temp1    = randn(TT+100,1)*3;
%     temp2    = randn(TT+100,1)*2;
%     temp3    = randn(TT+100,1);
%     temp4    = randn(TT+100,1);
%     for t = 3:TT+100;
%         temp1(t,1) = 0.8*temp1(t-1,1)+0.17*temp1(t-2,1)+temp3(t,1)-0.5*temp3(t-1,1);
%         temp2(t,1) = 0.5*temp2(t-1,1)+0.45*temp2(t-2,1)+temp4(t,1);
%     end
%     data     = [temp1 3*temp1+temp2];
%     data     = data(101:end,:)-ones(TT,1)*mean(data(101:end,:));
% end
% 
% [TT temp]  = size(data);
% JJ         = 10;
% 
% % estimate statistics of interest
% temp     = cov(data);
% var_i    = temp(1,1);
% var_s    = temp(2,2);
% corr_is  = temp(2,1)/sqrt(var_i*var_s); %this is the correlation coefficient we are interested in.
% ratio_is = sqrt(var_i/var_s);
% 
% % construct the residuals, i.e., h(x_t;theta_hat)
% 
% h_resid  = XXX;
% 
% % obtain estimate for the 2x2 D matrix
% 
% D        = [XXX,XXX;XXX,XXX];
% 
% % obtain estimate for the sigma_0 matrix (the variance-covariance matrix)
% % using the truncated and the Newey-West kernel
% 
% sigma_0_tr    = (h_resid*h_resid')/TT;
% sigma_0_nw    = (h_resid*h_resid')/TT;
% 
% for j = 1:JJ
%     temp = h_resid(:,j+1:end)*(h_resid(:,1:end-j)')/TT;
%     sigma_0_tr = sigma_0_tr + (temp+temp')*XXX;
%     sigma_0_nw = sigma_0_nw + (temp+temp')*XXX;
% end
% 
% % obtain an estimate for sigma_0 using VARHAC
% 
% 
% % varhac has the following structure: 
% %
% % varhac(datin,it1,it2,k1,k2,imax,ilag,imodel,imean)
% 
% % datin:      input matrix where each row corresponds to a different observation
% % it1:        first observation to be used
% % it2:        last  observation to be used
% % k1:         first column to be used
% % k2:         last  column to be used
% % imax:       maximum lag order considered  (if imax = 0, then no
% %             correction for serial correlation will be made)
% % ilag:       if equal to 1, then all elements enter with the same lag
% %             if equal to 2, then the own lag can enter with a different lag
% %             if equal to 3, then only the own lag enters
% % imodel:     if equal to 1, then AIC is used
% %             if equal to 2, then SCHWARTZ is used
% %             if equal to 3, then a fixed lag order equal to imax is used
% % imean:      if equal to 1, then the mean will be subtracted from each series
% 
% it1 = 1;      
% it2 = TT;     
% k1 = 1;
% k2 = 2;
% imax = XXX;
% ilag = 1;   % this is the easiest choice to begin with (limits the number of VARs considered)
% imodel = 3; % this is the easiest choice to begin with (limits the number of VARs considered)
% imean  = 1; % the analysis is based on the hypothesis that the mean is zero. This is not exactly true 
%             % even though the elements of BB  add up to zero
%             
% 
% sigma_0_vh=varhac(h_resid',it1,it2,k1,k2,imax,ilag,imodel,imean);
% 
% % calculate the V matrix V = inv(D*inv(sigma_0)*D');
% 
% V_tr = inv((D/sigma_0_tr)*D');
% V_nw = inv((D/sigma_0_nw)*D');
% V_vh = inv((D/sigma_0_vh)*D');
% 
% % calculate standard error corr_is
% 
% SE_tr = sqrt(V_tr(2,2)/TT);
% SE_nw = sqrt(V_nw(2,2)/TT);
% SE_vh = sqrt(V_vh(2,2)/TT);
% 
% disp([SE_tr SE_nw SE_vh])
% disp([V_tr(2,2) V_nw(2,2) V_vh(2,2)])
% 
% 
