function [base]=stoch_simul_MMB(base)

% Copyright (C) 2001-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

cd(base.setpath(base.models(base.epsilon),:)); % MODELBASE: change directory to the specific model folder to be solved

global M_ options_ oo_ it_

options_old = options_;
options_.order = 1;     % ADDED BY MODELBASE TEAM
if options_.linear
    options_.order = 1;
end
if options_.order == 1
    options_.replic = 1;
elseif options_.order == 3
    options_.k_order_solver = 1;
end

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

TeX = options_.TeX;

iter_ = max(options_.periods,1);
if M_.exo_nbr > 0
    oo_.exo_simul= ones(iter_ + M_.maximum_lag + M_.maximum_lead,1) * oo_.exo_steady_state';
end

%check_model;  % commented out by Sebastian, otherwise Matlab breaks down
%under Dynare 4.2 (not 4.1)

%warning off;
[oo_.dr, info] = resol(oo_.steady_state,0); % solve

if info(1)
disp(' ');
disp('NO SOLUTION FOUND');
disp(' ');

      if strcmp(base.innos(1,:),'all_shocks')  % this is the case if only one model and then all shocks for this model have been chosen
        shocks= M_.exo_names(M_.exo_names_orig_ord,:);  
        inv_lgx_orig_ord_(M_.exo_names_orig_ord)=(1:M_.exo_nbr)';
        base.innos = shocks; % put all shocks in the choice vector for the IRFs
        base.namesshocks = shocks; % put the right shock names for correct plots
      end      
      base.info(base.models(base.epsilon)) = 1;
      for p=1:size(base.innos,1)
          base.pos_shock(p,base.models(base.epsilon))=0;
      end
     base.info(base.models(base.epsilon)) = 1;
      for p=1:size(base.innos,1)
          base.pos_shock(p,base.models(base.epsilon))=0;
      end
      
else
     base.info(base.models(base.epsilon)) = info(1);
      
     
      
%Theoretical ACFs and Variances
      
        nvar  = length(oo_.dr.order_var);         
        ivar  = transpose(1:nvar);   
        options_.irf = base.horizon; % horizon for ACFs
        cd('..');
        [Gamma_y,stationary_vars] = th_autocovariances(oo_.dr,ivar,M_,options_,1); 
        
        oo_.var = Gamma_y{1}; % Variances
        base.VAR.(num2str(deblank(base.names(base.models(base.epsilon),:)))) =Gamma_y{1};
        base.VARendo_names.(num2str(deblank(base.names(base.models(base.epsilon),:))))=M_.endo_names;
        
        
        R=[];
        for i=1:base.horizon
           oo_.autocorr{i}=Gamma_y{i+1};
           R= [ R, diag(oo_.autocorr{i}) ];           
        end
        
        if base.option(1) == 1; % If ACF are selected...
            base.AUTR.(num2str(deblank(base.names(base.models(base.epsilon),:))))(:,:)=[ones(size(R,1),1),R]; 
            base.AUTendo_names.(num2str(deblank(base.names(base.models(base.epsilon),:))))(:,:)=M_.endo_names(ivar,:);
        else 
            base.AUTR.(num2str(deblank(base.names(base.models(base.epsilon),:)))) = [];
            base.AUTendo_names.(num2str(deblank(base.names(base.models(base.epsilon),:)))) = [];
        end
        
        
% Impulse response functions
        
        if base.option(2)==1   % iIf IRF are selected...
        options_.irf = base.horizon; % horizon for IRFs
        shocks = M_.exo_names(M_.exo_names_orig_ord,:);  % put shocks in the right order for Dynare
        inv_lgx_orig_ord_(M_.exo_names_orig_ord)=(1:M_.exo_nbr)'; % save the order
        
        if base.innos(1,:) == 'allshocks'  % this is the case if only one model and then all shocks for this model have been chosen
            
            % set up a menu to choose specific shocks
            shockmenu = [num2str('menu(''Choose the shocks you would like to pick. Press the button "continue" when finished.'', ''All ') num2str(size(shocks,1)) num2str(' shocks''')];
            for j=1:size(shocks,1)
                shockmenu = [shockmenu num2str(', ') num2str('''') shocks(j,:) num2str('''')]; 
            end
            shockmenu= [shockmenu ', ''Continue''' num2str(');')];
            shock = 1; % initialize shock
            chosenshocks = [];
            if size(shocks,1)>40 % display a warning if there is a large number of shocks
                menu('A T T E N T I O N : The model features more than 40 shocks. Choosing all shocks might cause a breakdown of Matlab.','Continue');
            end
            disp(' ');
            disp('Selected shocks of the model:');
            while shock<size(shocks,1)+2 && size(chosenshocks,2)<size(shocks,1)  % menu disappears if the button "continue" is pushed or once all available shocks have been chosen.
                shock = eval(shockmenu);
                if shock == 1
                    disp('all shocks');
                    chosenshocks=[1:size(shocks,1)];
                elseif shock<(size(shocks,1)+2) && isempty(find(chosenshocks==(shock-1), 1))  % the first shows that the "continue" button has not been pushed, the second checks for double chosen models
                    shock = shock-1; % Neccessary as the "all shocks" button is set first. Therefore the button for shock one produces "shock=2". --> Substract 1
                    disp(deblank(shocks(shock,:)));
                    chosenshocks=[chosenshocks; shock];
                end
            end
            disp(' ');
            shocks = shocks(chosenshocks,:);
            base.innos = shocks; % put all chosen shocks in the choice vector for the IRFs
            base.namesshocks1 = shocks; % put the right shock names for correct plots
            
          
            % replace the names for interest_ and fiscal_ by mon.pol.shock and fiscal pol. shock (this is only for the legend of the graphs)
            if size(base.namesshocks1,2)<size(base.namesshocks,2)  % check if long names for monetary and fiscal shock fit in the names vector
                nblanks=size(base.namesshocks,2)-size(base.namesshocks1,2);  % adjust vector length
                for j= 1:size(base.namesshocks1,1)
                   blankvector(j,:) = blanks(nblanks);
                end
                base.namesshocks1 = [base.namesshocks1 blankvector];
            end
            if loc(base.namesshocks1,'interest_')~=0
                base.namesshocks1(loc(base.namesshocks1,'interest_'),:)=char('Mon. Pol. Shock      ');  %put nice names for 'interest_' and 'fiscal_'
            end
            if loc(base.namesshocks1,'fiscal_')~=0
                base.namesshocks1(loc(base.namesshocks1,'fiscal_'),:)  =char('Fiscal Pol. Shock    ');
            end
            base.namesshocks = base.namesshocks1; clear base.namesshocks1;
            
            %Choose if you want to plot all variables or interest rate, inflation, outputgap and
            %output
            base.option(4) = menu('Do you want to plot interest rate, inflation, outputgap and output or all variables?', 'Selected variables','All variables');
            disp(' ');
            if base.option(4) ==1
                disp('You decided to plot selected variables.');
                disp(' ');
            else disp('You decided to plot all variables.');
                disp(' ');
            end
        else
            base.option(4)=1;
        end
        if base.option(3)==1 %Several innovations are shocked contemporaneously
            cd('..');
            SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord)=M_.Sigma_e+1e-14*eye(M_.exo_nbr);
            cs = zeros(size(SS,1),1); % this line produces in the end IRF that are independent of any covariance structure. Shocks are one unit shocks. (as in Küster, Wieland)
            for p=1:size(base.innos,1)
                ii=loc(M_.exo_names(inv_lgx_orig_ord_,:),base.innos(p,:)); %Position of the shock
                cd(base.setpath(base.models(base.epsilon),:));
                base.pos_shock(p,base.models(base.epsilon))=ii;
                if ii==0
                disp(['No ' deblank(num2str(base.namesshocks(p,:))) ' is available for Model: ' num2str(base.names(base.models(base.epsilon),:))]);
                end;
                if base.variabledim(base.models(base.epsilon)) ==1
                    cs(ii,1)=1;
                end
                if base.variabledim(base.models(base.epsilon)) ==2
                    cs(ii,1)=1/100; % in case that models are written in percent/100 terms, shocks are 0.01 shocks
                end
            end;
            %Compute the IRFs
            R=irf(oo_.dr,cs(M_.exo_names_orig_ord,1), options_.irf, options_.drop, options_.replic, options_.order);
            base.IRF.(num2str(deblank(base.names(base.models(base.epsilon),:))))(:,:,1) = [zeros(size(R,1),1),R];
            base.IRFendo_names.(num2str(deblank(base.names(base.models(base.epsilon),:))))(:,:)=M_.endo_names;
        else
           for p=1:size(base.innos,1)
                cd('..');
                ii=loc(M_.exo_names(inv_lgx_orig_ord_,:),base.innos(p,:)); %Position of the shock
                cd(base.setpath(base.models(base.epsilon),:));
                base.pos_shock(p,base.models(base.epsilon))=ii;
                if ii==0
                    disp(['No ' deblank(num2str(base.namesshocks(p,:))) ' is available for Model: ' num2str(base.names(base.models(base.epsilon),:))]);
    %                 base.IRF.(num2str(deblank(base.names(base.models(base.epsilon),:))))(:,:,p) = [];
    %                 base.IRFlgy_.(num2str(deblank(base.names(base.models(base.epsilon),:))))(:,:) = [];
                else
                    % Computing the IRFs
                    SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord)=M_.Sigma_e+1e-14*eye(M_.exo_nbr);
                    %cs = transpose(chol(SS)); % this line should be taken later as an option to produce IRFs that are dependent on the covariances where the shock is one standard deviation (as in Dynare)
                    if base.variabledim(base.models(base.epsilon)) ==1
                        cs = eye(size(SS,1)); % this line produces in the end IRF that are independent of any covariance structure. Shocks are one unit shocks. (as in Küster, Wieland)
                    end
                    if base.variabledim(base.models(base.epsilon)) ==2
                        cs = eye(size(SS,1))*(1/100); % in case that models are written in percent/100 terms, shocks are 0.01 shocks
                    end
                    R=irf(oo_.dr,cs(M_.exo_names_orig_ord,ii), options_.irf, options_.drop, options_.replic, options_.order);
                    base.IRF.(num2str(deblank(base.names(base.models(base.epsilon),:))))(:,:,p) = [zeros(size(R,1),1),R];
                    base.IRFendo_names.(num2str(deblank(base.names(base.models(base.epsilon),:))))(:,:)=M_.endo_names;
                end;
            end;
         end;
      else
          base.IRF.(num2str(deblank(base.names(base.models(base.epsilon),:)))) = [];
          base.IRFendo_names.(num2str(deblank(base.names(base.models(base.epsilon),:))))= [];
      end
end
end
      
      
      
      
      

% if ~options_.noprint
%     disp(' ')
%     disp('MODEL SUMMARY')
%     disp(' ')
%     disp(['  Number of variables:         ' int2str(M_.endo_nbr)])
%     disp(['  Number of stochastic shocks: ' int2str(M_.exo_nbr)])
%     disp(['  Number of state variables:   ' ...
%           int2str(length(find(oo_.dr.kstate(:,2) <= M_.maximum_lag+1)))])
%     disp(['  Number of jumpers:           ' ...
%           int2str(length(find(oo_.dr.kstate(:,2) == M_.maximum_lag+2)))])
%     disp(['  Number of static variables:  ' int2str(oo_.dr.nstatic)])
%     my_title='MATRIX OF COVARIANCE OF EXOGENOUS SHOCKS';
%     labels = deblank(M_.exo_names);
%     headers = strvcat('Variables',labels);
%     lh = size(labels,2)+2;
%     dyntable(my_title,headers,labels,M_.Sigma_e,lh,10,6);
%     disp(' ')
%     if options_.order <= 2
%         disp_dr(oo_.dr,options_.order,var_list);
%     end
% end

% if options_.periods == 0 && options_.nomoments == 0
%    disp_th_moments(oo_.dr,var_list); 
% elseif options_.periods ~= 0
%     if options_.periods < options_.drop
%         disp(['STOCH_SIMUL error: The horizon of simulation is shorter' ...
%               ' than the number of observations to be DROPed'])
%         options_ =options_old;
%         return
%     end
%     oo_.endo_simul = simult(repmat(oo_.dr.ys,1,M_.maximum_lag),oo_.dr);
%     dyn2vec;
%     if options_.nomoments == 0
%         disp_moments(oo_.endo_simul,var_list);
%     end
% end



% if options_.irf 
%     if size(var_list,1) == 0
%         var_list = M_.endo_names(1:M_.orig_endo_nbr, :);
%         if TeX
%             var_listTeX = M_.endo_names_tex(1:M_.orig_endo_nbr, :);
%         end
%     end
% 
%     n = size(var_list,1);
%     ivar=zeros(n,1);
%     if TeX
%         var_listTeX = [];
%     end
%     for i=1:n
%         i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
%         if isempty(i_tmp)
%             error (['One of the specified variables does not exist']) ;
%         else
%             ivar(i) = i_tmp;
%             if TeX
%                 var_listTeX = strvcat(var_listTeX,deblank(M_.endo_names_tex(i_tmp,:)));
%             end
%         end
%     end
% 
%     if TeX
%         fidTeX = fopen([M_.fname '_IRF.TeX'],'w');
%         fprintf(fidTeX,'%% TeX eps-loader file generated by stoch_simul.m (Dynare).\n');
%         fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
%         fprintf(fidTeX,' \n');
%     end
%     olditer = iter_;% Est-ce vraiment utile ? Il y a la même ligne dans irf... 
%     SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord)=M_.Sigma_e+1e-14*eye(M_.exo_nbr);
%     cs = transpose(chol(SS));
%     tit(M_.exo_names_orig_ord,:) = M_.exo_names;
%     if TeX
%         titTeX(M_.exo_names_orig_ord,:) = M_.exo_names_tex;
%     end
%     for i=1:M_.exo_nbr
%         if SS(i,i) > 1e-13
%             y=irf(oo_.dr,cs(M_.exo_names_orig_ord,i), options_.irf, options_.drop, ...
%                   options_.replic, options_.order);
%             if options_.relative_irf
%                 y = 100*y/cs(i,i); 
%             end
%             irfs   = [];
%             mylist = [];
%             if TeX
%                 mylistTeX = [];
%             end
%             for j = 1:n
%                 assignin('base',[deblank(M_.endo_names(ivar(j),:)) '_' deblank(M_.exo_names(i,:))],...
%                          y(ivar(j),:)');
%                 eval(['oo_.irfs.' deblank(M_.endo_names(ivar(j),:)) '_' ...
%                       deblank(M_.exo_names(i,:)) ' = y(ivar(j),:);']); 
%                 if max(y(ivar(j),:)) - min(y(ivar(j),:)) > 1e-10
%                     irfs  = cat(1,irfs,y(ivar(j),:));
%                     mylist = strvcat(mylist,deblank(var_list(j,:)));
%                     if TeX
%                         mylistTeX = strvcat(mylistTeX,deblank(var_listTeX(j,:)));
%                     end
%                 end
%             end
%             if options_.nograph == 0
%                 number_of_plots_to_draw = size(irfs,1);
%                 [nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
%                 if nbplt == 0
%                 elseif nbplt == 1
%                     if options_.relative_irf
%                         hh = figure('Name',['Relative response to' ...
%                                             ' orthogonalized shock to ' tit(i,:)]);
%                     else
%                         hh = figure('Name',['Orthogonalized shock to' ...
%                                             ' ' tit(i,:)]);
%                     end
%                     for j = 1:number_of_plots_to_draw
%                         subplot(nr,nc,j);
%                         plot(1:options_.irf,transpose(irfs(j,:)),'-k','linewidth',1);
%                         hold on
%                         plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
%                         hold off
%                         xlim([1 options_.irf]);
%                         title(deblank(mylist(j,:)),'Interpreter','none');
%                     end
%                     eval(['print -depsc2 ' M_.fname '_IRF_' deblank(tit(i,:)) '.eps']);
%                     if ~exist('OCTAVE_VERSION')
%                         eval(['print -dpdf ' M_.fname  '_IRF_' deblank(tit(i,:))]);
%                         saveas(hh,[M_.fname  '_IRF_' deblank(tit(i,:)) '.fig']);
%                     end
%                     if TeX
%                         fprintf(fidTeX,'\\begin{figure}[H]\n');
%                         for j = 1:number_of_plots_to_draw
%                             fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{$%s$}\n'],deblank(mylist(j,:)),deblank(mylistTeX(j,:)));
%                         end
%                         fprintf(fidTeX,'\\centering \n');
%                         fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_IRF_%s}\n',M_.fname,deblank(tit(i,:)));
%                         fprintf(fidTeX,'\\caption{Impulse response functions (orthogonalized shock to $%s$).}',titTeX(i,:));
%                         fprintf(fidTeX,'\\label{Fig:IRF:%s}\n',deblank(tit(i,:)));
%                         fprintf(fidTeX,'\\end{figure}\n');
%                         fprintf(fidTeX,' \n');
%                     end
%                     %   close(hh)
%                 else
%                     for fig = 1:nbplt-1
%                         if options_.relative_irf == 1
%                             hh = figure('Name',['Relative response to orthogonalized shock' ...
%                                                 ' to ' tit(i,:) ' figure ' int2str(fig)]);
%                         else
%                             hh = figure('Name',['Orthogonalized shock to ' tit(i,:) ...
%                                                 ' figure ' int2str(fig)]);
%                         end
%                         for plt = 1:nstar
%                             subplot(nr,nc,plt);
%                             plot(1:options_.irf,transpose(irfs((fig-1)*nstar+plt,:)),'-k','linewidth',1);
%                             hold on
%                             plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
%                             hold off
%                             xlim([1 options_.irf]);
%                             title(deblank(mylist((fig-1)*nstar+plt,:)),'Interpreter','none');
%                         end
%                         eval(['print -depsc2 ' M_.fname '_IRF_' deblank(tit(i,:)) int2str(fig) '.eps']);
%                         if ~exist('OCTAVE_VERSION')
%                             eval(['print -dpdf ' M_.fname  '_IRF_' deblank(tit(i,:)) int2str(fig)]);
%                             saveas(hh,[M_.fname  '_IRF_' deblank(tit(i,:)) int2str(fig) '.fig']);
%                         end
%                         if TeX
%                             fprintf(fidTeX,'\\begin{figure}[H]\n');
%                             for j = 1:nstar
%                                 fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{$%s$}\n'],deblank(mylist((fig-1)*nstar+j,:)),deblank(mylistTeX((fig-1)*nstar+j,:)));
%                             end
%                             fprintf(fidTeX,'\\centering \n');
%                             fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_IRF_%s%s}\n',M_.fname,deblank(tit(i,:)),int2str(fig));
%                             if options_.relative_irf
%                                 fprintf(fidTeX,['\\caption{Relative impulse response' ...
%                                                 ' functions (orthogonalized shock to $%s$).}'],deblank(titTeX(i,:)));
%                             else
%                                 fprintf(fidTeX,['\\caption{Impulse response functions' ...
%                                                 ' (orthogonalized shock to $%s$).}'],deblank(titTeX(i,:)));
%                             end
%                             fprintf(fidTeX,'\\label{Fig:BayesianIRF:%s:%s}\n',deblank(tit(i,:)),int2str(fig));
%                             fprintf(fidTeX,'\\end{figure}\n');
%                             fprintf(fidTeX,' \n');
%                         end
%                         %                                       close(hh);
%                     end
%                     hh = figure('Name',['Orthogonalized shock to ' tit(i,:) ' figure ' int2str(nbplt) '.']);
%                     m = 0; 
%                     for plt = 1:number_of_plots_to_draw-(nbplt-1)*nstar;
%                         m = m+1;
%                         subplot(lr,lc,m);
%                         plot(1:options_.irf,transpose(irfs((nbplt-1)*nstar+plt,:)),'-k','linewidth',1);
%                         hold on
%                         plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
%                         hold off
%                         xlim([1 options_.irf]);
%                         title(deblank(mylist((nbplt-1)*nstar+plt,:)),'Interpreter','none');
%                     end
%                     eval(['print -depsc2 ' M_.fname '_IRF_' deblank(tit(i,:)) int2str(nbplt) '.eps']);
%                     if ~exist('OCTAVE_VERSION')
%                         eval(['print -dpdf ' M_.fname  '_IRF_' deblank(tit(i,:)) int2str(nbplt)]);
%                         saveas(hh,[M_.fname  '_IRF_' deblank(tit(i,:)) int2str(nbplt) '.fig']);
%                     end
%                     if TeX
%                         fprintf(fidTeX,'\\begin{figure}[H]\n');
%                         for j = 1:m
%                             fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{$%s$}\n'],deblank(mylist((nbplt-1)*nstar+j,:)),deblank(mylistTeX((nbplt-1)*nstar+j,:)));
%                         end
%                         fprintf(fidTeX,'\\centering \n');
%                         fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_IRF_%s%s}\n',M_.fname,deblank(tit(i,:)),int2str(nbplt));
%                         if options_.relative_irf
%                             fprintf(fidTeX,['\\caption{Relative impulse response functions' ...
%                                             ' (orthogonalized shock to $%s$).}'],deblank(titTeX(i,:)));
%                         else
%                             fprintf(fidTeX,['\\caption{Impulse response functions' ...
%                                             ' (orthogonalized shock to $%s$).}'],deblank(titTeX(i,:)));
%                         end
%                         fprintf(fidTeX,'\\label{Fig:IRF:%s:%s}\n',deblank(tit(i,:)),int2str(nbplt));
%                         fprintf(fidTeX,'\\end{figure}\n');
%                         fprintf(fidTeX,' \n');
%                     end
%                     %                           close(hh);
%                 end
%             end
%         end
%         iter_ = olditer;
%         if TeX
%             fprintf(fidTeX,' \n');
%             fprintf(fidTeX,'%% End Of TeX file. \n');
%             fclose(fidTeX);
%         end
%     end
% end
% 
% % if options_.SpectralDensity == 1
% %     [omega,f] = UnivariateSpectralDensity(oo_.dr,var_list);
% % end
% 
% 
% options_ = options_old;
