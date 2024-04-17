%These are notes on the use of code for computing optimal monetary
%policy that was prepared by Andy Levin (Federal Reserve Board). This code
%may be used freely. However, if you use it, please include an acknowledgement,
%'the code for this is taken from  Andrew Levin, Lopez-Salido, J.D., 2004.
%"Optimal Monetary Policy with Endogenous Capital Accumulation", manuscript,
%Federal Reserve Board, and Andrew Levin, Onatski, A., Williams, J., 
%Williams, N., 2005. "Monetary Policy under Uncertainty in Microfounded 
%Macroeconometric Models." In: NBER Macroeconomics Annual 2005, 
%Gertler, M., Rogoff, K., eds. Cambridge, MA: MIT Press.
%The Levin code may be downloaded from
%http://faculty.wcas.northwestern.edu/~lchrist/d16/d1606/Get_Ramsey.zip
%Dynare (another software package that can be freely used, with
%attribution) may be downloaded from
%http://faculty.wcas.northwestern.edu/~lchrist/d16/d1606/dyn_mat_v3_064.zip


%Instructions for running the optimal policy code
%
%Make sure that dynare and Andy Levin's code are in your MATLAB path. For 
%example, in my computer I enter this at the MATLAB prompt:
%
%       addpath('c:\projects\solve\dynare')
%       addpath('c:\projects\solve\levin')
%
%In case dynare and Levin's code are in the same path where your programs
%reside, then there is no need to bother with the addpath commands.
%
%Suppose you place the basic monetary model with the exogenous policy rule in model.mod
%To study the properties of the model with the exogenous policy rule, simply type 
%
%   dynare model
%
%at the MATLAB command prompt.
%
%To get model.mod ready for optimal policy, you have to do several things.
%You must have two var statements, one with the
%variables to be optimized in the Ramsey problem, and one with the rest.
%(One exception is that Util and Welf must be included in the endogenous variables list and Andy
%Levin's code, sensibly, does not optimize with respect to these variables.)
%
%The var statement with the variables to be optimized in the 
%Ramsey problem must be preceded by the comment:
%
%         // Endogenous variables
%
%The var statement with the other 'endogenous' variables must be preceded
%by
%         // Exogenous variables
%
%In practice, you would include the name of your economic shocks in the
%latter var statement.
%
%You must have, among the equations of the model, an equation that 
%specifies an exogenous policy, preceded by the statement:
%
%        // Monetary Policy Rule
%
%You must apply Levin's translation code to model.mod.
%This requires making sure the model.mod file is set up correctly. You will
%need to have two additional equations. One is exactly the following
%
%     Welf = Util + nbeta*Welf(+1);
%
%These exact names MUST be used. The variables, Welf and Util must be included
%in the list of endogenous variables. The parameter, nbeta, must be set and
%included in the list of parameters. nbeta is the Ramsey problem discount
%rate which you might just want to equate to the private agents' 
%discount rate. 
%The second equation you must include is one that defines Util
%
%    Util = .....   ;
%
%where the definition of your period utility should be placed where the dots
%are. Now you are ready to apply Levin's translation code.
%At the Matlab prompt, type
%
%           infilename  = 'model';
%           outfilename = 'model_OUT';
%           get_ramsey
%
%this writes a new dynare file, model_OUT.mod, in which the monetary policy
%rule is replaced by the first order conditions associated with the 
%Lagrangian representation of the Ramsey problem.
%
%Once the translation is complete, things need to be done to the new
%program, model_OUT.mod, so that it can be run by Dynare.
%
%1. three var statements must be merged - the ones for endogenous and
%exogenous variables, and the one for the multipliers.
%
%2. you need to do a little (not much!) work to compute the steady state of 
%the Ramsey problem. The call to the steady state program in model.mod 
%computes the steady state in the version of the model with the exogenous 
%monetary policy rule. In the optimal policy problem, the policy rule has 
%been dropped and replaced by the Ramsey first order conditions. With two 
%exceptions, get_ramsey computes the first order conditions with respect to 
%all the variables you listed in the var statement after the 
%// Endogenous variables command in model_OUT.mod. The two exceptions are 
%Util and Welf. Getramsey does not differentiate with respect 
%to these variables, because they do not represent binding restrictions 
%on the level of utility that can be achieved in the Ramsey problem (i.e.,
%we know that their multipliers would be exactly zero).
%Suppose the number of variables with which you differentiate is N (N=4 in 
%modelans_OUT.mod, which is what is output by getramsey when getramsey is 
%applied to modelans.mod, my answer to homework #9). 
%So, application of getramsey to your model_OUT.mod program replaces the 
%monetary policy rule by N first order conditions. In addition, there are 
%N-1 multipliers. The steady state value of these multipliers must be 
%determined, in addition to the steady state value of the 
%N endogenous variables. So, the steady state problem is different in the
%Ramsey problem than it is in the model with exogenous policy.
%The good news is that getramsey takes care of the most difficult part of 
%this steady state problem.
%
%The strategy for finding the steady state of the Ramsey problem is as
%follows. You specify a conjectured value for the 
%steady state inflation rate, say x. Then, you compute the values of the N-1 
%endogenous variables that solve the steady state
%equilibrium conditions. For this step, you can simply use the steady state 
%program that you use for the case when monetary policy was exogenous). 
%Then, the steady state multipliers are found which
%solve the N first order conditions of the Ramsey problem. In general, you 
%will not be able to solve this system, because you only have N-1 multipliers. 
%The way the problem is solved, is you adjust the value of x until the N-1 
%multipliers you compute solve all N Ramsey first order conditions.
%Then you are done. In our example, life is simple. The optimal inflation
%rate in the Ramsey problem is unity. (Many problems have this form, where you
%can easily conjecture the Ramsey steady state inflation rate. To see
%a series of these problems, have a look at the examples in 
%  http://faculty.wcas.northwestern.edu/~lchrist/d16/d1606/ramsey.htm
%A research paper which uses getramsey, and where it is easy to conjecture
%the steady state can be found in 
%http://www.faculty.econ.northwestern.edu/faculty/christiano/research/ECB/m
%irage.htm).
%
%The only painful part of the Ramsey steady state computation is the 
%computation of the N-1 multipliers conditional on the steady state values 
%of the endogenous variables. Fortunately, getramsey.m writes a program that
%does this for you! To see how this works, have a look at the program, ssnew.
%
%The program which computes the steady state of the N-1 equilibrium
%conditions is called sstate.m. The program which computes the steady state
%of the multipliers is Levin's endogenously generated code, model_out_lmss.m
%
%make sure ssnew appears after the parameters statement with the multipliers
%
%
%3.%model_out_lmss.m must be adjusted in two ways:
%
%a. comment out the global statements at the start of the code.
%
%b. it is a good idea to remove the following commands:
% 
% if errcheck > 1e-08, 
%   disp('Warning: steady states of lagrange multipliers cannot be accurately determined'); 
%   disp(['         errcheck = ',num2str(errcheck,12)]); 
% end; 
% 
%this check need not be put anywhere else because it appears at the end of ssnew. 
%leaving the check in model_out_lmss.m is a nuisance because 
%this program is evaluated repeatedly if you need to try out different values of x
%in the process of computing a steady
%state. Until the steady state is found, the check would generate a warning.
