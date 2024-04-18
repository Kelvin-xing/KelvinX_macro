Tinit = 100;       % periods to initialize distribution of assets
T     = 50;        % periods to simulate allocations
Nf    = 50000;     % # workers

randn('state', 100);
rand('state',  100);

asim     = zeros(Nf,T);
xsim     = zeros(Nf,T); 
eesim    = zeros(Nf,T);

Pcum      = Pw*triu(ones(size(Pw)));   % cumulative ergodic distributon

unif = nodeunif(Nf,eps^(1/2),1-eps^(1/2));
unif = unif(randperm(Nf));

ppi       = [0; cumsum(Pwerg)];
[n,bin]   = histc(unif,ppi);

esim      = zeros(Nf,T+Tinit);
esim(:,1) = bin;

for t = 2 : Tinit + T
   
  unif = unif(randperm(Nf));
  Pnew = Pw(esim(:,t-1),:);
  Pcum = [zeros(Nf,1), cumsum(Pnew,2)];
  esim(:,t) = ((repmat(unif,1,2)<Pcum(:,2:end)).*(repmat(unif,1,2)>Pcum(:,1:end-1)))*(1:2)';
  
end


for t = 1 : Tinit+T

    if t==1
         
        e      = esim(:,t);
        state  = [ones(Nf,1), e]; 
        a      = state(:,1);
        x      = funeval(cxw,fspacew,state);  x = min(max(x, smin(1)),smax(1));
                
    else

        e      = esim(:,t);                  % "exit" is random, so ok to use incumbent's e as draw from ergodic
        state  = [x, e];
        a      = state(:,1); 
        x      = funeval(cxw,fspacew,state);  x = min(max(x, smin(1)),smax(1));
        
        
          if t > Tinit
             
              asim(:,t-Tinit)    = a;
              xsim(:,t-Tinit)    = x;
              eesim(:,t-Tinit)   = e;
          end
    end
end

a    = asim';
x    = xsim';
e    = eesim';

fprintf('Asset Supplied Workers       = %9.3f \n',   mean(a(:)));
