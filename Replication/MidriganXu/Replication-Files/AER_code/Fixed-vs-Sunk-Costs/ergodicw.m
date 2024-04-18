T = 1000;  % periods we wait for weights to converge

N = 5001;

agrid = nodeunif(N, smin(1), smax(1));

edges = [agrid(1); (agrid(2:end-1) + agrid(1:end-2))/2; agrid(end)];

state = gridmake(agrid, (1:2)');

Nw = size(state, 1); 
nw = 1/Nw*ones(Nw,1);

% Decision rules

x   = funeval(cxw, fspacew, state); x = min(max(x, smin(1)),smax(1));  % don't allow outside bounds 
       
% Sort savings into bins

[~, xbin] = histc(x, edges);     % bin by asset
xbin      = [xbin; (1:1:N)'];    % include all

% Iterate over measures

for t = 1 : T
    
nwold = nw;

Pwprime     = Pw(state(:,2),:);
Pwprime     = bsxfun(@times, nw, Pwprime); 

Pwprime = [Pwprime; zeros(N,2)];    % make sure all guys included

Pwnew = zeros(N, 2); 

for i = 1:2
    
   Pwnew(:,i) = accumarray(xbin,Pwprime(:,i));    % add all guys in a given bin
    
end

nw = reshape(Pwnew, N*2, 1); 

%fprintf('%4i %6.2e \n',[t, norm(nw-nwold)]);

if norm(nw-nwold) < 1e-10, break, end

end


Aw = nw'*state(:,1); 

fprintf('\n');

fprintf('Asset Supplied Workers       = %9.3f \n',   Aw);