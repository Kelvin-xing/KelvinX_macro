
switch 'twodim'
    
    case 'onedim'

ecutoff = zeros(size(grid{1}));


for i = 1 : length(grid{1});
    
    state = gridmake(grid{1}(i), (1:1:k)');
    [~, vcm, vct] = valfunc4(ct,fspacet,state);
    
    if any(vcm>vct)
        
    ecutoff(i) = min(state(vcm>vct,2));
    
    else 
        
        ecutoff(i) = smax(2);
    end

end

    case 'onedime'

acutoff = zeros(size(grid{2}));


for i = 1 : length(grid{2});
    
    state = gridmake(nodeunif(1000,sminm(1)+kappa,smax(1)),grid{2}(i));
    [~, vcm, vct] = valfunc4(ct,fspacet,state);
    
    if any(vcm>vct)
        
    acutoff(i) = min(state(vcm>vct,1));
    
    else 
        
        acutoff(i) = smax(2);
    end

end


    case 'twodim'

state = gridmake(nodeunif(300,smin(1),1), (1:1:k)');

[v, vcm, vct] = valfunc4(ct,fspacet,state);   % evaluate worker's decision
figure(3)
scatter(state(vcm>=vct,1),state(vcm>=vct,2))
hold on
scatter(state(vcm<vct,1),state(vcm<vct,2), 'r')

end