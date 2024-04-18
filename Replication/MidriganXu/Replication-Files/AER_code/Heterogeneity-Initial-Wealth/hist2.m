% Modified version of something written by  Laszlo Balkay
% this one allows weights
% also at this point keep second vector the bigger of the two

function histmat  = hist2(x, y, xedges, yedges,w)

[xn, xbin] = histc(x,xedges);
[yn, ybin] = histc(y,yedges);

xnbin = length(xedges);
ynbin = length(yedges);

if xnbin >= ynbin
    xy = ybin*(xnbin) + xbin;
      indexshift =  xnbin; 
else
    xy = xbin*(ynbin) + ybin;
      indexshift =  ynbin; 
end

xyuni = unique(xy);
[hstres,bin] = histc(xy,xyuni);
idx=cumsum(hstres);


junk=sortrows([w,bin],2);
w=junk(:,1);

for i=1:length(xyuni)
    if i>1
    hstres(i)=sum(w(idx(i-1)+1:idx(i)));    
    else
    hstres(i)=sum(w(1:idx(i))); 
    end
    
end

histmat = zeros(ynbin,xnbin);

histmat(xyuni-indexshift) = hstres;
histmat = histmat';