function shade(start,finish,colorstr)

%Holding figure
hold on

%Creating y-coordinates for shading (same for every recession)
current_ax = axis;
y = [current_ax(3) current_ax(4) current_ax(4) current_ax(3)];

%Creating x-coordinates for shading (different for every recession)
for i = 1:length(start),
    x = [start(i) start(i) finish(i) finish(i)];
    fill(x,y,colorstr)
end

%Visibility
set(findobj(gca,'Type','Patch'),'FaceAlpha',0.5)
set(findobj(gca,'Type','Patch'),'EdgeColor','none')
set(gca,'Layer','top')

end