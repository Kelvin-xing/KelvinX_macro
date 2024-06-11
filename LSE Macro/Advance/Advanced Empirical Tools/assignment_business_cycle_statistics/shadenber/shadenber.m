function shadenber

%Loading NBER recession periods
load nberdates

%Selecting which recessions to include
current_ax = axis;
index1     = find(finish > current_ax(1),1,'first');            %#ok<NODEF>
index2     = find(start  < current_ax(2),1,'last' );            %#ok<NODEF>

%Truncating if recession is only partly included
if start (index1) < current_ax(1), start (index1) = current_ax(1); end
if finish(index2) > current_ax(2), finish(index2) = current_ax(2); end

%Shading
colorstr = [159 182 205]/256;
shade(start(index1:index2),finish(index1:index2),colorstr)

end