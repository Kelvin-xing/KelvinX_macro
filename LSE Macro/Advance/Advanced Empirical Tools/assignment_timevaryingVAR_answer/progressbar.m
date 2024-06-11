function progressbar(fractiondone)

% You have to initialize the progressbar with "progressbar(0)". After that,
% you can update the progressbar with "progressbar(fractiondone)". Finally,
% you can close the progressbar with "progressbar(1)".
%
% Author: Joris de Wind
% Date:   July 31, 2009

%% Argument checking
%--------------------------------------------------------------------------

if ~exist('fractiondone','var')                          %Required argument
    return
end

%% Declaring persistent variables
%--------------------------------------------------------------------------

persistent progfig progpatch starttime lastupdate

%% Getting the old progfig handle or creating a new (empty) progfig handle
%--------------------------------------------------------------------------

try
    dummy = get(progfig,'UserData');                            %#ok<NASGU>
    if fractiondone == 0
        progfig = [];
    end
catch                                                            %#ok<CTCH>
    progfig = [];
end

%% Setting-up the new (empty) progfig handle
%--------------------------------------------------------------------------

if isempty(progfig) && fractiondone == 0
    
    progfig = figure(...
        'Units',            'normalized',...
        'Position',         [0.375 0.49 0.25 0.02],...
        'NumberTitle',      'off',...
        'Resize',           'off',...
        'MenuBar',          'none',...
        'BackingStore',     'off',...
        'CloseRequestFcn',  @closeBar);
    progaxes = axes(...
        'Position',         [0.02 0.15 0.96 0.70],...
        'XLim',             [0 1],...
        'YLim',             [0 1],...
        'Box',              'on',...
        'ytick',            [],...
        'xtick',            []);                                %#ok<NASGU>
    progpatch = patch(...
        'XData',            [0 0 0 0],...
        'YData',            [0 0 1 1],...
        'EraseMode',        'none',...
        'FaceColor',        [.1 1 .1]);
    
    lastupdate = clock - 1;                             %Just to initialize
    
    starttime = clock;
            
end

%% Clearing
%--------------------------------------------------------------------------

if isempty(progfig) && fractiondone ~= 0                %True if user closed the progressbar
    delete(progfig)
    clear progfig progpatch starttime lastupdate
    return    
end

%% Updating
%--------------------------------------------------------------------------

if etime(clock,lastupdate) < 0.01 && ~(fractiondone == 1)       %Don't update too often
    return
end

set(progpatch,'XData',[0 fractiondone fractiondone 0])          %Update the progpatch

if fractiondone == 0
    titlebarstr = ' 0%';
else
    runtime     = etime(clock,starttime);
    timeleft    = runtime/fractiondone - runtime;
    timeleftstr = sec2timestr(timeleft);
    titlebarstr = sprintf('%2d%%    %s remaining',floor(100*fractiondone),timeleftstr);
end
set(progfig,'Name',titlebarstr)                                 %Update the title

drawnow

if fractiondone == 1                                            %Close if finished
   delete(progfig)
   clear progfig progpatch starttime lastupdate
   return
end

lastupdate = clock;                                             %Record time of this update

%% Function sec2timestr 

function timestr = sec2timestr(sec)

% Convert seconds to other units
d = floor(sec/86400);               %Days
sec = sec - d*86400;
h = floor(sec/3600);                %Hours
sec = sec - h*3600;
m = floor(sec/60);                  %Minutes
sec = sec - m*60;
s = floor(sec);                     %Seconds

% Create time string
if d > 0
    if d > 9
        timestr = sprintf('%d day',d);
    else
        timestr = sprintf('%d day, %d hr',d,h);
    end
elseif h > 0
    if h > 9
        timestr = sprintf('%d hr',h);
    else
        timestr = sprintf('%d hr, %d min',h,m);
    end
elseif m > 0
    if m > 9
        timestr = sprintf('%d min',m);
    else
        timestr = sprintf('%d min, %d sec',m,s);
    end
else
    timestr = sprintf('%d sec',s);
end

end

%% Function closebar

function closeBar(src,evnt)              %#ok<INUSD> %Required input arguments, implicitly passed by Matlab
 
selection = questdlg('Stop?','Stop?','Yes','No','Yes');

switch selection
  case 'Yes'
    delete(gcf)                                      %Delete current handle     
  case 'No'
    return
end

end

end