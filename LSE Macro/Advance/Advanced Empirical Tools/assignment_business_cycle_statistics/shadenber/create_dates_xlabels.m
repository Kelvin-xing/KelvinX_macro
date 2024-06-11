% program to create the "weird" dates used in shade nber

Year_first = 1947;
Year_last  = 2013;

LL = (Year_last-Year_first);

dateX = [ Year_first 1; ...
          Year_first 2; ...
          Year_first 3; ...
          Year_first 4];

for i = 1:LL
    dateX = [ dateX; ...
             Year_first+i 1;
             Year_first+i 2;
             Year_first+i 3;
             Year_first+i 4];
end

dateX = datenum([dateX(:,1),dateX(:,2),ones(length(dateX),1),zeros(length(dateX),3)]);

save dateX