%enter starting month and year
startYear=2013;
startMonth=2; %number
endYear=2013;
endMonth=3;

numev=2; %minimum number of events to keep a set, any sets with fewer than numev will be thrown out
minxcorr=0.7; %minimum cross correlation with stack, events with a lower xcorr with the stack will be thrown out (note some events can have lower xcorrs with the stack than 0.7 because they were combined with an entire month of data)

%make vector of starts of each month in datenum format
monthvec=monthvector(startYear,startMonth,endYear,endMonth);

for i=1:length(monthvec)-1
    mon=monthvec(i);
%change back to home directory if not already there
cd Users\geo-user\Documents\KateA\MultCodes\HTWmultsearch
   
%change to the new directory
cd(datestr(mon,'yymm'))
    fprintf([datestr(mon,'yymm'),'\n'])
%fix multtempl
name=dir('*MULTTEMPL.EHZ*');
if ~isempty(name)
load(name.name)
fixMulttempl(multtempl,numev,minxcorr)
end

end