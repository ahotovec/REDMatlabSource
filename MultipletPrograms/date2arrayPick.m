function  [err,outfile,ind,S] = date2arrayPick(pickDate)
% This function extracts the array pick file matching most closely in the
% time the given input 6x1 date vector pickDate.
%
% USAGES
% [err,outfile,ind,S] = date2arrayPick(pickDate)
%
% INPUT
% pickDate: An 6x1 date vector input a row vector.
%
% OUTPUT
% err:      The time error in seconds between pickDate and pick date in:
%           S(1,ind).recStartTime. 
% outfile:  A string giving the array pick file, of the form
%           ARRAYPICKS*.mat
% ind:      The column of the MxN coral structure array stored in the
%           outfile that matches the pickDate.
% S:        The MxN coral structure whose field recStartTime at column ind
%           matches input pickDate most closely.
%-----------------------------------------------------------------------
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 27.May.2010
%
% Testing phase.
%-----------------------------------------------------------------------

boy     = pickDate;
boy(2:3)= 1;
boy(4:6)= 0;
day     = 1+floor(timediff(pickDate,boy)./(3600*24));
day     = sprintf('%.3i.%s',day,'mat');
year    = sprintf('%i',boy(1));

flist   = getFiles('ARRAYPICKS',year,day);

load(flist{1});

raDates = arraypicks.pTime;

[m,ind] = sort(abs(timediff( [pickDate,raDates] )),2,'ascend');

ind     = ind(2)-1;
err     = timediff( pickDate, raDates(:,ind) );
S       = arraypicks.coralstruc;
outfile = flist{1};

return;