function [xcorFunc,xcorLags]=coralTemplateXcorr(template,dataStream)
%code written by Justin Sweet for using a template to pull out signals from
%a longer time series (used to be named LFExcorr.m)

%INPUTS
%template: a coralstructure of the template event (1)
%dataStream: a coralstructure of the event to compare the template against

%OUTPUTS
%xcorFunc: the normalized cross correlation between template and dataStream at each point in time 
%xcorLags: the time lags, in seconds, from the start of the dataStream

%modifed from original version by Kate Allstadt 14-Mar-2013 to output event times in seconds from
%start of dataStream instead of from start of template and to use xcorr 
%instead of conv for speed

%%%%%BEGIN CC CODE%%%%
y=dataStream.data;
x=template.data;

M=length(y);
N=length(x);

%OLD METHOD
%xy=conv(x(N:-1:1),y); %this gives xcorr, but backwards...why does it work?
%xy=xy(N:end-N+1); %xy should have M elements

xy=xcorr(y,x); %xcorr is faster than conv, if length of y is >2^20 and you have old version of matlab (older than 2011), need to download updated version of fftfilt.m from here http://www.mathworks.com/support/bugreports/640328
xy=xy(M:end-N+1); %make it the right length

x2=sum(x.*x);

csy2=[0; cumsum(y.*y)];

I=1:M-N+1;
y2=csy2(N+I) - csy2(I);

xynorm=xy./sqrt(x2*y2);

indx=1:length(xy);
time=indx*dataStream.recSampInt;
%abs_time=timeadd(dataStream.recStartTime,time); %convert relative time to absolute time
%abs_time2=datenum(abs_time(:,:)');

xcorFunc=xynorm;
xcorLags=time;

return