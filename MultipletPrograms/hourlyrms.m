function [rms datevec]=hourlyrms(S)
%calculates hourly rms value
%input 
% S is a coral structure of data you want to calculate rms of

totalhrs=floor((S(1).recNumData*S(1).recSampInt)/3600);
%preallocate
rms=zeros(totalhrs,1);
datevec=datenum(S(1).recStartTime')+1/48:1/24:datenum(S(1).recStartTime')+(1/24)*(totalhrs-1)+1/48; %put datevec at center of time period of calculation

for i=1:totalhrs
    for j=1:size(S,2)
        samps=3600/S(j).recSampInt; %figure out number of samples in each hour-long chunk
        datachunk=S(j).data((i-1)*samps+1:i*samps);        
        rms(i,j)=sqrt(mean(datachunk.^2));
    
    end
end

end
