function K=coralKurtosis(D,trange)
%compute the kurtosis of signals in a coral structure, corrected for bias
%a normal distribution has a kurtosis close to 3, those with more outliers
%are much larger than 3, with fewer outliers, less than 3.

%INPUTS
%D coral structures of signals (can be vector of multiple signals)
%trange = time range of window to compute K over, in seconds. ex [2 4] computes K for 2 to 4 seconds
%Insert 0 to compute K for entire signal

%OUTPUTS
%K = Kurtosis of the signal in each coral structure input

K=zeros(size(D));

for i=1:size(D,1)
    for j=1:size(D,2)
        tvec=0:D(i,j).recSampInt:D(i,j).recSampInt*(D(i,j).recNumData-1);
        dat=D(i,j).data;
        
        if sum(trange)~=0 %if a time window was specified
            dat=dat(tvec>=trange(1) & tvec<=trange(2));
        end
        K(i,j)=kurtosis(dat,0);
    end
end
end