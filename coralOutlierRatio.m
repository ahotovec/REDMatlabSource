function OP=coralOutlierRatio(D,trange)
%determine the ratio of points taht are technically outliers

%INPUTS
%D coral structures of signals (can be vector of multiple signals)
%trange = time range of window to compute OP over, in seconds. ex [2 4] computes OP for 2 to 4 seconds
%Insert 0 to compute OP for entire signal

%OUTPUTS
%OP = percentage of absolute value of data points that are outliers

OP=zeros(size(D,1),size(D,2),1);

for i=1:size(D,1)
    for j=1:size(D,2)
        tvec=0:D(i,j).recSampInt:D(i,j).recSampInt*(D(i,j).recNumData-1);
        dat=D(i,j).data;
        
        if sum(trange)~=0 %if a time window was specified
            dat=dat(tvec>=trange(1) & tvec<=trange(2));
        end
        dat=abs(dat);
        z=(dat-median(dat))./mad(dat);
        OP(i,j)=sum(z>4.45)./length(z);
    end
end
end