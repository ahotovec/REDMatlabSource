function Fd=coralDominantFreq(D,trange)
%compute the dominant frequency using formula from Douma and Snieder 2006
%INPUTS
%D coral structures of signals (can be vector of multiple signals)
%trange = time range of window to compute Fd over, in seconds. ex [2 4] computes Fd for 2 to 4 seconds
%Insert 0 to compute Fd for entire signal

%OUTPUTS
%Fd = dominant frequency, in Hz

Fd=zeros(size(D,1),size(D,2),1);


for i=1:size(D,1)
    for j=1:size(D,2)
        tvec=0:D(i,j).recSampInt:D(i,j).recSampInt*(D(i,j).recNumData-1)';
        vel=D(i,j).data;
        acc=diff(D(i,j).data)./D(i,j).recSampInt;
        
        if sum(trange)~=0 %if a time window was specified
            tvec=tvec(tvec>=trange(1) & tvec<=trange(2));
            vel=vel(tvec>=trange(1) & tvec<=trange(2));
            acc=acc(tvec>=trange(1) & tvec<=trange(2));
            Fd(i,j)=(sqrt(trapz(tvec(2:end),acc(2:end).^2)./trapz(tvec(2:end),vel(2:end).^2)))./(2*pi);
        else
            Fd(i,j)=(sqrt(trapz(tvec(2:end),acc.^2)./trapz(tvec(2:end),vel(2:end).^2)))./(2*pi);
            
        end
    end
end

end