function [pTime,sL] = pickerRatio(D, sWindow, lWindow, rCutoff, nInterval, test)
% USAGE:
%       [pTime] = pickerRatio(D, sWindow, lWindow, rCutoff, nInterval, test)
%
% pickerRatio.m: function to use an sta/lta ratio on a hilbert transformed
%                waveform to determine picks
%
% INPUTS:
%   D: coral structure for picking (one channel at a time)
%   sWindow: short term window for picker in seconds (0.8 s is suggested)
%   lWindow: long term window for picker in seconds (7 s is suggested)
%   rCutoff: ratio cutoff for declaring event. Positive numbers are cutoffs
%         with respect to the actual ratio.  Negative numbers are cutoffs
%         with respect to standard deviations from the mean of the noise.
%         (2.5 or -2 is suggested)
%         ex.
%           2 means that the ratio lWindow/sWindow is two
%          -2 fits a normal distribution to the ratios of the time series,
%             (assumes that noise is normally distributed) and picks the
%             cutoff to be where the ratios are 2 standard deviations above
%             the mean
%   nInterval: Number of seconds between events (10 s is suggested)
%   test: Boolean, plot out the triggers? 1 = yes, 0 = no
%
% OUTPUT:
%   pTime: Cell array of pick times in coral format
%
% EDIT:
% added sL as an output (Joshua D Carmichael, 30.Oct.2011)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test = 1;

x = abs((hilbert(D.data)));
%x = (D.data).^(2);
S = round((sWindow/D.recSampInt));
L = round((lWindow/D.recSampInt));
% Calculate ratio (from Creager)
sL = zeros(size(x));  % initialize output vector
 if length(x) > S+L;
   X       = cumsum(x);       % cumulative sum 
   k       = L+2:length(X)-S+1; % index vector
   sL(k)   = [X(k+S-1)-X(k-1)] ./ [X(k-1)-X(k-L-1)] * L/S;
   sL(L+1) = (X(L+S)-X(L))/X(L) *L/S; % special case for first point in time series
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate stats on ratio
muhat=mean(sL);
sigmahat=std(sL);

% ! CHANGE TO WES CODE !
% [muhat, sigmahat] = normfit(sL);
% ! CHANGE TO WES CODE !

% Find ratio
if rCutoff < 1
    rCut = muhat + abs(rCutoff)*sigmahat;
else
    rCut = rCutoff;
end

% Find events above ratio
r = find(sL >= rCut);

% Get triggers
eTrace = zeros(length(sL),1);

eTrace(r) = 1;
for i = 1 : length(r)
        % Cut out appropriate window
        if (r(i)+nInterval*1/D.recSampInt < length(eTrace))
%             tSamp = eTrace(r(i):r(i)+nInterval*1/D.recSampInt);
%             tSampV = sL(r(i):r(i)+nInterval*1/D.recSampInt);
            tSamp = eTrace(r(i):round(r(i)+nInterval*1/D.recSampInt));
            tSampV = sL(r(i):round(r(i)+nInterval*1/D.recSampInt));
        else
            tSamp = eTrace(r(i):end);
            tSampV = sL(r(i):end);
        end
        k = find(tSamp == 1);
        % if multiple triggers
        if (length(k) > 1)
            kk = find(tSampV~=max(tSampV));
            tSamp(kk) = 0;
            ind = r(i)+kk;
            eTrace(ind) = 0;
        end
        if (r(i)+nInterval*1/D.recSampInt < length(eTrace))

            eTrace(r(i):round(r(i)+nInterval*1/D.recSampInt)) = tSamp;
        else
            if(length(eTrace(r(i):end))==length(tSamp));
                eTrace(r(i):end) = tSamp;
            end;
            if(length(eTrace(r(i):end))~=length(tSamp));
                minL=min(length(eTrace(r(i):end)),length(tSamp));
                eTrace(end-minL+1:end) = tSamp(1:minL);
            end;
        end
end

disp(sprintf('Number of events: %d', length(find(eTrace ==1))));
% Get absolute picktimes
rr = find(eTrace == 1);
pTime = cell(length(rr),1);
for i = 1 : length(rr)
    pTime{i} = timeadd(D.recStartTime, rr(i)*D.recSampInt);
end

if length(rr) < 1
    disp('No triggers');
    pTime{1} = 0;
    return;
end

% Testing plot
if test ==1
    figure; subplot(3,1,1);
    coralPlot(D);
    tvec = [0:length(eTrace)-1]*D.recSampInt;
    plot(tvec, (x/max(x) * 0.5), '-b');
    for i = 1 : length(rr)
        plot(tvec(rr(i)), 0, 'or');
    end
    xlim([0 tvec(end)]);
    ylabel('Counts');
    ylim([-0.5 0.5]);
    legend('Waveform', 'Hilbert', 'Trigger', 'Location', 'NorthEast');
    
    subplot(3,1,2);
    plot(tvec, sL, '-k');
    hold on;
    for i = 1 : length(rr)
        plot(tvec(rr(i)), rCut, 'or');
    end
    xlim([0 tvec(end)]);
    xlabel('Time, s');
    ylabel('STA/LTA');
    axis tight
    legend('STA/LTA', 'Trigger', 'Location', 'NorthEast');
    subplot(3,1,3);
    title('Histogram of sta/lta values');
    h = hist(sL, [0:0.1:max(sL)]);
    hist(sL, [0:0.1:max(sL)]);
    hold on;
%     nn = (1/(sigmahat*sqrt(2*pi)))*exp((-0.5*((muhat+sigmahat*randn(length(sL),1))-muhat).^2)./(muhat.^2));
    plot([rCut rCut], [0 max(h)], '-r', 'LineWidth', 3);
    xlabel('STA/LTA Ratio');
    ylabel('Occurences');
    text(rCut+0.2, mean([0 max(h)]), 'cutoff', 'HorizontalAlignment', 'center', 'Rotation', 90, 'Color', 'r');
    %text(rCut+0.2, mean([0 max(h)]));
    axis tight;
end

