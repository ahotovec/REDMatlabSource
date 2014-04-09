function [tau,varargout] = coralSignalWidth(S)
% Computes the duration of the data in MxN input coral structure S, using
% common definition initial given (Gabor 1946) and also used in Seismology
% (i.e., Meier and Lee, Geophyics, Vol. 74, 2009).
% Required Fields: data, recSampInt, recNumData
%
%
% USAGES
% [tau]         = coralSignalWidth(S);
% [tau,phi]     = coralSignalWidth(S);
%
% INPUT
% S:    The MxN coral structure array with time series in field data.
%        
% OUTPUT
% tau:  An MxN array of effective signal width values in seconds
% phi:  An MxN array of effective frequency width values in hertz
%-----------------------------------------------------------------------
% Latest Edit: 29.April.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------
[M,N]    = size(S);

%Pad data so that columns of S have same data length
Sc       = mat2cell(S,M,ones(N,1));
S        = cellfun(@coralPad,Sc,'UniformOutput',false);
S        = cell2mat(S);

tau = zeros(M,N);
phi = tau;

%loop over columns (record sections)
for k = 1:N
    
    ndata   =[S(1,k).recNumData];
    %sample rate at all receivers
    srate   = [S(:,k).recSampInt];
    %logical test to see if srate is uniform. If isloop is non-zero, then
    %enter loop below. Otherwise, vectorize calculations
    isloop  = max(diff(srate));
    
    if(isloop > eps)
        
        %if sample rate not uniform, loop over rows as well
        for n = 1:M
            %get time sample index to weight data for width computation
            t       = 0:srate(n):ndata(k)*srate(n)-srate(n);
            t       = t(:);
            w       = abs(S(n,k).data);
            w1      = norm(t.*w).^2/(norm(w).^2); %first term
            w2      = norm(sqrt(t).*w).^2/(norm(w).^2); %second term
            tdur    = 2*pi*( w1 - w2.^2 ); %get signal width
            tdur    = sqrt(tdur); %square root to get pulse width
            tau(n,k)= tdur;
            
            %put fft stuff in here
            if(nargout>=2)
                
                %get time sample index to weight data for width computation
                f       = (1./srate(n))*(0:1/ndata(n):1-1/ndata(k));
                f       = f(:);
                w       = abs(fft(S(n,k).data));
                
                if(rem(ndata(n),2)<1),
                    
                    f       = f(1:ndata(n)/2);
                    w       = w(1:ndata(n)/2);
                    f(1)    = 0;
                    
                else
                    
                    %f should retain it's DC value
                    f       = f(1:(ndata(n)-1)/2);
                    %w should use the first sample over Nyquist as 0 hertz
                    w(1)    = w((ndata(n)+1)/2);
                    w       = w(1:(ndata(n)-1)/2);
                    
                end;
                
                w1      = norm(f.*w).^2/(norm(w).^2); %first term
                w2      = norm(sqrt(f).*w).^2/(norm(w).^2); %second term
                fdur    = 2*pi*( w1 - w2.^2 ); %get signal width
                fdur    = sqrt(fdur); %square root to get pulse width
                phi(n,k)= fdur;
      
            end;

        end;
    
    else
        
        %then all sample rates are uniform
        srate   = srate(1);
        %compute the time series matrix
        w           = abs([S(:,k).data]);
        t           = 0:srate:ndata*srate-srate;
        t           = repmat(t(:),1,M);
        
        [temp,n1]   = normColumns(t.*w);
        [temp,d1]   = normColumns(w);
        [temp,n2]   = normColumns(sqrt(t).*w);
        
        tdur        = 2*pi*( (n1./d1).^2 - (n2./d1).^4 );
        
        %square root to get pulse width
        tau(:,k)    = sqrt(tdur);
        
        %put fft stuff in here
        if(nargout>=2)
            
            %then all sample rates are uniform
            srate       = srate(1);
            %compute the time series matrix
            w           = abs(fft([S(:,k).data]));
            f           = (1./srate)*(0:1/ndata:1-1/ndata);
            f           = repmat(f(:),1,M);
            
            if(rem(ndata(1),2)<1),
                
                f       = f(1:ndata(1)/2,:);
                w       = w(1:ndata(1)/2,:);
                f(1,:)  = 0;
                
            else
                
                %f should retain it's DC value
                f       = f(1:(ndata(1)-1)/2,:);
                %w should use the first sample over Nyquist as 0 hertz
                w(1,:)  = w((ndata(1)+1)/2,:);
                w       = w(1:(ndata(1)-1)/2,:);
                
            end;

            [temp,n1]   = normColumns(f.*w);
            [temp,d1]   = normColumns(w);
            [temp,n2]   = normColumns(sqrt(f).*w);

            fdur        = 2*pi*( (n1./d1).^2 - (n2./d1).^4 );

            %square root to get pulse width
            phi(:,k)    = sqrt(fdur);

        end;

    end;
    
end;

if(nargout>=2), varargout{1} = phi;     end;

return;