function [cc, t, indx, ccLag, ccMax, ierr] = coralCrossCorr(data, opt)
%   coralCrossCorr    cross correlate data in coral structure
% USAGE: [cc, t, indx, CLag, maxC, ierr] = coralCrossCorr(data, opt);
%
% Cross correlate data contained in the data structure
% Calculations are done in the frequency domain because it is fast
% but are identical to cross correlation in the time domain.  Cross correlation 
% results are the same as the matlab routine xcorr using the 'coeff' option and the mlag
% option.  All time series must have the same sample interval and same number of samples.
%
% See coral for explanation of data structure
% opt is an optional structure with optional fields: 
% opt.mlag (only keep lags from -opt.mlag:opt.mlag so the
%          cross correlograms have 2*mlag+1 elements.
% there are several options for determining which pairs of seismograms to cross correlate
%
% opt.ndiag (integer) provides an option to calculate correlograms for some, but not all pairs
%           of seismograms. calculate only pairs within opt.ndiag of each other.
%           for example, opt.ndiag=1 calculates only the auto-correlations, 2 calculates
%           cross correlations for all adjacent seismograms
%
% opt.refseis (integer) is index to one of the seismograms which is cross correlated 
%           against the rest.  This overrides opt.ndiag
%           
%
% opt.indx  nx2 matrix of indices pointing to desired pairs of seismograms to correlate 
%           (same as output parameter indx) This overrides opt.ndiag and opt.refseis
%
% OUTPUT: 
% 
%  cc(2*opt.mlag+1,N) contains correlations as column vectors ordered by:
%  11, 12, 13, ... 1ndiag, 22, 23, 24, ...2ndaig+1, ...
%  t   is a time lag vector (s)
%  indx is a (2xN) matrix containing integer pointers to the ith and jth seismogram for 
%         each of the N cross correlograms
%  ccLag = sparse(ndata,ndata) matrix containing the lags(s) corresponding to the best correlations
%  ccMax = sparce(ndata,ndata) matrix containing the maximum correlation for each pair of seismograms
%  ierr is a vector of zeros if there are no errors, 
%       = sample intervals for each seismogram if they are not all the same
%       = number of samples in each seismogram if they are not all the same
%
% K. Creager  kcc@ess.washington.edu   10/18/2004   last mod: 2/28/2004

% initialize outputs
ndata= length(data);
cc   = [];
t    = [];
indx = [];
ccLag= [];
ccMax= [];
ierr = zeros(ndata,1);

% Check that sample intervals are all the same to within a tolerance
sampInt =  [data.recSampInt]';
if (max(sampInt) - min(sampInt)) / sampInt(1) > .001;
  disp('WARNING: sample intervals must be the same to run coralCrossCorr')
  ierr = sampInt;
  return
end

% Check that the number of samples are the same for each seismogram
numData=zeros(ndata,1);
for k=1:ndata
  numData(k) = length(data(k).data);
end
if max(numData) ~= min(numData)
  disp('WARNING: each seismogram must have the same number of samples in coralCrossCorr')
  ierr = numData;
  return
end


% set defalut values for mlag and ndiag
mlag    = numData(1) - 1;
ndiag   = ndata;
refseis = [];
indx    = [];
flds      = {'mlag' , 'ndiag', 'refseis', 'indx'};
% if values are passed in throught 'opt' change them from their defaults
if nargin >1;
  if isstruct(opt);
    for k=1:length(flds);
      fld = flds{k};
      if any(strcmp(fields(opt),fld));
        eval(sprintf('%s=opt.%s;',fld,fld));
      end
    end
  end
end
ndiag = min(ndiag,ndata);   % ndiag can not be bigger than the number of data

% calculate indx matrix (this defines all the pairs of seismograms to 
% cross correlate and their output order
if length(indx)>0;         % this indx matrix was already supplied in opt.indx, don't need to create it
elseif length(refseis)>0;  % calculate cross correlations against a reference seismogram
  if length(refseis)==1 & refseis(1)== max(min(round(refseis(1)),ndata),1); % check to see whether refseis is an integer between 1 and ndata
     indx=zeros(2,ndata); 
     indx(1,:)=refseis; 
     indx(2,:) = [1:ndata];
  else
    disp('WARNING: opt.refseis must be an integer between 1 and ndata in coralCrossCorr')
    ierr = opt.refseis;
    return
  end
else
  ncc= ndiag*(2*ndata-ndiag+1)/2;                 % number of cross correlations to calculate
  indx = zeros(2,ncc);
  k=0;
  for i=1:ndata; 
    for j=i:min(i+ndiag-1,ndata); 
      k=k+1;
      indx(:,k) = [i;j];
    end 
  end
end 
ncc  =  size(indx,2);
indx(3,:) = (indx(1,:)-1)*ndata + indx(2,:); 


% pad data with zeros to double the number of samples to avoid wrap around
% then pad with more zeros to get to the next power of 2 for efficient fft calculation
n  = numData(1);                                 % number of samples per seismogram
n2 = 2^ceil(log(2*n)/log(2));                    % number of samples in fft 

ncut     = n2/2 - mlag - 1;
mlagInd = [ncut+2:n2-ncut];                      % index to extract -mlag:mlag from cross correlations
t       = [-mlag:1:mlag]'*sampInt(1);            % time axis (centered at 0 lag)

% initialize arrays
cc    = zeros(2*mlag+1,ncc); 
ccLag = sparse(ndata,ndata);
ccMax = sparse(ndata,ndata);
norm0 = zeros(ndata,1);
ff    = zeros(n2,ndata);   
for idata=1:ndata                                % calculate fourier transforms
  ff(:,idata) = fft(data(idata).data,n2);
end


% loop over desired pairs of seismograms and calculate cross-correlation in the frequency
% domain, then apply the inverse fft to move to the time domain.  time shift cross correlogram 
% so zero lag is at center of vectors and keep just the part from -mlag to mlag
for k=1:ncc;
  i=indx(1,k);
  j=indx(2,k);
  tmp = real( ifft( ff(:,i).*conj(ff(:,j)),n2 ) ) ;
  if i == j; norm0(i)=sqrt(tmp(1)); end;
  tmp = fftshift(tmp);
  cc(:,k)=tmp(mlagInd);
end

% calculate normalizations if not already done
for i=1:ndata
  if norm0(i)==0;
    tmp = real( ifft( ff(:,i).*conj(ff(:,i)),n2 ) ) ;
    norm0(i)=sqrt(tmp(1)); 
  end
end


% go back through each time series to normalize them such that each autocorrelation
% equals one at zero lag.  This is the same as using the 'coeff' option in the matlab
% supplied cross correlation program xcorr (which is much slower than this routine).
for k=1:ncc;
  i=indx(1,k);
  j=indx(2,k);
  cc(:,k)=cc(:,k)/(norm0(i)*norm0(j));
  [tmpMax,tmpInd]=max(cc(:,k));
  ccLag(i,j)=t(tmpInd);  ccLag(j,i)=t(tmpInd);
  ccMax(i,j)=tmpMax;     ccMax(j,i)=tmpMax;
end