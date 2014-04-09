function [S] = arithmStruct(S,refField,opt,W)
% arithmStruct uses arithmetic function handle opt to compute an operation
% on S.refField using W. S must be 1x1. Results are extended to MxN
% structure arrays using arrayfun with @arithmStruct.
%
% USAGES
% [S] = arithmStruct(S,refField,opt,W);
%
% INPUT
% S:        A 1x1 structure array
% refField: A field in S
% opt:      A function handle, such as @plus and @times
% W:        The second required argument of opt, i.e. if opt is @times,
%           then W and S.(refField) are multiplied together.
%
% OUTPUT
% S:    The input coral structure with field refField changed by W.
%
% Example: Change the polarity of a time series in field data.
% x         = zeros(100,1);
% x(49:51)  = 1;
% A.data    = x;
% A.channel = 'X';
% [A] = arithmStruct(A,'data',@times,-1);
% plot(x,'k'); hold on; plot(A.data,'r');
%-----------------------------------------------------------------------
% Latest Edit: 11.Mar.2011
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

temp            = opt(W,S.(refField));
S.(refField)    = temp;

return;