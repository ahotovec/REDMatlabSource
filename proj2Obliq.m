function [P,Ehs,Esh,Ia] = proj2Obliq(H,S)
%proj2Obliq computes oblique projectors based upon subspaces span(H) and
%span(S). The oblique projector that projects onto span(H) has zero
%projection onto span(S). Likewise, the projector that projects onto
%span(S) has zero projection onto span(H).
%The singular values and angles between subspaces are optionally returned
%too.
%
% USAGE
% [P,Ehs,Esh] = proj2Obliq(H,S)
%
% INPUT
% H:    An NxM matrix (N>M) whose columns span the projector range
% S:    An NxT matrix (N>T) whose columns span the null space of the
%       projector. M + T <= N.
%
% OUTPUT
% P:    The projector onto span([H,S]) where P = Ehs + Ehs.
% Ehs:  The projector onto span(H) with zero projection onto span(S).
% Esh:  The projector onto span(S) with zero projection onto span(H).
% Ia:   The projector onto the space perpendicular to span([H,S]);
%

[n,m]   = size(H);
[nr,t]  = size(S);

if(m+t > n)
    
    error('Range + null dimension exceeds column number in first matrix');
    
end;

if( nr ~= n)
    
    error('Range and null space incompatible in row number of matrices'); 
    
end;

R       = [[H'*H, H'*S];[S'*H, S'*S]];
P       = ([H,S]*(R\eye(size(R)))*[H,S]');
Ehs     = ([H,zeros(size(S))]*(R\eye(size(R)))*[H,S]');
Esh     = ([zeros(size(H)),S]*(R\eye(size(R)))*[H,S]');

%Ia likely incorrect. Need to check.
Ia      = eye(length(H)) - Ehs - Esh;

return;