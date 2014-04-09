function [S] = coralExtract(S,refField,fname,varargin)
% Extracts a PxN structure array from an MxN structure array S by
% extracting those rows of S with refField 'fname'. This function uses
% dynamic field names.
%
% USAGES
% [S] = coralExtract(S, refField, fname);
%
% INPUT
% S:        An Mx1 or MxN coral structure
% refField: The reference field over which to match fields between the
%           structures S and Sref. Must be input as string.
% fname:    The name of refField to extract
% type:     An optional input. If a boolean 0 or false is input,
%           coralExtract does not assume a record section is input.
%
% OUTPUT
% S:    The output Px1 or PxN coral structure containing only rows with
%       refField fname.
%
% Example:
% [S] = coralExtract(S, 'staChannel', 'ELZ');
% [S] = coralExtract(S, 'staChannel', 'ELZ',0); %not a record section
%-----------------------------------------------------------------------
% Latest Edit: 12.Jan.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

%type is true for record sections, false otherwise.
isRecSection = true;

if(not(isempty(varargin)))
    isRecSection = varargin{1};
end;

%a record section is being input
if(isRecSection)

    ind     = regexp( {S(:,1).(refField)}, fname );
    ind     = logical( cellfun(@length,ind) );
    S       = S(ind,:);

%general case
else

    [M,N]   = size(S);
    ind     = regexp( {S.(refField)}, fname );
    ind     = logical( cellfun(@length,ind) );

    % if input is a column vector, force output to be a column vector;
    if N==1
        
        S = S(ind(:));
        
    else

        ind     = reshape(ind,[M,N]);
        S       = S(ind);
        
    end;

end;

return;