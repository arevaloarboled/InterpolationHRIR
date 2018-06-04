function h = extractHRIRs(aed, lr,source)
% EXTRACTHRIRS extract the HRIRs speficied by the rows of aed array and
% returns them in the mxnx2 h matrix, where m is the number of taps, n the
% number of hrirs, and left and right channels are stored in the third 
% dimension.
%
% SYNOPSIS: h = extractHRIRs(aed,source)
%
% INPUT aed: an array of locations with amplitude, elevation, and distance
%            for columns. Amplitude and elevation are given in degrees.
%       source: What database to use. Default: ['PKU&IOA']
%
% OUTPUT h: a matrix of hrirs
%
% REMARKS 
%
% SEE ALSO 
%
% AUTHOR    : Julian Villegas
% $DATE     : 23-Mar-2017 14:20:14 $
% $Revision : 1.00 $
% DEVELOPED : 9.1.0.441655 (R2016b)
% FILENAME  : extractHRIRs.m
switch nargin
    case 0
        error('not enough arguments');
    case 1
        lr='r';
        source = 'PKU&IOA';
    case 2
        source = 'PKU&IOA';
end

switch source
    case 'PKU&IOA'
        m = 1024;% PKU&IOA has 1024 taps per HRIR
        %%%path = '/Users/julian/Documents/hrtf~/equalizedHrir/';
        path = 'equalizedHrir/';
        n = size(aed,1);
        h = zeros(m,n,2);
        for i=1:n
            h(:,i,:) = readEqualizedHrir(path,aed(i,3),aed(i,2),aed(i,1),lr);
        end
end

end
% ===== EOF ====== [extractHRIRs.m] ======
