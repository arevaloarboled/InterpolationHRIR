function h = interpolateMinPhase(hrirs,ws)
% interpolateMinPhase linearly interpolates the HRIRs specified in the matrix 
% hrirs with the weights ws. Hrirs should be in minimum-phase form
%  
% SYNOPSIS: h = interpolateMinPhase(hrirs,ws)  
% 
% INPUT hrirs: a mxnx2 matrix of HRIRs where m is the number of taps, 
%              n is the number of HRIRs, the third dimension correspons to 
%              the the left (1) and right channel (2) of a given HRIR.
%       ws: a vector of weigths of length n used for the interpolation
%
% OUTPUT h: an mx2 matrix with the interpolated HRIRs  
% 
% REMARKS  
% 
% SEE ALSO 
% 
% AUTHOR    : Julian Villegas 
% $DATE     : 23-Mar-2017 14:00:24 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.1.0.441655 (R2016b) 
% FILENAME  : interpolateHRTFs.m 
%
n = size(hrirs,2);

if n ~= length(ws)
    error('the number of weights and minimum-phase do not match')
end

for i=1:n
    hrirs(:,i,:) = hrirs(:,i,:) * ws(i);
end
h = sum(hrirs,2);
h = [h(:,:,1), h(:,:,2)];
end 
% ===== EOF ====== [interpolateHRTFs.m] ======  
