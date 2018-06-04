function s = spatialize(y,p)
% SPATIALIZE binaurally spatialize the vector y so its apparent location
% becomes p
%
% SYNOPSIS: s = spatialize(y,p)
%
% INPUT y: an audio vector sampled at 2^16 Hz
%       p: a location in spherical coordinates (azimuth, elevation, distance)
%
% OUTPUT s: a binaural version of y
%
% REMARKS There seems to be a bug with the implementation of DeLauney
% triangulation in Matlab, especifically, the function pointLocation()
% returns NaN for some locations, to circumvent that, the location of such
% points is randomly changed in the vicinity of (-0.5, 0.5) in each axis.
%
% SEE ALSO 
%
% AUTHOR    : Julian Villegas
% $DATE     : 26-Mar-2017 16:25:13 $
% $Revision : 1.00 $
% DEVELOPED : 9.2.0.538062 (R2017a)
% FILENAME  : spatialize.m
load('Delaunay_PKU&IOA_HRTF.mat')
sr = 2^16;
source_pos = sphToRect(p);
idx = pointLocation(dt,source_pos);
while isnan(idx)
    % dunno why Matlab fails
    tmpP = [source_pos(1)+rand()-0.5,...
            source_pos(2)+rand()-0.5,...
            source_pos(3)+rand()-0.5];
    idx = pointLocation(dt,tmpP);
    disp('**')
end
c = dt.ConnectivityList(idx,:);
c = dt.Points(c,:);
ws = cartesianToBarycentric(dt,idx,source_pos);
[a,e,d] = cart2sph(c(:,1),c(:,2),c(:,3));
%[a,e,d] = rad_to_deg(a,e,d);
a=rad2deg(a);
e=rad2deg(e);
aed = [a,e,d];

if sum(ws>1) || isnan(sum(ws))
    ws = ones(1,length(ws))/length(ws);
    disp('**');
end
hrirs = extractHRIRs(aed);
itds = zeros(4,3);
for i=1:4
    [itds(i,1),itds(i,2),itds(i,3)] = computeITD(hrirs(:,i,:),sr);
    itds(i,:) = round(itds(i,:)*sr);
end
mPhase = minPhaseHRIR(hrirs);

h = interpolateMinPhase(mPhase,ws);
t = round(interpolateITDs(itds,ws));

% convolve with sound
y = [conv(y,h(:,1)), conv(y,h(:,2))];

% apply ITD
s = [[zeros(t(2),1); y(:,1); zeros(t(3),1)], ...
     [zeros(t(3),1); y(:,2); zeros(t(2),1)]];
end
% ===== EOF ====== [spatialize.m] ======
