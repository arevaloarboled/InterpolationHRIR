function s = spatialInterpol(y,p,m,lr)
load('Delaunay_PKU&IOA_HRTF.mat')
sr = 2^16;
source_pos = sphToRect(p);
if m==-1
    h=readEqualizedHrir('equalizedHrir',p(3),p(2),p(1),lr);
    s = [conv(y,h(:,1)), conv(y,h(:,2))];
    return;
end
%delete point in the mesh;
%%>>>>>>>>>>>>>>>>>>>>>>>>
if m==1
    idx = pointLocation(dt,source_pos);
    ctmp = dt.ConnectivityList(idx,:);
    c = dt.Points(ctmp,:);
    [a,e,d] = cart2sph(c(:,1),c(:,2),c(:,3));
    a=rad2deg(wrapTo2Pi(a));
    e=rad2deg(e);
    aed = [a,e,d];
    for i=1:4    
        if all(abs(aed(i,:) - p)<=1e-10)
            dt.Points(ctmp(i),:)=[];
            break;
        end
    end
end
%%<<<<<<<<<<<<<<<<<<<<<<<<
%start to search other points
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
a=rad2deg(wrapTo2Pi(a));
e=rad2deg(e);
aed = [a,e,d];

if sum(ws>1) || isnan(sum(ws))
    ws = ones(1,length(ws))/length(ws);
    disp('**');
end
hrirs = extractHRIRs(aed,lr);
%%%itds = zeros(4,3);
for i=1:4
    [itds(i,1),itds(i,2),itds(i,3)] = computeITD(hrirs(:,i,:),sr);
    itds(i,:) = round(itds(i,:)*sr);
end
mPhase=zeros(1024,4,2);
for i=1:2
    mPhase(:,:,i)=minPhaseHRIR(hrirs(:,:,i));
end
%mPhase = minPhaseHRIR(hrirs);

h = interpolateMinPhase(mPhase,ws);
%%%t = round(interpolateITDs(itds,ws));

% convolve with sound
%y = [conv(y,h(:,1)), conv(y,h(:,2))];
s = [conv(y,h(:,1)), conv(y,h(:,2))];

% apply ITD
%s = [[zeros(t(2),1); y(:,1); zeros(t(3),1)], ...
%     [zeros(t(3),1); y(:,2); zeros(t(2),1)]];
end
% ===== EOF ====== [spatialize.m] ======
