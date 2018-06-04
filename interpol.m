function s=interpol(p,m)
load('Delaunay_PKU&IOA_HRTF.mat')
source_pos = sphToRect(p);
if m==0
   idx = pointLocation(dt,source_pos); 
else
   idx=nan;
end
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
%d=rad2deg(d);
aed = [a,e,d];
if sum(ws>1) || isnan(sum(ws))
    ws = ones(1,length(ws))/length(ws);
    disp('**');
end
hrirs = extractHRIRs(aed);
mPhase=zeros(1024,4,2);
for i=1:2
    mPhase(:,:,i)=minPhaseHRIR(hrirs(:,:,i));
end
%mPhase = minPhaseHRIR(hrirs);
s = interpolateMinPhase(mPhase,ws);
end