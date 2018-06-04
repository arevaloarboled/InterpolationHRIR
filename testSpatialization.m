function testSpatialization(method, s)
%
% This script tests Qu's HRIR and HRTF data.
%
% By J. Villegas, 2016
% University of Aizu
f = filesep;
sr = 65536;

switch nargin
    case 0
        method = 'hrir';
        s = rand(floor(sr/4),1) * 2 - 1;
    case 1
        s = rand(floor(sr/4),1) * 2 - 1;
    case 2
        [s, tsr] = audioread(s);
        if tsr~=sr
            s = resample(s,sr,tsr);
        end
end

distance = (2:16) * 10;
elevation = (-4:1:9) * 10;
mAzi = {(0:71) * 5 % <= 50
    (0:35) * 10 % == 60
    (0:23) * 15 %== 70
    (0:11) * 30 % == 80
    0}; % measured azimuths depend on the elevation


if strcmpi(method, 'hrir')
    filepath = ['.' f 'equalizedHrir' f];
else
    filepath = ['.' f 'equalizedHrtf' f];
end

for i=1:length(distance)
    dis = distance(i);
    for j = 1:length(elevation)
        ele = elevation(j);
        switch  ele
            case 90
                a = mAzi{5};
            case 80
                a = mAzi{4};
            case 70
                a = mAzi{3};
            case 60
                a = mAzi{2};
            otherwise
                a = mAzi{1};
        end
        for k = 1:length(a)
            azi = a(k);
            disp(['dist: ' num2str(dis)...
                ' -- ele: ' num2str(ele) ...
                ' -- azi: ' num2str(azi)])
            disp('Press any key or mouse button over the Figure to continue')
            if strcmpi(method, 'hrir')
                [h] = readEqualizedHrir(filepath, dis, ele, azi);
                L = conv(s, h(:,1));
                R = conv(s, h(:,2));
            else
                [h] = readEqualizedHrtf(filepath, dis, ele, azi);
                L = conv(s, real(ifft(h(:,1))));
                R = conv(s, real(ifft(h(:,2))));
            end
            spatialized = [L R];
            sound(spatialized,sr);
            java.lang.Runtime.getRuntime.freeMemory;
            waitforbuttonpress;
        end
    end
end
end
