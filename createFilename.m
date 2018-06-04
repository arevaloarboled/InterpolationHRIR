function name = createFilename(di,el,az)
% CREATEFILENAME creates the file name of the HRIR defined by a,e,d.
%
% SYNOPSIS: name = createFilename(a,e,d)
%
% INPUT     az: azimuth in degrees
%           el: elevation in degrees
%           di: distance in cm
%
% OUTPUT name: the string with the name
%
% REMARKS based on the homonymous script provided by PKU&IOA
%
% SEE ALSO 
%
% AUTHOR    : Julian Villegas
% $DATE     : 23-Mar-2017 14:20:14 $
% $Revision : 1.00 $
% DEVELOPED : 9.1.0.441655 (R2016b)
% FILENAME  : createFilename.m
extension = '.dat';

az = round(az);
di = round(di);
if di < 100
    d = ['0' num2str(di)];
else
    d = num2str(di);
end

if el < 0
    e = ['m' num2str(el * -1)];
else
    e = ['0' num2str(el)];
    if el == 0
        e = ['0' e];
    end
end

if az < 100
    a = ['0' num2str(az)];
    if az < 10
        a = ['0' a];
    end
    
else
    a = num2str(az);
end
name = ['D' d '_E' e '_A' a extension];
end
