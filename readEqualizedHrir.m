function hrir = readEqualizedHrir(filepath,d,e,a,lr)
% readEqualizedHrir reads the diffused-field equalized version of PKU&IOA
% database.
%
% SYNOPSIS: hrir = readEqualizedHrir(filepath,a,e,d)
%
% INPUT filepath: path to the directory where the equalized hrir are stored.
%       a: Azimuth in degress, counterclock-wise
%       e: elevation in degrees
%       d: distance in cm
%       lr: which pinna to use? 'l' for large pinna on the left, and 'r' for
%       small pinna on the right
%
% OUTPUT hrir: a 1024x2 matrix with the corresponding hrir
%
% REMARKS
% Distance can only be 20 30 40 50 75 100 130 160
% Elevation from -40 to 90 in step of 10
% Azimulth from 0 to 355 in step of 5 (elev <= 50),
%          from 0 to 350 in step of 10 (elev == 60),
%          from 0 to 345 in step of 15 (elev == 70),
%          from 0 to 330 in step of 30(elev == 80), and
%          0 (elev == 90)
%
% Version 2 includes the flag for pinna size
%
% AUTHOR    : Julian Villegas
% $DATE     : 23-Mar-2017 14:20:14 $
% $Revision : 2.00 $
% DEVELOPED : 9.1.0.441655 (R2016b)
% FILENAME  : readEqualizedHrir.m
%
if nargin ~= 5
    lr='r';
end

f = filesep;
if a == 360
    a = 0;
end

flipa = 360 - a;
if (flipa == 360)
    flipa = 0;
end

if e == 90
    a = 0;
end

filename = [filepath f createFilename(d,e,a)];
p = fopen(filename,'r');
filename = [filepath f createFilename(d,e,flipa)];
q = fopen(filename,'r');
if p == -1 || q == -1
    error(['Unable to open ',filename]);
else
    hrir1 = fread(p, 'single');
    hrir1 = reshape(hrir1, 1024, 2);
    hrir2 = fread(q, 'single');
    hrir2 = reshape(hrir2, 1024, 2);
    if (strcmpi('L',lr))
        hrir = [hrir1(:,1), hrir2(:,1)];
    elseif (strcmpi('R',lr))
        hrir = [hrir2(:,2), hrir1(:,2)];
    end
    fclose(p);
    fclose(q);
end
end