function [el,er,sdl,sdr,y,soHs]=comparingInter(p,lr)
switch nargin
    case 1
        lr='r';        
end
N=1024;
sr=2^16;
binRes = sr/N;
minFreq = 70;
maxFreq = 16000;
%maxFreq = 20000;
t = [1:1:1024];
f=440;
x = sin(2*pi*f/sr*t);
y=spatialInterpol(x,p,1,lr);
soHs=spatialInterpol(x,p,0,lr);
%y=interpol(p,1);
%soHs=interpol(p,0);
%%soHs=extractHRIRs(p);
%%soHs=soHs(:,1,:);
%%soHs=readEqualizedHrir('equalizedHrir',p(1),p(1),p(1),'r');
minBin = ceil(minFreq/binRes);
maxBin = floor(maxFreq/binRes);
mag1 = 20*log10(abs(fft(y)));
mag2 =  20*log10(abs(fft(soHs)));
mag1 = mag1(minBin:maxBin,:);
mag2 = mag2(minBin:maxBin,:);
el = mag2(:,1)-mag1(:,1);
er = mag2(:,2)-mag1(:,2);
sdl = (mean(el.^2))^0.5; % SD of the left channel
sdr = (mean(er.^2))^0.5; % SD o the right channel
end
% ===== EOF ====== [comparingInter.m] ======
