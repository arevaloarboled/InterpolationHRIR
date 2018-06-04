function comparingMethods()
% This function compares the HRTF of a desired location with the HRTF
% obtained by using traditional VBAP, Equalizing filters, and Transaural
% audio. The sound y at the ears of the listener is given by the expression
% y = Hxhat, where H = [H_ll H_rl; H_lr H_rr] and xhat is the signal at the
% loudspeakers. H_xy refers to the transfer function from the loudspeaker x
% to the ear y.
%
% SYNOPSIS:
%
% REMARKS
%
% SEE ALSO
%
% AUTHOR    : Julian Villegas
% $DATE     : 12-Apr-2017 09:08:21 $
% $Revision : 1.00 $
% DEVELOPED : 9.2.0.538062 (R2017a)
% FILENAME  : comparingMethods.m
db = 'PKU&IOA'; %Qu's database
sr = 2^16; % Qu's database sampling rate
N=1024;
halfN = N/2;
binRes = sr/N;
minFreq = 70;
maxFreq = 16000;
%maxFreq = 20000;
minBin = ceil(minFreq/binRes);
maxBin = floor(maxFreq/binRes);

% Define loudspeakers
spAzi = [0,30,110,250,330];
%distance from the center of the head to the loudspeakers in the lab
spDis = 160;
spEle = 0;
m = SPK_defineSpeakers(spAzi,spEle,spDis);
spHs = extractHRIRs(m(1:end-1,:)',db); % retrieve speaker impulse responses

elEl = -40:10:90;
diDi = 160;
if maxFreq> 16000
    f_all = fopen('HRTFAzimuthOnlycorrelations.csv','w','n','UTF-8');
else
    f_all = fopen('HRTFComparisons70_16K.csv','w','n','UTF-8');
end
fprintf(f_all,'Azimuth, Elevation, Distance, Method, rmseL, rmseR\n');

cm = colormap();
for j=1:length(elEl)
    if elEl(j) <= 50
        azAngle = 5;
    elseif elEl(j)==60
        azAngle = 10;
    elseif elEl(j) == 70
        azAngle = 15;
    elseif elEl(j) == 80
        azAngle = 30;
    elseif elEl(j) == 90
        azAngle = 360;
    end
    % define angles for azimuth
    azAz = 0:azAngle:180;
    fePlots = zeros(length(azAz),3,maxBin-minBin+1,2);
    for i = 1:length(azAz)
        %%%% Equalizing filters %%%%%
        source_loc = [azAz(i),elEl(j),diDi];
        soHs = extractHRIRs(source_loc);
        % find the azimuthal spatialization
        idx = SPK_findSpeakers(m,source_loc(1));
        oriRMS = reshape(rms(soHs),[1,2]);
        hLen = length(spHs);
        corHs = SPK_getCorrectingFilters(spHs,idx, reshape(soHs,[hLen,2]));
        corHs = corHs(1:hLen,:); % 1024 is sufficient
        
        yA = computeCrossTalk(idx,corHs,spHs);
        if sum(m(1,idx) > 90)==2 && sum(m(1,idx) < 270)== 2
            % we're at the back of the listener
            yA = fliplr(yA);
        end
        
        %find elevation spatialization
        hLen = length(spHs);
        idx = find(m(1,:) > 180 & m(1,:) < 360);
        leftHs = SPK_getCorrectingFilters(spHs,idx, reshape(soHs,[hLen,2]));
        leftHs = leftHs(1:hLen,:); % 1024 is sufficient
        yL = computeCrossTalk(idx,leftHs,spHs);
        
        idx = find(m(1,:) > 0 & m(1,:) < 180);
        rightHs = SPK_getCorrectingFilters(spHs,idx, reshape(soHs,[hLen,2]));
        rightHs = rightHs(1:hLen,:); % 1024 is sufficient
        yR = computeCrossTalk(idx,rightHs,spHs);
        
        yE = yL+yR;
        
        y = SPK_crossFade(yA(1:1024,:), yE(1:1024,:), elEl(j));
        
        y = y.*oriRMS./rms(y);
        
        mag1 = 20*log10(abs(fft(y)));
        mag2 =  20*log10(abs(fft(soHs)));
        mag1 = mag1(minBin:maxBin,:);
        mag2 = mag2(minBin:maxBin,:);
        error1 = mag2(:,1)-mag1(:,1);
        error2 = mag2(:,2)-mag1(:,2);
        
        fePlots(i,1,:,1) = error1;
        fePlots(i,1,:,2) = error2;
        
        eqFrmsel = (mean(error1.^2))^0.5;
        eqFrmser = (mean(error2.^2))^0.5;
        fprintf(f_all,'%d,%d,%d,%s,%.10f,%.10f\n',...
            azAz(i),elEl(j),diDi,'Equalizer',eqFrmsel,eqFrmser);
        
        plotHRTFs(sr,N,idx,soHs,spHs,y,azAz,elEl,diDi,i,j,'EqFilter');
        
        %         %%%% Applying only HRTF %%%%%
        %         y = computeCrossTalk(idx,soHs,spHs);
        %         if sum(m(1,idx) > 90)==2 && sum(m(1,idx) < 270)== 2
        %             % we're at the back of the listener
        %             y = fliplr(y);
        %         end
        %         mag1 = 20*log10(abs(fft(y(1:halfN+1,:))));%abs(fft(y(1:1024,:))).^2;
        %         mag2 =  20*log10(abs(fft(soHs(1:halfN+1,:))));%abs(fft(soHs(:,1,:))).^2;
        %         onlyHRTFcoysol = 1;%(mag1(:,1),mag2(:,1),'type','Spearman');
        %         onlyHRTFcoysor = 1;%corr(mag1(:,2),mag2(:,2),'type','Spearman');
        %
        % %         error1 = mag2(:,1)/max(mag2(:)) - mag1(:,1)/max(mag1(:));
        % %         error2 = mag2(:,2)/max(mag2(:)) - mag1(:,2)/max(mag1(:));
        %         error1 = mag2(:,1) - mag1(:,1);
        %         error2 = mag2(:,2) - mag1(:,2);
        %         fePlots(i,2,:,1) = error1;
        %         fePlots(i,2,:,2) = error2;
        %
        %         onlyHRTFrmsel = (mean(error1.^2))^0.5;
        %         onlyHRTFrmser = (mean(error2.^2))^0.5;
        %
        %         fprintf(f_all,'%d,%d,%d,%s,%.3f,%.3f,%.10f,%.10f\n',...
        %             azAz(i),elEl(j),diDi,'onlyHRTF',onlyHRTFcoysol,onlyHRTFcoysor,onlyHRTFrmsel,onlyHRTFrmser);
        %
        %         plotHRTFs(sr,N,idx,soHs,spHs,y,azAz,elEl,diDi,i,j,'HRTFOnly');
        
        %%%% Applying VBAP %%%%%
        if elEl(j) == 0
            % find the azimuthal spatialization
            idx = SPK_findSpeakers(m,source_loc(1));
            ls_groups = findLsPairs(spAzi);
            layoutInvMtx = invertLsMtx(spAzi, ls_groups);
            gains2D = vbap(azAz(i), ls_groups, layoutInvMtx); % compute vbap gains
            x = zeros(1024,2);
            x(1,:) = [1 1];
            y=computeCrossTalk(idx,x,[gains2D gains2D(1)].* spHs);
            
            mag1 = 20*log10(abs(fft(y)));
            mag2 =  20*log10(abs(fft(soHs)));
            mag1 = mag1(minBin:maxBin,:);
            mag2 = mag2(minBin:maxBin,:);
            
            error1 = mag2(:,1) - mag1(:,1);
            error2 = mag2(:,2) - mag1(:,2);
            
            fePlots(i,3,:,1) = error1;
            fePlots(i,3,:,2) = error2;
            VBAPrmsel = (mean(error1.^2))^0.5;
            VBAPrmser = (mean(error2.^2))^0.5;
            fprintf(f_all,'%d,%d,%d,%s,%.10f,%.10f\n',...
                azAz(i),elEl(j),diDi,'VBAP',VBAPrmsel,VBAPrmser);
            
            plotHRTFs(sr,N,idx,soHs,spHs,y,azAz,elEl,diDi,i,j,'VBAP');
        end
        
        %%%%% Applying Cross-talk cancelation %%%%%
        % As if only the L-R surround speakers were used for spatialization
        idx = [3;4];
        % Following the nomenclature of  W. Garner "Transaural 3-D audio",
        % Fig.4
        if length(idx) == 1
            yl = spHs(:,idx(1),1);
            yr = spHs(:,idx(1),2);
        else
            G = conv(minPhaseHRIR(spHs(:,idx(1),2)),minPhaseHRIR(spHs(:,idx(2),1)))-...
                conv(minPhaseHRIR(spHs(:,idx(1),1)),minPhaseHRIR(spHs(:,idx(2),2)));
            G = SPK_invert(G,25);
            G = minPhaseHRIR(G);% 1/H
            xl = conv(minPhaseHRIR(soHs(:,1)),minPhaseHRIR(spHs(:,idx(1),2)))-...
                conv(minPhaseHRIR(soHs(:,2)),minPhaseHRIR(spHs(:,idx(2),2)));
            xl = conv(xl,G);
            xr = conv(minPhaseHRIR(soHs(:,2)),minPhaseHRIR(spHs(:,idx(2),1)))-...
                conv(minPhaseHRIR(soHs(:,1)),minPhaseHRIR(spHs(:,idx(1),1)));
            xr = conv(xr,G);
            
            yl = conv(xl,spHs(:,idx(2),1)) + ...
                conv(xr,spHs(:,idx(1),1));
            yl = yl(1:1024);
            yr = conv(xl,spHs(:,idx(2),2)) + ...
                conv(xr,spHs(:,idx(1),2));
            yr = yr(1:1024);
        end
        y = [yl, yr];
        
        mag1 = 20*log10(abs(fft(y)));
        mag2 =  20*log10(abs(fft(soHs)));
        mag1 = mag1(minBin:maxBin,:);
        mag2 = mag2(minBin:maxBin,:);
        
        error1 = mag2(:,1) - mag1(:,1);
        error2 = mag2(:,2) - mag1(:,2);
        
        fePlots(i,4,:,1) = error1;
        fePlots(i,4,:,2) = error2;
        
        CTCrmsel = (mean(error1.^2))^0.5;
        CTCrmser = (mean(error2.^2))^0.5;
        
        fprintf(f_all,'%d,%d,%d,%s,%.10f,%.10f\n',...
            azAz(i),elEl(j),diDi,'CTC',CTCrmsel,CTCrmser);
        
        plotHRTFs(sr,N,idx,soHs(:,1,:),spHs,y,azAz,elEl,diDi,i,j,'CTC');
    end
    % let's plot the error magnitude
    plotErrors(sr,N,fePlots,1,halfN,halfN,azAz,...
        ['Error_Elevation_' num2str(elEl(j)) '_Equalizing_filters'])
    %plotErrors(sr,N,fePlots,2,maxBin,halfN,azAz,'Only_HRTF')
    if elEl(j) ==0
        plotErrors(sr,N,fePlots,3,halfN,halfN,azAz,...
            ['Error_Elevation_' num2str(elEl(j)) '_VBAP'])
    end
    plotErrors(sr,N,fePlots,4,halfN,halfN,azAz,...
        ['Error_Elevation_' num2str(elEl(j)) '_Cross-talk_cancelation'])
end
fclose(f_all);


colormap(cm);
end

function plotHRTFs(sr,N,idx,soHs,spHs,y,azAz,elEl,diDi,i,j,method)
subplot(1,2,1);
plotmaglogf(soHs(:,1), sr, N, 'b', '-') % desired HRTF
hold on
% left signal at ear
plotmaglogf(y(:,1), sr, N, 'k', '-')
if elEl(j) == 0
    plotmaglogf(spHs(:,idx(1),1), sr, N, 'g', ':') % left signal, speaker 1
    if length(idx) == 1
        plotmaglogf(spHs(:,idx(1),1), sr, N, 'g', ':') % left signal, speaker 1
    else
        plotmaglogf(spHs(:,idx(2),1), sr, N, 'm', ':') % left signal, speaker 2
    end
    legend('target', 'spk1', 'spk2', 'atEar','Location','southwest')
else
    legend('target', 'atEar','Location','southwest')
end
grid on
title('Contralateral')
hold off

subplot(1,2,2);
plotmaglogf(soHs(:,2), sr, N, 'r', '-')
hold on
% right signal at ear
plotmaglogf(y(:,2), sr, N, 'k', '-')
if elEl(j) == 0
    plotmaglogf(spHs(:,idx(1),2), sr, N, 'g', ':') % right signal, speaker 1
    if length(idx) == 1
        plotmaglogf(spHs(:,idx(1),2), sr, N, 'g', ':') % right signal, speaker 1
    else
        plotmaglogf(spHs(:,idx(2),2), sr, N, 'm', ':') % right signal, speaker 2
    end
    legend('target', 'spk1', 'spk2', 'atEar','Location','southwest')
else
    legend('target', 'atEar','Location','southwest')
end
grid on
title('Ipsilateral')
hold off
name = sprintf('A%dE%dD%d_%s',azAz(i),elEl(j),diDi,method);

set(gcf,'NextPlot','add');
axes;
h = title(name,'Interpreter', 'none');%
set(gca,'Visible','off');
set(h,'Visible','on');
set(h,'Position',get(h,'Position')+[0 .03 0]);  % move up slightly
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',14)
set(fig,'PaperPositionMode','auto');
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[20 20 800 400]);
print(['./plots/' name '.pdf'],'-dpdf')
end

function plotErrors(sr,N,fePlots,idx,minBin,maxBin,halfN,azAz,plotName)
% fePlots=abs(fePlots);
% topC = max(fePlots(:));
% freq = (minBin : maxBin) * sr / N;
% myYticks = round(linspace(minBin,maxBin,7));
% myYLabels = round(linspace(freq(1),freq(end),7)/1000,1);
%
% degrees = sprintf('%c', char(176));
% colormap(1-gray)
% subplot(1,2,1);
% tp = reshape(fePlots(:,idx,:,1),[length(azAz),maxBin])';
% imagesc(tp)
% title('Contralateral')
% if length(azAz)>1
%     xticks(linspace(1,length(azAz),7))
%     xticklabels(linspace(azAz(1),azAz(end),7))
%     yticks(myYticks)
%     yticklabels(myYLabels)
% end
% set(gca,'YDir','normal')
% ylabel('Freq./kHz');
% ylim([minBin,maxBin])
% xlabel(['Azimuth/' degrees]);
% c = colorbar('southoutside');
% c.Label.String = 'dB';
% caxis manual
% caxis([0 topC]);
%
% subplot(1,2,2);
% tp = reshape(fePlots(:,idx,:,2),[length(azAz),maxBin])';
% imagesc(tp)
% title('Ipsilateral')
% if length(azAz)>1
%     xticks(linspace(1,length(azAz),7))
%     xticklabels(linspace(azAz(1),azAz(end),7))
%     yticks(myYticks)
%     yticklabels(myYLabels)
% end
% set(gca,'YDir','normal')
% ylabel('Freq./kHz');
% ylim([0,maxBin])
% xlabel(['Azimuth/' degrees]);
% c = colorbar('southoutside');
% c.Label.String = 'dB';
% caxis manual
% caxis([0 topC]);
%
% name = plotName;
%
% set(gcf,'NextPlot','add');
% axes;
% h = title(name,'Interpreter', 'none');%
% set(gca,'Visible','off');
% set(h,'Visible','on');
% set(h,'Position',get(h,'Position')+[0 .04 0]);  % move up slightly
% fig=gcf;
% set(findall(fig,'-property','FontSize'),'FontSize',14)
% set(fig,'PaperPositionMode','auto');
% set(fig,'PaperOrientation','landscape');
% set(fig,'Position',[20 20 800 400]);
% print(['./plots/errorAnalysis/' name '.pdf'],'-dpdf')
end
% ===== EOF ====== [comparingMethods.m] ======
