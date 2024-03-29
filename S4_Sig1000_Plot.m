%% ************************************************************************
%***** Nortek Signature 1000 Processing - Step 4                      *****
%*****                                                                *****
%***** S4_Sig1000_Plot.m                                              *****
%*****                                                                *****
%***** A script to plot ADCP velocity data as checkerboard plots and  *****
%***** mean velocity profiles. Additionaly turbidity data is plotted. *****
%*****                                                                *****
%***** Input: .mat files generated by S3_Sig1000_QC.m                 *****
%***** Output: .figs and .png from plots                              *****
%*****                                                                *****
%*****                                                                *****
%***** Prepared by J. Tiede, Januar 2018                              *****
%**************************************************************************

%% Clean workspace
clc; clear; close all;

%% Add tools and functions to MATLAB path.
addpath('D:\Studium\MasterThesis\70_Matlab\ADCP\NortekSig1000\Functions')
addpath('D:\Studium\MasterThesis\70_Matlab\ADCP\NortekSig1000\Tools')

%% Define ADCP base folder.
basefolder = 'D:\LuFI\Projekte\LuggagePoint_JournalPaper\60_Data\Nortek_061117\',

%% Folder where files live.
DataFolder = [basefolder '05_QCed\'];
if ~isdir(DataFolder)
  errorMessage = sprintf('Error: Folder does not exist\n%s', DataFolder);
  uiwait(warndlg(errorMessage));
  return;
end

%% User specified variables
fs = 8 % Sampling frequency
Ravt = 5; % Period in minutes for averaging Reynolds stresses
Vavt = 10; % Period in minutes for averaging horizontal velocity
npoints = fs * 60 * Ravt; % 8 Hz * 60 s * Ravt min
vpoints = fs * 60 * Vavt % 8 Hz * 60s * Vavt min

%% Loop over all files/folders
p = 1
%% Define Data folder.
datafolder = [basefolder, '05_QCed\'];

%% Load quality controlled data
load([datafolder,'PNCU' num2str(p) '_p_qc.mat'])

%% Average of Speed
for i = 1: size(ADCP.Speed,2)
    spe = ADCP.Speed(:,i);
    M  = size(spe, 1) - mod(size(spe,1), vpoints);
    y  = reshape(spe(1:M), vpoints, []);
    ADCP.Speed_Avg(:,i) = transpose(sum(y, 1) / vpoints);
end
for i = 1: size(ADCP.Time,2)
    tim = ADCP.Time(:,i);
    M  = size(tim, 1) - mod(size(tim,1), vpoints);
    y  = reshape(tim(1:M), vpoints, []);
    ADCP.Time_avg(:,i) = transpose(sum(y, 1) / vpoints);
end

%% Vertical mean of speed
ADCP.Speed_VerMean = nanmean(ADCP.Speed_Avg,2)
ADCP.Speed_VerMean = smooth(ADCP.Speed_VerMean)

%% Plot mean speed
[img,ax] = figur([16.5 12],{'Calibri',10,'normal','latex'},'on')
grid on
xlim([ADCP.Time_avg(1) ADCP.Time_avg(end)])
% title('Velocity Magnitude depth averaged - Flood Stream 06.11.17')
datetick('keeplimits')
ylabel('Velocity [m/s]')
xlabel('Time [HH:MM]')

plot(ADCP.Time,ADCP.Speed,'Color',[0.7 0.7 0.7])
hold on
plot(ADCP.Time_avg,ADCP.Speed_VerMean,'b','Linewidth',1.5)

% 
save2word('VEL_mean.doc')
%% Pcolor Plot Speed Amplitude Correlation
[img,ax] = figur([16.5 14.5],{'Calibri',10,'normal','latex'},'on')
subplot(3,1,1)
    PS = pcolor(ADCP.Time,ADCP.Range,ADCP.Speed')
    set(PS, 'EdgeColor', 'none');
    title('Horizontal Velocity')
    cb = colorbar
    ylabel(cb, 'Velocity [m/s]')
    caxis([0 1])
    ylim([0.6 12.5])
    set(gca,'YDir','normal')
    ylabel('m above bed [m]')
    datetick('keeplimits')
    set(gca,'layer','top')

subplot(3,1,2)
PC = pcolor(ADCP.Time,ADCP.Range,ADCP.c1')
    set(PC, 'EdgeColor', 'none');
    title('Correlation')
    cb = colorbar
    ylabel(cb, 'Correlation [%]')
    caxis([0 100])
    ylim([0.6 12.5])
    set(gca,'YDir','normal')
    ylabel('m above bed [m]')
    datetick('keeplimits')
    set(gca,'layer','top')

    subplot(3,1,3)
PA = pcolor(ADCP.Time,ADCP.Range,ADCP.a1')
    set(PA, 'EdgeColor', 'none');
    title('Amplitude')
    cb = colorbar
    ylabel(cb, 'Amplitude [dB]')
    caxis([0 85])
    ylim([0.6 12.5])
    set(gca,'YDir','normal')
    ylabel('m above bed [m]')
    datetick('keeplimits')
    set(gca,'layer','top')    
    
% save2word('Speed_Cor_Amp.doc')
%% Device Orientation
[img,ax] = figur([16.5 14.5],{'Calibri',10,'normal','latex'},'on')
subplot(3,1,1)
plot(ADCP.Time,ADCP.phi1)
xlim([ADCP.Time(1) ADCP.Time(end)])
datetick('keeplimits')
title('Heading')
grid on
ylabel('Orientation [�]')
hold on

Head = mean(ADCP.phi1)

subplot(3,1,2)
plot(ADCP.Time,ADCP.phi2)
xlim([ADCP.Time(1) ADCP.Time(end)])
datetick('keeplimits')
title('Pitch')
grid on
ylabel('Tilt [�]')
Pitch = mean(ADCP.phi2)

subplot(3,1,3)
plot(ADCP.Time,ADCP.phi3)
xlim([ADCP.Time(1) ADCP.Time(end)])
datetick('keeplimits')
title('Roll')
grid on
ylabel('Tilt [�]')
xlabel('Time [HH:MM]')
Roll = mean(ADCP.phi3)

% save2word('Orientation.doc')