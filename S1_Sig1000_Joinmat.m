%% ************************************************************************
%***** Nortek Signature 1000 Processing - Step 1                      *****
%*****                                                                *****
%***** S1_Sig1000_Joinmat.m                                           *****
%*****                                                                *****
%***** A script for reading in ADCP .mat files, concatenating         *****
%***** relevant parameters and saving as one .mat file.               *****
%*****                                                                *****
%***** Input: .mat files generated by Signature Deployment Software   *****
%***** Output: .mat files                                             *****
%*****                                                                *****
%*****                                                                *****
%***** Prepared by J. Tiede, Januar 2018                              *****
%**************************************************************************

%% Clean workspace
clc; clear; close all;

%% Add tools and functions to MATLAB path.
addpath('C:\Users\j_tie\OneDrive\LuFI\Projekte\LuggagePoint_JournalPaper\70_Matlab\ADCP\NortekSig1000\Functions')
addpath('C:\Users\j_tie\OneDrive\LuFI\Projekte\LuggagePoint_JournalPaper\70_Matlab\tools')

%% Define ADCP base folder.
basefolder = 'C:\Users\j_tie\OneDrive\LuFI\Projekte\LuggagePoint_JournalPaper\60_Data\Nortek_061117\'

%% Folder where files live.
DataFolder = [basefolder '02_Converted\'];
if ~isdir(DataFolder)
  errorMessage = sprintf('Error: Folder does not exist\n%s', DataFolder);
  uiwait(warndlg(errorMessage));
  return;
end

%% List Folders
ListFolders = dir(DataFolder);
ListFolders(1:2) = []

%% Folder where files are saved to
savepath = [basefolder '03_Joined\'];

%% Loop over all files/folders
p = 1
%% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile([DataFolder, num2str(p) '\']);
ListFiles = dir(filePattern);
ListFiles(1:2) = [];

  baseFileName = ListFiles(1).name;
  fullFileName = fullfile([ListFolders(p).folder '\' num2str(p) '\', baseFileName]);
  fprintf(1, 'Now reading %s\n', fullFileName);
  data = load(fullFileName);  


%% Define Range
ADCP_CFG.Cs = double(data.Config.Burst_CellSize);
ADCP_CFG.Bd = double(data.Config.Burst_BlankingDistance);
ADCP_CFG.Cn = double(data.Config.Burst_NCells);


%% Concatenate Variables
fnumber = length(ListFiles);
ADCP.Speed = data.Data.Burst_VelSpeed;
ADCP.VelDirection = data.Data.Burst_VelDirection;
ADCP.b1 = data.Data.Burst_VelBeam1;
ADCP.b2 = data.Data.Burst_VelBeam2;
ADCP.b3 = data.Data.Burst_VelBeam3;
ADCP.b4 = data.Data.Burst_VelBeam4;
ADCP.b5 = data.Data.IBurst_VelBeam5;
ADCP.VelN = data.Data.Burst_VelNorth;
ADCP.VelE = data.Data.Burst_VelEast;
ADCP.VelX = data.Data.Burst_VelX;
ADCP.VelY = data.Data.Burst_VelY;
ADCP.VelUp1 = data.Data.Burst_VelUp1;
% ADCP.VelUp2 = data.Data.Burst_VelUp2;
ADCP.phi1 = data.Data.Burst_Heading  ;
ADCP.phi2 = data.Data.Burst_Pitch;
ADCP.phi3 = data.Data.Burst_Roll;
ADCP.Time = data.Data.Burst_Time;  
% ADCP.AmbVel = data.Data.Burst_AmbiguityVel;
ADCP.Press = data.Data.Burst_Pressure;
ADCP.a1 = data.Data.Burst_AmpBeam1;
ADCP.a2 = data.Data.Burst_AmpBeam2;
ADCP.a3 = data.Data.Burst_AmpBeam3;
ADCP.a4 = data.Data.Burst_AmpBeam4;
ADCP.a5 = data.Data.IBurst_AmpBeam5;
ADCP.c1 = data.Data.Burst_CorBeam1;
ADCP.c2 = data.Data.Burst_CorBeam2;
ADCP.c3 = data.Data.Burst_CorBeam3;
ADCP.c4 = data.Data.Burst_CorBeam4;
ADCP.c5 = data.Data.IBurst_CorBeam5;
% ADCP.Temp = data.Data.Burst_Temperature;
% ADCP.SoundV = data.Data.Burst_Soundspeed;
ADCP.Range = ADCP_CFG.Bd+0.3:ADCP_CFG.Cs:ADCP_CFG.Cs.*(ADCP_CFG.Cn)+0.3;
clearvars data

%%
for k = 2 : length(ListFiles)
  baseFileName = ListFiles(k).name;
  fullFileName = fullfile([ListFolders(p).folder '\' num2str(p) '\', baseFileName]);
  fprintf(1, 'Now reading %s\n', fullFileName);
  data = load(fullFileName);  


%% Define Range
ADCP_CFG.Cs = double(data.Config.Burst_CellSize);
ADCP_CFG.Bd = double(data.Config.Burst_BlankingDistance);
ADCP_CFG.Cn = double(data.Config.Burst_NCells);

%% Clear variables


%% Concatenate Variables
fprintf(1, 'Now concatenating... %s\n');
disp(' ')
fnumber = length(ListFiles);

ADCP.Speed = [ADCP.Speed; data.Data.Burst_VelSpeed];
ADCP.VelDirection = [ADCP.VelDirection; data.Data.Burst_VelDirection];
ADCP.b1 = [ADCP.b1; data.Data.Burst_VelBeam1];
ADCP.b2 = [ADCP.b2; data.Data.Burst_VelBeam2];
ADCP.b3 = [ADCP.b3; data.Data.Burst_VelBeam3];
ADCP.b4 = [ADCP.b4; data.Data.Burst_VelBeam4];
ADCP.b5 = [ADCP.b5; data.Data.IBurst_VelBeam5];
ADCP.VelN = [ADCP.VelN; data.Data.Burst_VelNorth];
ADCP.VelE = [ADCP.VelE; data.Data.Burst_VelEast];
ADCP.VelX = [ADCP.VelX; data.Data.Burst_VelX];
ADCP.VelY = [ADCP.VelY; data.Data.Burst_VelY];
ADCP.VelUp1 = [ADCP.VelUp1; data.Data.Burst_VelUp1];
% ADCP.VelUp2 = [ADCP.VelUp2; data.Data.Burst_VelUp2];
ADCP.phi1 = [ADCP.phi1; data.Data.Burst_Heading];
ADCP.phi2 = [ADCP.phi2; data.Data.Burst_Pitch];
ADCP.phi3 = [ADCP.phi3; data.Data.Burst_Roll];
ADCP.Time = [ADCP.Time; data.Data.Burst_Time];
% ADCP.AmbVel = [ADCP.AmbVel; data.Data.Burst_AmbiguityVel];
ADCP.Press = [ADCP.Press; data.Data.Burst_Pressure];
ADCP.a1 = [ADCP.a1; data.Data.Burst_AmpBeam1];
ADCP.a2 = [ADCP.a2; data.Data.Burst_AmpBeam2];
ADCP.a3 = [ADCP.a3; data.Data.Burst_AmpBeam3];
ADCP.a4 = [ADCP.a4; data.Data.Burst_AmpBeam4];
ADCP.a5 = [ADCP.a5; data.Data.IBurst_AmpBeam5];
ADCP.c1 = [ADCP.c1; data.Data.Burst_CorBeam1];
ADCP.c2 = [ADCP.c2; data.Data.Burst_CorBeam2];
ADCP.c3 = [ADCP.c3; data.Data.Burst_CorBeam3];
ADCP.c4 = [ADCP.c4; data.Data.Burst_CorBeam4];
ADCP.c5 = [ADCP.c5; data.Data.IBurst_CorBeam5];
% ADCP.Temp = [ADCP.Temp; data.Data.Burst_Temperature];
% ADCP.SoundV = [ADCP.SoundV; data.Data.Burst_Soundspeed];

clear data
end
%% Add Start and End time of measurement to ADCP structure
ADCP.Begin = datestr(ADCP.Time(1));
ADCP.End = datestr(ADCP.Time(end));

%% Clear unneeded variables.
clearvars data 

%% Save to savepath as .mat
fprintf(1, 'Now saving file PNCU%s.mat...', num2str(p));
disp(' ')

save([savepath 'PNCU' num2str(p) '.mat'],'ADCP','-v7.3');

%% Clean
clear ADCP data


fprintf(1, 'Done! Saved %s file(s).', num2str(p));