%% All user input should be completed in this top section

%%% Change this variable to wherever your TRUST processing directory is
trustDir='/Users/spencerwaddle/Documents/PPE/Scanning/20210301_Donahue_241665';
mrid = '20191102Donahue_239005';
hct = 0.42; %%% we should have an hct value for everybody
trustflip=1;

dataDir=[trustDir,'/Acquired'];

%%% USER NOTES:
%%% This program will identify the TRUST file based on a few key words in
%%% the trust source file. A few key words include SOURCE, and TRUST_VEIN,
%%% as well as the user input mrid. I suggest naming your trust files
%%% similarly to the exmample in the Data directory, like this:
%%% Donahue_139953_WIP_SOURCE_-_MJD_TRUST_VEIN.PAR
%%% To avoid any issues with the program identifying the TRUST file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(trustDir))

%% TRUST Processing
cd(trustDir);

[trust_T2,HbAA_Yv,HbF_Yv,Bovine_Yv,HbS_Yv] = Process_TRUST_SA_v3(mrid,hct,'both',trustDir,dataDir,trustflip); %% Options are both, first, second
trust_reply = questdlg('Does the decay curve look okay?','Saving TRUST Results','Yes','No','Yes');
close all;

%% Write to file

% userinp=input(['Write processing file for ',mrid,'? y or n', newline],'s');
userinp='y';

if userinp=='y'
    dataFile=fopen([dataDir,'/',mrid,'_TRUSTresults.txt'],'w');
    try
        fwrite(dataFile,['TRUST T2(1), ',num2str(trust_T2(1)),', TRUST T2(2), ',num2str(trust_T2(2)),',  Mean T2:, ',num2str(mean(trust_T2)),newline])
        fwrite(dataFile,['HbAA Yv(1), ',num2str(HbAA_Yv(1)),', HbAA Yv(2), ',num2str(HbAA_Yv(2)),',  Mean HbAA YV:, ',num2str(mean(HbAA_Yv)),newline])
        fwrite(dataFile,['HbF Yv(1), ',num2str(HbF_Yv(1)),', HbF Yv(2), ',num2str(HbF_Yv(2)),',  Mean HbF YV:, ',num2str(mean(HbF_Yv)),newline])
        fwrite(dataFile,['Bovine Yv(1), ',num2str(Bovine_Yv(1)),', Bovine Yv(2), ',num2str(Bovine_Yv(2)),',  Mean Bovine YV:, ',num2str(mean(Bovine_Yv)),newline])
        fwrite(dataFile,['HbS Yv(1), ',num2str(HbS_Yv(1)),', HbS Yv(2), ',num2str(HbS_Yv(2)),',  Mean HbS YV:, ',num2str(mean(HbS_Yv)),newline])
        
    end
    fclose(dataFile);
    disp('File was written')
else
    disp('File was NOT written')
end