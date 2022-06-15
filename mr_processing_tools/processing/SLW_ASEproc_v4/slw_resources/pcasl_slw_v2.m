function [ doneflag,m0file ] = pcasl_slw_v2( pcasldir,bet2_path, dcm2nii64_path )
%%% SLW wrote this to make a single .m file that processes CBF data. This
%%% uses the simlified model.
%This program uses single compartment model to calculate CBF from
% an input ASL scan. Includes slice-time correction
% Grabs pcasl, M0, LD, and PLD automatically
% Should only be 1 pCASL parrec in pcaslDir.
% This function also saves the CBF nii and a colage of CBF images in
% pcasldir
% also saves a pcasl_readme.txt with details about processing in pcasldir
% Needs SOURCE*pCASL, and _PLDxxxx and _LDxxxx to identify pld and ld
% For example: Donahue_241665_06_02_17.21.33_(WIP_SOURCE_-_MJD_pCASL_PLD1650_LD1525).PAR
% This program does intensity thresholding on m0 to get mask, not BET
% TODO: Parse PAR file to get out voxel information so I can save CBF .nii

%%%%% all from eqn 1 in Recommended implementation... Alsop 2015
slicetime_increment=30; %% 30 ms between each slice acquisition. used for slice time correction
lambda=0.9; %mL/g partition coefficient
T1blood=1650; %% T1 blood at 3T
alpha=0.85; % pcasl labeling efficienct
imgsperwidth=5; % How many images per width will there be in the plot


dirout_asl=dir([pcasldir,'/*SOURCE*pCASL*.PAR']);

dirout_m0=dir([pcasldir,'/*pCASL*M0*.PAR']);

aslfile=dirout_asl(1).name;
m0file=dirout_m0(1).name;

% % %%%%%%%% grab voxel information from PAR file. Not finished. Assuming
% % [3.44 3.44 7]
% % pause
% % partemp=fopen([pcasldir,'/',aslfile],'r');
% % while 1
% %     fgetl(partemp)
% %     pause
% % end
% % fclose(partemp);
% % pause
% % 
% % voxel_size=[]
% % %%%%%%%%

%%%%%%%%%%%%%% Auto grab ld and pld from filename
for i=1:length(aslfile(1:end-4))
    if strcmp(aslfile(i:i+3),'_PLD')
        pldind=i;
    end
    if strcmp(aslfile(i:i+2),'_LD')
        ldind=i;
    end
end
pld=aslfile(pldind+4:pldind+7);
pld=str2num(pld);
if isempty(pld)
    pld=aslfile(pldind+4:pldind+6);
    pld=str2num(pld);
end
ld=aslfile(ldind+3:ldind+6);
ld=str2num(ld);
if isempty(ld)
    ld=aslfile(ldind+3:ldind+5);
    ld=str2num(ld);
end
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Get scan resolution

%%%%%%%%%%%

asl=readrec_V4_2([pcasldir,'/',aslfile]);
m0=readrec_V4_2([pcasldir,'/',m0file]);

%%% Get asl mask.
system([dcm2nii64_path,' -i n -e n -d n -g y -4 y -v n -p n -f y -o ',pcasldir,' ',[pcasldir,'/',m0file]])
m0nii=[pcasldir,'/',m0file(1:end-4),'.nii.gz'];
system([bet2_path,' ',m0nii,' ',[pcasldir,'/tempm0bet'], ' -m -f 0.2']);
m0mask=load_untouch_nii([pcasldir,'/tempm0bet_mask.nii.gz']);
m0mask=double(m0mask.img);
m0mask=rot90(m0mask);
delete([pcasldir,'/tempm0bet_mask.nii.gz']);delete([pcasldir,'/tempm0bet.nii.gz']);

aslsize=size(asl);
dyns=aslsize(5);
slices=aslsize(3);
cbfbig=zeros(aslsize([1:3 5]));
slicetime=0;
for slice=1:slices
    for dyn=1:dyns
        subtimg=asl(:,:,slice,1,dyn,1,1,1)-asl(:,:,slice,1,dyn,1,1,2);
        numer=6e6*lambda*subtimg*exp((pld+slicetime)/T1blood); %%6e6 converts from mL/g/ms to mL/100g/min
        denom=2*alpha*T1blood*m0(:,:,slice)*(1-exp(-ld/T1blood));
        cbfbig(:,:,slice,dyn)=numer./denom;
    end
    slicetime=slicetime+slicetime_increment;
end

cbf=mean(cbfbig,4).*m0mask;
cbf(cbf<0)=0;
cbf(cbf>200)=0;

cbfmax=max(cbf(:));
maxval=10*round(0.5*cbfmax/10);

centervar=(imgsperwidth-ceil(size(cbf,3)/imgsperwidth))*1/5*1/2;%% shifts images up a bit to eliminate white space
mcount=0;
ncount=0;
clf

for i=1:size(cbf,3)
    subplot('Position',[ncount/imgsperwidth mcount/imgsperwidth+centervar 1/imgsperwidth 1/imgsperwidth])
    imagesc(cbf(:,:,i))
    caxis([0 maxval])
    colormap jet
    axis square
    axis off
    ncount=ncount+1;
    if ncount==imgsperwidth
        ncount=0;
        mcount=mcount+1;
    end
end
set(gcf,'Position',[100 100 750 750]);
saveas(gcf,[pcasldir,'/pCASLimg.png']);
close all

tempnii=make_nii(flipud(rot90(cbf,3)),[3.44 3.44 7]);

[~,newfilename]=fileparts(aslfile);
save_nii(tempnii,[pcasldir,'/CBF.nii.gz']);

dataFile=fopen([pcasldir,'/',newfilename,'_readme.txt'],'w');
fwrite(dataFile,[aslfile,newline,'slicetime_increment: ',num2str(slicetime_increment),newline,'partition coefficient: ', num2str(lambda), ...
    newline,'T1blood: ',num2str(T1blood),newline,'labeling efficiency: ', num2str(alpha), ...
    newline,'PLD: ',num2str(pld),newline,'LD: ',num2str(ld), newline, 'caxis (mL/100g/min): ', ...
    num2str([0 maxval])])
fclose(dataFile);

doneflag='pCASL finished processing';
end

