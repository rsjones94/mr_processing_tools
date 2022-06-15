%%%% You will probably need to run matlab from terminal, so that it can
%%%% have access to all of the fsl tools and functions
function [ master_doneflag ] = ASEmaster_Bioclinica(pt_name)

% addpath('/Users/spencerwaddle/Documents/PPE/Scanning/CodeTools')
% addpath('/Users/spencerwaddle/Documents/PPE/Scripts')
%standaloneDir='/Users/spencerwaddle/Documents/SLWtools/SLW_ASEproc_v4'; %Directory with all ASE processing scripts in it
standaloneDir=mfilename('fullpath');
addpath(genpath(standaloneDir))

slw_resourcesDir=[standaloneDir,'/slw_resources']; %%% directory for slw_resources

%%

%datadir1='/Users/spencerwaddle/Documents/SLWtools/SLW_ASEproc_v4/DataDir';
datadir1=[standaloneDir,'/DataDir'];

% patstruct={'20210803_Donahue_242703', '20210804_Donahue_243000', '20210805_Donahue_242991', ...
%     '20210805_Donahue_243038', '20210809_Donahue_243013', '20210809_Donahue_243015', ...
%     '20210813_Donahue_242993', '20210823_Donahue_243037', '20210823_Donahue_243062', ...
%     '20210824_Donahue_243034'}; %%% Names of all patient directories to process.

patstruct={pt_name}; %%% Names of all patient directories to process.

warning('off','all')

sliceGrab=0; %%% 0 or 1. This won't do anything if threeD is 0
slice2grab=6;
echo2Exclude=0; %%% 0 or 1. Should be 0 if only one echo was collected

TR=4400; %ms
TE1=64; %ms
TE2=107;
hctUncorr=0.42; %% Uncorrected hematocrit. Estimated Value 0.42

%%%% This is the ordering of tau-shiftings in the ASE acquisition. This
%%%% has been standardized. Changing this should be a very carefully made
%%%% decision. You will ruin the ASE processing if this is changed and the
%%%% ASE acquisition isn't altered accordingly
tauvec=[19 20 6.5 8.5 1 16.5 1.5 14.5 3.5 13 9 19.5 17 2.5 ...
  4.5 12.5 2 9.5 6 14 7.5 5.5 12 8 11 0.5 16 10.5 13.5 5 ...
  10 0 18.5 17.5 3 11.5 7 18 15.5 15 4];


% ASEresolution=[3.44,3.44,3]; %%% Gotta include slice gap in 3 dimension size
ASEresolution=[3.44/2,3.44/2,3]; %%% Gotta include slice gap in 3 dimension size

%%% We found that tau>20 has nonlinearity that is bad for ASE, So you can
%%% exclude them here.
excludedyns_1echo=[]; %% Should be the same as exclude dyns, but only the first echo indices
excludedyns=[];%% Exclude [42:46 88:92] to remove nonlinearity at end of scan
                    %%% You won't need this if you keep tau <= 20ms

%%%%%%%%%%%%%%%%%%%%%%
mni_brain=[slw_resourcesDir,'/atlases/MNI152_T1_2mm_brain.nii.gz'];
MNIdir=[standaloneDir,'/MNI_valGrab_sky/MNI_valGrab_masks'];
dcm2nii64_path=[slw_resourcesDir,'/dcm2nii64'];
dcm2niix_path=[slw_resourcesDir,'/dcm2niix'];
flirt_path=[slw_resourcesDir,'/FSLTOOLS/flirt'];
bet_path=[slw_resourcesDir,'/FSLTOOLS/bet'];
reorient_path=[slw_resourcesDir,'/FSLTOOLS/fslreorient2std'];
%%%%

for participantInd=1:length(patstruct)

    %%% Define parameters
%     patdir='/Users/spencerwaddle/Documents/PPE/Scanning/20210701_Donahue_242810'; %%% patdir should have an Acquired folder with all scans in it
    patdir=[datadir1,'/', patstruct{participantInd}];
    [~,patid]=system(['basename ',patdir]);patid=patid(1:end-1);

    mkdir([patdir,'/2dCoreg'])
    mkdir([patdir,'/Results'])

    %%% Find B0 file
    B0=dir([patdir,'/Acquired/*B0*.PAR']);
    try
        B0=[patdir,'/Acquired/',B0(1).name];
    end
    if isempty(B0)
        B0=3;
    end
    %%%%%%%%%% USE B0=3 Tesla
    B0=3;
    %%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ASE
    
    ASEcheck_v7_fn( patdir, hctUncorr,sliceGrab,slice2grab,echo2Exclude,tauvec,TR,TE1,TE2,excludedyns_1echo,excludedyns,B0,dcm2nii64_path,flirt_path,ASEresolution)
    disp('finished ASEcheck')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Begin TRUST processing.
    
    mkdir([patdir,'/Acquired/TRUSTstore_movetoacquired'])
    source = dir([patdir,'/Acquired/*SOURCE*TRUST*.PAR']);
    
    %%%% Move files away to process each one individually. This is the
    %%%% easiest way to do it, because the TRUST script is only designed to
    %%%% process one file at a time
    for i=1:length(source)
        movefile([patdir,'/Acquired/',source(i).name],[patdir,'/Acquired/TRUSTstore_movetoacquired'])
        movefile([patdir,'/Acquired/',source(i).name(1:end-4),'.REC'],[patdir,'/Acquired/TRUSTstore_movetoacquired'])
    end
    
    for i=1:length(source)
        %%% Move file in for processing
        movefile([patdir,'/Acquired/TRUSTstore_movetoacquired/',source(i).name],[patdir,'/Acquired'])
        movefile([patdir,'/Acquired/TRUSTstore_movetoacquired/',source(i).name(1:end-4),'.REC'],[patdir,'/Acquired'])
        
        Master_SA_v2_fn(patdir,hctUncorr)
        try
            movefile([patdir,'/Acquired/TRUSTresults.txt'],[patdir,'/Results/TRUSTresults',num2str(i),'.txt'])
        end
        
        %%%% Move file back out to process next trust file
        movefile([patdir,'/Acquired/',source(i).name],[patdir,'/Acquired/TRUSTstore_movetoacquired'])
        movefile([patdir,'/Acquired/',source(i).name(1:end-4),'.REC'],[patdir,'/Acquired/TRUSTstore_movetoacquired'])
    end
        
    %%%% Move files back for storage
    for i=1:length(source)
        movefile([patdir,'/Acquired/TRUSTstore_movetoacquired/',source(i).name],[patdir,'/Acquired'])
        movefile([patdir,'/Acquired/TRUSTstore_movetoacquired/',source(i).name(1:end-4),'.REC'],[patdir,'/Acquired'])
    end
    rmdir([patdir,'/Acquired/TRUSTstore_movetoacquired'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end TRUST processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin pCASL Processing
    aslfiles=dir([patdir,'/Acquired/*SOURCE*pCASL*.PAR']);
    m0files=dir([patdir,'/Acquired/*pCASL*M0*.PAR']);
    if length(aslfiles) ~= length(m0files)
        error('there arent the same number of asl and m0 files. There is a probably a problem, resolve this to continue processing')
    end
    %%%% Move ASL files out of acquired. That is how we process them
    mkdir([patdir,'/Acquired/ASLstore_moveToAcquired'])
    for i=1:length(aslfiles)
        movefile([patdir,'/Acquired/',aslfiles(i).name],[patdir,'/Acquired/ASLstore_moveToAcquired/',aslfiles(i).name])
        movefile([patdir,'/Acquired/',aslfiles(i).name(1:end-4),'.REC'],[patdir,'/Acquired/ASLstore_moveToAcquired/',aslfiles(i).name(1:end-4),'.REC'])
        movefile([patdir,'/Acquired/',m0files(i).name],[patdir,'/Acquired/ASLstore_moveToAcquired/',m0files(i).name])
        movefile([patdir,'/Acquired/',m0files(i).name(1:end-4),'.REC'],[patdir,'/Acquired/ASLstore_moveToAcquired/',m0files(i).name(1:end-4),'.REC'])
    end
    for i=1:length(aslfiles)
        aslname=aslfiles(i).name(1:end-4);
        movefile([patdir,'/Acquired/ASLstore_moveToAcquired/',aslfiles(i).name],[patdir,'/Acquired/',aslfiles(i).name])
        movefile([patdir,'/Acquired/ASLstore_moveToAcquired/',aslname,'.REC'],[patdir,'/Acquired/',aslname,'.REC'])
        movefile([patdir,'/Acquired/ASLstore_moveToAcquired/',m0files(i).name],[patdir,'/Acquired/',m0files(i).name])
        movefile([patdir,'/Acquired/ASLstore_moveToAcquired/',m0files(i).name(1:end-4),'.REC'],[patdir,'/Acquired/',m0files(i).name(1:end-4),'.REC'])
        %%% Process pCASL. This will get confused if there are 2 pCASL scans in
        %%% patdir
        [~,M0file]=pcasl_slw_v2([patdir,'/Acquired'],bet_path, dcm2nii64_path);
        M0file=[patdir,'/Acquired/',M0file];
        %%%%%%

        %%%%%%%%%%% move output of processing to results
        try
        movefile([patdir,'/Acquired/',aslname,'_readme.txt'],[patdir,'/Results/',aslname,'_readme.txt'])
        end
        try
        movefile([patdir,'/Acquired/pCASLimg.png'],[patdir,'/Results/',aslname,'pCASLimg.png'])
        end
        try
        movefile([patdir,'/Acquired/CBF.nii.gz'],[patdir,'/Results/',aslname,'_CBF.nii.gz'])
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% Begin Coregistration
        coregdir=[patdir,'/Results/coreg_junk'];
        coregResults=[patdir,'/Results/coregResults'];
        T1file=dir([patdir,'/Acquired/*3DT1*.PAR']);
        if isempty(T1file)
            T1file=dir([patdir,'/Acquired/*3D_T1*.PAR']);
        end
        T1file=T1file.name; T1file=[patdir,'/Acquired/',T1file];
        
        %%%%%%%%%%%%% CBF to MNI
        %%% M0 to T1
        ASE_coreg_slw_v2( M0file, [coregdir,'/M02T1.nii.gz'], T1file, [coregdir,'/M02T1.omat'], 'save', dcm2nii64_path, reorient_path, flirt_path, '-dof 6', bet_path, '-f 0.2','-B', 1, 1 )
        %%% T1 to MNI
        ASE_coreg_slw_v2( T1file, [coregdir,'/T12MNI.nii.gz'], mni_brain, [coregdir,'/T12MNI.omat'], 'save', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '-B','', 1, 0 )
        %%% CBF to T1
        ASE_coreg_slw_v2( [patdir,'/Results/',aslname,'_CBF.nii.gz'], [coregdir,'/',aslname,'_CBF2T1.nii.gz'], T1file, [coregdir,'/M02T1.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        %%% CBF to MNI
        ASE_coreg_slw_v2( [coregdir,'/',aslname,'_CBF2T1.nii.gz'], [coregResults,'/',aslname,'_CBF2MNI.nii.gz'], mni_brain, [coregdir,'/T12MNI.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        % % % % % % % % MNI_valGrab_v2(inFile, roidir, fxnHandles, csvname)
        MNI_valGrab_v2([coregResults,'/',aslname,'_CBF2MNI.nii.gz'],MNIdir,{@mean},[coregResults,'/',patid,'_ProcResults.csv'],patid,'CBF',aslname,'mean');
        MNI_valGrab_v2([coregResults,'/',aslname,'_CBF2MNI.nii.gz'],MNIdir,{@std},[coregResults,'/',patid,'_ProcResults.csv'],patid,'CBF',aslname,'STD');
        movefile([patdir,'/Acquired/',aslfiles(i).name],[patdir,'/Acquired/ASLstore_moveToAcquired/',aslfiles(i).name])
        movefile([patdir,'/Acquired/',aslname,'.REC'],[patdir,'/Acquired/ASLstore_moveToAcquired/',aslname,'.REC'])
        movefile([patdir,'/Acquired/',m0files(i).name],[patdir,'/Acquired/ASLstore_moveToAcquired/',m0files(i).name])
        movefile([patdir,'/Acquired/',m0files(i).name(1:end-4),'.REC'],[patdir,'/Acquired/ASLstore_moveToAcquired/',m0files(i).name(1:end-4),'.REC'])
    end


    %%% Move ASL files back to acquired
    for i=1:length(aslfiles)
        movefile([patdir,'/Acquired/ASLstore_moveToAcquired/',aslfiles(i).name],[patdir,'/Acquired/',aslfiles(i).name])
        movefile([patdir,'/Acquired/ASLstore_moveToAcquired/',aslfiles(i).name(1:end-4),'.REC'],[patdir,'/Acquired/',aslfiles(i).name(1:end-4),'.REC'])
        movefile([patdir,'/Acquired/ASLstore_moveToAcquired/',m0files(i).name],[patdir,'/Acquired/',m0files(i).name])
        movefile([patdir,'/Acquired/ASLstore_moveToAcquired/',m0files(i).name(1:end-4),'.REC'],[patdir,'/Acquired/',m0files(i).name(1:end-4),'.REC'])
    end
    
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End pCASL processing
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% OEF 2 MNI
    coregdir=[patdir,'/Results/coreg_junk'];
    coregResults=[patdir,'/Results/coregResults'];
    T1file=dir([patdir,'/Acquired/*3DT1*.PAR']);
    if isempty(T1file)
        T1file=dir([patdir,'/Acquired/*3D_T1*.PAR']);
    end
    T1file=T1file.name; T1file=[patdir,'/Acquired/',T1file];
    coregdir=[patdir,'/Results/coreg_junk'];
    coregResults=[patdir,'/Results/coregResults'];
    T1file=dir([patdir,'/Acquired/*3DT1*.PAR']);
    if isempty(T1file)
        T1file=dir([patdir,'/Acquired/*3D_T1*.PAR']);
    end
    T1file=T1file.name; T1file=[patdir,'/Acquired/',T1file];
    
    %%%%%% Multislice first
    names=dir([patdir,'/Acquired/*ms_ASE*.PAR']);
    system('. /usr/local/fsl/etc/fslconf/fsl.sh')
    %%%% Find which 3d matrices we will coregister each 2d acquisition to
    [ ase2dstruct, ase3dstruct, simvec ] = Find2dSimVec_v1( patdir );
    %%%%%%%%%%%%%%%%%%%%%
    for i=1:length(names)
        try
            rmdir(coregdir,'s')
        end
        tempname=names(i).name;
        scanname=tempname(1:end-4);
        %%% ASE to T1
        ASE_coreg_slw_v2( [patdir,'/Acquired/',tempname], [coregdir,'/ASEn2T1.nii.gz'], T1file, [coregdir,'/ASEn2T1.omat'], 'save', dcm2nii64_path, reorient_path, flirt_path, '-dof 6', bet_path, '-f 0.2','-B', 1, 1 )
        %%% T1 to MNI
        ASE_coreg_slw_v2( T1file, [coregdir,'/T12MNI.nii.gz'], mni_brain, [coregdir,'/T12MNI.omat'], 'save', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '-B','', 1, 0 )
        %%% R2prime to T1
        ASE_coreg_slw_v2( [patdir,'/Results/',scanname,'_R2prime.nii.gz'], [coregdir,'/R2prime2T1.nii.gz'], T1file, [coregdir,'/ASEn2T1.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        %%% R2prime to MNI
        ASE_coreg_slw_v2( [coregdir,'/R2prime2T1.nii.gz'], [coregResults,'/',scanname,'_R2prime_coreg.nii.gz'], mni_brain, [coregdir,'/T12MNI.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        %%% rOEF to T1
        ASE_coreg_slw_v2( [patdir,'/Results/',scanname,'_rOEF.nii.gz'], [coregdir,'/rOEF2T1.nii.gz'], T1file, [coregdir,'/ASEn2T1.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        %%% rOEF to MNI
        ASE_coreg_slw_v2( [coregdir,'/rOEF2T1.nii.gz'], [coregResults,'/',scanname,'_rOEF_coreg.nii.gz'], mni_brain, [coregdir,'/T12MNI.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        %%% rvCBV to T1
        ASE_coreg_slw_v2( [patdir,'/Results/',scanname,'_rvCBV.nii.gz'], [coregdir,'/rvCBV2T1.nii.gz'], T1file, [coregdir,'/ASEn2T1.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        %%% rvCBV to MNI
        ASE_coreg_slw_v2( [coregdir,'/rvCBV2T1.nii.gz'], [coregResults,'/',scanname,'_rvCBV_coreg.nii.gz'], mni_brain, [coregdir,'/T12MNI.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        %%% Rsquared to T1
        ASE_coreg_slw_v2( [patdir,'/Results/',scanname,'_Rsquared.nii.gz'], [coregdir,'/Rsquared2T1.nii.gz'], T1file, [coregdir,'/ASEn2T1.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        %%% Rsquared to MNI
        ASE_coreg_slw_v2( [coregdir,'/Rsquared2T1.nii.gz'], [coregResults,'/',scanname,'_Rsquared_coreg.nii.gz'], mni_brain, [coregdir,'/T12MNI.omat'], 'use', dcm2nii64_path, reorient_path, flirt_path, '', bet_path, '','', 0, 0 )
        %%%% Grab regional mean values
        MNI_valGrab_v2([coregResults,'/',scanname,'_R2prime_coreg.nii.gz'],MNIdir,{@mean},[coregResults,'/',patid,'_ProcResults.csv'],patid,'R2prime',scanname,'mean');
        MNI_valGrab_v2([coregResults,'/',scanname,'_rOEF_coreg.nii.gz'],MNIdir,{@mean},[coregResults,'/',patid,'_ProcResults.csv'],patid,'rOEF',scanname,'mean');
        MNI_valGrab_v2([coregResults,'/',scanname,'_rvCBV_coreg.nii.gz'],MNIdir,{@mean},[coregResults,'/',patid,'_ProcResults.csv'],patid,'rvCBV',scanname,'mean');
        MNI_valGrab_v2([coregResults,'/',scanname,'_Rsquared_coreg.nii.gz'],MNIdir,{@mean},[coregResults,'/',patid,'_ProcResults.csv'],patid,'Rsquared',scanname,'mean');
        %%% Save omats for later if we are going to use them to coregister 2d
        %%% Otherwise they get overwritten in each coregistration
        for j=1:length(simvec)
            if strcmp([scanname,'.PAR'], ase3dstruct(simvec(j)).name)
                copyfile([coregdir,'/ASEn2T1.omat'],[patdir,'/2dCoreg/',scanname,'ASE2T1.omat'])
                copyfile([coregdir,'/T12MNI.omat'],[patdir,'/2dCoreg/',scanname,'T12MNI.omat'])          
            end
        end
    end
    
    %%%%% Single slice second
    % % % % % % [ ase2dstruct, ase3dstruct, simvec ] = Find2dSimVec_v1( patdir );
    for i=1:length(ase2dstruct)
        ASE2dcoreg=coreg2dto3d_v1(ase2dstruct(i).name, ase3dstruct(simvec(i)).name, T1file, mni_brain,patdir,ASEresolution,dcm2nii64_path, reorient_path,flirt_path,bet_path);
        tempname=ase2dstruct(i).name(1:end-4);
        MNI_valGrab_v2([patdir,'/2dCoreg/',tempname,'_R2prime_3dcoreg.nii.gz'],MNIdir,{@mean},[coregResults,'/',patid,'_ProcResults.csv'],patid,'R2prime',ase2dstruct(i).name,'mean');
        MNI_valGrab_v2([patdir,'/2dCoreg/',tempname,'_rOEF_3dcoreg.nii.gz'],MNIdir,{@mean},[coregResults,'/',patid,'_ProcResults.csv'],patid,'rOEF',ase2dstruct(i).name,'mean');
        MNI_valGrab_v2([patdir,'/2dCoreg/',tempname,'_rvCBV_3dcoreg.nii.gz'],MNIdir,{@mean},[coregResults,'/',patid,'_ProcResults.csv'],patid,'rvCBV',ase2dstruct(i).name,'mean');
        MNI_valGrab_v2([patdir,'/2dCoreg/',tempname,'_Rsquared_3dcoreg.nii.gz'],MNIdir,{@mean},[coregResults,'/',patid,'_ProcResults.csv'],patid,'Rsquared',ase2dstruct(i).name,'mean');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Coregistration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
master_doneflag='ASE (master) has finished';