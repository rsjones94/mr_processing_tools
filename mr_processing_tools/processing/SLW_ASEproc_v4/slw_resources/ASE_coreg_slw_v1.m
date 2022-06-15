function [ coreg_doneflag ] = ASE_coreg_slw_v1( input, output_nii, reference, omat, omat_status, dcm2nii64_path, flirt_path, flirt_options, bet2_path, bet_input_options,bet_ref_options, bet_input, bet_reference )
%input the image we are coregistering. Will convert a .PAR to a .nii if necessary
%output_nii the name of the coregistered output. intermediate coregistrations will be saved in the same directory
%reference the image we are coregistering to, will convert a .PAR to a .nii if necessary
%omat the name of the omat, or 'none' if no omat is involved
%omat_status, either 'none', 'use', 'save', whether the input omat name will be used in the coregistration, or saved for use later
%dcm2nii64_path is the path to the file format converter
%flirt_path, path to flirt binary
% flirt_options can be an empty string, in which case flirt defaults will be used flirt_options may be -dof, for example
% bet2_path is the path to bet2
% bet_options can include -f, for example, or it can be '' for defaults
% bet2_input and bet2_reference are binary, 0 or 1. 0 If we should skip brain extraction for that volume
% will not re-bet if a brain extracted file exists

%%%%%%% Put \ before ( and ), otherwise, program can't process\
input=nameParenthesisPreslasher(input);
reference=nameParenthesisPreslasher(reference);
output_nii=nameParenthesisPreslasher(output_nii);
%%%%%%%%%%%%%%%%%5


[~,outdir]=system(['dirname ', output_nii]); outdir=outdir(1:end-1);
coregdir=[outdir,'/coreg_steps']; mkdir(coregdir);
tempdir=[coregdir,'/tempqieudncmcm'];

%%%%%%%%%% convert input and reference images if needed from .par to .nii.gz
if strcmpi(input(end-3:end),'.par')
    [~,basename1]=fileparts(input);
    [~,exists]=system(['find ',coregdir,' -maxdepth 1 -name ',basename1,'.nii.gz']);
    if isempty(exists)
        mkdir(tempdir)
        system([dcm2nii64_path,' -i n -e n -d n -g y -4 y -v n -p n -f y -o ',tempdir,' ',input])
        outputname=dir([tempdir,'/*.nii.gz']);
        outputname=outputname(1);
        [~,outputname]=fileparts(outputname.name);
        system(['mv ',tempdir,'/',outputname,'.gz ',coregdir,'/',basename1,'.nii.gz'])
        rmdir(tempdir,'s')
    end
    input=[coregdir,'/',basename1,'.nii.gz'];
end
reference
if strcmpi(reference(end-3:end),'.par')
    [~,basename1]=fileparts(reference);
    [~,exists]=system(['find ',coregdir,' -maxdepth 1 -name ',basename1,'.nii.gz']);
    if isempty(exists)
        mkdir(tempdir)
        system([dcm2nii64_path,' -i n -e n -d n -g y -4 y -v n -p n -f y -o ',tempdir,' ',reference])
        outputname=dir([tempdir,'/*.nii.gz']);
        outputname=outputname(1);
        [~,outputname]=fileparts(outputname.name);
        system(['mv ',tempdir,'/',outputname,'.gz ',coregdir,'/',basename1,'.nii.gz'])
        rmdir(tempdir,'s')
    end
    reference=[coregdir,'/',basename1,'.nii.gz'];
end
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% gzip input and reference if they aren't already.
% % if strcmpi(input(end-3:end),'.nii')
% %     system(['gzip ',input]);
% %     input=[input,'.gz'];
% % end
% % if strcmpi(reference(end-3:end),'.nii')
% %     system(['gzip ',reference]);
% %     input=[reference,'.gz'];
% % end

%%%%%%%%%%%%%%%% Brain extraction
if bet_input    
    [~,basename,ext]=fileparts(input);
    if strcmpi(basename(end-3:end),'.nii')
        basename=basename(1:end-4);
        ext=['.nii.gz'];
    end
    bet_name=[coregdir,'/',basename,'_BET',ext];
    [~,exists]=system(['find ',coregdir,' -maxdepth 1 -name ',basename,'_BET',ext]);
    if isempty(exists)
        system([bet2_path,' ',input,' ',bet_name,' ',bet_input_options]); %% This triggers FSLOUTPUT error
    end
    input=bet_name;
end

if bet_reference
    [~,basename,ext]=fileparts(reference);
    if strcmpi(basename(end-3:end),'.nii')
        basename=basename(1:end-4);
        ext=['.nii.gz'];
    end
    bet_name=[coregdir,'/',basename,'_BET',ext];
    [~,exists]=system(['find ',coregdir,' -maxdepth 1 -name ',basename,'_BET',ext]);
    if isempty(exists)
        system([bet2_path,' ',reference,' ',bet_name,' ',bet_ref_options]);
    end
    reference=bet_name;
end
%%%%%%%%%%%%%%%%%%%

%%%%%%% Do coregistration
if strcmp(omat_status,'none')
    if ~strcmp(omat,'none')
        error('omat_status was set to none, but omat was not input as none. There may be an error. Set omat to none if not')
    end
    system([flirt_path,' ',flirt_options,' -in ',input,' -ref ',reference,' -out ',output_nii])
    
elseif strcmp(omat_status,'use')
    system([flirt_path,' ',flirt_options,' -in ',input,' -ref ',reference,' -applyxfm -init ',omat,' -out ',output_nii])
elseif strcmp(omat_status,'save')
    system([flirt_path,' ',flirt_options,' -in ',input,' -ref ',reference,' -omat ',omat,' -out ',output_nii])
else
    error('omat_status can only be none, use, or save')
end
%%%%%%%%%%%%%% finished coregistration

coreg_doneflag='coregistration has finished';
end

