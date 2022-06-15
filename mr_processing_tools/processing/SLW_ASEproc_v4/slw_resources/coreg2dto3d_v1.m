function [ finflag ] = coreg2dto3d_v1(ase2dname,ase3dname,t1file,mni_brain,patdir,ASEresolution,dcm2nii64_path, reorient_path,flirt_path,bet2_path)
% This is the function that coregisters 2D images to 3D images
coregstepsdir=[patdir,'/2dCoreg/coreg_steps'];
mkdir(coregstepsdir)

ase2t1omat=[patdir,'/2dCoreg/',ase3dname(1:end-4),'ASE2T1.omat'];
t12mniomat=[patdir,'/2dCoreg/',ase3dname(1:end-4),'T12MNI.omat'];

disp(ase2dname)

%%%% Figure out how many slices we need to make ase2d
ase3d=readrec_V4_2([patdir,'/Acquired/',ase3dname]);
numslices=size(ase3d,3);
%%%%%%%%%%%%%%%

%%%% Make non-coregistered 3d versions of each 2d image
files2d=dir([patdir,'/Results/*',ase2dname(1:end-4),'*.nii.gz']);
mask2dname=[coregstepsdir,'/',files2d(1).name(1:end-15),'_2dmask_2dto3d.nii'];
%%%% Turn each 2d functional image into a 3d image, with 2d as center slice
for i=1:length(files2d)
    newfile3d=[coregstepsdir,'/',files2d(i).name(1:end-7),'_2dto3d.nii'];
    funcmat=load_untouch_nii([patdir,'/Results/',files2d(i).name]);
    newmat=zeros([size(funcmat.img), numslices]);
%     newmat=zeros([size(funcmat.img), numslices])/0;
    newmat(:,:,ceil(end/2))=funcmat.img;
    
    tempnii=make_nii(newmat,ASEresolution);
    save_nii(tempnii,newfile3d);
    
end

%%%%% This will be used to pick out which voxels in 3d space to use.
mask2d=zeros([size(funcmat.img), numslices]);
mask2d(:,:,ceil(end/2))=1; %%% This will be used to pick which voxels in 3d space are most like the ones from 2d
tempnii=make_nii(mask2d,ASEresolution);
save_nii(tempnii,mask2dname)

%%%% Make coregistered 3d versions of each 2d image
files2d=dir([patdir,'/2dCoreg/coreg_steps/*',ase2dname(1:end-4),'*.nii']);
%%%% Turn each 2d functional image into a 3d image, with 2d as center slice
for i=1:length(files2d)
    newfile3d=[patdir,'/2dCoreg/',files2d(i).name(1:end-11),'_3dcoreg.nii'];
    ase2t1_tempfile=[coregstepsdir,'/',files2d(i).name(1:end-11),'_ase2t1.nii'];
    input1=nameParenthesisPreslasher([coregstepsdir,'/',files2d(i).name]);
    output1=nameParenthesisPreslasher(ase2t1_tempfile);
    %%%% Originally did not have nearest neighbor, and instead ran the
    %%%% weighting correction below. Using nearest neighbor fixes this
    %%%% problem more simply though.
    ASE_coreg_slw_v2( input1, output1, nameParenthesisPreslasher(t1file), nameParenthesisPreslasher(ase2t1omat), 'use', dcm2nii64_path, reorient_path, flirt_path, ' -interp nearestneighbour ', bet2_path, '','', 0, 0 )
    input1=output1;
    output1=nameParenthesisPreslasher(newfile3d);
    ASE_coreg_slw_v2( input1, output1, nameParenthesisPreslasher(mni_brain), nameParenthesisPreslasher(t12mniomat), 'use', dcm2nii64_path, reorient_path, flirt_path, ' -interp nearestneighbour ', bet2_path, '','', 0, 0 )
end

%%%% Make sure pseudo3d volume is only 6mm thick, similar to original ASE

%%%% Not each voxel in the psudo3d space corresponds to one voxel in the
%%%% original 2d image, so I go through and pick which pseudo3d voxels to
%%%% use, and correct them based on the 2dmask_3dcoreg file
mask3dcoreg=load_untouch_nii([patdir,'/2dCoreg/',ase2dname(1:end-4),'_2dmask_3dcoreg.nii.gz']);
voxelsize_3d=mask3dcoreg.hdr.dime.pixdim(2:4);
mask3dcoreg=mask3dcoreg.img;
ssmask3dcoreg=zeros(size(mask3dcoreg));
sscorrectmap=zeros(size(mask3dcoreg));
for i=1:size(mask3dcoreg,1)
    for j=1:size(mask3dcoreg,2)
        maxval=0;
        maxind=1;
        for k=1:size(mask3dcoreg,3)
            if mask3dcoreg(i,j,k)>maxval
                maxval=mask3dcoreg(i,j,k);
                maxind=k;
            end
        end
        ssmask3dcoreg(i,j,maxind)=1; %%%% get rid of the below 6 lines to do a probability correction instead of a nearest neighbor fix
        try %%% Putting these in trys because it is possible the volume will hit the edge of the space, and I want to avoid indexing error
            ssmask3dcoreg(i,j,maxind+1)=1;
        end
        try
            ssmask3dcoreg(i,j,maxind-1)=1;
        end
        sscorrectmap(i,j,maxind)=maxval;
    end
end

%%%%% Correct and resave each pseudo3d image
files2d=dir([patdir,'/2dCoreg/*',ase2dname(1:end-4),'*3dcoreg.nii.gz']);
%%%% Turn each 2d functional image into a 3d image, with 2d as center slice
for i=1:length(files2d)
    filename1=files2d(i).name;
%     if ~contains(filename1,'2dmask')
        tempmat=load_untouch_nii(filename1);
        tempmat=tempmat.img;
        tempmat=tempmat.*ssmask3dcoreg;
%         tempmat=tempmat./sscorrectmap; %%% Correct for combining 2d voxels with 0 voxels
        tempmat=flipud(tempmat);
        tempnii=make_nii(tempmat,voxelsize_3d);
        delete([patdir,'/2dCoreg/',filename1])
        save_nii(tempnii,[patdir,'/2dCoreg/',filename1(1:end-3)])
        gzip([patdir,'/2dCoreg/',filename1(1:end-3)])
        delete([patdir,'/2dCoreg/',filename1(1:end-3)])
        
%     end
end


finflag='Finished 2d to 3d';

end

