
matname='/Users/spencerwaddle/Documents/SLWtools/ASEcheck/MNI_valGrab_masks/MNI_2mm_WM_R.nii.gz';
matout='/Users/spencerwaddle/Documents/SLWtools/ASEcheck/MNI_valGrab_masks/MNI_2mm_WM_R_eroded.nii.gz';

mat=load_untouch_nii(matname);
mat=mat.img;

SE=strel('sphere',1);
mat=imerode(mat,SE);

% for i=1:size(mat,3)
%     tempmat=mat(:,:,i);
%     SE=strel(100);
%     mat(:,:,i)=imerode(tempmat,SE);
% end


tempnii=make_nii(flipud(mat),[2 2 2]);
save_nii(tempnii,matout)