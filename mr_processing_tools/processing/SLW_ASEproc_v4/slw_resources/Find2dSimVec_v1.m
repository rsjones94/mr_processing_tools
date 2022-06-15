function [ ase2dstruct, ase3dstruct, simvec ] = Find2dSimVec_v1( patdir )
%%% SLW wrote this to find which 3D image a 2D image of interest is most
%%% similar to. This helps in coregistering 2D to 3D, because I use the
%%% .omat files from the similar 3D image to coregister the 2D image.
    cd([patdir,'/Acquired'])

    %%%%%%%%%%%% Find which 3d file to coregister each 2d file to
    %%%%%%%%%%%%

    %%%% Find similar 3D files for ss ASE
    ase3dstruct=dir([patdir,'/Acquired/*ms_ASE*.PAR']);
    %%% grab masks from mid slice of each 3d image. This will be compared to
    %%% masks from 2d images to see if any major shifting occured between scans
    ase3dmaskmat=[];
    for i = 1:length(ase3dstruct)
        asemat=readrec_V4_2(ase3dstruct(i).name);
        firstslice=asemat(:,:,ceil(end/2),1,1);
        firstslicemaxnorm=firstslice/max(firstslice(:));
        brainmask=firstslicemaxnorm>0.2;
        ase3dmaskmat(:,:,i) = brainmask;
    end
    %%% Do 2d ASE first, VASOASE next, the process will be different, I think
    %%% Load 2D ase scans, and see which multi-slice scan it is most similar
    %%% to. That is the scan we will coreg each 2D scan to.
    ase2dss=dir([patdir,'/Acquired/*ss_ASE*.PAR']);
    simvecss=[]; % vector of most similar 3d mat for each 2d mat
    for i=1:length(ase2dss)
        asemat=readrec_V4_2(ase2dss(i).name);
        firstslice=asemat(:,:,1,1,1);
        firstslicemaxnorm=firstslice/max(firstslice(:));
        brainmask=firstslicemaxnorm>0.2;
        simvectemp=[]; %%% stores similarities for each 2d vs 3d mat
        for j=1:size(ase3dmaskmat,3)
            diffmat=abs(brainmask-ase3dmaskmat(:,:,j));
            simvectemp(j)=sum(diffmat(:));
        end
        min3dmats=find(simvectemp==min(simvectemp));
        simvecss(i)=min3dmats(end);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Find Similar 3D files for VASOASE
    ase3dstruct=dir([patdir,'/Acquired/*ms_ASE*.PAR']);
    %%% grab masks from mid slice of each 3d image. This will be compared to
    %%% masks from 2d images to see if any major shifting occured between scans
    ase3dmaskmat=[];
    for i = 1:length(ase3dstruct)
        asemat=readrec_V4_2(ase3dstruct(i).name);
        firstslice=asemat(:,:,ceil(end/2),1,1);
        firstslicemaxnorm=firstslice/max(firstslice(:));
        csfmask=firstslicemaxnorm>0.7; %%% Doing CSF because I think it is the thing I can idenfity most reliably
        %%%% between normal ASE and VASOASE
        ase3dmaskmat(:,:,i) = csfmask;
    end

    %%% VASOASE 
    %%% Load 2D ase scans, and see which multi-slice scan it is most similar
    %%% to. That is the scan we will coreg each 2D scan to.
    ase2dvaso=dir([patdir,'/Acquired/*VASOASE*.PAR']);
    simvecvaso=[]; % vector of most similar 3d mat for each 2d mat
    for i=1:length(ase2dvaso)
        asemat=readrec_V4_2(ase2dvaso(i).name);
        firstslice=asemat(:,:,1,1,1);
        firstslicemaxnorm=firstslice/max(firstslice(:));
        brainmask=firstslicemaxnorm>0.75;
        simvectemp=[]; %%% stores similarities for each 2d vs 3d mat
        for j=1:size(ase3dmaskmat,3)
            diffmat=abs(brainmask-ase3dmaskmat(:,:,j));
            simvectemp(j)=sum(diffmat(:));
        end
        min3dmats=find(simvectemp==min(simvectemp));
        simvecvaso(i)=min3dmats(end);
    end

    %%%% simvec gives the index of ase3d that each ase2d image should be
    %%%% coregistered to.
    ase2dstruct=[dir([patdir,'/Acquired/*ss_ASE*.PAR']);dir([patdir,'/Acquired/*VASOASE*.PAR'])];
    simvec=[simvecss simvecvaso];

end

