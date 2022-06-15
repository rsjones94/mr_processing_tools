function [ ASE_doneflag ] = ASEcheck_v7_fn_fornifti(patdir,hctUncorr,sliceGrab,slice2grab,echo2Exclude,tauvec,TR,TE1,TE2,excludedyns_1echo,excludedyns,B0,ASEresolution,RFon,vasoase)
% SLW wrote this to fit ASE data.
% RSAJ modified this to use NiFTis as input and be a part of a larger
% pipeline

%%%%%%%%%%%%%%%%%%%
% cd([patdir,'/Acquired'])
% add all subdirectories to path, c. Benjamin Rinauto
here = mfilename('fullpath');
[path, ~, ~] = fileparts(here);
addpath(genpath(path));

slashes=find(patdir=='/');
patid=patdir(slashes(end)+1:end);

resultsdir=[patdir,'/Results'];
mkdir(resultsdir)

hct=hctUncorr*0.85; %%% microvascular to macrovascular hematocrit ratio is ~0.85
dChi0=0.18;
gamma=267.5; % Gamma is in rad because in "MRI susceptometry: image-based measurement
% of absolute susceptibility of MR contrast agents and human blood." They
% were measuring phase in SWI images (eqn 1

%%%% Grab ASE files

% % %  load_untouch_nii
% taunift1=load_untouch_nii('Jordan_237757ASE_tau0x10then10x2step0.5to20x2_1_mc.nii.gz');
% taumat1=taunift1.img;
% taumat1=rot90(taumat1);

dataFile=fopen([resultsdir,'/',patid,'_PROC.csv'],'w');

fwrite(dataFile,['scanname,B0, HCTuncorrected, sliceGrab,']);
fwrite(dataFile,['slice2grab, size(taumat), tauvec, TR, TE1, TE2, excludedyns_1echo, excludedyns ',newline]);

names=dir([patdir,'/Acquired/*ASE*']);
for namecount=1:length(names)
    tempname=names(namecount).name;
    if strcmpi(tempname(end-6:end),'.nii.gz')
        
        %%% Taumat inds are 
        %%% 1, 2, 3,      4,      5
        %%% x, y, slices, echoes, dynamics
        loadname=[patdir,'/Acquired/',tempname];
        %taumat=readrec_V4_2([patdir,'/Acquired/',tempname]);
        %taumat=readrec_V4_2(loadname);
        origtaumat=niftiread(loadname);
        origtaumatsize=size(origtaumat);
        % unpack the interleaved 4th dimension
        nechoes = 2;
        ndyns = origtaumatsize(4) / nechoes;
        
        taumatsize = [origtaumatsize(1:3) nechoes ndyns];
        
        taumat = zeros(taumatsize);
        
        echo1inds = 1:2:origtaumatsize(4);
        echo2inds = 2:2:origtaumatsize(4);
        
        for i=1:length(echo1inds)
            taumat(:,:,:,1,i)=origtaumat(:,:,:,echo1inds(i));
        end
        for i=1:length(echo2inds)
            taumat(:,:,:,2,i)=origtaumat(:,:,:,echo2inds(i));
        end
            
        
        if taumatsize(3)>1
           threeD=1;
        else
            threeD=0;
        end
        nechoes=taumatsize(4);
        nslices=taumatsize(3);
        ndyn=taumatsize(end);
        
        %%%%%%%%%% SLice grab
        if threeD && sliceGrab
            taumat=taumat(:,:,slice2grab,:,:);
        end
        %%%%%%%%%% Exclude second echo, if set to do so
        if echo2Exclude && nechoes==2
            taumat=taumat(:,:,:,1,:);
            nechoes=1;
            warning('echo 2 is being excluded')
        end
        %%%%%%%%%%%
        
        %%%% Sort tau inds
        [tauvec_sort,sortinds]=sort(tauvec);
        size(taumat)
        taumat=taumat(:,:,:,:,sortinds);
        
        
        %%%% smooth collected ASE images
        %%% Change kernel to define how image is smoothed
        kernel=ones(3);
        kernel=kernel/sum(kernel(:));
        depadLen=(size(kernel,1)-1)/2;
        Imatsmooth=zeros(size(taumat));
        for i=1:size(taumat,3)
            for j=1:nechoes
                for k=1:ndyn
                    tempimg=taumat(:,:,i,j,k);
                    tempimg2=conv2(tempimg,kernel);
                    Imatsmooth(:,:,i,j,k)=tempimg2(depadLen+1:end-depadLen,depadLen+1:end-depadLen);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Adjust cuts to brain mask this was experimentally determined,
        %%% and does not need to be perfect. However, if it is really
        %%% wrong, it could cut parts off of the brain.
        if vasoase
            maxcut=0.99;
            mincut=0.2;
            wmmaxcut=0.94;
            wmmincut=0.8;
        else
            maxcut=0.8;
            mincut=0.2;
            wmmaxcut=0.48;
            wmmincut=0.42;
            if ~threeD || sliceGrab
                wmmaxcut=0.55;
                wmmincut=0.45;
            end
        end
        
%         %%%% Override
%         maxcut=0.99;
%         mincut=0.01;
%         wmmaxcut=0.94;
%         wmmincut=0.8;
%         %%%% END Override
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%
        %%%% Generate brain mask. Note that the mask is generated from the first
        %%%% tau timepoint
        if threeD && ~sliceGrab
            tempimg=Imatsmooth(:,:,:,1);
            maxval=max(tempimg(:));
            roimask=(Imatsmooth(:,:,:,1) < maxval*maxcut) & (Imatsmooth(:,:,:,1) > maxval*mincut);% 3d
            wmmask=(Imatsmooth(:,:,:,1) < maxval*wmmaxcut) & (Imatsmooth(:,:,:,1) > maxval*wmmincut);
        else
            tempimg=Imatsmooth(:,:,1);
            maxval=max(tempimg(:));
            roimask=(Imatsmooth(:,:,1) < maxval*maxcut) & (Imatsmooth(:,:,1) > maxval*mincut); %2d
            wmmask=(Imatsmooth(:,:,1) < maxval*wmmaxcut) & (Imatsmooth(:,:,1) > maxval*wmmincut);
        end
        Imat=Imatsmooth.*roimask;
        
% %         %%%%%%%%%%%%%%%%%%%%%
% %         size(roimask)
% %         imagesc(roimask(:,:,6))
% %         pause
% %         %%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%% Grab B0 map
        
        if isnumeric(B0) & B0~=3
            error('B0 can only be set to 3, if no B0 map is provided')
        end
        B0map=roimask.*3;
        
        %%%%%%%%%%%%%%%%%%%%%% End grab B0
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Full Model Fit! CALCULATE!!! This is where all the time is
        %%%%%% spent
        if threeD && ~sliceGrab
            rOEFmap=zeros(size(taumat(:,:,:,1)));
            rvCBVmap=zeros(size(taumat(:,:,:,1)));
            R2primemap=zeros(size(taumat(:,:,:,1)));
            Rsquaredmap=zeros(size(taumat(:,:,:,1)));
        else
            rOEFmap=zeros(size(taumat(:,:,1)));
            rvCBVmap=zeros(size(taumat(:,:,1)));
            R2primemap=zeros(size(taumat(:,:,1)));
            Rsquaredmap=zeros(size(taumat(:,:,1)));
        end

        findmat=find(roimask);
        maxlen=length(findmat);
        count=0; prog=0; progstore=0; tvec=[];taudata=[];
        for i=1:nechoes
            taudata=[taudata tauvec_sort];
        end
        taudata=taudata/1000; % Must be in seconds. Tau values for first echo, then second echo
        TEvec=[TE1*ones(double(nechoes>=1)*size(tauvec_sort)) ...
                TE2*ones(double(nechoes>=2)*size(tauvec_sort))];
        tic
        %%% Full model fit
        for i=findmat'
            count=count+1;
            Ivec=[];
            [a,b,c]=ind2sub(size(roimask),i); %%% Having c here is fine for 2d because c will always be 1

            for vecind=1:nechoes
                Ivec=[Ivec; squeeze(Imat(a,b,c,vecind,:))];
            end

            %%% Full Yablonskiy fit
            if nechoes==1
                includedyns=1:length(tauvec); 
                includedyns(excludedyns_1echo)=[];
            elseif nechoes==2
                includedyns=1:2*length(tauvec); 
                includedyns(excludedyns)=[];
            else
                error('echoes should be 1 or 2')
            end
            
            %%%%%%%%%% SMOOTH DATA. This makes processing like, 5x
            %%%%%%%%%% longer, and doesn't improve the results at all.
            %%%%%%%%%% Also, if you don't have the filtering package, some
            %%%%%%%%%% of the options in line_filter_slw might not be
            %%%%%%%%%% available to you.
            longtaucheck=taudata>=0.01;
            TE1check=TEvec==TE1;
            Ivec(longtaucheck & TE1check)=line_filter_slw(Ivec(longtaucheck & TE1check)','fft');
            Ivec(longtaucheck & ~TE1check)=line_filter_slw(Ivec(longtaucheck & ~TE1check)','fft');
            Ivec(~longtaucheck & TE1check)=line_filter_slw(Ivec(~longtaucheck & TE1check)','smooth');
            %%%%%%%%%%%% END SMOOTH DATA
            
            for rmnegind=find(Ivec<=0)'
                if rmnegind==1
                    Ivec(rmnegind)=Ivec(2);
                elseif rmnegind==length(Ivec)
                    Ivec(rmnegind)=Ivec(rmnegind-1);
                else
                    Ivec(rmnegind)=Ivec(rmnegind-1);
                end
            end
            
            [ lambda,OEF,otherR,R2prime, Rsquared ] = ASE_YablonskiyFit_v4(taudata(includedyns),TEvec(includedyns),hct,RFon,Ivec(includedyns),B0map(a,b,c)); %%% SV hct corr above
            R2primemap(a,b,c)=R2prime;
            rvCBVmap(a,b,c)=lambda;
            rOEFmap(a,b,c)=OEF;
            Rsquaredmap(a,b,c)=Rsquared;

            prog=floor((count/maxlen)*100);
            if prog~=progstore
                tvec=[tvec toc];
                tic
                t=mean(tvec);
                tsec=floor(t*(100-prog));
                disp(['Progress: ',num2str(prog),'%  Approx. ',num2str(floor(tsec/60)),' minutes and ' ...
                    ,num2str(mod(tsec,60)),' seconds remaining.'])
                progstore=prog;     
            end
        end

        %%%%% Remove bad data
        removeInds=zeros(size(rvCBVmap));
        removeInds=removeInds|rvCBVmap<0|isnan(rvCBVmap)|rOEFmap>1|isinf(rvCBVmap);
        removeInds=removeInds|isinf(rOEFmap)|isnan(rOEFmap)|rOEFmap<0|rOEFmap>1;
%         removeInds=removeInds|R2primemap>9;
        rvCBVmap(removeInds)=0;
        rOEFmap(removeInds)=0;
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% Plot figures
        sliceView=ceil(size(rOEFmap,3)/2);
        clf

        subplot(2,2,1)
        imagesc(rOEFmap(:,:,sliceView))
        title(['OEF (decimal). Mean = ',num2str(mean(rOEFmap(rOEFmap>0)))])
        axis off
        axis square
        colorbar
        colormap jet
        set(gca,'fontsize',12)
        caxis([0 0.4])

        subplot(2,2,2)
        imagesc(rvCBVmap(:,:,sliceView))
        title(['vCBV (decimal). Mean = ',num2str(mean(rvCBVmap(rvCBVmap>0)))])
        axis off
        axis square
        colorbar
        colormap jet
        set(gca,'fontsize',12)
        caxis([0 0.30])

        subplot(2,2,3)
        imagesc(R2primemap(:,:,sliceView))
        title(['R2prime Mean = ',num2str(mean(R2primemap(R2primemap>0)))])
        axis off
        axis square
        colorbar
        colormap jet
        set(gca,'fontsize',12)
        caxis([0 15])

        subplot(2,2,4)
        imagesc(Rsquaredmap(:,:,sliceView))
        title(['Rsquared Mean = ',num2str(mean(Rsquaredmap(Rsquaredmap>0)))])
        axis off
        axis square
        colorbar
        colormap jet
        set(gca,'fontsize',12)
        caxis([0 1])
        
        scanname=tempname(1:end-7);
        saveas(gcf,[resultsdir,'/',scanname,'.png'])
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%% Save .nii
        tempnii=make_nii(rOEFmap,ASEresolution);
        save_nii(tempnii,[resultsdir,'/',scanname,'_rOEF.nii.gz'])
        tempnii=make_nii(rvCBVmap,ASEresolution);
        save_nii(tempnii,[resultsdir,'/',scanname,'_rvCBV.nii.gz'])
        tempnii=make_nii(R2primemap,ASEresolution);
        save_nii(tempnii,[resultsdir,'/',scanname,'_R2prime.nii.gz'])
        tempnii=make_nii(Rsquaredmap,ASEresolution);
        save_nii(tempnii,[resultsdir,'/',scanname,'_Rsquared.nii.gz'])
        
%         save_avw_slw(rot90(rOEFmap,3),[resultsdir,'/',scanname,'_rOEF.nii.gz'],'d',[3.44 3.44 6]);
        %%% Voxel size is [3.44 3.44 6];
        
        %%%%%%%%%%%%%%%%%%%%


%         
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         %%%%% WM values
% %         tempOEF=rOEFmap.*wmmask;
% %         OEFwm= mean(tempOEF(tempOEF>0));
% %         
% %         tempvCBV=rvCBVmap.*wmmask;
% %         vCBVwm= mean(tempvCBV(tempvCBV>0));
% %         
% %         tempR2prime=R2primemap.*wmmask;
% %         R2primewm= mean(tempR2prime(tempR2prime>0));
% %         
% %         tempRsquared=Rsquaredmap.*wmmask;
% %         Rsquaredwm= mean(tempRsquared(tempRsquared>0));
% %         
        fwrite(dataFile,[scanname,',',num2str(B0),',',num2str(hctUncorr),',']);
        fwrite(dataFile,[num2str(sliceGrab),',',num2str(slice2grab),',',num2str(size(taumat)),',',num2str(tauvec),',',num2str(TR),',']);
        fwrite(dataFile,[num2str(TE1),',',num2str(TE2),',',num2str(excludedyns_1echo),',',num2str(excludedyns),newline]);
       

    end
end
% fclose(dataFile);

%%%%%%%%%%%%%%%%%%%%%
% % %% SAVE!!! maps as niftis
% % 
% % %%% 3.44,3.44,4 is voxel size, with the 1mm slice gap added to the z
% % %%% dimension
% % rOEFmaptemp=fliplr(rot90(rOEFmap));
% % tempnii=make_nii(rOEFmaptemp,[220/64,220/64,4]);
% % save_nii(tempnii,[patdir,'/Pat3_OEF_RFoff.nii.gz'])

ASE_doneflag='ASE has finished';
quit
end

