function [trust_T2,HbAA_Yv,HbF_Yv,Bovine_Yv,HbS_Yv] = Process_TRUST_SA_v3(mrid,hct,trust_op,trustDir,dataDir,trustflip)

addpath(genpath(trustDir))

% Determine release type
source = dir([dataDir,'/*SOURCE*TRUST*']);
if isempty(source)
    r53 = dir([dataDir,'/*R53*TRUST*']);
    if isempty(r53)
        trust.version = 'old';
    else
        trust.version = 'new';
    end
else
    trust.version = 'new';
end


% Extract TRUST inputs.
trust.whichavg = trust_op;
    switch trust.whichavg
        case 'first'
            trust.scan=[1];
        case 'second'
            trust.scan=[2];
        case 'both'
            trust.scan=[1 2];
    end
    
% trust_T2=[0 0];
% Bovine_Yv=[0 0];

% Loop over number of scans to process
for numscan = trust.scan
    switch trust.version
        case 'old'
            % Input TRUST Data
            trust.list = dir([dataDir,'/*TRUST_VEIN*REC*']);
            if isempty(trust.list)
                trust_T2(numscan)=''; 
                Bovine_Yv(numscan)='';
                Hbf_Yv(numscan)='';
                HbAA_Yv(numscan)='';
                HbS_Yv(numscan)='';
                continue; 
            end;
            [trust.data, handles.Header] = readrec_V4_2([dataDir,'/',trust.list(numscan).name]);
            [trust.nRows,trust.nColumns,trust.nSlices,trust.nDynamics,trust.nEchoes] = size(trust.data);

            
            % Gather Parameters
            trust.eTE = [0 40 80 160]; %Effective TEs in ms
            trust.nTE = length(trust.eTE); %Number of eTEs
            trust.Measures = trust.nEchoes/(trust.nTE*2); %Number of Averages

            %Order of Control & Label Images
            trust.nLabel = [1:2:trust.nEchoes];
            trust.nControl = [2:2:trust.nEchoes];

            % Average Across Measurements
            counter = 0;
            for i = 1:trust.nTE
                for j = 1:trust.Measures
                    counter = counter + 1;
                    interm_label(:,:,:,:,j) = trust.data(:,:,:,:,trust.nLabel(counter));
                    interm_control(:,:,:,:,j) = trust.data(:,:,:,:,trust.nControl(counter));
                end
                trust.avgLabel(:,:,:,:,i) = sum(interm_label(:,:,:,:,:),5)/trust.Measures;
                trust.avgControl(:,:,:,:,i) = sum(interm_control(:,:,:,:,:),5)/trust.Measures;
            end
        case 'new'
            % Input TRUST Data
            disp('here')
            trust.list = dir([dataDir,'/*SOURCE*TRUST_VEIN*REC*']);
               
            if isempty(trust.list)
                trust_T2(numscan)=0; 
                Bovine_Yv(numscan)=0;
                Hbf_Yv(numscan)=0;
                HbAA_Yv(numscan)=0;
                HbS_Yv(numscan)=0;
                continue; 
            end;
            
            [trust.data, handles.Header] = readrec_V4_2([dataDir,'/',trust.list(1).name]);
            [trust.nRows,handles.trust.nColumns,trust.nSlices,trust.nDynamics,trust.nEchoes,~,~,trust.nRepeats] = size(trust.data);
            % Gather Parameters
            trust.eTE = [0 40 80 160]; %Effective TEs in ms
            trust.nTE = length(trust.eTE); %Number of eTEs
            trust.Measures = trust.nEchoes/(trust.nTE*2); %Number of Averages
            %Order of Control & Label Images
            trust.nLabel = [2:2:trust.nEchoes];
            trust.nControl = [1:2:trust.nEchoes];

            % Average Across Measurements
            counter = 0;
            for i = 1:trust.nTE
                for j = 1:trust.Measures
                    counter = counter + 1;
                    interm_label(:,:,:,:,j) = trust.data(:,:,:,:,trust.nLabel(counter),:,:,numscan);
                    interm_control(:,:,:,:,j) = trust.data(:,:,:,:,trust.nControl(counter),:,:,numscan);
                end
                trust.avgLabel(:,:,:,:,i) = sum(interm_label(:,:,:,:,:),5)/trust.Measures;
                trust.avgControl(:,:,:,:,i) = sum(interm_control(:,:,:,:,:),5)/trust.Measures;
            end
    end
    % Subtract Label from Control
    trust.Diff = trust.avgControl - trust.avgLabel;
    
    diffmat=trust.Diff(:,:,1,1,1,1);
    meandiff=mean(diffmat(:));
    
    if meandiff<0
        trust.Diff=-trust.Diff;
    end
    
%     if trustflip
%         trust.Diff=-trust.Diff;
%     end

    %%%%% Creat SSS pre-mask Automated
    xlen=size(trust.avgLabel,1);
    ylen=size(trust.avgLabel,2);
    maskmat=zeros(xlen,ylen);
    maskmat(round(ylen/2):end,round(xlen/4):round(3*xlen/4))=1;
    trust.BW = maskmat;
    %%%%%%%%%%%%%%%%
    
%     %%%%%%%%%% Create trust pre-mask manual
%     imagesc(trust.Diff(:,:,1,1,1,1))
%     colormap gray
%     trust.BW = roipoly();

    % Isolate Sagittal Sinus with SSS pre-mask
    for i = 1:trust.nTE
        trust.SagSinus =squeeze(trust.Diff(:,:,:,:,i)).*trust.BW;
        temp = sort(trust.SagSinus(:),1,'descend');
        trust.MeanSagSinus(i) = mean(temp(1:4));
    end

    % Blood T1
    trust.T1Blood = 1./(0.52*hct+0.38); %s - Jordan et al.

    % Fitting Parameters (tau CPMG = 10 ms) - Lu et al.
    trust.a1 = -13.5; trust.a2 = 80.2; trust.a3 = -75.9;
    trust.b1 = -0.5; trust.b2 = 3.4;
    trust.c1 = 247.4;

    % Fit to T2-Decay Curve
    X = trust.eTE/1000;
    Y = trust.MeanSagSinus;

    params = lsqnonlin(@fit_simp,[0.100 1000000],[],[],[],X,Y);

    trust.C = params(1);
    trust.S0 = params(2);
    % Compute Blood T2
    trust.T2Blood = 1/((1/trust.T1Blood)-trust.C);
    trust_T2(numscan)=trust.T2Blood;
    HbF_Yv(numscan)=lu_model(trust.T2Blood,hct,10,2);
    Bovine_Yv(numscan)=lu_model(trust.T2Blood,hct,10,1);
    HbS_Yv(numscan)=wood_model(trust.T2Blood,hct,10,2);
    HbAA_Yv(numscan)=wood_model(trust.T2Blood,hct,10,1);
end

% mjd fix interpolation 17 Feb 2020
% interp_eTE=0:1:160;
% 
% s_interp=trust.MeanSagSinus(1)*exp(-(interp_eTE/1000)/mean(trust_T2)); % calculate the curve
% figure; scatter(trust.eTE,trust.MeanSagSinus,'ro');hold on; plot(interp_eTE,s_interp,'k--')

% % % % old way
% % % %s=trust.MeanSagSinus(1)*exp(-(trust.eTE/1000)/mean(trust_T2));
% % % % figure; scatter(trust.eTE,trust.MeanSagSinus,'k*');hold on; plot(trust.eTE,s,'r-')

