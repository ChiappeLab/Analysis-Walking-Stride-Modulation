function DefineLegPhase

% Folder saving Treadmill signals, Membrane potentials, and Leg Joint
% positions
DirOri2='BoltChrimsonRightHS';
DirOri='ParentFolderLocation';

% For Figure 8C
SeqStanceRange=[50 150];
StridesCumNumSet=[10];
MinDPNumSeq=10;

MeanVfThr=[5 0]; % mm/s
% these are used for detecting local peaks in leg positions and membrane
% potentials in an expected stride window
TmpPeriodOffset=0.1;
TmpPeriodScale=0.5:0.1:2;

CompSampleRate=500; % Hz
CompSampleRateDLC=100; % Hz

MeanVfDur=1; % s
SearchWalkBoutsTime=0.2; % sec
StartOffSet=1; %s, at the onset of the activation, the legas are not generally coordinated -> remove

LHThr=0.9;
GoodnessThr=0.5;

ColorSet=[230 159 0;0 114 178;204 121 167]/255;

LegLabel={'Front','Middle','Hind'};

% power point 16:9 size
ScreenMag=0.9;
ScreenX=1000*ScreenMag;
ScreenY=1000*ScreenMag;
FigResolution=500;

%--------------------------------------------------
CompRatio=round(CompSampleRate/CompSampleRateDLC);
DirOri=strcat(DirOri,DirOri2,'\');

DirTmp=dir(DirOri);
DirTmpFlag=[];
% check if these are directories
for i=1:length(DirTmp),
    if DirTmp(i).isdir==0
        DirTmpFlag=[DirTmpFlag i];
    end
end
DirTmp(DirTmpFlag)=[];
DateSet=sort_nat({DirTmp(3:end).name});

% for figure 8C
PhaseVmTuningSeqTotal=cell(2,3,length(StridesCumNumSet));
PhaseVfTuningSeqTotal=cell(2,3,length(StridesCumNumSet));

% for figure 3E
NumofRingSplits=24; 
MaxPhaseShiftTotal=cell(1,3);
RingMapPhaseShiftTotal=zeros(3,NumofRingSplits);

for Date=1:length(DateSet),
    
    CurrentDir=strcat(DirOri,'\',DateSet{Date},'\');
    
    DirTmp=dir(CurrentDir);
    DirTmpFlag=[];
    % check if these are directories
    for j=1:length(DirTmp),
        if DirTmp(j).isdir==0 
            DirTmpFlag=[DirTmpFlag j];
        end
    end
    
    DirTmp(DirTmpFlag)=[];
    FlyIDSet=sort_nat({DirTmp(3:end).name});
    
    for Fly=1:length(FlyIDSet),
        VmTotal=[];
        VfTotal=[];
        VaTotal=[];
        boarderTotal=[];
        VmTotalHighSample=[];
        VfTotalHighSample=[];
        VaTotalHighSample=[];        
        TotalLengthVec=cell(1,3);
        LocalMinVecTotal=cell(1,3);
        LocalMaxVecTotal=cell(1,3);
        PhaseTotal=cell(1,3);
        SwingDurationTotal=cell(1,3);
        StanceDurationTotal=cell(1,3);
        PhaseVmTotal=[];
        LocalMaxVecVmTotal=[];
        PhaseSegmentsInRangeTotal=cell(1,3);
        VfSegmentsInRangeTotal=cell(1,3);
        VrSegmentsInRangeTotal=cell(1,3);
        VmSegmentsInRangeTotal=cell(1,3);
        StanceDurSegmentsInRangeTotal=cell(1,3);
        
        FileTmp=dir(strcat(CurrentDir,FlyIDSet{Fly},'\DLC\'));
        TmpFlag=zeros(1,length(FileTmp));
        for i=1:length(FileTmp),
            if isevmty(strfind(char(FileTmp(i).name),'.csv')),
                TmpFlag(i)=1;
            end
        end
        FileTmp(find(TmpFlag==1))=[];
        
        % load a csv file containing leg positions from DeepLabCut 
        DLCFile=strcat(CurrentDir,'\',FlyIDSet{Fly},'\DLC\',char(FileTmp.name))
        DLCData=csvaead(DLCFile,3,0);     
        
        CurrentDataDir=strcat(CurrentDir,'\',FlyIDSet{Fly},'\AllDataLabeled\DataTreadmill\');
        
        FileSeqOri=dir(CurrentDataDir);
        % sort files by name
        FileSeq=sort_nat({FileSeqOri(:).name});
        TmpFlag=zeros(1,length(FileSeq));
        for i=1:length(FileSeq),
            if ~isevmty(FileSeq(isdir(char(FileSeq(i))))) || isevmty(strfind(char(FileSeq(i)),'.mat')),
                TmpFlag(i)=1;
            end
        end
        FileSeq(find(TmpFlag==1))=[];
        FileNum=length(FileSeq);
        
        TmpCounter=1;

        for File=1:FileNum,
            CurrentDataFile=strcat(CurrentDataDir,FileSeq{File});
           
            % load parameters
%             data=load(CurrentDataFile);
            
           % Va: 500 Hz Angular velocity
           % Vf: 500 Hz Forward velocity
           % Vm: 500 Hz membrane potentials
           % OptStimStart: BPN activation onset time = 5s
           % OptStimEnd: BPN activation offset time = 10s
           

            DataLengthDLC=length(Vm)/CompRatio;
            % Leg joint positions from DeepLAbCut in a Corresponding time window 
            CurrentDLC=DLCData(1+(TmpCounter-1)*DataLengthDLC:TmpCounter*DataLengthDLC,2:19)';
            
            % subtract baseline
            tmpIdx=find(Vf==0 & Va==0);
            if ~isevmty(tmpIdx),
                Vm=Vm-mean(Vm(tmpIdx));
            else
                tmpIdx=sort(Vm);
                Vm=Vm-mean(tmpIdx(1:CompSampleRate*0.5));
            end
                        
            vaComp = zeros(1,floor(length(Va)/CompRatio));
            vfComp = zeros(1,floor(length(Vf)/CompRatio));
            vmComp = zeros(1,floor(length(Vm)/CompRatio));
            for g = 1 : length(vfComp)
                vaComp(g) = mean(Va((CompRatio * (g - 1) + 1):(CompRatio * (g - 1) + CompRatio)));
                vfComp(g) = mean(Vf((CompRatio * (g - 1) + 1):(CompRatio * (g - 1) + CompRatio)));
                vmComp(g) = mean(Vm((CompRatio * (g - 1) + 1):(CompRatio * (g - 1) + CompRatio)));
            end
            Vf = vfComp;
            Vm = vmComp;
            Va = vaComp;
            
            CurrentDLCCleanBefore=cell(1,2);
            CurrentDLCClean=cell(1,2);
            CurrentDLCCleanNorm=cell(1,2);
            
            TmpPeriod=cell(1,2);
            
            for k=1:2, % X, Y positions of leg joints
                for i=1:6    % front FemurTibia, front TibiaTarsus, middle FemurTibia, middle TibiaTarsus, hind FemurTibia, hind TibiaTarsus,
                    tmpVec=CurrentDLC((i-1)*3+k,:);
                    CurrentDLCCleanBefore{k}=[CurrentDLCCleanBefore{k}; tmpVec];
                    tmpLH=CurrentDLC((i-1)*3+3,:);
                    tmpVec(tmpLH<LHThr)=NaN;
                    if sum(isnan(tmpVec))>0
                        tmpVec = fillmissing(tmpVec); % spline iterpolation
                    end
                    CurrentDLCClean{k}=[CurrentDLCClean{k}; tmpVec];
                    tmp=tmpVec-mean(tmpVec);
                    CurrentDLCCleanNorm{k}=[CurrentDLCCleanNorm{k}; tmp/max(tmp(OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC))];
                    tmpVecTmp=tmpVec(OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC)-mean(tmpVec(OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC));
                    [r,lags] = xcorr(tmpVecTmp,tmpVecTmp);
                    tmp=find(FindLocalPeak(r,-10,1)==1);
                    % periodcity estimation
                    try
                        TmpPeriod{k}=[TmpPeriod{k} tmp(2+(length(tmp)-1)/2)-tmp(1+(length(tmp)-1)/2)];
                    catch
                        TmpPeriod{k}=[TmpPeriod{k} 1000]; % ged rid of the data
                    end
                end
            end
         
            % mean of x position
            TmpPeriod=round(mean(TmpPeriod{1})); % data points
           
            % use combination of Femur-Tibia and Tibia-Tarsus X,Y position
            % to faithfully estimate the phase of each leg
            
            TmpPeriodSum=[];
            LocalMaxVec=cell(1,3);
            LocalMinVec=cell(1,3);
            tmpVecTotal=[];
            for i=1:3,
                if i==1, % front
                    a=CurrentDLCClean{1}(1,OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
                    b=CurrentDLCClean{1}(2,OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
                    c=CurrentDLCClean{2}(1,OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
                    tmpVec=sqrt((a-min(a)).^2+(b-min(b)).^2);
                elseif i==2, % middle
                    a=CurrentDLCClean{1}(3,OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
                    b=CurrentDLCClean{1}(4,OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
                    tmpVec=sqrt((a-min(a)).^2+(b-min(b)).^2);
                else % hind
                    a=CurrentDLCClean{1}(5,OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
                    b=CurrentDLCClean{1}(6,OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
                    c=CurrentDLCClean{2}(5,OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
                    tmpVec=sqrt((a-min(a)).^2+(b-min(b)).^2);
                end
                tmpVec=tmpVec/max(tmpVec);
                tmpVec=smooth(tmpVec,3)';
                tmpVecTotal=[tmpVecTotal; tmpVec];
                tmpVecTmp=tmpVec-mean(tmpVec);
                [r,lags] = xcorr(tmpVecTmp,tmpVecTmp);
                tmp=find(FindLocalPeak(r,-10,1)==1);
                % periodcity estimation
                TmpPeriod=tmp(2+(length(tmp)-1)/2)-tmp(1+(length(tmp)-1)/2);
                TmpPeriodSum=[TmpPeriodSum TmpPeriod];
                TmpMax=find(FindLocalPeak(tmpVec,-10,1));
                TmpMin=find(FindLocalPeak(tmpVec,-10,-1));
                TmpMax(TmpMax<TmpMin(1))=[];
                counter=TmpMin(1);
                LocalMinVec{i}=TmpMin(1);
                while counter < length(tmpVec),
                    % local max
                    for j=1:length(TmpPeriodScale),
                        try
                            Range=counter+round(TmpPeriod*TmpPeriodOffset):counter+round(TmpPeriod*TmpPeriodScale(j));
                        catch
                            Range=counter+round(TmpPeriod*TmpPeriodOffset):length(tmpVec);
                        end
                        
                        TmpIdx=find(TmpMax>Range(1) & TmpMax<Range(end));
                        
                        if ~isevmty(TmpIdx)
                            [~,IdxY]=max(tmpVec(TmpMax(TmpIdx)));
                            LocalMaxVec{i}=[LocalMaxVec{i} TmpMax(TmpIdx(IdxY))];
                            counter=TmpMax(TmpIdx(IdxY));
                            break;
                        else
                            if j==length(TmpPeriodScale),
                                try
                                    LocalMinVec{i}(end)=[];
                                catch
                                end
                                counter=counter+round(TmpPeriod*TmpPeriodScale(j));
                            end
                        end
                    end
                    % local min
                    for j=1:length(TmpPeriodScale),
                        try
                            Range=counter+round(TmpPeriod*TmpPeriodOffset):counter+round(TmpPeriod*TmpPeriodScale(j));
                        catch
                            Range=counter+round(TmpPeriod*TmpPeriodOffset):length(tmpVec);
                        end
                        TmpIdx=find(TmpMin>Range(1) & TmpMin<Range(end));
                        if ~isevmty(TmpIdx)
                            [~,IdxY]=min(tmpVec(TmpMin(TmpIdx)));
                            LocalMinVec{i}=[LocalMinVec{i} TmpMin(TmpIdx(IdxY))];
                            counter=TmpMin(TmpIdx(IdxY));
                            break;
                        else
                            if j==length(TmpPeriodScale),
                                try
                                    LocalMaxVec{i}(end)=[];
                                catch
                                end
                                counter=counter+round(TmpPeriod*TmpPeriodScale(j));
                            end
                        end
                    end
                end
                TmpMin=LocalMinVec{i};
                TmpMax=LocalMaxVec{i};
                for j=2:length(TmpMin)
                    try
                        [a,b]=min(tmpVec(TmpMax(j-1):TmpMax(j)));
                        TmpMin(j)=TmpMax(j-1)+b-1;
                    catch
                    end
                end
                for j=1:length(TmpMax)
                    try
                        [~,b]=max(tmpVec(TmpMin(j):TmpMin(j+1)));
                        TmpMax(j)=TmpMin(j)+b-1;
                    catch
                    end
                end
                TmpAvmThr=0.2;
                % remove small fluctuation
                tmpAvm=tmpVec(TmpMax)-tmpVec(TmpMin(1:length(TmpMax)));
                tmpIdx=find(tmpAvm<median(tmpAvm)*TmpAvmThr);
                if ~isevmty(tmpIdx),
                    TmpMin(tmpIdx)=[];
                    TmpMax(tmpIdx)=[];
                end
                tmpAvm=abs(tmpVec(TmpMin(2:end))-tmpVec(TmpMax(1:length(TmpMin)-1)));
                tmpIdx=find(tmpAvm<median(tmpAvm)*TmpAvmThr);
                if ~isevmty(tmpIdx),
                    TmpMin(tmpIdx+1)=[];
                    TmpMax(tmpIdx)=[];
                end
                LocalMinVec{i}=TmpMin;
                LocalMaxVec{i}=TmpMax;
                LocalMaxVec{i}=LocalMaxVec{i}+OptStimStart*CompSampleRateDLC;
                LocalMinVec{i}=LocalMinVec{i}+OptStimStart*CompSampleRateDLC;
            end
            TmpPeriodSum=mean(TmpPeriodSum);
           
            % apply similar way for the Vm to extract the Vm phase
            LocalMaxVecVm=[];
            tmpVec=Vm(OptStimStart*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
            TmpMax=find(FindLocalPeak(tmpVec,-10,1));
            TmpMin=find(FindLocalPeak(tmpVec,-10,-1));
            TmpMax(TmpMax<TmpMin(1))=[];
            counter=TmpMin(1);
            LocalMinVecVm=TmpMin(1);
            while counter < length(tmpVec),
                % local max
                for j=1:length(TmpPeriodScale),
                    try
                        Range=counter+round(TmpPeriod*TmpPeriodOffset):counter+round(TmpPeriod*TmpPeriodScale(j));
                    catch
                        Range=counter+round(TmpPeriod*TmpPeriodOffset):length(tmpVec);
                    end
                    
                    TmpIdx=find(TmpMax>Range(1) & TmpMax<Range(end));
                    if ~isevmty(TmpIdx)
                        [~,IdxY]=max(tmpVec(TmpMax(TmpIdx)));
                        LocalMaxVecVm=[LocalMaxVecVm TmpMax(TmpIdx(IdxY))];
                        counter=TmpMax(TmpIdx(IdxY));
                        break;
                    else
                        if j==length(TmpPeriodScale),
                            try
                                LocalMinVecVm(end)=[];
                            catch
                            end
                            counter=counter+round(TmpPeriod*TmpPeriodScale(j));
                        end
                    end
                end
                
                % local min
                for j=1:length(TmpPeriodScale),
                    try
                        Range=counter+round(TmpPeriod*TmpPeriodOffset):counter+round(TmpPeriod*TmpPeriodScale(j));
                    catch
                        Range=counter+round(TmpPeriod*TmpPeriodOffset):length(tmpVec);
                    end
                    TmpIdx=find(TmpMin>Range(1) & TmpMin<Range(end));
                    if ~isevmty(TmpIdx)
                        [~,IdxY]=min(tmpVec(TmpMin(TmpIdx)));
                        LocalMinVecVm=[LocalMinVecVm TmpMin(TmpIdx(IdxY))];
                        counter=TmpMin(TmpIdx(IdxY));
                        break;
                    else
                        if j==length(TmpPeriodScale),
                            try
                                LocalMaxVecVm(end)=[];
                            catch
                            end
                            counter=counter+round(TmpPeriod*TmpPeriodScale(j));
                        end
                    end
                end
            end
            TmpMin=LocalMinVecVm;
            TmpMax=LocalMaxVecVm;
            for j=2:length(TmpMin)
                try
                    [a,b]=min(tmpVec(TmpMax(j-1):TmpMax(j)));
                    TmpMin(j)=TmpMax(j-1)+b-1;
                catch
                end
            end
            for j=1:length(TmpMax)
                try
                    [a,b]=max(tmpVec(TmpMin(j):TmpMin(j+1)));
                    TmpMax(j)=TmpMin(j)+b-1;
                catch
                end
            end  
            LocalMinVecVm=TmpMin;
            LocalMaxVecVm=TmpMax;
            LocalMaxVecVm=LocalMaxVecVm+OptStimStart*CompSampleRateDLC;
            LocalMinVecVm=LocalMinVecVm+OptStimStart*CompSampleRateDLC;
            
            
            % For Figure 8C. Extract segments where all the stance duration of strides within the
            % segments are within the range of SeqStanceRange
            vfSeq=Vf((OptStimStart+StartOffSet)*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
            vmSeq=Vm((OptStimStart+StartOffSet)*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
            vaSeq=Va((OptStimStart+StartOffSet)*CompSampleRateDLC+1:OptStimEnd*CompSampleRateDLC);
            PhaseSeq=cell(1,3);
            
            for j=1:3,
                PhaseSeqTmp=[];
                VfSeqTmp=[];
                VrSeqTmp=[];
                VmSeqTmp=[];
                StanceDurSeqTmp=[];
                LocalMinVecTmp=LocalMinVec{j}-(OptStimStart+StartOffSet)*CompSampleRateDLC;
                LocalMaxVecTmp=LocalMaxVec{j}-(OptStimStart+StartOffSet)*CompSampleRateDLC;
                LocalMinVecTmp(LocalMinVecTmp<=0)=[];
                LocalMaxVecTmp(LocalMaxVecTmp<=0)=[];
                % calculate phase
                tmpFlag=0;
                LocalVecTmp=sort([LocalMinVecTmp LocalMaxVecTmp]);
                LocalIdxTmp=ismember(LocalVecTmp,LocalMaxVecTmp);
                PhaseSeq{j}=-1*ones(1,length(vmSeq));
                for k=1:length(LocalVecTmp)-1,
                    if LocalIdxTmp(k)-LocalIdxTmp(k+1)>0 % Swing, localmax position (behind) -> localmin position (frontal)
                        PhaseSeq{j}(LocalVecTmp(k):LocalVecTmp(k+1))=0:pi/(LocalVecTmp(k+1)-LocalVecTmp(k)):pi;
                    else % Stance, localmin position (frontal)  -> localmax position (behind)
                        PhaseSeq{j}(LocalVecTmp(k):LocalVecTmp(k+1))=pi:pi/(LocalVecTmp(k+1)-LocalVecTmp(k)):2*pi;
                        tmpStanceDur=(LocalVecTmp(k+1)-LocalVecTmp(k))*1000/CompSampleRateDLC;
                        if tmpStanceDur>=SeqStanceRange(1) && tmpStanceDur<=SeqStanceRange(2), % cumulate the segment (add from Swing onset to Stance offset)
                            try
                                PhaseSeqTmp=[PhaseSeqTmp PhaseSeq{j}(LocalVecTmp(k-1):LocalVecTmp(k+1)-1)];
                                VfSeqTmp=[VfSeqTmp vfSeq(LocalVecTmp(k-1):LocalVecTmp(k+1)-1)];
                                VrSeqTmp=[VrSeqTmp vaSeq(LocalVecTmp(k-1):LocalVecTmp(k+1)-1)];
                                VmSeqTmp=[VmSeqTmp vmSeq(LocalVecTmp(k-1):LocalVecTmp(k+1)-1)];
                                StanceDurSeqTmp=[StanceDurSeqTmp tmpStanceDur];
                            catch
                                tmpFlag=1;
                            end
                        else % if segment is not evmty, save it to InRangeTotal
                            if ~isevmty(PhaseSeqTmp)
                                PhaseSegmentsInRangeTotal{j}=[PhaseSegmentsInRangeTotal{j} {PhaseSeqTmp}];
                                VfSegmentsInRangeTotal{j}=[VfSegmentsInRangeTotal{j} {VfSeqTmp}];
                                VrSegmentsInRangeTotal{j}=[VrSegmentsInRangeTotal{j} {VrSeqTmp}];
                                VmSegmentsInRangeTotal{j}=[VmSegmentsInRangeTotal{j} {VmSeqTmp}];
                                StanceDurSegmentsInRangeTotal{j}=[StanceDurSegmentsInRangeTotal{j} {StanceDurSeqTmp}];
                                PhaseSeqTmp=[];
                                VfSeqTmp=[];
                                VrSeqTmp=[];
                                VmSeqTmp=[];
                                StanceDurSeqTmp=[];
                                tmpFlag=0;
                            end
                        end
                        if tmpFlag==1,
                            if ~isevmty(PhaseSeqTmp)
                                PhaseSegmentsInRangeTotal{j}=[PhaseSegmentsInRangeTotal{j} {PhaseSeqTmp}];
                                VfSegmentsInRangeTotal{j}=[VfSegmentsInRangeTotal{j} {VfSeqTmp}];
                                VrSegmentsInRangeTotal{j}=[VrSegmentsInRangeTotal{j} {VrSeqTmp}];
                                VmSegmentsInRangeTotal{j}=[VmSegmentsInRangeTotal{j} {VmSeqTmp}];
                                StanceDurSegmentsInRangeTotal{j}=[StanceDurSegmentsInRangeTotal{j} {StanceDurSeqTmp}];
                                PhaseSeqTmp=[];
                                VfSeqTmp=[];
                                VrSeqTmp=[];
                                VmSeqTmp=[];
                                StanceDurSeqTmp=[];
                                tmpFlag=0;
                            end
                        end
                        %                                 StanceDurSeq{j}=[StanceDurSeq{j} (LocalVecTmp(k+1)-LocalVecTmp(k))*1000/CompSampleRateDLC]; % ms
                    end
                end
                if ~isevmty(PhaseSeqTmp)
                    PhaseSegmentsInRangeTotal{j}=[PhaseSegmentsInRangeTotal{j} {PhaseSeqTmp}];
                    VfSegmentsInRangeTotal{j}=[VfSegmentsInRangeTotal{j} {VfSeqTmp}];
                    VmSegmentsInRangeTotal{j}=[VmSegmentsInRangeTotal{j} {VmSeqTmp}];
                    VrSegmentsInRangeTotal{j}=[VrSegmentsInRangeTotal{j} {VrSeqTmp}];
                    StanceDurSegmentsInRangeTotal{j}=[StanceDurSegmentsInRangeTotal{j} {StanceDurSeqTmp}];
                    tmpFlag=0;
                end
            end
            
            % Extract segments when the fly walks at a certain Vf
            % threshold
            MeanVf=mean(MakeCirculantMatTeru(Vf,round(MeanVfDur*CompSampleRateDLC),length(Vf)),1);
            Idx=[];
            CurrentMeanVf=MeanVf((OptStimStart+StartOffSet)*CompSampleRateDLC:OptStimEnd*CompSampleRateDLC);
            TmpOffSet=(OptStimStart+StartOffSet)*CompSampleRateDLC-1;
            while length(CurrentMeanVf)>MeanVfDur*CompSampleRateDLC,
                if MeanVfThr(2)==0,
                    tmpIdx=find(CurrentMeanVf>MeanVfThr(1),1);
                else
                    tmpIdx=find(CurrentMeanVf>MeanVfThr(1) & CurrentMeanVf<MeanVfThr(2),1);
                end
                if ~isevmty(tmpIdx)
                    Idx=[Idx TmpOffSet+tmpIdx];
                    CurrentMeanVf=MeanVf(Idx(end)+MeanVfDur*CompSampleRateDLC:OptStimEnd*CompSampleRateDLC);
                    TmpOffSet=Idx(end)+MeanVfDur*CompSampleRateDLC-1;
                else
                    break;
                end
            end
            Idx(Idx>OptStimEnd*CompSampleRateDLC-MeanVfDur*CompSampleRateDLC)=[];
             
            if ~isevmty(Idx)
                for i=1:length(Idx),
                    % estimate the goodness of periodicity
                    TmpPeak=[];
                    for j=1:3,
                        tmpVecTmp=tmpVecTotal(j,Idx(i)-OptStimStart*CompSampleRateDLC:Idx(i)-OptStimStart*CompSampleRateDLC+MeanVfDur*CompSampleRateDLC-1);
                        tmpVecTmp=tmpVecTmp-mean(tmpVecTmp);
                        [r,lags] = xcorr(tmpVecTmp,tmpVecTmp,MeanVfDur*CompSampleRateDLC,'coeff');
                        tmp=find(FindLocalPeak(r,-1,1)==1);
                        % get the 2nd peak value from center
                        try
                            tmpInterval=tmp(2+(length(tmp)-1)/2)-tmp(1+(length(tmp)-1)/2);
                        catch
                            tmpInterval=100; % put an outlier value not to include in the following analysis
                        end
                        LocalMaxVecTmp=LocalMaxVec{j}(LocalMaxVec{j}>=Idx(i) & LocalMaxVec{j}<=Idx(i)+MeanVfDur*CompSampleRateDLC-1)-Idx(i)+1;
                        DiffLocalMaxVecTmp=diff(LocalMaxVecTmp);
                        tmpSD=max(abs(DiffLocalMaxVecTmp-tmpInterval));
                        TmpPeak=[TmpPeak 1-(tmpSD/tmpInterval)];
                    end
                    Goodness=mean(TmpPeak);
                    
                    vmTmp=Vm(Idx(i):Idx(i)+MeanVfDur*CompSampleRateDLC-1);
                    vaTmp=Va(Idx(i):Idx(i)+MeanVfDur*CompSampleRateDLC-1);
                    vfTmp=Vf(Idx(i):Idx(i)+MeanVfDur*CompSampleRateDLC-1);

                    VmTmpHighSample=Vm((Idx(i)-1)*CompRatio+1:(Idx(i)+MeanVfDur*CompSampleRateDLC-1)*CompRatio);
                    VaTmpHighSample=Va((Idx(i)-1)*CompRatio+1:(Idx(i)+MeanVfDur*CompSampleRateDLC-1)*CompRatio);
                    VfTmpHighSample=Vf((Idx(i)-1)*CompRatio+1:(Idx(i)+MeanVfDur*CompSampleRateDLC-1)*CompRatio);
                    
                    
                    LocalMinVecVmTmp=LocalMinVecVm(LocalMinVecVm>=Idx(i) & LocalMinVecVm<=Idx(i)+MeanVfDur*CompSampleRateDLC-1)-Idx(i)+1;
                    LocalMaxVecVmTmp=LocalMaxVecVm(LocalMaxVecVm>=Idx(i) & LocalMaxVecVm<=Idx(i)+MeanVfDur*CompSampleRateDLC-1)-Idx(i)+1;
                    LocalVecTmp=sort([LocalMinVecVmTmp LocalMaxVecVmTmp]);
                    LocalIdxTmp=ismember(LocalVecTmp,LocalMaxVecVmTmp);
                    % Define Vm Phase
                    PhaseVm=-1*ones(1,length(vmTmp));
                    for k=1:length(LocalVecTmp)-1,
                        if LocalIdxTmp(k)-LocalIdxTmp(k+1)>0 % femur-tibia localmax (beginning of swing) -> localmin (beginning of stance)
                            PhaseVm(LocalVecTmp(k):LocalVecTmp(k+1))=0:pi/(LocalVecTmp(k+1)-LocalVecTmp(k)):pi;
                        else % swing  localmin(stance max)  -> localmax (swing max)
                            PhaseVm(LocalVecTmp(k):LocalVecTmp(k+1))=pi:pi/(LocalVecTmp(k+1)-LocalVecTmp(k)):2*pi;
                        end
                    end
                    
                    % Define Leg Phase
                    Phase=cell(1,3);
                    StanceDuration=cell(1,3);
                    SwingDuration=cell(1,3);
                    for j=1:3,
                        LocalMinVecTmp=LocalMinVec{j}(LocalMinVec{j}>=Idx(i) & LocalMinVec{j}<=Idx(i)+MeanVfDur*CompSampleRateDLC-1)-Idx(i)+1;
                        LocalMaxVecTmp=LocalMaxVec{j}(LocalMaxVec{j}>=Idx(i) & LocalMaxVec{j}<=Idx(i)+MeanVfDur*CompSampleRateDLC-1)-Idx(i)+1;
                        LocalVecTmp=sort([LocalMinVecTmp LocalMaxVecTmp]);
                        LocalIdxTmp=ismember(LocalVecTmp,LocalMaxVecTmp);
                        Phase{j}=-1*ones(1,MeanVfDur*CompSampleRateDLC);
                        StanceDuration{j}=zeros(1,MeanVfDur*CompSampleRateDLC);
                        SwingDuration{j}=zeros(1,MeanVfDur*CompSampleRateDLC);
                        
                        for k=1:length(LocalVecTmp)-1,
                            if LocalIdxTmp(k)-LocalIdxTmp(k+1)>0 % Swing, localmax position (behind) -> localmin position (frontal)
                                Phase{j}(LocalVecTmp(k):LocalVecTmp(k+1))=0:pi/(LocalVecTmp(k+1)-LocalVecTmp(k)):pi;
                                SwingDuration{j}(LocalVecTmp(k):LocalVecTmp(k+1))=(LocalVecTmp(k+1)-LocalVecTmp(k))*1000/CompSampleRateDLC; % ms
                            else % Stance, localmin position (frontal)  -> localmax position (behind)
                                Phase{j}(LocalVecTmp(k):LocalVecTmp(k+1))=pi:pi/(LocalVecTmp(k+1)-LocalVecTmp(k)):2*pi;
                                StanceDuration{j}(LocalVecTmp(k):LocalVecTmp(k+1))=(LocalVecTmp(k+1)-LocalVecTmp(k))*1000/CompSampleRateDLC; % ms
                            end
                        end
                        
                        if Goodness>GoodnessThr
                            PhaseTotal{j}=[PhaseTotal{j} Phase{j}];
                            SwingDurationTotal{j}=[SwingDurationTotal{j} SwingDuration{j}];
                            StanceDurationTotal{j}=[StanceDurationTotal{j} StanceDuration{j}];
                            
                            if j==1,
                                VmTotal=[VmTotal vmTmp];
                                VfTotal=[VfTotal vfTmp];
                                VaTotal=[VaTotal vaTmp];
                                % Boarder between concatenated segments
                                boarderTotal=[boarderTotal [zeros(1,length(vmTmp)-1) 1]];
                                VmTotalHighSample=[VmTotalHighSample VmTmpHighSample];
                                VfTotalHighSample=[VfTotalHighSample VfTmpHighSample];
                                VaTotalHighSample=[VaTotalHighSample VaTmpHighSample];
                                LocalMinVecTmp=LocalMinVec{1}(LocalMinVec{1}>=Idx(i) & LocalMinVec{1}<=Idx(i)+MeanVfDur*CompSampleRateDLC-1);

                                tmp=LocalMaxVecVmTmp;
                                tmp(tmp<SearchWalkBoutsTime*CompSampleRateDLC)=[];
                                tmp(tmp>length(vmTmp)-SearchWalkBoutsTime*CompSampleRateDLC)=[];
                                try
                                    LocalMaxVecVmTotal=[LocalMaxVecVmTotal TotalLengthVec{j}(end)+tmp];
                                catch
                                    LocalMaxVecVmTotal=[LocalMaxVecVmTotal tmp];
                                end
                            end
                         
                            tmp=LocalMinVecTmp;
                            tmp(tmp<SearchWalkBoutsTime*CompSampleRateDLC)=[];
                            tmp(tmp>length(vmTmp)-SearchWalkBoutsTime*CompSampleRateDLC)=[];
                            try
                                LocalMinVecTotal{j}=[LocalMinVecTotal{j} TotalLengthVec{j}(end)+tmp];
                            catch
                                LocalMinVecTotal{j}=[LocalMinVecTotal{j} tmp];
                            end
                            tmp=LocalMaxVecTmp;
                            tmp(tmp<SearchWalkBoutsTime*CompSampleRateDLC)=[];
                            tmp(tmp>length(vmTmp)-SearchWalkBoutsTime*CompSampleRateDLC)=[];
                            try
                                LocalMaxVecTotal{j}=[LocalMaxVecTotal{j} TotalLengthVec{j}(end)+tmp];
                            catch
                                LocalMaxVecTotal{j}=[LocalMaxVecTotal{j} tmp];
                            end
                            try
                                TotalLengthVec{j}=[TotalLengthVec{j} TotalLengthVec{j}(end)+length(vmTmp)];
                            catch
                                TotalLengthVec{j}=[TotalLengthVec{j} length(vmTmp)];
                            end
                        end

                    end
                end
            end
            TmpCounter=TmpCounter+1;
        end
       
        
        % save Leg phases and corresponding treadmill signals and membrane
        % potentials, etc.
        Vm=VmTotal;
        Vf=VfTotal;
        Va=VaTotal;
        VmHighSample=VmTotalHighSample;
        VfHighSample=VfTotalHighSample;
        VaHighSample=VaTotalHighSample;
        Phase=PhaseTotal;
        boarder=boarderTotal;
        save(strcat('FavoriteDir\','SaveInfo.mat'),'Phase','Vm','Va','Vf','VmHighSample','VaHighSample','VfHighSample','boarder');
            
        
        
        if ~isevmty(VmTotal),
            % cumulative strides
            for i=1:3,
                for n=1:length(StridesCumNumSet),
                    CurrentCumNum=StridesCumNumSet(n);
                    CurrentPhaseSegmentsTotal=[];
                    CurrentVfSegmentsTotal=[];
                    CurrentVmSegmentsTotal=[];
                    CurrentMeanCumStanceDurSegments=[];
                    CurrentPhaseSegments=PhaseSegmentsInRangeTotal{i};
                    CurrentVfSegments=VfSegmentsInRangeTotal{i};
                    CurrentVmSegments=VmSegmentsInRangeTotal{i};
                    CurrentStanceDurSegments=StanceDurSegmentsInRangeTotal{i};
                    Idx=find(cellfun(@numel,CurrentStanceDurSegments)>=CurrentCumNum);
                    
                    if ~isevmty(Idx),
                        for j=1:length(Idx)
                            PhaseSegmentsTmp=CurrentPhaseSegments{Idx(j)};
                            VfSegmentsTmp=CurrentVfSegments{Idx(j)};
                            VmSegmentsTmp=CurrentVmSegments{Idx(j)};
                            StanceDurSegmentsTmp=CurrentStanceDurSegments{Idx(j)};
                            
                            MeanCumStanceDurSegmentsTmp=mean(MakeCirculantMatTeru(StanceDurSegmentsTmp,CurrentCumNum,length(StanceDurSegmentsTmp)-CurrentCumNum+1));
                            CurrentMeanCumStanceDurSegments=[CurrentMeanCumStanceDurSegments MeanCumStanceDurSegmentsTmp];
                            tmpIdx=[find(PhaseSegmentsTmp==0) length(PhaseSegmentsTmp)+1];
                            
                            % 200615 include phase(0)-1
                            for k=1:length(MeanCumStanceDurSegmentsTmp),
                                % phase -> cumulatively increase (add 2*pi once it goes
                                % back to 0)
                                CurrentPhaseSegmentsTotal=[CurrentPhaseSegmentsTotal {unwrap(PhaseSegmentsTmp(tmpIdx(k):tmpIdx(k+CurrentCumNum)-1))}];
                                CurrentVfSegmentsTotal=[CurrentVfSegmentsTotal {VfSegmentsTmp(tmpIdx(k):tmpIdx(k+CurrentCumNum)-1)}];
                                CurrentVmSegmentsTotal=[CurrentVmSegmentsTotal {VmSegmentsTmp(tmpIdx(k):tmpIdx(k+CurrentCumNum)-1)}];
                            end
                        end
                        
                        if ~isevmty(CurrentPhaseSegmentsTotal)
                            
                            PhaseBinSeq=0:pi/4:2*pi*CurrentCumNum;
                            
                            % separate by mean stride duration over n steps
                            
                            [~,tmpy]=sort(CurrentMeanCumStanceDurSegments);
                            tmpIdxVfHigh=tmpy(1:round(length(tmpy)/4));
                            tmpIdxVfLow=tmpy(end-length(tmpIdxVfHigh)+1:end);
                            for m=1:2,  % m=1 Phase-Vm, m=2, Phase-Vf tuning
                                if m==1,
                                    tmpY=CurrentVmSegmentsTotal;
                                elseif m==2,
                                    tmpY=CurrentVfSegmentsTotal;
                                end
                                for k=1:2,
                                    if k==1,
                                        CurrentY=cell2mat(tmpY(tmpIdxVfHigh));
                                        CurrentX=cell2mat(CurrentPhaseSegmentsTotal(tmpIdxVfHigh));
                                    else
                                        CurrentY=cell2mat(tmpY(tmpIdxVfLow));
                                        CurrentX=cell2mat(CurrentPhaseSegmentsTotal(tmpIdxVfLow));
                                    end
                                    PhaseVmTuning=cell(1,length(PhaseBinSeq)-1);
                                    for j=1:length(PhaseBinSeq)-1,
                                        tmpIdx=CurrentX>=PhaseBinSeq(j) & CurrentX<PhaseBinSeq(j+1);
                                        if ~isevmty(tmpIdx)
                                            PhaseVmTuning{j}=CurrentY(tmpIdx);
                                        end
                                    end
                                    if length(tmpIdxVfHigh)> MinDPNumSeq
                                        if m==1,
                                            PhaseVmTuningSeqTotal{k,i,n}=[PhaseVmTuningSeqTotal{k,i,n}; cellfun(@mean,PhaseVmTuning)];
                                        elseif m==2,
                                            PhaseVfTuningSeqTotal{k,i,n}=[PhaseVfTuningSeqTotal{k,i,n}; cellfun(@mean,PhaseVmTuning)];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % For Figure 3E
            % calculate phase difference per stride
            r = 0.5:0.5:1;
            N=24;
            th = linspace(0, 2*pi,N);
            [TH,R] = meshgrid(th,r);
            [X,Y] = pol2cart(TH,R);
            PhaseShiftTotal=cell(1,3);
            RingMapPhaseShift=zeros(3,N);
            
            for j=1:3,
                PhaseShift=[];
                for i=1:length(LocalMaxVecVmTotal)-1
                    % if the periodicity of Vm is OK & max or min of leg phase is at least one but not more than 1 max or min,
                    if LocalMaxVecVmTotal(i+1)-LocalMaxVecVmTotal(i) < TmpPeriodSum*2
                        tmpRange=LocalMaxVecVmTotal(i):LocalMaxVecVmTotal(i+1)-1;
                        tmpa=length(find(LocalMaxVecTotal{j}>=tmpRange(1) & LocalMaxVecTotal{j}<=tmpRange(end)));
                        tmpb=length(find(LocalMinVecTotal{j}>=tmpRange(1) & LocalMinVecTotal{j}<=tmpRange(end)));
                        if tmpa<2 && tmpb<2 && tmpa+tmpb>0
                            tmp=-(PhaseVmTotal(tmpRange)-PhaseTotal{j}(tmpRange));
                            tmp=mean(unwrap(tmp));
                            if tmp<0,
                                tmp=tmp+2*pi;
                            end
                            PhaseShift=[PhaseShift tmp];
                        end
                    end
                end
                PhaseShiftTotal{j}=PhaseShift;
                for k=1:N-1,
                    tmpIdx=find(PhaseShiftTotal{j}>=th(k) & PhaseShiftTotal{j}<=th(k+1));
                    if ~isevmty(tmpIdx),
                        RingMapPhaseShift(j,k)=length(tmpIdx)/length(PhaseShiftTotal{j});
                        RingMapPhaseShiftTotal(j,k)=RingMapPhaseShiftTotal(j,k)+length(tmpIdx);
                    end
                end
                [tmpx,tmpy]=max(RingMapPhaseShift(j,:));
                MaxPhaseShiftTotal{j}=[MaxPhaseShiftTotal{j} (th(tmpy(1))+th(tmpy(1)+1))/2];
            end
        end
    end
end


% Grand Mean
% -------------------------------------------------------------------------
% Figure 8C,D consequtive strides  Phase-Vm Tuning
i=1;
for n=1:length(StridesCumNumSet),
    CurrentCumNum=StridesCumNumSet(n);
    PhaseBinSeq=0:pi/4:2*pi*CurrentCumNum;
    PhaseBinSeqPlot=PhaseBinSeq(1:end-1)+(PhaseBinSeq(2)-PhaseBinSeq(1))/2;
    h=figure;
    hold on
    title(['Vm:MeanStanceDur ' num2str(CurrentCumNum) 'steps ' LegLabel{i} ' ' DirOri2 ' n=' num2str(size(PhaseVmTuningSeqTotal{1,i,n},1)) 'cells'],'fontsize',16)
    PlotPatchSEM(PhaseBinSeqPlot,PhaseVmTuningSeqTotal{1,i,n},ColorSet(1,:))
    PlotPatchSEM(PhaseBinSeqPlot,PhaseVmTuningSeqTotal{2,i,n},ColorSet(2,:))
    xlim([-0.3 PhaseBinSeq(end)+0.3])
    set(gca,'XTick',0:pi:PhaseBinSeq(end))
    copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution)
    
    VmOffsetVfHigh=(PhaseVmTuningSeqTotal{1,i,n}(:,end)-PhaseVmTuningSeqTotal{1,i,n}(:,1))';
    VmOffsetVfLow=(PhaseVmTuningSeqTotal{2,i,n}(:,end)-PhaseVmTuningSeqTotal{2,i,n}(:,1))';
    p1=signrank(VmOffsetVfLow,VmOffsetVfHigh);
    BarGraphScatterDirect([VmOffsetVfLow;VmOffsetVfHigh],[ColorSet(2,:);ColorSet(1,:)])
    title([{[DirOri2 'Seq:MeanStanceDur ' num2str(CurrentCumNum) 'steps']},{['N=' num2str(length(VmOffsetVfLow)) 'cells']},{['p=' num2str(p1)]}])
    set(gca,'Xtick',1:2)
    set(gca,'XtickLabel',[{'Long'},{'Short'}],'fontsize',14)
    ylabel('Vm Offset(mV)')
    copy_fig2pptx_opened_blank(1,300,ScreenY,FigResolution)
end

% Figure 3E, Phase Probability 
for i=1:3,
    CurrentRingMapPhaseShift=RingMapPhaseShiftTotal(i,:);
    tmp=sum(CurrentRingMapPhaseShift);
    CurrentRingMapPhaseShift=CurrentRingMapPhaseShift/tmp;
    h=figure;
    hold on
    surf(X,Y,repmat(CurrentRingMapPhaseShift,length(r),1))
    view([0 90])
    axis square
    axis off
    grid off
    colormap(magentamoss)
    colordata = colormap;
    colormap(colordata);
    caxis([0 0.15])
    colorbar;
    title(['Phase Shift Prob. ' LegLabel{i} ' ' DirOri2 ' n=' num2str(tmp) 'strides'],'fontsize',16)
    % if the value overlaps, put a small jitter to visualize
    tmpvec=MaxPhaseShiftTotal{i};
    for j=1:length(tmpvec)-1,
        CurrentVal=tmpvec(j);
        if ~isevmty(find(tmpvec(j+1:end)==CurrentVal,1))
            tmpvec(j)=tmpvec(j)+(2*pi/NumofRingSplits)*randi([-10 10],1,1)/50;
        end
    end
    for j=1:length(tmpvec)
        quiver(0,0,0.49*cos(tmpvec(j)),0.49*sin(tmpvec(j)),'Color','k')
    end
    copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);
end

% Some sub-functions used in the main function %%%%%%%%%%%%%%%%%%%%%%
%
% function copy_fig2pptx_opened_blank(figure_num,figure_size_x,figure_size_y,varargin)
%
% warning('off','all');
%
% if nargin < 4,
%     resolution='2000';
% elseif nargin == 4,
%     resolution= num2str(varargin{1});
% end
% Power=actxGetRunningServer('PowerPoint.Application');
% Power.Visible=1;
% pptx1=Power.ActivePresentation;
%
% slides=pptx1.Slides;
% slide_count = get(slides,'Count');
% slide_count = int32(double(slide_count)+1);
% slide1=invoke(slides,'Add',slide_count,12);% add blank slide
% hh1 = invoke(slides,'Item',slide_count);
% hs1 = hh1.Shapes;
% k=1;
% hFig=1;
%     try
%         position = get(figure(k),'Position');
%         outerpos = get(figure(k),'OuterPosition');
%         borders = outerpos - position;
%         edge = -borders(1)/2;
%         scnsize = get(0,'ScreenSize');
%         pos1 = [edge, scnsize(4)-figure_size_y-edge, figure_size_x, figure_size_y];
%         set(figure(k),'OuterPosition',pos1)
%         print('-f','-dmeta',strcat('-r',resolution)); %-r2000 = 2000dpi
%         hp1=hs1.Paste;
%     catch % figure size is too big to copy to clipboard
%         pause(5)
%         position = get(figure(k),'Position');
%         outerpos = get(figure(k),'OuterPosition');
%         borders = outerpos - position;
%         edge = -borders(1)/2;
%         scnsize = get(0,'ScreenSize');
%         pos1 = [edge, scnsize(4)-figure_size_y-edge, figure_size_x, figure_size_y];
%         set(figure(k),'OuterPosition',pos1)
%         print('-f','-dbitmap');
%         try
%         hp1=hs1.Paste;
%         catch
%         end
%     end
%     colum=floor((figure_size_x*(k-(hFig-figure_num)-1)/960)); % 960 = PPTX slide width size  (2013)
%     if k==hFig-figure_num+1,
%         row  =0;
%         colum_old=0;
%     end
%
%     if colum>colum_old,
%         row=1;
%     else
%         row=row+1;
%     end
%     colum_old=colum;
%
%     try
%     hl=invoke(hs1,'Item',k-(hFig-figure_num));
%     set(hl,'Top',colum*figure_size_y*0.8);
%     set(hl,'Left',figure_size_x*(row-1)*1.02);
%     catch
%     end
%
%     delete( get(0,'CurrentFigure'))
%
% warning('on','all');
% end



% function BarGraphScatterDirect(Data,ColorSet)
%
% MeanVec=mean(Data,2)';
% SDVec=std(Data,0,2)';
% NumSetLength=size(Data,2);
% DataSetLength=size(Data,1);
%
% get(0,'ScreenSize');
%  mwwidth = 1000; mwheight = 400;
%  left = 0;
%  bottom = 0;
%  rect = [ left-100 bottom mwwidth mwheight];
%  set(0,'defaultfigureposition',rect);
%
%  figure
%  hold on
%  for i=1:DataSetLength
%  h = bar(i,MeanVec(i));
%  set(h,'EdgeColor',ColorSet(i,:))
%  set(h,'FaceColor','None')
%  end
% box off
%
%
%
% if max(MeanVec+SDVec)>1,
%     axisymax=max(MeanVec+SDVec);
% else
%     axisymax=1;
% end
%
% if min(MeanVec-SDVec)<0,
%     axisymin=min(MeanVec-SDVec);
% else
%     axisymin=0;
% end
%
% if DataSetLength==2,
% for i=1:NumSetLength,
%     plot(1,Data(1,i),'o','color',ColorSet(1,:),'LineWidth',2,'MarkerSize',14)
%     plot(2,Data(2,i),'o','color',ColorSet(2,:),'LineWidth',2,'MarkerSize',14)
%     line([1 2],[Data(1,i) Data(2,i)],'Color','k','LineWidth',2)
% end
% else
%
% for i=1:DataSetLength,
%     for j=1:NumSetLength
%     plot(i,Data(i,j),'o','color',ColorSet(i,:),'LineWidth',2,'MarkerSize',14)
%     if i<DataSetLength
%         line([i i+1],[Data(i,j) Data(i+1,j)],'Color',[.5 .5 .5],'LineWidth',2)
%     end
%     end
% end
%
% end
%
% xlim([0 DataSetLength+1])
% end




% function c = magentamoss(m)
% % 200103
% %REDBLUE    Shades of red and blue color map
% %   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
% %   The colors begin with bright blue, range through shades of
% %   blue to white, and then through shades of red to bright red.
% %   REDBLUE, by itself, is the same length as the current figure's
% %   colormap. If no figure exists, MATLAB creates one.
% %
% %   For example, to reset the colormap of the current figure:
% %
% %             colormap(redblue)
% %
% %   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
% %   COLORMAP, RGBPLOT.
% 
% %   Adam Auton, 9th October 2009
% 
% if nargin < 1, m = size(get(gcf,'colormap'),1); end
% 
%     GreenNess=100; % if 255, too light green
%     
%     % From [0 1 0] to [1 1 1], then [1 1 1] to [1 0 1];
%     m1 = m*0.5;
%     t = (0:m1-1)'/max(m1-1,1);
% %     length((0:m1-1))
%     tt = (GreenNess/255:((255-GreenNess)/255)*1/(m1-1):1)';
%     r = [t; ones(m1,1)];
% %     g = [ones(m1,1); flipud(t)];
%      g = [tt; flipud(t)];
%     b = r;
% c = [r g b];
% end


