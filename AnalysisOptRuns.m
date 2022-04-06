function AnalysisOptRuns

% Folder saving Treadmill signals, Membrane potentials, and Leg Phases
DirOri2='BoltChrimsonRightHS';
DirOri='ParentFolderLocation';

VaThr=50; % deg/s for local Vf max/min Trig plot (Figure 5E)
VfFiltThr=1; % mm/s

PhaseBin=0:pi/4:2*pi;
CurrentButFiltRange=[5 0];

MinDataNumTrigAve=5;

CompSampleRate=500; % Hz
CompSampleRateDLC=100; % Hz

SearchWalkBoutsTime=0.2; % sec

ColorSet=[230 159 0;0 114 178;204 121 167]/255;
LegLabel={'Front','Middle','Hind'};
%--------------------------------------------------------
DirOri=strcat(DirOri,DirOri2,'\');

PhaseBinPlot=PhaseBin(1:end-1)+(PhaseBin(2)-PhaseBin(1))/2;

if CurrentButFiltRange(1)==0,
    [bfilt,afilt]=butter(1,2*CurrentButFiltRange(2)/CompSampleRateDLC,'low');
elseif CurrentButFiltRange(2)==0,
    [bfilt,afilt]=butter(1,2*CurrentButFiltRange(1)/CompSampleRateDLC,'high');
else
    [bfilt,afilt]=butter(1,[2*CurrentButFiltRange(1)/CompSampleRateDLC 2*CurrentButFiltRange(2)/CompSampleRateDLC]);
end

if CurrentButFiltRange(1)==0,
    [bhs,ahs]=butter(1,2*CurrentButFiltRange(2)/CompSampleRate,'low');
elseif CurrentButFiltRange(2)==0,
    [bhs,ahs]=butter(1,2*CurrentButFiltRange(1)/CompSampleRate,'high');
else
    [bhs,ahs]=butter(1,[2*CurrentButFiltRange(1)/CompSampleRate 2*CurrentButFiltRange(2)/CompSampleRate]);
end



% power point 16:9 size
ScreenMag=0.9;
ScreenX=1000*ScreenMag;
ScreenY=1000*ScreenMag;
FigResolution=500;


DirTmp=dir(DirOri);
DirTmpFlag=[];
% check if these are directories
for i=1:length(DirTmp),
    if DirTmp(i).isdir==0
        DirTmpFlag=[DirTmpFlag i];
    end
end
DirTmp(DirTmpFlag)=[];
DateSet=sort_nat({DirTmp(3:end).name})

VmTotal=[];
VfTotal=[];
VaTotal=[];
DVmTotal=[];
BoarderTotal=[];

% For Figure 3F
PhaseVmTuningTotal=cell(3,3);

% For Figure 5B
FiltVfPosVaPosTrigMPVaSortedTotal=cell(3,3,3);
FiltVfPosVaPosTrigVaVaSortedTotal=cell(3,3,3);
FiltVfPosVaPosTrigVfVaSortedTotal=cell(3,3,3);
FiltVfPosVaPosTrigDMPVaSortedTotal=cell(3,3,3);

% For Figure 5C
VaRange=[0 50;150 200];
PhaseDVmTuningTotalVaRange=cell(size(VaRange,1),3);

% For Figure 8A,B
PhaseVmTuningStanceShortLongTotal=cell(2,3,2);
MeanVfStanceShortLongTotal=cell(2,3);
MeanStancePeriodStanceShortLongTotal=cell(2,3);


CurrentCellNum=1;
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
        
        % load parameters
%         load(strcat(CurrentDir,FlyIDSet{Fly},'\SavedInfo\SaveInfo.mat'));     
      %  Vm: 100Hz membrane potential
      %  Va: 100Hz angular velocity
      %  Vf: 100Hz forward velocity
      %  Phase{3}: 100Hz Front, Middle, Hind leg phases 
      %  boarder: boarder between concatenated segments
      
      %  VmHighSample: 500Hz membrane potential
      %  VaHighSample: 500Hz angular velocity
      %  VfHighSample: 500Hz fprward velocity       
      %  boarderHighSample: 500Hz boarder between concatenated segments
      
        %-- pool all the fly data --%
        VmTotal=[VmTotal Vm];
        VfTotal=[VfTotal Vf];
        VaTotal=[VaTotal Va];
        DVm=[0 diff(Vm)]*CompSampleRateDLC;
        DVmTotal=[DVmTotal DVm];
        BoarderTotal=[BoarderTotal boarder];
        DVmHighSample=[0 diff(VmHighSample)]*CompSampleRate;

        
        if ~isempty(Vm),
            vfFiltHighSample=filtfilt(bhs,ahs,VfHighSample);
            % get rid of where there is a boarder among patched segments
            TmpFlagBoarder=circshift(sum(MakeCirculantMatTeru(boarderHighSample,2*SearchWalkBoutsTime*CompSampleRate+1,length(boarderHighSample)),1),[0 SearchWalkBoutsTime*CompSampleRate]);
            DevPortion=0.5;
            
            % Figure 5B
            for k=1:3, % Front, Middle, Hind legs
                
                IdxLeftDrift=FindLocalPeak(vfFiltHighSample,VfFiltThr,1).*VaHighSample>VaThr;
                
                % get rid of points close to boarders among concatenated segments
                IdxLeftDrift(TmpFlagBoarder>0)=0;
                tmpIdx=find(IdxLeftDrift==1);
                tmpIdx(tmpIdx<SearchWalkBoutsTime*CompSampleRate)=[];
                tmpIdx(tmpIdx>length(VaHighSample)-SearchWalkBoutsTime*CompSampleRate)=[];
                
                % how much Va is compensated in the following 200 ms
                tmpDiffVa=zeros(1,length(tmpIdx));
                for j=1:length(tmpIdx),
                    tmpDiffVa(j)=VaHighSample(tmpIdx(j))-mean(VaHighSample(tmpIdx(j):tmpIdx(j)+SearchWalkBoutsTime*CompSampleRate));
                end
                [~,IdxY]=sort(tmpDiffVa);
                
                tmpNum=[];
                tmpMeanMP=cell(1,3);
                tmpMeanVa=cell(1,3);
                tmpMeanVf=cell(1,3);
                tmpMeanDMP=cell(1,3);

                for n=1:2, % n=1 top Idx, n=2 bottom idx
                    Idx=zeros(1,length(VmHighSample));
                    if n==1,
                        Idx(tmpIdx(IdxY(1:round(length(IdxY)*DevPortion))))=1;
                        tmpColor=ColorSet(k,:);
                    elseif n==2,
                        Idx(tmpIdx(IdxY(end-round(length(IdxY)*DevPortion)+1:end)))=1;
                        tmpColor=ColorSet(k,:)/2;
                    end
                    
                    TrigMeanMP=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
                    TrigSDMP=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
                    
                    TrigMeanVa=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
                    TrigSDVa=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
                    
                    TrigMeanVf=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
                    TrigSDVf=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
                    
                    TrigMeanDMP=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
                    TrigSDDMP=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
                    
                    KerMat=zeros(2*SearchWalkBoutsTime*CompSampleRate+1,length(VmHighSample));
                    for j=1:2*SearchWalkBoutsTime*CompSampleRate+1,
                        TmpShift=j-(1+SearchWalkBoutsTime*CompSampleRate);
                        KerMat(j,:)=circshift(Idx,[0 TmpShift]);
                    end
                    for j=1:2*SearchWalkBoutsTime*CompSampleRate+1,
                        TmpVec=KerMat(j,:).*VmHighSample;
                        TmpVec=TmpVec(KerMat(j,:)==1);
                        TrigMeanMP(j)=mean(TmpVec);
                        TrigSDMP(j)=std(TmpVec);
                        
                        TmpVec=KerMat(j,:).*VfHighSample;
                        TmpVec=TmpVec(KerMat(j,:)==1);
                        TrigMeanVf(j)=mean(TmpVec);
                        TrigSDVf(j)=std(TmpVec);
                        
                        TmpVec=KerMat(j,:).*VaHighSample;
                        TmpVec=TmpVec(KerMat(j,:)==1);
                        TrigMeanVa(j)=mean(TmpVec);
                        TrigSDVa(j)=std(TmpVec);
                        
                        TmpVec=KerMat(j,:).*DVmHighSample;
                        TmpVec=TmpVec(KerMat(j,:)==1);
                        TrigMeanDMP(j)=mean(TmpVec);
                        TrigSDDMP(j)=std(TmpVec);
                    end
                    
                    tmpNum=[tmpNum sum(Idx)];
                    tmpMeanMP{n}=TrigMeanMP;
                    tmpMeanVa{n}=TrigMeanVa;
                    tmpMeanVf{n}=TrigMeanVf;
                    tmpMeanDMP{n}=TrigMeanDMP;
                end

                if tmpNum(1) >= MinDataNumTrigAve && tmpNum(2) >= MinDataNumTrigAve,
                    for n=1:3,
                        FiltVfPosVaPosTrigMPVaSortedTotal{k,n}=[FiltVfPosVaPosTrigMPVaSortedTotal{k,n}; mean(tmpMeanMP{n},1)];
                        FiltVfPosVaPosTrigVaVaSortedTotal{k,n}=[FiltVfPosVaPosTrigVaVaSortedTotal{k,n}; mean(tmpMeanVa{n},1)];
                        FiltVfPosVaPosTrigVfVaSortedTotal{k,n}=[FiltVfPosVaPosTrigVfVaSortedTotal{k,n}; mean(tmpMeanVf{n},1)];
                        FiltVfPosVaPosTrigDMPVaSortedTotal{k,n}=[FiltVfPosVaPosTrigDMPVaSortedTotal{k,n}; mean(tmpMeanDMP{n},1)];
                    end
                end
                
                
                % Figure 5C: Phase-DVm Tuning at differnt Va range
                for n=1:length(VaRange),
                    CurrentVaRange=Va>VaRange(n,1) & Va<=VaRange(n,2);
                    tmpDMPVec=DVm(CurrentVaRange);
                    tmpPhaseVec=Phase{k}(CurrentVaRange);
                    PhaseDVmTuning=cell(1,length(PhaseBin)-1);
                    
                    for j=1:length(PhaseBin)-1,
                        tmptmpIdx=find(tmpPhaseVec>=PhaseBin(j) & tmpPhaseVec<PhaseBin(j+1));
                        if ~isempty(tmptmpIdx)
                            PhaseDVmTuning{j}=tmpDMPVec(tmptmpIdx);
                        else
                            PhaseDVmTuning{j}=0;
                        end
                    end
                    PhaseDVmTuningTotalVaRange{n,k}=[PhaseDVmTuningTotalVaRange{n,k}; cellfun(@mean,PhaseDVmTuning)];
                end
            end
            
            vfFiltTmp=filtfilt(bfilt,afilt,Vf);
            mpFiltTmp=filtfilt(bfilt,afilt,Vm);
            
            % Figure 3F, 4F
            for i=1:3
                for m=1:3, % 1: Vm, 2: Vf, 3:Vm (normalized per cell)
                    PhaseVmTuning=cell(1,length(PhaseBin)-1);
                    for j=1:length(PhaseBin)-1,
                        tmpIdx=find(Phase{i}>=PhaseBin(j) & Phase{i}<PhaseBin(j+1));
                        if ~isempty(tmpIdx)
                            if m==1 || m==3,
                                PhaseVmTuning{j}=mpFiltTmp(tmpIdx);
                            elseif m==2
                                PhaseVmTuning{j}=vfFiltTmp(tmpIdx);
                            end
                        end
                    end    
                        tmp=cellfun(@mean,PhaseVmTuning);
                        if m==3
                            tmp=tmp/max(abs(tmp));
                        end
                    PhaseVmTuningTotal{i,m}=[PhaseVmTuningTotal{i,m}; tmp];
                end
            end
            
            for i=1:3
                % devide into stance duration short and long
                [~,tmpy]=sort(StancePeriod{i});
                tmpIdxVfHigh=tmpy(1:round(length(tmpy)/4));
                tmpIdxVfLow=tmpy(end-length(tmpIdxVfHigh)+1:end);
                
                for m=1:2,  % m=1 Phase-Vm, m=2 Phase-Vf,
                    for k=1:2,
                        if k==1,
                            CurrentVf=cell2mat(VfStancePeriod{i}(tmpIdxVfHigh));
                            CurrentVm=cell2mat(VmStancePeriod{i}(tmpIdxVfHigh));
                            tmpColor=ColorSet(1,:);
                            CurrentStancePeriod=StancePeriod{i}(tmpIdxVfHigh);
                            CurrentPhase=cell2mat(PhaseStancePeriod{i}(tmpIdxVfHigh));
                        else
                            CurrentVf=cell2mat(VfStancePeriod{i}(tmpIdxVfLow));
                            CurrentVm=cell2mat(VmStancePeriod{i}(tmpIdxVfLow));
                            tmpColor=ColorSet(2,:);
                            CurrentStancePeriod=StancePeriod{i}(tmpIdxVfLow);
                            CurrentPhase=cell2mat(PhaseStancePeriod{i}(tmpIdxVfLow));
                        end
                        if m==1
                            MeanVfStanceShortLongTotal{k,i}=[MeanVfStanceShortLongTotal{k,i} mean(CurrentVf)];
                            MeanStancePeriodStanceShortLongTotal{k,i}=[MeanStancePeriodStanceShortLongTotal{k,i} mean(CurrentStancePeriod)];
                        end
                        PhaseVmTuning=cell(1,length(PhaseBin)-1);
                        for j=1:length(PhaseBin)-1,
                            tmpIdx=CurrentPhase>=PhaseBin(j) & CurrentPhase<PhaseBin(j+1);
                            if ~isempty(tmpIdx)
                                if m==1,
                                    PhaseVmTuning{j}=CurrentVm(tmpIdx);
                                elseif m==2,
                                    PhaseVmTuning{j}=CurrentVf(tmpIdx);
                                end
                            end
                        end
                        PhaseVmTuningStanceShortLongTotal{k,i,m}=[PhaseVmTuningStanceShortLongTotal{k,i,m}; cellfun(@mean,PhaseVmTuning)];
                    end
                end
            end
        end
        CurrentCellNum=CurrentCellNum+1;
    end
end

% Plotting Grand Mean-------------------------------------------------

TmpTimeHighSample=(-SearchWalkBoutsTime*CompSampleRate:SearchWalkBoutsTime*CompSampleRate)/CompSampleRate;

% Figure 5B
k=1; % Front leg

h=figure;
subplot1(4,1,'FontS',16,'Gap',[0 0.05]);
subplot1(1)
hold on
title([{['VfMaxVaPos VaSorted VaThr=' num2str(VaThr)  'deg/s:VfThr=' num2str(VfFiltThr) 'mm/s N=' num2str(size(FiltVfPosVaPosTrigMPVaSortedTotal{k,1},1)) 'cells']},{[DirOri2 'Population']}])
for n=1:2,
    if n==1,
        tmpColor=ColorSet(k,:);
    elseif n==2,
        tmpColor=ColorSet(k,:)/2;
    end
    subplot1(1)
    hold on
    PlotPatchSEM(TmpTimeHighSample,FiltVfPosVaPosTrigMPVaSortedTotal{k,n},tmpColor);
    subplot1(2)
    hold on
    PlotPatchSEM(TmpTimeHighSample,FiltVfPosVaPosTrigVaVaSortedTotal{k,n},tmpColor);
    subplot1(3)
    hold on
    PlotPatchSEM(TmpTimeHighSample,FiltVfPosVaPosTrigDMPVaSortedTotal{k,n},tmpColor);
    subplot1(4)
    hold on
    PlotPatchSEM(TmpTimeHighSample,FiltVfPosVaPosTrigVfVaSortedTotal{k,n},tmpColor);
end
subplot1(1)
hold on
ylabel('MP[mV]')
set(gca,'box','off')
xlim([-SearchWalkBoutsTime SearchWalkBoutsTime])
line([0 0],get(gca,'ylim'),'LineStyle','--','LineWidth',2,'color','k')
subplot1(2)
hold on
ylabel('Va[deg/s]')
xlim([-SearchWalkBoutsTime SearchWalkBoutsTime])
line([0 0],get(gca,'ylim'),'LineStyle','--','LineWidth',2,'color','k')
set(gca,'box','off')
subplot1(3)
hold on
ylabel('dVm/dt[mV]')
line([0 0],get(gca,'ylim'),'LineStyle','--','LineWidth',2,'color','k')
set(gca,'box','off')
subplot1(4)
hold on
ylabel('Vf[mm/s]')
xlim([-SearchWalkBoutsTime SearchWalkBoutsTime])
line([0 0],get(gca,'ylim'),'LineStyle','--','LineWidth',2,'color','k')
set(gca,'box','off')
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

% Figure 3F, 4F
for i=1:3, % Front, Middle, Hind Legs
    h=figure;
    hold on
    title(['Phase-VmFilt Tuning ' LegLabel{i} ' ' DirOri2 ' n=' num2str(size(PhaseVmTuningTotal{i,1},1)) 'cells'],'fontsize',16)
    PlotPatchSEM(PhaseBinPlot,PhaseVmTuningTotal{i,1},ColorSet(i,:))
    xlim([0 2*pi])
    set(gca,'XTick',[0 pi 2*pi])
    set(gca,'XTickLabel',{'0','pi','2pi'})
    copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);
    
    h=figure;
    hold on
    title(['Phase-VfFilt Tuning ' LegLabel{i} ' ' DirOri2 ' n=' num2str(size(PhaseVmTuningTotal{i,2},1)) 'cells'],'fontsize',16)
    PlotPatchSEM(PhaseBinPlot,PhaseVmTuningTotal{i,2},ColorSet(i,:))
    xlim([0 2*pi])
    set(gca,'XTick',[0 pi 2*pi])
    set(gca,'XTickLabel',{'0','pi','2pi'})
    copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);
    
    h=figure;
    hold on
    title(['Phase-VmFilt Tuning (Norm) ' LegLabel{i} ' ' DirOri2 ' n=' num2str(size(PhaseVmTuningTotal{i,3},1)) 'cells'],'fontsize',16)
    PlotPatchSEM(PhaseBinPlot,PhaseVmTuningTotal{i,3},ColorSet(i,:))
    xlim([0 2*pi])
    set(gca,'XTick',[0 pi 2*pi])
    set(gca,'XTickLabel',{'0','pi','2pi'})
    copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);
end

i=1; % Front Leg

% Figure 5C, D
h=figure;
hold on
title([{['DVm VaRange:' LegLabel{i} ' ' DirOri2 ' N=' num2str(size(PhaseDVmTuningTotalVaRange{1,i},1))]},{['VaRange=' num2str(VaRange(1,1)) '~'  num2str(VaRange(end,2)) 'deg/s VaThr=' num2str(VaThr) ' deg/s VfThr=' num2str(VfFiltThr) 'mm/s']}])
for n=1:size(VaRange,1),
    PlotPatchSD_MeanSD(PhaseBinPlot,MeanSkipZero(PhaseDVmTuningTotalVaRange{n,i}),SDSkipZero(PhaseDVmTuningTotalVaRange{n,i})/sqrt(size(PhaseDVmTuningTotalVaRange{n,i},1)),ColorSet(i,:)/n);
end
xlim([0 2*pi])
set(gca,'XTick',[0 pi 2*pi])
set(gca,'XTickLabel',{'0','pi','2pi'})
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

AUCtmp=zeros(size(VaRange,1),size(PhaseDVmTuningTotalVaRange{1,i},1));
for n=1:size(VaRange,1),
    for m=1:size(PhaseDVmTuningTotalVaRange{1,i},1),
        CurrentCurve=PhaseDVmTuningTotalVaRange{n,i}(m,:);
        AUCtmp(n,m)=mean(CurrentCurve);
    end
end
% compare only smallest and largest
ColorSetTmp=[ColorSet(i,:);ColorSet(i,:)/4;];
BarGraphScatterDirectWithConnection(AUCtmp([1 end],:),ColorSetTmp)
[p1,ff,stats]=signrank(AUCtmp(1,:),AUCtmp(end,:));
try
    zval=stats.zval;
catch
    zval=0;
end
hold on
title([{['AUC Phase-DVmTuning-VaRange ' strrep(DirOri2,'_','-')]},{['p=' num2str(p1)]},{['z=' num2str(zval)]}])
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

% Figure 8B
for i=1:3
    % Phase-Vm Tuning: devide into Stance Long vs Short
    for m=1:2, % 1: Vm, 2:Vf
        h=figure;
        hold on
        if m==1,
            tmplabel='Phase-Vm';
        else
            tmplabel='Phase-Vf';
        end
        title([tmplabel ' Tuning Stance Long vs Short ' LegLabel{i} ' ' DirOri2 ' n=' num2str(size(PhaseVmTuningStanceShortLongTotal{1,i,m},1)) 'cells'],'fontsize',16)
        
        PlotPatchSEM(PhaseBinPlot,PhaseVmTuningStanceShortLongTotal{1,i,m},ColorSet(i,:))
        PlotPatchSEM(PhaseBinPlot,PhaseVmTuningStanceShortLongTotal{2,i,m},ColorSet(i,:)/2)
        xlim([0 2*pi])
        set(gca,'XTick',[0 pi 2*pi])
        set(gca,'XTickLabel',{'0','pi','2pi'})
        copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution)
        
        if m==1,
            VmOffsetVfHigh=(PhaseVmTuningStanceShortLongTotal{1,i,m}(:,end)-PhaseVmTuningStanceShortLongTotal{1,i,m}(:,1))';
            VmOffsetVfLow=(PhaseVmTuningStanceShortLongTotal{2,i,m}(:,end)-PhaseVmTuningStanceShortLongTotal{2,i,m}(:,1))';
            
            [p1,ff,stats]=signrank(VmOffsetVfLow,VmOffsetVfHigh);
            try
                zval=stats.zval;
            catch
                zval=0;
            end
            BarGraphScatterDirect([VmOffsetVfLow;VmOffsetVfHigh],[ColorSet(i,:)/2;ColorSet(i,:)])
            title([{DirOri2 ' VmOffSet Stance Long vs Short'},{['N=' num2str(length(VmOffsetVfLow)) 'cells']},{['p=' num2str(p1)]},{['z=' num2str(zval)]}])
            set(gca,'Xtick',1:2)
            set(gca,'XtickLabel',[{'Long'},{'Short'}],'fontsize',14)
            ylabel('Vm Offset(mV)')
            copy_fig2pptx_opened_blank(1,300,ScreenY,FigResolution)
        end
    end
     
    % Figure 8A
    [p1,ff,stats]=signrank(MeanVfStanceShortLongTotal{1,i},MeanVfStanceShortLongTotal{2,i});
    try
        zval=stats.zval;
    catch
        zval=0;
    end
    figure;
    hold on;
    title([{DirOri2 ' MeanVf Stance Long vs Short'},{['N=' num2str(length(MeanVfStanceShortLongTotal{1,i})) 'cells']},{['p=' num2str(p1)]},{['z=' num2str(zval)]}])
    for n=1:length(MeanStancePeriodStanceShortLongTotal{1,i}),
        plot(MeanStancePeriodStanceShortLongTotal{1,i}(n),MeanVfStanceShortLongTotal{1,i}(n),'Color',ColorSet(1,:),'Marker','o','MarkerSize',14)
        plot(MeanStancePeriodStanceShortLongTotal{2,i}(n),MeanVfStanceShortLongTotal{2,i}(n),'Color',ColorSet(2,:),'Marker','o','MarkerSize',14)
        line([MeanStancePeriodStanceShortLongTotal{1,i}(n) MeanStancePeriodStanceShortLongTotal{2,i}(n)],[MeanVfStanceShortLongTotal{1,i}(n) MeanVfStanceShortLongTotal{2,i}(n)],'color','k')
    end
    copy_fig2pptx_opened_blank(1,900,ScreenY,FigResolution)
end


if CurrentButFiltRange(1)==0,
    [b,a]=butter(1,2*CurrentButFiltRange(2)/CompSampleRateDLC,'low');
elseif CurrentButFiltRange(2)==0,
    [b,a]=butter(1,2*CurrentButFiltRange(1)/CompSampleRateDLC,'high');
else
    [b,a]=butter(1,[2*CurrentButFiltRange(1)/CompSampleRateDLC 2*CurrentButFiltRange(2)/CompSampleRateDLC]);
end

vfFilt=filtfilt(bfilt,afilt,VfTotal);

% Figure 5E
% % the relation between dVm/dt at PreWindow (-200:0) and Drift
% % attenuation  at postWindow (0:200)
SearchWalkBoutsTimeTmp=0.2;
Idx=FindLocalPeak(vfFilt,VfFiltThr,1).*VaTotal>VaThr;
Idx(1:SearchWalkBoutsTimeTmp*CompSampleRateDLC)=0;
Idx(length(Idx)-SearchWalkBoutsTimeTmp*CompSampleRateDLC:end)=0;
Idx=find(Idx==1);

PreDVmDt=zeros(1,length(Idx));
PostDVa=zeros(1,length(Idx));
tmpFlag=zeros(1,length(Idx));
for j=1:length(Idx),
    % remove if it includes a boarder across concatenated segments
    tmpFlag(j)=sum(BoarderTotal(Idx(j)-SearchWalkBoutsTimeTmp*CompSampleRateDLC:Idx(j)+SearchWalkBoutsTimeTmp*CompSampleRateDLC));
    PreDVmDt(j)=mean(DVmTotal(Idx(j)-SearchWalkBoutsTimeTmp*CompSampleRateDLC:Idx(j)));
    PostDVa(j)=-mean(VaTotal(Idx(j):Idx(j)+SearchWalkBoutsTimeTmp*CompSampleRateDLC)-VaTotal(Idx(j)));
end

PreDVmDt(tmpFlag>0)=[];
PostDVa(tmpFlag>0)=[];

h=figure;
hold on
tmpx=PreDVmDt;
tmpxlabel='dVm/dt (mV/s)';
tmpy=PostDVa;
plot(tmpx,tmpy,'Color','k','LineStyle','o')
P = polyfit(tmpx,tmpy,1);

x = min(tmpx):(max(tmpx)-min(tmpx))/100:max(tmpx);

y = P(1)*x+P(2);

plot(x,y,'r','lineWidth',3)
xlabel(tmpxlabel)
ylabel('Drift attenuation [deg/s]')

x=get(gca,'xlim');
y=get(gca,'ylim');
[r,p]=corrcoef(tmpx,tmpy);
title([{['LocalVfMax' strrep(DirOri2,'_','-')]},{['R=' num2str(r(2,1),'%.2f') ':P=' num2str(p(2,1))]}])
set(gca,'box','off')
ylim([-100 100])
copy_fig2pptx_opened_blank(1,1000,1000,200);

% Some sub-functions used in the main function %%%%%%%%%%%%%%%%%%%%%%

%     function [OutVec,OutVec2] = FindLocalPeak(Vec,Thr,Sign)
%
%         % Find local peaks in a Vector (Vec) larger than a Threshold (Thr).
%         % Sign: 1 = local max, -1= local min
%
%         OutVec=zeros(1,length(Vec));
%
%         if length(Thr)==2,
%             if Sign>0
%                 OutVec=Vec>Thr(1) & Vec<Thr(2);
%             else
%                 OutVec=Vec>-Thr(2) & Vec<-Thr(1);
%             end
%         else
%             if Sign>0
%                 OutVec(Vec>Thr)=1;
%             else
%                 OutVec(Vec<-Thr)=1;
%             end
%         end
%
%         diffVec=[diff(Vec) 0];
%         if Sign>0
%             tmpIdx=diffVec.*circshift(diffVec,[0 1])<0 & diffVec<0;
%         else
%             tmpIdx=diffVec.*circshift(diffVec,[0 1])<0 & diffVec>0;
%         end
%         OutVec=OutVec.*tmpIdx;
%         OutVec2=diffVec.*circshift(diffVec,[0 1]);
%     end
% end

%
%     function PlotPatchSD_MeanSD(X,Mean,SD,varargin) % usually X=time, Y=data set, varargin = color
%
%         edge_transparency=0;%
%         face_transparency=0.2;%
%
%         X(isnan(Mean))=[];
%         Mean(isnan(Mean))=[];
%         SD(isnan(SD))=[];
%
%
%         MeanY=Mean;
%         StdY=SD;
%
%         err_x=[X fliplr(X)];
%         err_y=[MeanY+StdY fliplr(MeanY-StdY)];
%
%         if nargin > 3,
%             ColorSet=varargin{1};
%             plot(X,MeanY,'color',ColorSet,'LineWidth',2)
%             hold on
%             patch(err_x,err_y,ColorSet,'EdgeColor',ColorSet,'EdgeAlpha',edge_transparency,'FaceAlpha',face_transparency)
%         else
%             plot(X,MeanY,'k','LineWidth',2)
%             hold on
%             patch(err_x,err_y,'k','EdgeColor','k','EdgeAlpha',edge_transparency,'FaceAlpha',face_transparency)
%         end
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

% function PlotPatchSEM(X,Y,varargin) % usually X=time, Y=data set, varargin = color
% 
% edge_transparency=0;%
% face_transparency=0.2;%
% 
% MeanY=mean(Y,1);
% StdY=std(Y,0,1)/sqrt(size(Y,1));
% 
% err_x=[X fliplr(X)];
% err_y=[MeanY+StdY fliplr(MeanY-StdY)];
% 
% if nargin > 2,
%     ColorSet=varargin{1};
%     plot(X,MeanY,'color',ColorSet,'LineWidth',2)
%     hold on
%     patch(err_x,err_y,ColorSet,'EdgeColor',ColorSet,'EdgeAlpha',edge_transparency,'FaceAlpha',face_transparency)
% else
%     plot(X,MeanY,'k','LineWidth',2)
%     hold on
%     patch(err_x,err_y,'k','EdgeColor','k','EdgeAlpha',edge_transparency,'FaceAlpha',face_transparency)
% end
% end

