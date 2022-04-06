function OptogeneticHSActivation

% Folder saving Treadmill signals and Leg's position and phases
DirOri2='R27B03ADVT048588DBD_ReachR_CL';
DirOri='ParentFolderLocation';

% for plotting
TimeCut=[-0.2 0.2]; %sec

ResponseLatencyBaselineWindow=-0.1; % s consider -0.1 - 0 s as the baseline
SaveInfo.AveTime=0.5; % s for comparing the ipsi leg X position between pre and post
SearchWalkBoutsTime=0.5;

PreVfThr=5; % mm/s
CompSampleRate=100; % Hz

SwStDurLimit=500; % ms % remove too long stance durations
%--------------------------------
DirOri=strcat(DirOri,DirOri2,'\');

% power point 16:9 size
ScreenMag=0.9;
ScreenX=1000*ScreenMag;
ScreenY=1000*ScreenMag;
FigResolution=500;

SearchWalkBoutsTimeDP=SearchWalkBoutsTime*CompSampleRate;
TmpTime=(-SearchWalkBoutsTimeDP:SearchWalkBoutsTimeDP)/CompSampleRate;

ColorSet=[0.5 0.5 0.5;0.5 0.5 0.5;0.5 0.5 0.5;0.5 0.5 0.5;0.5 0.5 0.5];
%----------------------------------------
StimOnsetTrigIpsiLegXposTotal=[];
StimOnsetTrigVaNSTotal=[];

StimOnsetTrigVaNSSwingAtStimTotal=[];
StimOnsetTrigVaNSStanceAtStimTotal=[];

StanceDurationStimAtStanceTotal=[];
StanceDurationStimAtStanceShortPreStanceTotal=[];
StanceDurationStimAtStanceLongPreStanceTotal=[];

VaNSLatencyTotal=[];
VaNSLatencySwingAtStimTotal=[]; % ms
VaNSLatencyStanceAtStimTotal=[]; % ms

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
        
        % load Parameters
        %  load(strcat(CurrentDir,FlyIDSet{Fly},'\',DLCDir,'\','SaveInfoCombined.mat'));
        %  Va: 100Hz angular velocity
        %  VaNS: 100Hz non-smoothed Va (to precisely determine the response
        %  latency)
        %  Vf:100Hz forward velocity
        %  IpsiFrontLegXpos: X position of the ipsilateral front leg's tip
        %  StSwVec: Swing(=0)/Stance(=1) Chart
        %  OptStimEvents: =1 during optogenetic activation (50 ms/stim)
                
        StimOnsetIdx=zeros(1,length(vrCurrent));
        StimOnsetIdx([0 diff(OptStimEvents)]==1)=1;
        StimOnsetIdxDP=find(StimOnsetIdx==1);
        
        tmpIpsiLegXposMat=[];
        tmpVaMatNS=[];
        tmpVfMat=[];
        for j=1:length(StimOnsetIdxDP)
            tmpIpsiLegXposMat=[tmpIpsiLegXposMat; IpsiFrontLegXpos(StimOnsetIdxDP(j)-SearchWalkBoutsTimeDP:StimOnsetIdxDP(j)+SearchWalkBoutsTimeDP)];
            tmpVfMat=[tmpVfMat; Vf(StimOnsetIdxDP(j)-SearchWalkBoutsTimeDP:StimOnsetIdxDP(j)+SearchWalkBoutsTimeDP)];
            tmpVaMatNS=[tmpVaMatNS; VaNS(StimOnsetIdxDP(j)-SearchWalkBoutsTimeDP:StimOnsetIdxDP(j)+SearchWalkBoutsTimeDP)];
        end
        
        MidPoint=1+(size(tmpIpsiLegXposMat,2)-1)/2;
        % find stim events when the fly walked > 5mm/s before stim
        MeanPreVfStim=mean(tmpVfMat(:,1:(size(tmpVfMat,2)-1)/2),2)';
        PreVfIdx=find(MeanPreVfStim>PreVfThr);
        
        DiffStSwVec=[0 diff(StSwVec)];
        
        % -2, -1, 0, +1, +2
        StanceDurationStimAtStance=[];
        
        if ~isempty(PreVfIdx),
            for i=1:length(PreVfIdx),
                CurrentDP=StimOnsetIdxDP(PreVfIdx(i));
                % Ipsi front leg phase at the stim
                CurrentSwSt=StSwVec(CurrentDP);
                tmp=find(DiffStSwVec==1);
                PreStanceOnset=tmp(tmp<=CurrentDP);
                PostStanceOnset=tmp(tmp>CurrentDP);
                tmp=find(DiffStSwVec==-1);
                PreSwingOnset=tmp(tmp<=CurrentDP);
                PostSwingOnset=tmp(tmp>CurrentDP);
                if CurrentSwSt % stance at stim
                    try
                        StanceOnsetSeq=[PreStanceOnset(end-2:end) PostStanceOnset(1:2)];
                        SwingOnsetSeq=[PreSwingOnset(end-1:end) PostSwingOnset(1:3)];
                        StanceDuration=(SwingOnsetSeq-StanceOnsetSeq)*1000/CompSampleRate;
                        if max(StanceDuration)<SwStDurLimit % remove too large value
                            StanceDurationStimAtStance=[StanceDurationStimAtStance; StanceDuration];
                        end
                    catch
                    end
                end
            end
        end
        
        PreVfIdxVec=zeros(1,length(Vf));
        PreVfIdxVec(StimOnsetIdxDP(PreVfIdx))=1;
        tmpIdx=StimOnsetIdx.*PreVfIdxVec;
        
        % TrigAve
        Info.TrigMeanIpsiLegXpos=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        Info.TrigSDIpsiLegXpos=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        Info.TrigMeanVaNS=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        Info.TrigSDVaNS=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        
        KerMat=zeros(2*SearchWalkBoutsTime*CompSampleRate+1,length(Va));
        for j=1:2*SearchWalkBoutsTime*CompSampleRate+1,
            TmpShift=j-(1+SearchWalkBoutsTime*CompSampleRate);
            KerMat(j,:)=circshift(tmpIdx,[0 TmpShift]);
        end
        % cut start and end edge
        KerMat(:,1:SearchWalkBoutsTime*CompSampleRate)=0;
        KerMat(:,end-SearchWalkBoutsTime*CompSampleRate:end)=0;
        for j=1:2*SearchWalkBoutsTime*CompSampleRate+1,
            % front leg's tip
            TmpVec=KerMat(j,:).*IpsiFrontLegXpos;
            TmpVec=TmpVec(KerMat(j,:)==1);
            Info.TrigMeanIpsiLegXpos(j)=mean(TmpVec);
            Info.TrigSDIpsiLegXpos(j)=std(TmpVec);
            
            TmpVec=KerMat(j,:).*VaNS;
            TmpVec=TmpVec(KerMat(j,:)==1);
            Info.TrigMeanVaNS(j)=mean(TmpVec);
            Info.TrigSDVaNS(j)=std(TmpVec);
        end
        
        if sum(tmpIdx)>=MinDataNum
            % For figure 6K, Response Latency
            MidPoint=1+(length(Info.TrigMeanVaNS)-1)/2;
            tmpStart=find(TmpTime>=ResponseLatencyBaselineWindow,1);
            tmpMean=mean(Info.TrigMeanVaNS(tmpStart:MidPoint-1));
            tmpSD=std(Info.TrigMeanVaNS(tmpStart:MidPoint-1));
            tmpThr=tmpMean-2*tmpSD;
            tmpIdx=find(Info.TrigMeanVaNS(MidPoint+1:end)<tmpThr,1);
            if ~isempty(tmpIdx)
                VaNSLatencyTotal=[VaNSLatencyTotal tmpIdx*1000/CompSampleRate]; % ms
            end
            StimOnsetTrigIpsiLegXposTotal=[StimOnsetTrigIpsiLegXposTotal; Info.TrigMeanIpsiLegXpos-mean(Info.TrigMeanIpsiLegXpos(tmpStart:MidPoint-1))];
            StimOnsetTrigVaNSTotal=[StimOnsetTrigVaNSTotal; Info.TrigMeanVaNS];
        end
        
        
        % Stim at Stance
        tmpIdx=StimOnsetIdx.*PreVfIdxVec.*StSwVec;
        
        % TrigAve
        Info.TrigMeanIpsiLegXpos=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        Info.TrigSDIpsiLegXpos=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        Info.TrigMeanVaNS=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        Info.TrigSDVaNS=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        
        KerMat=zeros(2*SearchWalkBoutsTime*CompSampleRate+1,length(Va));
        for j=1:2*SearchWalkBoutsTime*CompSampleRate+1,
            TmpShift=j-(1+SearchWalkBoutsTime*CompSampleRate);
            KerMat(j,:)=circshift(tmpIdx,[0 TmpShift]);
        end
        % cut start and end edge
        KerMat(:,1:SearchWalkBoutsTime*CompSampleRate)=0;
        KerMat(:,end-SearchWalkBoutsTime*CompSampleRate:end)=0;
        for j=1:2*SearchWalkBoutsTime*CompSampleRate+1,
            % front leg's tip
            TmpVec=KerMat(j,:).*IpsiFrontLegXpos;
            TmpVec=TmpVec(KerMat(j,:)==1);
            Info.TrigMeanIpsiLegXpos(j)=mean(TmpVec);
            Info.TrigSDIpsiLegXpos(j)=std(TmpVec);
            
            TmpVec=KerMat(j,:).*VaNS;
            TmpVec=TmpVec(KerMat(j,:)==1);
            Info.TrigMeanVaNS(j)=mean(TmpVec);
            Info.TrigSDVaNS(j)=std(TmpVec);
        end
        
        tmpIdxTmp=tmpIdx;
        Info.TrigMeanIpsiLegXposTmp=Info.TrigMeanIpsiLegXpos;
        Info.TrigMeanVaNSTmp=Info.TrigMeanVaNS;
        
        % Stim at Swing
        tmpIdx=StimOnsetIdx.*PreVfIdxVec.*(~StSwVec);
        
        % TrigAve
        Info.TrigMeanIpsiLegXpos=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        Info.TrigSDIpsiLegXpos=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        Info.TrigMeanVaNS=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        Info.TrigSDVaNS=zeros(1,2*SearchWalkBoutsTime*CompSampleRate+1);
        
        KerMat=zeros(2*SearchWalkBoutsTime*CompSampleRate+1,length(Va));
        for j=1:2*SearchWalkBoutsTime*CompSampleRate+1,
            TmpShift=j-(1+SearchWalkBoutsTime*CompSampleRate);
            KerMat(j,:)=circshift(tmpIdx,[0 TmpShift]);
        end
        % cut start and end edge
        KerMat(:,1:SearchWalkBoutsTime*CompSampleRate)=0;
        KerMat(:,end-SearchWalkBoutsTime*CompSampleRate:end)=0;
        for j=1:2*SearchWalkBoutsTime*CompSampleRate+1,
            % front leg's tip
            TmpVec=KerMat(j,:).*IpsiFrontLegXpos;
            TmpVec=TmpVec(KerMat(j,:)==1);
            Info.TrigMeanIpsiLegXpos(j)=mean(TmpVec);
            Info.TrigSDIpsiLegXpos(j)=std(TmpVec);
            
            TmpVec=KerMat(j,:).*VaNS;
            TmpVec=TmpVec(KerMat(j,:)==1);
            Info.TrigMeanVaNS(j)=mean(TmpVec);
            Info.TrigSDVaNS(j)=std(TmpVec);
        end
        
        if sum(tmpIdx)>=MinDataNum && sum(tmpIdxTmp)>=MinDataNum
            % Response Latency
            MidPoint=1+(length(Info.TrigMeanVaNS)-1)/2;
            tmpStart=find(TmpTime>=ResponseLatencyBaselineWindow,1);
            BaselineWindow=tmpStart:MidPoint-1;
            tmpMean=mean(Info.TrigMeanVaNS(BaselineWindow));
            tmpSD=std(Info.TrigMeanVaNS(BaselineWindow));
            tmpThr=tmpMean-2*tmpSD;
            tmpIdx2=find(Info.TrigMeanVaNS(MidPoint+1:end)<tmpThr,1);
            MidPoint=1+(length(Info.TrigMeanVaNSTmp)-1)/2;
            tmpMean=mean(Info.TrigMeanVaNSTmp(BaselineWindow));
            tmpSD=std(Info.TrigMeanVaNSTmp(BaselineWindow));
            tmpThr=tmpMean-2*tmpSD;
            tmpIdx=find(Info.TrigMeanVaNSTmp(MidPoint+1:end)<tmpThr,1);
            
            if ~isempty(tmpIdx) && ~isempty(tmpIdx2)
                VaNSLatencySwingAtStimTotal=[VaNSLatencySwingAtStimTotal tmpIdx2*1000/CompSampleRate]; % ms
                VaNSLatencyStanceAtStimTotal=[VaNSLatencyStanceAtStimTotal tmpIdx*1000/CompSampleRate]; % ms
                StimOnsetTrigVaNSSwingAtStimTotal=[StimOnsetTrigVaNSSwingAtStimTotal; Info.TrigMeanVaNS];
                StimOnsetTrigVaNSStanceAtStimTotal=[StimOnsetTrigVaNSStanceAtStimTotal; Info.TrigMeanVaNSTmp];
            end
        end
        
        if ~isempty(StanceDurationStimAtStance)
            [h,L,MX,MED,bw]=violin(StanceDurationStimAtStance,{'-2','-1','0','1','2'});
            StanceDurationStimAtStanceTotal=[StanceDurationStimAtStanceTotal; MX];
            close;
            % further divide into if the preceding stance duration (at -1) is long or short
            tmpVec=StanceDurationStimAtStance(:,2)';
            [IdxX IdxY]=sort(tmpVec);
            StanceDurationStimAtStanceShortPreStanceTotal=[StanceDurationStimAtStanceShortPreStanceTotal; mean(StanceDurationStimAtStance(IdxY(1:round(IdxY/2)),:),1)];
            StanceDurationStimAtStanceLongPreStanceTotal=[StanceDurationStimAtStanceLongPreStanceTotal; mean(StanceDurationStimAtStance(IdxY(round(IdxY/2)+1:end),:),1)];
        end
    end
end

% Figure 6I
h=figure;
hold on
title([strrep(DirOri2,'_','-')])
tmpMean=mean(mean(StimOnsetTrigIpsiLegXposTotal(:,1:(size(StimOnsetTrigIpsiLegXposTotal,2)-1)/2),1),2);
PlotPatchSEM(TmpTime,StimOnsetTrigIpsiLegXposTotal,[0 1 1]);
set(gca,'box','off')
line(get(gca,'xlim'),[tmpMean tmpMean],'LineStyle','--','LineWidth',2,'color','k')
line([0 0],get(gca,'ylim'),'LineStyle','--','LineWidth',2,'color','g')
xlabel('Time since stim onset [s]')
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);


% Figure 6J
SaveInfo.TmpVaPre=mean(StimOnsetTrigIpsiLegXposTotal(:,MidPoint-SaveInfo.AveTime*CompSampleRate+1:MidPoint),2);
SaveInfo.TmpPost=mean(StimOnsetTrigIpsiLegXposTotal(:,MidPoint+1:MidPoint+SaveInfo.AveTime*CompSampleRate),2);
[p1,ff,stats]=signrank(SaveInfo.TmpVaPre',SaveInfo.TmpPost');
try
    zval=stats.zval;
catch
    zval=0;
end
[num2str(mean(SaveInfo.TmpVaPre)) ' ' num2str(std(SaveInfo.TmpVaPre)/sqrt(length(SaveInfo.TmpVaPre))) 'vs. '  num2str(mean(SaveInfo.TmpPost)) ' ' num2str(std(SaveInfo.TmpPost)/sqrt(length(SaveInfo.TmpPost)))]

BarGraphScatterDirect([SaveInfo.TmpVaPre SaveInfo.TmpPost]',[0.5 0.5 0.5;1 0 0])
hold on
title([{['IpsiLegX:' strrep(DirOri2,'_','-') ' n=' num2str(length(SaveInfo.TmpVaPre))]},{['pre:post=' num2str(SaveInfo.AveTime) 's']},{['p=' num2str(p1)]},{['z=' num2str(zval)]}])
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

BarGraphScatterDirect((SaveInfo.TmpPost-SaveInfo.TmpVaPre)',[0.5 0.5 0.5])
hold on
title([{['IpsiLegX:' strrep(DirOri2,'_','-') ' n=' num2str(length(SaveInfo.TmpVaPre))]},{['pre:post=' num2str(SaveInfo.AveTime) 's']},{['p=' num2str(p1)]},{['z=' num2str(zval)]}])
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

tmpStart=find(TmpTime>=TimeCut(1),1);
tmpEnd=find(TmpTime<=TimeCut(2),1,'last');
tmpWindow=tmpStart:tmpEnd;
TmpTimeCut=TmpTime(tmpWindow);
StimOnsetTrigVaNSStanceAtStimTotalCut=StimOnsetTrigVaNSStanceAtStimTotal(:,tmpWindow);
StimOnsetTrigVaNSSwingAtStimTotalCut=StimOnsetTrigVaNSSwingAtStimTotal(:,tmpWindow);

% Figure 6I
h=figure;
hold on
BaseSwing=repmat(mean(StimOnsetTrigVaNSSwingAtStimTotalCut(:,1:(size(StimOnsetTrigVaNSSwingAtStimTotalCut,2)-1)/2),2),1,size(StimOnsetTrigVaNSSwingAtStimTotal,2));
BaseStance=repmat(mean(StimOnsetTrigVaNSStanceAtStimTotalCut(:,1:(size(StimOnsetTrigVaNSStanceAtStimTotalCut,2)-1)/2),2),1,size(StimOnsetTrigVaNSStanceAtStimTotal,2));
PlotPatchSEM(TmpTime,StimOnsetTrigVaNSSwingAtStimTotal-BaseSwing,[1 0 1]);
PlotPatchSEM(TmpTime,StimOnsetTrigVaNSStanceAtStimTotal-BaseStance,[0 1 0]);
ylabel('Va[deg/s]')
set(gca,'box','off')
xlim([-SearchWalkBoutsTime SearchWalkBoutsTime])
line([0 0],get(gca,'ylim'),'LineStyle','--','LineWidth',2,'color','k')
xlabel('Time since stim onset [s]')
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

% Figure 6I, Inlet
h=figure;
hold on
title([strrep(DirOri2,'_','-') ' IpsLeg=Sw(m) vs. St(g) at Stim'])
BaseSwing=repmat(mean(StimOnsetTrigVaNSSwingAtStimTotalCut(:,1:(size(StimOnsetTrigVaNSSwingAtStimTotalCut,2)-1)/2),2),1,size(StimOnsetTrigVaNSSwingAtStimTotalCut,2));
BaseStance=repmat(mean(StimOnsetTrigVaNSStanceAtStimTotalCut(:,1:(size(StimOnsetTrigVaNSStanceAtStimTotalCut,2)-1)/2),2),1,size(StimOnsetTrigVaNSStanceAtStimTotalCut,2));
PlotPatchSEM(TmpTimeCut,StimOnsetTrigVaNSSwingAtStimTotalCut-BaseSwing,[1 0 1]);
PlotPatchSEM(TmpTimeCut,StimOnsetTrigVaNSStanceAtStimTotalCut-BaseStance,[0 1 0]);
ylabel('Va[deg/s]')
set(gca,'box','off')
line([0 0],get(gca,'ylim'),'LineStyle','--','LineWidth',2,'color','k')
xlabel('Time since stim onset [s]')
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

% Figure 6L
Data=StanceDurationStimAtStanceTotal';
BarGraphScatterDirectNoLine(Data,ColorSet)
hold on
h=gcf;
[DataSetLength,NumSetLength]=size(Data);
tmpColor=cool;
for i=1:DataSetLength,
    for j=1:NumSetLength
        plot(i,Data(i,j),'o','MarkerFaceColor',tmpColor(round(64*j/NumSetLength),:),'MarkerEdgeColor','none','LineWidth',2,'MarkerSize',10)
        if i<DataSetLength
            line([i i+1],[Data(i,j) Data(i+1,j)],'Color',tmpColor(round(64*j/NumSetLength),:),'LineWidth',2)
        end
    end
end
title([{['Stim at Stance ' strrep(DirOri2,'_','-')]},{['n=' num2str(size(StanceDurationStimAtStanceTotal,1))]}])
plot(1:size(StanceDurationStimAtStanceTotal,2),mean(StanceDurationStimAtStanceTotal,1),'r','LineWidth',5)
set(gca,'XTick',1:size(StanceDurationStimAtStanceTotal,2))
set(gca,'XTickLabel',{'-2','-1','0','1','2'})
xlabel('Stance No.')
ylabel('Stance Duration [ms]')
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

% Figure 6L
% compare stance duraion at -1 vs. 0
tmpx=StanceDurationStimAtStanceTotal(:,2)';
tmpy=StanceDurationStimAtStanceTotal(:,3)';
[p,hh,stats] =signrank(tmpx,tmpy);
try
    zval=stats.zval;
catch
    zval=0;
end
BarGraphScatterDirect([tmpx;tmpy],ColorSet)
hold on
ylabel('Stance Duration [ms]')
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'-1','0'})
title([{['StimAtSt (N=' num2str(length(tmpx)) ')'] } {['p=' num2str(round(p*100000)/100000) '(signrank)']} {['z=' num2str(zval)]}])
copy_fig2pptx_opened_blank(1,500,ScreenY,FigResolution);

% Figure 6M
% compare stance duraion at -1 vs. 0
tmpx=StanceDurationStimAtStanceShortPreStanceTotal(:,2)';
tmpy=StanceDurationStimAtStanceShortPreStanceTotal(:,3)';
[p,hh,stats] =signrank(tmpx,tmpy);
try
    zval=stats.zval;
catch
    zval=0;
end
BarGraphScatterDirect([tmpx;tmpy],ColorSet)
h=gcf;
hold on
ylabel('Stance Duration [ms]')
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'-1','0'})
title([{['ShortPre (N=' num2str(length(tmpx)) ')'] } {['p=' num2str(round(p*100000)/100000) '(signrank)']} {['z=' num2str(zval)]}])
copy_fig2pptx_opened_blank(1,500,ScreenY,FigResolution);


% Figure 6M
% compare stance duraion at -1 vs. 0
tmpx=StanceDurationStimAtStanceLongPreStanceTotal(:,2)';
tmpy=StanceDurationStimAtStanceLongPreStanceTotal(:,3)';
[p,hh,stats] =signrank(tmpx,tmpy);
try
    zval=stats.zval;
catch
    zval=0;
end
BarGraphScatterDirect([tmpx;tmpy],ColorSet)
h=gcf;
hold on
ylabel('Stance Duration [ms]')
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'-1','0'})
title([{['LongPre (N=' num2str(length(tmpx)) ')'] } {['p=' num2str(round(p*100000)/100000) '(signrank)']} {['z=' num2str(zval)]}])
copy_fig2pptx_opened_blank(1,500,ScreenY,FigResolution);

% Figure 6K
tmpx=VaNSLatencySwingAtStimTotal;
tmpy=VaNSLatencyStanceAtStimTotal;
[p,hh,stats] =signrank(tmpx,tmpy);
try
    zval=stats.zval;
catch
    zval=0;
end
BarGraphScatterDirect([tmpx; tmpy],[1 0 1;0 1 0])
hold on
title([{['Va Latency StimAt Sw vs. St']} {['n=' num2str(length(VaNSLatencySwingAtStimTotal))]} {strrep(DirOri2,'_','-')} {['p=' num2str(round(p*100000)/100000) '(signrank)']} {['z=' num2str(zval)]}])
ylabel('Latency [ms]')
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{'Sw','St'})
copy_fig2pptx_opened_blank(1,500,ScreenY,FigResolution);


% Some sub-functions used in the main function %%%%%%%%%%%%%%%%%%%%%%

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
