
function PassiveLegMotion

% Folder saving membrane potentials
DirOri2='LegServoContra';
DirOri='ParentFolderLocation';
% Folder saving videos (tracking fictive swing/stance)
VideoDirOri='VideoFolderLocation';

TimeCutRange=[4 11]; % for plotting
PhaseBin=0:pi/4:2*pi;

ColorSet=[230 159 0;0 114 178;204 121 167]/255;
StimStart=5; % s
StimDur=5; % s
%---------------------------------
PhaseBinPlot=PhaseBin(1:end-1)+(PhaseBin(2)-PhaseBin(1))/2;

DirOri=strcat(DirOri,DirOri2,'\');

% power point 16:9 size
ScreenMag=0.9;
ScreenX=1000*ScreenMag;
ScreenY=1000*ScreenMag;
FigResolution=100;

CompSampleRate=500; % Hz
CompSampleRateVideo=100;

CurrentButFiltRange=[5 0];
if CurrentButFiltRange(1)==0,
    [bfilt,afilt]=butter(1,2*CurrentButFiltRange(2)/CompSampleRate,'low');
elseif CurrentButFiltRange(2)==0,
    [bfilt,afilt]=butter(1,2*CurrentButFiltRange(1)/CompSampleRate,'high');
else
    [bfilt,afilt]=butter(1,[2*CurrentButFiltRange(1)/CompSampleRate 2*CurrentButFiltRange(2)/CompSampleRate]);
end


CompRatio=round(CompSampleRate/CompSampleRateVideo);

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

TmpCounter=1;

PhaseVmFiltTuningTotalNorm=[];
VmOsciAmpTotal=[];
MeanVmOsciAmp=[];
VmOsciAmpTotalBase=[];
MeanVmOsciAmpBase=[];

rng(1)
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
        
        FlyLabelCurrent=[strrep(DirOri2,'_','-') '-' strrep(DateSet{Date},'_','-') '-' strrep(FlyIDSet{Fly},'_','-')];
        CurrentDataDir=strcat(CurrentDir,'\',FlyIDSet{Fly},'\AllDataLabeled\DataTreadmill\');
        
        % load a video file (side-view video to determine fictive swing/stance phase)
        Video=mmread(strcat(VideoDirOri,'\',DateSet{Date},'\',FlyIDSet{Fly},'_MergedVideo\MergedVideo.avi'));
        
        FileSeqOri=dir(CurrentDataDir);
        FileSeq=sort_nat({FileSeqOri(:).name});
        TmpFlag=zeros(1,length(FileSeq));
        for i=1:length(FileSeq),
            if ~isempty(FileSeq(isdir(char(FileSeq(i))))) || isempty(strfind(char(FileSeq(i)),'.mat')),
                TmpFlag(i)=1;
            end
        end
        FileSeq(find(TmpFlag==1))=[];
        FileNum=length(FileSeq);
        
        VmTotal=[];
        VmFiltTotal=[];
        PhaseTotal=[];
        
        StanceOnsetVecTotal=[];
        StanceOffsetVecTotal=[];
        SwingOnsetVecTotal=[];
        SwingOffsetVecTotal=[];
        
        VmOsciAmp=[];
        VmOsciAmpBase=[];
        for File=1:FileNum,
            CurrentDataFile=strcat(CurrentDataDir,FileSeq{File});
            
            % load parameters
            % data=load(CurrentDataFile);
            % Vm: 500Hz membrane potentials
            % VmOri: 10Kz data for estimating Vm oscillation amplitude within a fictive stride
            
            VmFilt=filtfilt(bfilt,afilt,Vm);
            % subtract small Vm values for 1 s
            VmTmp=sort(Vm);
            base=mean(VmTmp(1:CompSampleRate));
            Vm=Vm-base;
            
            
            VmComp = zeros(1,floor(length(Vm)/CompRatio));
            VmCompFilt = zeros(1,floor(length(VmFilt)/CompRatio));
            for g = 1 : length(VmComp)
                VmComp(g) = mean(Vm((CompRatio * (g - 1) + 1):(CompRatio * (g - 1) + CompRatio)));
                VmCompFilt(g) = mean(VmFilt((CompRatio * (g - 1) + 1):(CompRatio * (g - 1) + CompRatio)));
            end
            
            VmTotal=[VmTotal; VmComp];
            VmFiltTotal=[VmFiltTotal; VmCompFilt];
            
            VideoLength=length(Vm)/(CompSampleRate/CompSampleRateVideo);
            % extract fictive Stance/Swing from the video
            VideoCurrent=Video(:,:,(File-1)*VideoLength+1:File*VideoLength);
            % only pick pixels where has a lot of variance (top quarter)
            VideoCurrent=imresize(VideoCurrent,0.25);
            
            
            SDPixel=std(VideoCurrent,0,3);
            [~,idx]=sort(SDPixel(:),'descend');
            [tmpx,tmpy]=ind2sub(size(SDPixel),idx);
            
            
            VideoMotionVec=[];
            for i=1:round(length(tmpx)/4),
                tmpVec=squeeze(VideoCurrent(tmpx(i),tmpy(i),:))';
                VideoMotionVec=[VideoMotionVec; abs([0 diff(tmpVec)])];
            end
            
            % detect fictive stance and swing segments
            SwingOnsetVec=[];
            
            MeanVideoMotionVec=mean(VideoMotionVec,1);
            StanceOnsetVec=[];
            StanceOffsetVec=[];
            for i=1:StimDur*5,
                [~,tmpIdx]=min(MeanVideoMotionVec((StimStart+(i-1)*0.2)*CompSampleRateVideo:(StimStart+(i-1)*0.2+0.05)*CompSampleRateVideo));
                StanceOnsetVec=[StanceOnsetVec round(tmpIdx(1)+(StimStart+(i-1)*0.2)*CompSampleRateVideo)-1];
                [~,tmpIdx]=min(MeanVideoMotionVec((StimStart+(i-1)*0.2+0.1)*CompSampleRateVideo:(StimStart+(i-1)*0.2+0.15)*CompSampleRateVideo));
                StanceOffsetVec=[StanceOffsetVec round(tmpIdx(1)+(StimStart+(i-1)*0.2+0.1)*CompSampleRateVideo)-1];
            end
            
            SwingOnsetVec=[];
            SwingOffsetVec=[];
            for i=1:StimDur*5,
                [~,tmpIdx]=min(MeanVideoMotionVec((StimStart+(i-1)*0.2+0.1)*CompSampleRateVideo:(StimStart+(i-1)*0.2+0.15)*CompSampleRateVideo));
                SwingOnsetVec=[SwingOnsetVec round(tmpIdx(1)+(StimStart+(i-1)*0.2+0.1)*CompSampleRateVideo)-1];
                [~,tmpIdx]=min(MeanVideoMotionVec((StimStart+(i-1)*0.2+0.2)*CompSampleRateVideo:(StimStart+(i-1)*0.2+0.25)*CompSampleRateVideo));
                SwingOffsetVec=[SwingOffsetVec round(tmpIdx(1)+(StimStart+(i-1)*0.2+0.2)*CompSampleRateVideo)-1];
            end
            
            StanceOnsetVecTotal=[StanceOnsetVecTotal; StanceOnsetVec];
            StanceOffsetVecTotal=[StanceOffsetVecTotal; StanceOffsetVec];
            SwingOnsetVecTotal=[SwingOnsetVecTotal; SwingOnsetVec];
            SwingOffsetVecTotal=[SwingOffsetVecTotal; SwingOffsetVec];
            
            
            VmOriRatio=10000/CompSampleRateVideo;
            for i=1:length(StanceOnsetVec)-1,
                CurrentMPtmp=VmOri(StanceOnsetVec(i)*VmOriRatio:StanceOnsetVec(i+1)*VmOriRatio);
                VmOsciAmp=[VmOsciAmp max(CurrentMPtmp)-min(CurrentMPtmp)];
                tmpIdx=randi([1 StimStart*10000-(StanceOnsetVec(i+1)*VmOriRatio-StanceOnsetVec(i)*VmOriRatio+1)]);
                CurrentMPtmpBase=VmOri(tmpIdx:tmpIdx+StanceOnsetVec(i+1)*VmOriRatio-StanceOnsetVec(i)*VmOriRatio);
                VmOsciAmpBase=[VmOsciAmpBase max(CurrentMPtmpBase)-min(CurrentMPtmpBase)];
            end
            
            Phase=-1*ones(1,length(VmComp));
            
            for i=1:length(SwingOnsetVec)
                Phase(SwingOnsetVec(i):SwingOffsetVec(i))=0:pi/(SwingOffsetVec(i)-SwingOnsetVec(i)):pi;
            end
            for i=1:length(StanceOnsetVec)
                Phase(StanceOnsetVec(i):StanceOffsetVec(i))=pi:pi/(StanceOffsetVec(i)-StanceOnsetVec(i)):2*pi;
            end
            
            PhaseTotal=[PhaseTotal; Phase];
            
        end
        MeanVmOsciAmp=[MeanVmOsciAmp mean(VmOsciAmp)];
        VmOsciAmpTotal=[VmOsciAmpTotal VmOsciAmp];
        MeanVmOsciAmpBase=[MeanVmOsciAmpBase mean(VmOsciAmpBase)];
        VmOsciAmpTotalBase=[VmOsciAmpTotalBase VmOsciAmpBase];
        
        mpTmp=reshape(VmTotal',1,size(VmTotal,1)*size(VmTotal,2));
        phaseTmp=reshape(PhaseTotal',1,size(PhaseTotal,1)*size(PhaseTotal,2));
        mpFiltTmp=reshape(VmFiltTotal',1,size(VmFiltTotal,1)*size(VmFiltTotal,2));
        
        PhaseVmTuning=cell(1,length(PhaseBin)-1);
        for j=1:length(PhaseBin)-1,
            tmpIdx=find(phaseTmp>=PhaseBin(j) & phaseTmp<PhaseBin(j+1));
            if ~isempty(tmpIdx)
                PhaseVmTuning{j}=mpTmp(tmpIdx);
            else
                PhaseVmTuning{j}=0;
            end
        end
        
        PhaseVmTuning=cell(1,length(PhaseBin)-1);
        for j=1:length(PhaseBin)-1,
            tmpIdx=find(phaseTmp>=PhaseBin(j) & phaseTmp<PhaseBin(j+1));
            if ~isempty(tmpIdx)
                PhaseVmTuning{j}=mpFiltTmp(tmpIdx);
            else
                PhaseVmTuning{j}=0;
            end
        end
        tmp=cellfun(@mean,PhaseVmTuning);
        PhaseVmFiltTuningTotalNorm=[PhaseVmFiltTuningTotalNorm; tmp/max(abs(tmp))];
        
        
        Time=(1:size(VmTotal,2))/CompSampleRateVideo;
        TimeCut=(TimeCutRange(1)*CompSampleRateVideo+1:TimeCutRange(2)*CompSampleRateVideo)/CompSampleRateVideo;
        VmTotalCut=VmTotal(:,TimeCutRange(1)*CompSampleRateVideo+1:TimeCutRange(2)*CompSampleRateVideo);
        
        % Figure 4D
        figure;
        hold on;
        title([FlyLabelCurrent])
        plot(TimeCut,mean(VmTotalCut,1),'k','LineWidth',4);
        xlabel('Time (s)')
        ylabel('Vm (mV)')
        ylim([min(mean(VmTotalCut,1)) max(mean(VmTotalCut,1))+0.1])
        y=get(gca,'ylim');
        ymax=y(2)-0.05;
        for i=1:size(StanceOnsetVecTotal,2),
            line([Time(round(mean(StanceOnsetVecTotal(:,i)))) Time(round(mean(StanceOffsetVecTotal(:,i))))],[ymax ymax],'Color',ColorSet(1,:),'LineWidth',4);
            line([Time(round(mean(SwingOnsetVecTotal(:,i)))) Time(round(mean(SwingOffsetVecTotal(:,i))))],[ymax ymax],'Color',ColorSet(2,:),'LineWidth',4);
        end
        copy_fig2pptx_opened_blank(1,1500,1000,FigResolution);
        
        % Only onset
        TmpCounter=TmpCounter-1;
    end
end

% GrandMean

% Figure 4F
h=figure;
hold on;
title(['Phase-VmFiltNorm Tuning N=' num2str(size(PhaseVmFiltTuningTotalNorm ,1)) 'cells'],'fontsize',16)
PlotPatchSD_MeanSD(PhaseBinPlot,MeanSkipZero(PhaseVmFiltTuningTotalNorm),MeanSkipZero(PhaseVmFiltTuningTotalNorm)/sqrt(size(PhaseVmFiltTuningTotalNorm,1)));
xlim([-0.3 PhaseBin(end)+0.3])
set(gca,'XTick',[0 pi 2*pi])
set(gca,'XTickLabel',{'0','pi','2pi'})
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

% Figure 4E
figure
hold on
title(['Osci Amp Stance-Stance (mV) n=' num2str(length(VmOsciAmpTotal))],'fontsize',16)
hist(VmOsciAmpTotal,[0:0.1:4])
xlim([0 4])
ylim([0 200])
set(gca,'box','off')
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

disp([num2str(mean(MeanVmOsciAmp)) '+-' num2str(std(MeanVmOsciAmp)/sqrt(length(MeanVmOsciAmp))) ' n=' num2str(length(MeanVmOsciAmp))])

figure
hold on
title(['Base Osci Amp Stance-Stance (mV) n=' num2str(length(VmOsciAmpTotalBase))],'fontsize',16)
hist(VmOsciAmpTotalBase,[0:0.1:4])
xlim([0 4])
ylim([0 200])
set(gca,'box','off')
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

disp([num2str(mean(MeanVmOsciAmpBase)) '+-' num2str(std(MeanVmOsciAmpBase)/sqrt(length(MeanVmOsciAmpBase))) ' n=' num2str(length(MeanVmOsciAmpBase))])


[p,ff,stats] = signrank(MeanVmOsciAmp,MeanVmOsciAmpBase);
BarGraphScatterDirectWithConnection([MeanVmOsciAmp;MeanVmOsciAmpBase],[1 0 0;0.5 0.5 0.5])
try
    zval=stats.zval;
catch
    zval=0;
end
ylabel('Mean Amp [mV]')
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'St','Base'})
title([{['N=' num2str(length(MeanVmOsciAmp))]} {['p=' num2str(p) ':z=' num2str(zval) '(ranksum)']}])
copy_fig2pptx_opened_blank(1,300,ScreenY,FigResolution)


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


