function ChemogeneticHSSilencing

% Folder saving Treadmill signals and Membrane potentials
DirOri2='HS-Ort-1mM_10ms';
DirOri='ParentFolderLocation';

% time window in angular velocity pre and post injection
SaveInfo.AveTime=2; %sec

PreVfDur=1; % s

VfRange=[0.1 1;1 50];
VmSet=[1 2 3 4 5 6 7];
VfRangeMinTrialNum=10;

ColorSet=[230 159 0;0 114 178;204 121 167]/255;

% power point 16:9 size
ScreenMag=1;
ScreenX=700*ScreenMag;
ScreenY=1500*ScreenMag;
FigResolution=50;

CompSampleRate=500; % Hz
TrigTimeWindow=5; % s
%-------------------------------------------------------------------
VmSetTmp=[];
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

TotalLowVfOnsetTrigVm=[];
TotalLowVfOnsetTrigVr=[];
TotalLowVfOnsetTrigVf=[];
TotalHighVfOnsetTrigVm=[];
TotalHighVfOnsetTrigVr=[];
TotalHighVfOnsetTrigVf=[];

tmpTime=(-TrigTimeWindow*CompSampleRate:TrigTimeWindow*CompSampleRate)/CompSampleRate;

TmpCounter=1;
for Date=1:length(DateSet),
    
    CurrentDir=strcat(DirOri,'\',DateSet{Date},'\');
    
    DirTmp=dir(CurrentDir);
    DirTmpFlag=[];
    % check if these are directories
    for j=1:length(DirTmp),
        if DirTmp(j).isdir==0,
            DirTmpFlag=[DirTmpFlag j];
        end
    end
    
    DirTmp(DirTmpFlag)=[];
    FlyIDSet=sort_nat({DirTmp(3:end).name});
    
    for Fly=1:length(FlyIDSet),
        OnsetTrigVm=[];
        OnsetTrigVr=[];
        OnsetTrigVf=[];   
        
        % folder where Treadmill signals and Membrane potentials are saved
        CurrentDataDir=strcat(CurrentDir,'\',FlyIDSet{Fly},'\AllDataLabeled\DataTreadmill\');
        
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

        for File=1:FileNum,
            
            CurrentDataFile=strcat(CurrentDataDir,FileSeq{File});
            
            % load parameters
            % load(CurrentDataFile);
            % Va: Angular velocity  
            % Vf: Forward velocity
            % Vm: Membrane potentials
            % OnsetVec: Histamine Injection Events (second)

            if OnsetVec(1)~=0,
                for i=1:length(OnsetVec),
                    try
                        OnsetTrigVm=[OnsetTrigVm; Vm(round((OnsetVec(i)-TrigTimeWindow)*CompSampleRate):round((OnsetVec(i)-TrigTimeWindow)*CompSampleRate)+2*TrigTimeWindow*CompSampleRate)];
                        OnsetTrigVr=[OnsetTrigVr; Va(round((OnsetVec(i)-TrigTimeWindow)*CompSampleRate):round((OnsetVec(i)-TrigTimeWindow)*CompSampleRate)+2*TrigTimeWindow*CompSampleRate)];
                        OnsetTrigVf=[OnsetTrigVf; Vf(round((OnsetVec(i)-TrigTimeWindow)*CompSampleRate):round((OnsetVec(i)-TrigTimeWindow)*CompSampleRate)+2*TrigTimeWindow*CompSampleRate)];
                    catch
                    end           
                end
            end
        end
        
        OnsetTime=round(TrigTimeWindow*CompSampleRate)+1;
        tmpPreVf=mean(OnsetTrigVf(:,OnsetTime-PreVfDur*CompSampleRate:OnsetTime),2);
        tmpIdxLow=tmpPreVf>VfRange(1,1) & tmpPreVf<=VfRange(1,2);
        tmpIdxHigh=tmpPreVf>VfRange(2,1) & tmpPreVf<=VfRange(2,2);

        if sum(tmpIdxLow)>0 && sum(tmpIdxHigh)>0
            LowVfOnsetTrigVm=OnsetTrigVm(tmpIdxLow,:);
            LowVfOnsetTrigVr=OnsetTrigVr(tmpIdxLow,:);
            LowVfOnsetTrigVf=OnsetTrigVf(tmpIdxLow,:);
            HighVfOnsetTrigVm=OnsetTrigVm(tmpIdxHigh,:);
            HighVfOnsetTrigVr=OnsetTrigVr(tmpIdxHigh,:);
            HighVfOnsetTrigVf=OnsetTrigVf(tmpIdxHigh,:);

            if size(LowVfOnsetTrigVr,1)>VfRangeMinTrialNum && size(HighVfOnsetTrigVr,1)>VfRangeMinTrialNum        
                TotalLowVfOnsetTrigVm=[TotalLowVfOnsetTrigVm; mean(LowVfOnsetTrigVm,1)];
                TotalLowVfOnsetTrigVr=[TotalLowVfOnsetTrigVr; mean(LowVfOnsetTrigVr,1)];
                TotalLowVfOnsetTrigVf=[TotalLowVfOnsetTrigVf; mean(LowVfOnsetTrigVf,1)];
                TotalHighVfOnsetTrigVm=[TotalHighVfOnsetTrigVm; mean(HighVfOnsetTrigVm,1)];
                TotalHighVfOnsetTrigVr=[TotalHighVfOnsetTrigVr; mean(HighVfOnsetTrigVr,1)];
                TotalHighVfOnsetTrigVf=[TotalHighVfOnsetTrigVf; mean(HighVfOnsetTrigVf,1)];
                if ~isempty(find(VmSet==TmpCounter, 1))
                    VmSetTmp=[VmSetTmp size(TotalLowVfOnsetTrigVm,1)];
                end

            end   
        end
        TmpCounter=TmpCounter+1;
    end
end



MidPoint=round(size(TotalLowVfOnsetTrigVr,2)/2);

% Figure 1E, F 
SaveInfo.TmpLowVfPre=mean(TotalLowVfOnsetTrigVr(:,MidPoint-SaveInfo.AveTime*CompSampleRate+1:MidPoint),2);
SaveInfo.TmpLowVfPost=mean(TotalLowVfOnsetTrigVr(:,MidPoint+1:MidPoint+SaveInfo.AveTime*CompSampleRate),2);
SaveInfo.TmpHighVfPre=mean(TotalHighVfOnsetTrigVr(:,MidPoint-SaveInfo.AveTime*CompSampleRate+1:MidPoint),2);
SaveInfo.TmpHighVfPost=mean(TotalHighVfOnsetTrigVr(:,MidPoint+1:MidPoint+SaveInfo.AveTime*CompSampleRate),2);
[p]=signrank(SaveInfo.TmpLowVfPost-SaveInfo.TmpLowVfPre,SaveInfo.TmpHighVfPost-SaveInfo.TmpHighVfPre)
BarGraphScatterDirect([SaveInfo.TmpLowVfPost-SaveInfo.TmpLowVfPre SaveInfo.TmpHighVfPost-SaveInfo.TmpHighVfPre]',[ColorSet(2,:);ColorSet(1,:)])
title([{[strrep(DirOri2,'_','-') ' VfLowHigh']},{['N=' num2str(size(TotalLowVfOnsetTrigVr,1)) 'flies']},{['p=' num2str(p)]}])
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

% Figure 1E, F 
h=figure;
hold on
subplot1(3,1,'FontS',16,'Gap',[0.05 0]);
subplot1(1)
title(['Total LowHighVf:' strrep(DirOri2,'_','-') ' N=' num2str(size(TotalLowVfOnsetTrigVm,1)) 'flies'])
PlotPatchSEM(tmpTime,TotalLowVfOnsetTrigVm(VmSetTmp,:),ColorSet(2,:))
PlotPatchSEM(tmpTime,TotalHighVfOnsetTrigVm(VmSetTmp,:),ColorSet(1,:))
xlim([tmpTime(1) tmpTime(end)]);
xlabel('time [s]')
ylabel('Vm')
line([0 0],get(gca,'Ylim'),'color','k','LineStyle','--')
set(gca,'box','off')
subplot1(2)
hold on
title(['Total LowHighVf: Vr ' strrep(DirOri2,'_','-') ' N=' num2str(size(TotalLowVfOnsetTrigVr,1)) 'flies'])
PlotPatchSEM(tmpTime,TotalLowVfOnsetTrigVr,ColorSet(2,:))
PlotPatchSEM(tmpTime,TotalHighVfOnsetTrigVr,ColorSet(1,:))
xlim([tmpTime(1) tmpTime(end)]);
xlabel('time [s]')
ylabel('Vr')
line(get(gca,'Xlim'),[0 0],'color','k','LineStyle','--')
line([0 0],get(gca,'Ylim'),'color','k','LineStyle','--')
set(gca,'box','off')
subplot1(3)
hold on
title(['Total LowHighVf: Vf ' strrep(DirOri2,'_','-') ' N=' num2str(size(TotalLowVfOnsetTrigVf,1)) 'flies'])
PlotPatchSEM(tmpTime,TotalLowVfOnsetTrigVf,ColorSet(2,:))
PlotPatchSEM(tmpTime,TotalHighVfOnsetTrigVf,ColorSet(1,:))
xlim([tmpTime(1) tmpTime(end)]);
xlabel('time [s]')
ylabel('Vf')
line(get(gca,'Xlim'),[0 0],'color','k','LineStyle','--')
line([0 0],get(gca,'Ylim'),'color','k','LineStyle','--')
set(gca,'box','off')
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);


% Some sub-functions used in the main function %%%%%%%%%%%%%%%%%%%%%%

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