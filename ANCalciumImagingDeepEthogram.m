function ANCalciumImagingDeepEthogram

% Folder saving treadmill signals and membrane potentials
DirOri2='R46A02VT023823-sytGC6s';
DirOri='ParentFolderLocation';

DEThrDur=0.5; % s
%-----------------------------------------------
ColorSet=[1 0 1;0.5 0.5 0.5];

DirOri=strcat(DirOri,DirOri2,'\');

% power point 16:9 size
ScreenMag=0.8;
ScreenX=1300*ScreenMag;
ScreenY=1000*ScreenMag;
FigResolution=100;

VideoRate=100; % Frames/s

WalkCaTotal=[];
GroomCaTotal=[];

WalkCaTotalDP=[];
GroomCaTotalDP=[];

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
        
        CaTotal=[];
        WalkTotal=[];
        GroomTotal=[];
        
        % Search for csv file from DeepEthogram
        CurrentDataDirDE=strcat(CurrentDir,'\',FlyIDSet{Fly},'\Deepethogram\');
        FileSeqOri=dir(CurrentDataDirDE);
        FileSeqDE=sort_nat({FileSeqOri(:).name});
        TmpFlag=zeros(1,length(FileSeqDE));
        for i=1:length(FileSeqDE),
            if ~isempty(FileSeqDE(isdir(char(FileSeqDE(i))))) || isempty(strfind(char(FileSeqDE(i)),'.csv')),
                TmpFlag(i)=1;
            end
        end
        FileSeqDE(find(TmpFlag==1))=[];
        
        for File=1:FileNum,
            CurrentDataFile=strcat(CurrentDataDir,FileSeq{File});
            % Load Imaging data
            % Ca: Calcium signal
            % CompSampleRate: Frame Rate
            load(CurrentDataFile);
            CaTotal=[CaTotal Ca];
           
            CurrentDataFileDE=strcat(CurrentDataDirDE,FileSeqDE{File});
           % load an output file from DeepEthogram
           % WalkVec: Walking events
           % GroomVec: Grooming events
            DEData=csvread(CurrentDataFileDE,1,0);
            WalkTotal=[WalkTotal WalkVec];
            GroomTotal=[GroomTotal GroomVec];
        end
        
        CompRatio=VideoRate/CompSampleRate;
        % find events more than thresholded duration
        % only include intiation of walking or grooming events
        
        % Normalize
        CaTotal=CaTotal/prctile(CaTotal,99.9);
        CaTotal(CaTotal>1)=nan;
        CaTotal=smooth(CaTotal,2)';
        
        WalkTotalFR=zeros(1,length(CaTotal));
        GroomTotalFR=zeros(1,length(CaTotal));
        
        for g=1:length(WalkTotalFR),
            WalkTotalFR(g)=median(WalkTotal(round(1+(g-1)*CompRatio):round(g*CompRatio)));
            GroomTotalFR(g)=median(GroomTotal(round(1+(g-1)*CompRatio):round(g*CompRatio)));
        end
        
        tmpIdx1=sum(MakeCirculantMatTeru(WalkTotalFR,round(CompSampleRate*DEThrDur),length(WalkTotalFR)),1);
        tmpIdx2=[0 diff(WalkTotalFR)];
        tmpIdxWalk=find(tmpIdx1==round(CompSampleRate*DEThrDur) & tmpIdx2==1);
        tmpIdxWalk(tmpIdxWalk>length(WalkTotalFR)-round(CompSampleRate*DEThrDur))=[];
        
        tmpIdx1=sum(MakeCirculantMatTeru(GroomTotalFR,round(CompSampleRate*DEThrDur),length(GroomTotalFR)),1);
        tmpIdx2=[0 diff(GroomTotalFR)];
        tmpIdxGroom=find(tmpIdx1==round(CompSampleRate*DEThrDur) & tmpIdx2==1);
        tmpIdxGroom(tmpIdxGroom>length(GroomTotalFR)-round(CompSampleRate*DEThrDur))=[];
        
        if ~isempty(tmpIdxWalk) && ~isempty(tmpIdxGroom)
            WalkCaMat=[];
            GroomCaMat=[];
            
            for i=1:length(tmpIdxWalk)
                WalkCaMat=[WalkCaMat; CaTotal(tmpIdxWalk(i):tmpIdxWalk(i)+round(CompSampleRate*DEThrDur)-1)];
            end
            for i=1:length(tmpIdxGroom)
                GroomCaMat=[GroomCaMat; CaTotal(tmpIdxGroom(i):tmpIdxGroom(i)+round(CompSampleRate*DEThrDur)-1)];
            end

            WalkCaTotal=[WalkCaTotal; mean(WalkCaMat,1)];
            GroomCaTotal=[GroomCaTotal; mean(GroomCaMat,1)];
            WalkCaTotalDP=[WalkCaTotalDP mean(WalkCaMat,2)'];
            GroomCaTotalDP=[GroomCaTotalDP mean(GroomCaMat,2)'];
        end
    end
end


% Figure 7I
MeanWalkCa=mean(WalkCaTotal,2);
GroomCa=mean(GroomCaTotal,2);

[p,hh,stats] =signrank(MeanWalkCa,GroomCa);
try
    zval=stats.zval;
catch
    zval=0;
end
BarGraphScatterDirectCellInput([{MeanWalkCa},{GroomCa}],[ColorSet(1,:);ColorSet(2,:)])
hold on
for i=1:length(MeanWalkCa)
    line([1 2],[MeanWalkCa(i) GroomCa(i)],'color','k')
end
ylabel('DF/F [Norm]')
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'Pre','Post'})
title([{['Walk vs. groom (N=' num2str(length(MeanWalkCa)) ')'] } {['p=' num2str(round(p*100000)/100000) '(signrank)']} {['z=' num2str(zval)]}])
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

% Figure 7H
xx = 0:0.05:1;
figure;
hold on
title(['Walk n=' num2str(length(WalkCaTotalDP))])
[heights,locations] = hist(WalkCaTotalDP,xx);
bar(locations,heights/sum(heights),'hist')
ylim([0 0.35])
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);

figure;
hold on
title(['Gloom n=' num2str(length(GroomCaTotalDP))])
[heights,locations] = hist(GroomCaTotalDP,xx);
bar(locations,heights/sum(heights),'hist')
ylim([0 0.35])
copy_fig2pptx_opened_blank(1,ScreenX,ScreenY,FigResolution);


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
