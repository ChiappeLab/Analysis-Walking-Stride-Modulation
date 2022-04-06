function ANCalciumImaging()

% Folder saving treadmill signals and membrane potentials
DirOri2='R46A02VT023823-sytGC6s';
DirOri='ParentFolderLocation';

% For Figure 7G
MinNumPoints=10;
TimeDelaySet=[666];
Vfbin=1;
Vabin=40;
VfStart=-4;
VfEnd=10;
VaStart=-400;
VaEnd=400;
Vacent = VaStart: Vabin : VaEnd;
Vfcent = VfStart: Vfbin : VfEnd;

FigResolution=200;
%-------------------------------------------------------------------

CRange=[0 1];

DirOri=strcat(DirOri,DirOri2,'\');

CaTotal=[];
VaTotal=[];
VfTotal=[];
CrossCovVrCaTotal=[];
CrossCovVfCaTotal=[];

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
        
        CaPerFly=[];
        VaPerFly=[];
        VfPerFly=[];
        for File=1:FileNum,
            
            CurrentDataFile=strcat(CurrentDataDir,FileSeq{File});
            % load parameters
            % data=load(CurrentDataFile);
            % CompSampleRate=data.FrameRate; % Hz
            
            % CompSampleRate: Calcium Imaging Rate
            % Ca: Calcium signal
            % Va: Angular Velocity (down-sampled to Calcium Imaging Rate)
            % Vf: Forward Velocity (down-sampled to Calcium Imaging Rate)
            
            CaPerFly=[CaPerFly Ca];
            VaPerFly=[VaPerFly Va];
            VfPerFly=[VfPerFly Vf];
            
        end
        
        CaPerFlyNorm=CaPerFly/prctile(CaPerFly,95);
        [cr, lags] = xcov(VaPerFly,CaPerFly,round(10*CompSampleRate),'coeff');
        [cf, lags] = xcov(VfPerFly,CaPerFly,round(10*CompSampleRate),'coeff');
        CrossCovVrCaTotal=[CrossCovVrCaTotal; cr];
        CrossCovVfCaTotal=[CrossCovVfCaTotal; cf];
        CaTotal=[CaTotal CaPerFlyNorm];
        VaTotal=[VaTotal VaPerFly];
        VfTotal=[VfTotal VfPerFly];
    end
end

% Figure 7F
[IdxX,IdxY]=max(mean(CrossCovVrCaTotal,1));
TimeMaxVr=1000*lags(IdxY)/CompSampleRate;
[IdxX,IdxY]=max(mean(CrossCovVfCaTotal,1));
TimeMaxVf=1000*lags(IdxY)/CompSampleRate;

h=figure;
hold on
title([{['Vr.Vf vs CaRatio Total ' DirOri2 ' N=' num2str(size(CrossCovVrCaTotal,1)) 'flies']},{['VrMax:' num2str(TimeMaxVr) 'ms VfMax:' num2str(TimeMaxVf) 'ms']}])
PlotPatchSEM(lags/CompSampleRate,CrossCovVrCaTotal,'r','LineWidth',3)
PlotPatchSEM(lags/CompSampleRate,CrossCovVfCaTotal,'b','LineWidth',2)
xlim([-10 10])
xlabel('lag [s]')
ylabel('Xcov')
line(get(gca,'Xlim'),[0 0],'color','k','LineStyle','--')
line([0 0],get(gca,'Ylim'),'color','k','LineStyle','--')
copy_fig2pptx_opened_blank(1,1000,1000,FigResolution);

% Figure 7G
for TimeDelayMode=1:length(TimeDelaySet), 
    TimeDelayOri=TimeDelaySet(TimeDelayMode);
    TimeDelay=round(TimeDelayOri*CompSampleRate/1000);
    CaCurrent=circshift(CaTotal,[0 -TimeDelay]);
    
    Ca = cell(length(Vfcent), length(Vacent));
    for x = 1 : length(Vfcent)
        tmpa=find(VfTotal >= Vfcent(1)+(x-1)*Vfbin);
        tmpb=find(VfTotal < Vfcent(1)+x*Vfbin);
        tmpx=tmpa(ismembc(tmpa,tmpb));
        for y = 1 : length(Vacent)
            tmpa=find(VaTotal >= Vacent(1)+(y-1)*Vabin);
            tmpb=find(VaTotal < Vacent(1)+y*Vabin);
            tmpy=tmpa(ismembc(tmpa,tmpb));
            tmp=tmpx(ismembc(tmpx,tmpy));
            if ~isempty(tmp)
                Ca{x,y}=CaCurrent(tmp);
            else
                Ca{x,y}=0;
            end
            
        end
    end
    CaM = cellfun(@mean,Ca);
    TmpNum = cellfun(@numel,Ca);
    MPM=CaM;
    MPM(TmpNum<MinNumPoints)=0;
    % Magenta-green
    MaskMat=MPM;
    MaskMat(MaskMat==0)=-1000;
    tmpy=size(magentamoss,1);
    MaskMat(MaskMat>=CRange(2)-(CRange(2)-CRange(1))/(tmpy-1))=CRange(2)-(CRange(2)-CRange(1))/(tmpy-1);
    MaskMat(MaskMat==-1000)=1000;
    
    imagesc(flipud(MaskMat));
    
    set(gca,'Ytick',1:2:length(Vfcent))
    set(gca,'YtickLabel',Vfcent(end):-2*Vfbin:Vfcent(1))
    set(gca,'Xtick',1:10:length(Vacent))
    set(gca,'XtickLabel',Vacent(1):10*Vabin:Vacent(end))
    set(gca,'box','off')
    
    colorbar
    title([{[num2str(TimeDelayOri) 'ms ' DirOri2 ' n=' num2str(size(CrossCovVrCaTotal,1)) 'flies']}])
    xlabel('Ang. Speed [deg/s]')
    ylabel('Forw. Speed [mm/s]')
    colormap magentamoss
    colordata = colormap;
    colordata(tmpy,:)=[0.5 0.5 0.5];
    colormap(colordata);
    set(gcf,'position',[50 100 900 700])
    caxis([0 1])
    copy_fig2pptx_opened_blank(1,1000,1000,FigResolution);
end

% Figure 7G
VfTuning=cellfun(@mean,cell2vec(Ca,0));
Center=find(Vfcent==0);
tmpVfTuning=VfTuning;
tmpVfTuning(Center)=[];
tmpVfcent=Vfcent;
tmpVfcent(Center)=[];
P = polyfit(tmpVfcent,tmpVfTuning,1);

figure
hold on
title([{[num2str(TimeDelayOri) 'ms ' DirOri2]}])
plot(Vfcent,VfTuning,'k','LineWidth',2)
plot(Vfcent,P(1)*Vfcent+P(2),'r','LineWidth',1)
xlim([Vfcent(1) Vfcent(end)])
line(get(gca,'xlim'),[VfTuning(Center) VfTuning(Center)],'LineStyle','--','LineWidth',1,'color','k')
line([0 0],get(gca,'ylim'),'LineStyle','--','LineWidth',1,'color','k')
xlabel('Vf [mm]')
ylabel('DCa [-]')
copy_fig2pptx_opened_blank(1,1000,1000,FigResolution);


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

