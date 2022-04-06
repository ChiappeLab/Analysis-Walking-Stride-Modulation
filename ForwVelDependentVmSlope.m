function ForwVelDependentVmSlope

% Folder saving treadmill signals and membrane potentials
DirOri2='HS';
% DirOri2='LegPanSensoryTNT';
% DirOri2='LegPanSensoryCtl';
DirOri='ParentFolderLocation';

VfRangeSet={[1 3;10 50]};
ColorSet=[0 0 1;1 0 0];

VfRange=0:15;
VaRange=-400:400;
VfBin=1;
VaBin=50;
MinDataNum=10;

FigResolution=200;

CompSampleRate=500; % Hz
CompSampleRateTmp=100; % further down-sample from 500 -> 100 Hz
%---------------------------------
BoutVfThr=VfRangeSet{1}(end,1);
DirOri=strcat(DirOri,DirOri2,'\');

vfcent = VfRange(1): VfBin : VfRange(end);
vacent = VaRange(1): VaBin : VaRange(end);

CompSampleRateOri=CompSampleRate; % 500 Hz
CompRatio=round(CompSampleRateOri/CompSampleRateTmp);

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
BoutsNumTotal=[];
SlopeTotal=cell(1,size(VfRangeSet{1},1));

VmMatTotal = cell(length(vfcent), length(vacent));


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
        
        % Load parameters
        % load(strcat(CurrentDir,FlyIDSet{Fly},'\ConcatenatedBouts\ConcBouts.mat'));
        % Va: 500 Hz Angular velocity
        % Vf: 500 Hz Forward velocity
        % Vm: 500 Hz membrane potentials
        % BoutsStart, BoutsEnd: vectors indicating the start and
        % end time of each walking bout
        
        VaTmp=[];
        VfTmp=[];
        VmTmp=[];
        BoutsFlag=[];
        CurrentBoutsNum=0;
        for Bouts=1:length(BoutsStart),
            CurrentVfBout=Vf(BoutsStart(Bouts):BoutsEnd(Bouts));
            CurrentvfComp = zeros(1,floor(length(CurrentVfBout)/CompRatio));
            for g = 1 : length(CurrentvfComp)
                CurrentvfComp(g) = mean(CurrentVfBout((CompRatio * (g - 1) + 1):(CompRatio * (g - 1) + CompRatio)));
            end
            if max(CurrentvfComp)>=BoutVfThr
                VaTmp=[VaTmp Va(BoutsStart(Bouts):BoutsEnd(Bouts))];
                VfTmp=[VfTmp Vf(BoutsStart(Bouts):BoutsEnd(Bouts))];
                VmTmp=[VmTmp Vm(BoutsStart(Bouts):BoutsEnd(Bouts))];
                BoutsFlag=[BoutsFlag Bouts];
                CurrentBoutsNum=CurrentBoutsNum+1;
            end
        end
        
        Va=VaTmp;
        Vf=VfTmp;
        Vm=VmTmp;
        
            vaComp = zeros(1,floor(length(Va)/CompRatio));
            vfComp = zeros(1,floor(length(Vf)/CompRatio));
            mpComp = zeros(1,floor(length(Vm)/CompRatio));
            for g = 1 : length(vaComp)
                vaComp(g) = mean(Va((CompRatio * (g - 1) + 1):(CompRatio * (g - 1) + CompRatio)));
                vfComp(g) = mean(Vf((CompRatio * (g - 1) + 1):(CompRatio * (g - 1) + CompRatio)));
                mpComp(g) = mean(Vm((CompRatio * (g - 1) + 1):(CompRatio * (g - 1) + CompRatio)));
            end
            
            Va=vaComp;
            Vf=vfComp;
            Vm=mpComp;
            
            VmMat = cell(length(vfcent), length(vacent));
            for x = 1 : length(vfcent)
                tmpa=find(Vf >= vfcent(1)+(x-1)*VfBin);
                tmpb=find(Vf < vfcent(1)+x*VfBin);
                tmpx=tmpa(ismembc(tmpa,tmpb));
                for y = 1 : length(vacent)
                    tmpa=find(Va >= vacent(1)+(y-1)*VaBin);
                    tmpb=find(Va < vacent(1)+y*VaBin);
                    tmpy=tmpa(ismembc(tmpa,tmpb));
                    tmp=tmpx(ismembc(tmpx,tmpy));
                    if ~isempty(tmp)
                        VmMat{x,y}=Vm(tmp);
                    end
                end
            end
            
            NumPerPixel=cellfun(@numel,VmMat);
            SumPerPixel=cellfun(@sum,VmMat);
            SumPerPixel(NumPerPixel<MinDataNum)=0;
            NumPerPixel(NumPerPixel<MinDataNum)=0;
            
            for j=1:length(VfRangeSet)
                
                VfRangeSetCurrent=VfRangeSet{j};
                
                % check if there is any points at the highest
                % Vf velocity
                for i=size(VfRangeSetCurrent,1):size(VfRangeSetCurrent,1),
                    tmpa=find(vfcent > VfRangeSetCurrent(i,1));
                    tmpb=find(vfcent <= VfRangeSetCurrent(i,2));
                    tmp=tmpa(ismembc(tmpa,tmpb));
                    SlopeCurrnet=nansum(SumPerPixel(tmp,:),1)./nansum(NumPerPixel(tmp,:),1);
                    SlopeCurrnetTmp=SlopeCurrnet;
                    SlopeCurrnetTmp(isnan(SlopeCurrnet))=[];
                end
                
                if ~isempty(SlopeCurrnetTmp)
                    BoutsNumTotal=[BoutsNumTotal CurrentBoutsNum];
                    for i=1:size(VfRangeSetCurrent,1),
                        tmpa=find(vfcent > VfRangeSetCurrent(i,1));
                        tmpb=find(vfcent <= VfRangeSetCurrent(i,2));
                        tmp=tmpa(ismembc(tmpa,tmpb));
                        SlopeCurrnet=nansum(SumPerPixel(tmp,:),1)./nansum(NumPerPixel(tmp,:),1);
                        SlopeTotal{i}=[SlopeTotal{i}; SlopeCurrnet];
                    end
                    VmMatTotal =cellfun(@(x,y) cat(2,x,y),VmMatTotal,SpeedInfo.VmMat,'UniformOutput',false);
                end
            end
        TmpCounter=TmpCounter+1;
    end
end

% Figure 1C, 8G
labeltmp=[];
h=figure;
hold on
title([DirOri2 ' VfBin:' num2str(VfBin) '-VaBin:' num2str(VaBin) '-TimeDelay' num2str(TimeDelayOri) 'ms-N=' num2str(size(SlopeTotal{1},1)) 'cells'])
for i=1:size(VfRangeSetCurrent,1),
    SlopeTotalCurrent=SlopeTotal{i};
    tmpSEM=[];
    for j=1:size(SlopeTotalCurrent,2)
        tmpBoutsNumTotal=BoutsNumTotal;
        tmpBoutsNumTotal(isnan(SlopeTotalCurrent(:,j)))=[];
        tmpVec=SlopeTotalCurrent(:,j);
        tmpVec(isnan(tmpVec))=[];
        tmpSEM=[tmpSEM sqrt(var(tmpVec,tmpBoutsNumTotal/sum(tmpBoutsNumTotal)))/sqrt(length(tmpBoutsNumTotal))];
        tmpSEM(isnan(tmpSEM))=[];
        SlopeTotalCurrent(:,j)=SlopeTotalCurrent(:,j).*BoutsNumTotal'/sum(tmpBoutsNumTotal);
    end
    tnmpvaCent=vacent;
    tnmpvaCent(sum(~isnan(SlopeTotalCurrent),1)==0)=[];
    SlopeTotalCurrent(:,sum(~isnan(SlopeTotalCurrent),1)==0)=[];
    PlotPatchSD_MeanSD(tnmpvaCent,nansum(SlopeTotalCurrent),tmpSEM,ColorSet(i,:));
    labeltmp=[labeltmp strcat(num2str(VfRangeSetCurrent(i,1)),'-',num2str(VfRangeSetCurrent(i,2)),'-')];
end
xlabel('Va (deg/s)')
ylabel('DVm (mV)')
legend(label)
copy_fig2pptx_opened_blank(1,1000,1000,FigResolution);

% bootstrap to estimate the ofsset and slope difference
Stats.SubSampleSize=0.01;
Stats.NumRepeat=1000;
Stats.PSlopetotal=zeros(size(VfRangeSetCurrent,1),Stats.NumRepeat);
Stats.POffsettotal=zeros(size(VfRangeSetCurrent,1),Stats.NumRepeat);
Stats.Rsqtotal=zeros(size(VfRangeSetCurrent,1),Stats.NumRepeat);
figure
hold on
for i=1:size(VfRangeSetCurrent,1),
    tmpa=find(vfcent > VfRangeSetCurrent(i,1));
    tmpb=find(vfcent <= VfRangeSetCurrent(i,2));
    tmp=tmpa(ismembc(tmpa,tmpb));
    tmp=VmMatTotal(tmp,:);
    DataPointsVecCurrnet=arrayfun(@(x) cat(2,tmp{:,x}),1:size(tmp,2),'UniformOutput',false);
    tmpNumDataPoints=cellfun(@numel,DataPointsVecCurrnet);
    tmpx=[];
    for j=1:length(tmpNumDataPoints)
        tmpx=[tmpx vacent(j)*ones(1,tmpNumDataPoints(j))];
    end
    tmpy=cell2mat(arrayfun(@(x) cat(2,DataPointsVecCurrnet{x,:}),1,'UniformOutput',false));
    
    for n=1:Stats.NumRepeat
        % use 1% of data
        tmpIdx=randi([1 length(tmpx)],1,round(Stats.SubSampleSize*length(tmpx)));
        p=polyfit(tmpx(tmpIdx),tmpy(tmpIdx),1);
        Stats.PSlopetotal(i,n)=p(1);
        Stats.POffsettotal(i,n)=p(2);
        yfit=polyval(p,tmpx(tmpIdx));
        Stats.Rsqtotal(i,n)=1-sum((tmpy(tmpIdx)-yfit).^2)/((length(tmpy(tmpIdx))-1)*var(tmpy(tmpIdx)));
        yfit=polyval(p,vacent);
        plot(vacent,yfit,'color',ColorSet(i,:))
    end
end
copy_fig2pptx_opened_blank(1,1000,1000,FigResolution);
save('Save the parameter Stats');




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

% 
% function PlotPatchSD_MeanSD(X,Mean,SD,varargin) % usually X=time, Y=data set, varargin = color
% 
% edge_transparency=0;%
% face_transparency=0.2;%
% 
% X(isnan(Mean))=[];
% Mean(isnan(Mean))=[];
% SD(isnan(SD))=[];
% 
% 
% MeanY=Mean;
% StdY=SD;
% 
% err_x=[X fliplr(X)];
% err_y=[MeanY+StdY fliplr(MeanY-StdY)];
% 
% if nargin > 3,
%     ColorSet=varargin{1};
%     plot(X,MeanY,'color',ColorSet,'LineWidth',2)
%     hold on
%     patch(err_x,err_y,ColorSet,'EdgeColor',ColorSet,'EdgeAlpha',edge_transparency,'FaceAlpha',face_transparency)
% else
%     plot(X,MeanY,'k','LineWidth',2)
%     hold on
%     patch(err_x,err_y,'k','EdgeColor','k','EdgeAlpha',edge_transparency,'FaceAlpha',face_transparency)
% end

