function PlotCoherence

% Folder saving treadmill signals and membrane potentials
DirOri2='Folder'; 
DirOri='ParentFolderLocation';

CompSampleRateDown=100; % further down-sample from 500 -> 100 Hz
FigResolution=200;
%-----------------------------------------------
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

TmpCounter=1;

MscohereVrVmTotalTotal=[];
MscohereVfVmTotalTotal=[];
MscohereVrVrTotalTotal=[];
MscohereVfVfTotalTotal=[];
MscohereVmVmTotalTotal=[];

for Date=1:length(DateSet),
    
    CurrentDir=strcat(DirOri,DateSet{Date},'\');
    
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
        
        MscohereVrVmTotal=[];
        MscohereVfVmTotal=[];
        
        MscohereVrVrTotal=[];
        MscohereVfVfTotal=[];
        MscohereVmVmTotal=[];
        
        VaConc=[];
        VfConc=[];
        VmConc=[];
        VaConcOri=[];
        VfConcOri=[];
        VmConcOri=[];
        tmpCounter=0;
        
        for File=1:FileNum,
  % Load parameters
      %  Vm: 100Hz membrane potential
      %  Va: 100Hz angular velocity
      %  Vf: 100Hz forward velocity
%             CurrentDataFile=strcat(CurrentDataDir,FileSeq{File});
%             data=load(CurrentDataFile);
 
            VfOri=Vf;
            VaOri=Va;
            VmOri=Vm;
            
            % take out when the fly is not moving
            tmpIdx=Vf+Va==0;
            Vf(tmpIdx)=[];
            Va(tmpIdx)=[];
            Vm(tmpIdx)=[];
            
            if strcmp(DirOri2,'BoltChrimsonRightHS') % Bolt activation per trial is short -> concatenate 3 trials
                % for Figure 2G: Opt-Run Protocol 
                tmpCounter=tmpCounter+1;
                if tmpCounter==3,
                    tmpCounter=0;
                    [Cxy,Fmscohere]= mscohere(VaConc,VmConc,[],[],256,CompSampleRateDown);
                    MscohereVrVmTotal=[MscohereVrVmTotal; Cxy'];
                    [Cxy,Fmscohere]= mscohere(VfConc,VmConc,[],[],256,CompSampleRateDown);
                    MscohereVfVmTotal=[MscohereVfVmTotal; Cxy'];
                    x=VaConcOri;
                    N = length(x);
                    xdft = fft(x);
                    xdft = xdft(1:N/2+1);
                    psdx = (1/(CompSampleRateDown*N)) * abs(xdft).^2;
                    psdx(2:end-1) = 2*psdx(2:end-1);
                    freq = 0:CompSampleRateDown/length(x):CompSampleRateDown/2;
                    MscohereVrVrTotal=[MscohereVrVrTotal; 10*log10(psdx)];
                    x=VfConcOri;
                    N = length(x);
                    xdft = fft(x);
                    xdft = xdft(1:N/2+1);
                    psdx = (1/(CompSampleRateDown*N)) * abs(xdft).^2;
                    psdx(2:end-1) = 2*psdx(2:end-1);
                    MscohereVfVfTotal=[MscohereVfVfTotal; 10*log10(psdx)];
                    x=VmConcOri;
                    N = length(x);
                    xdft = fft(x);
                    xdft = xdft(1:N/2+1);
                    psdx = (1/(CompSampleRateDown*N)) * abs(xdft).^2;
                    psdx(2:end-1) = 2*psdx(2:end-1);
                    MscohereVmVmTotal=[MscohereVmVmTotal; 10*log10(psdx)];
                    VaConc=[];
                    VfConc=[];
                    VmConc=[];
                    VaConcOri=[];
                    VfConcOri=[];
                    VmConcOri=[];
                else
                    VaConc=[VaConc Va];
                    VfConc=[VfConc Vf];
                    VmConc=[VmConc Vm];
                    VaConcOri=[VaConcOri VaOri];
                    VfConcOri=[VfConcOri VfOri];
                    VmConcOri=[VmConcOri VmOri];
                end
            else
                % for Figure 2D: Spontaneous-walking Protocol 
                [Cxy,Fmscohere]= mscohere(Va,Vm,[],[],256,CompSampleRateDown);
                MscohereVrVmTotal=[MscohereVrVmTotal; Cxy'];
                [Cxy,Fmscohere]= mscohere(Vf,Vm,[],[],256,CompSampleRateDown);
                MscohereVfVmTotal=[MscohereVfVmTotal; Cxy'];
                x=VaOri;
                N = length(x);
                xdft = fft(x);
                xdft = xdft(1:N/2+1);
                psdx = (1/(CompSampleRateDown*N)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                freq = 0:CompSampleRateDown/length(x):CompSampleRateDown/2;
                MscohereVrVrTotal=[MscohereVrVrTotal; 10*log10(psdx)];
                x=VfOri;
                N = length(x);
                xdft = fft(x);
                xdft = xdft(1:N/2+1);
                psdx = (1/(CompSampleRateDown*N)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                MscohereVfVfTotal=[MscohereVfVfTotal; 10*log10(psdx)];
                x=VmOri;
                N = length(x);
                xdft = fft(x);
                xdft = xdft(1:N/2+1);
                psdx = (1/(CompSampleRateDown*N)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                MscohereVmVmTotal=[MscohereVmVmTotal; 10*log10(psdx)];
            end
        end
       
        MscohereVrVmTotal(isnan(MscohereVrVmTotal(:,1)),:)=[];
        MscohereVfVmTotal(isnan(MscohereVfVmTotal(:,1)),:)=[];
        MscohereVrVrTotal(isnan(MscohereVrVrTotal(:,1)),:)=[];
        MscohereVfVfTotal(isnan(MscohereVfVfTotal(:,1)),:)=[];
        MscohereVmVmTotal(isnan(MscohereVmVmTotal(:,1)),:)=[];

        MscohereVrVmTotalTotal=[MscohereVrVmTotalTotal; mean(MscohereVrVmTotal,1)];
        FmscohereOri=Fmscohere;
        MscohereVfVmTotalTotal=[MscohereVfVmTotalTotal; mean(MscohereVfVmTotal,1)];
        try
            MscohereVrVrTotalTotal=[MscohereVrVrTotalTotal; mean(MscohereVrVrTotal,1)];
            freqOri=freq;
        catch
            MscohereVrVrTotalTotal=[MscohereVrVrTotalTotal; interp1(freq,mean(MscohereVrVrTotal,1),freqOri)];
        end
        try
            MscohereVfVfTotalTotal=[MscohereVfVfTotalTotal; mean(MscohereVfVfTotal,1)];
        catch
            MscohereVfVfTotalTotal=[MscohereVfVfTotalTotal; interp1(freq,mean(MscohereVfVfTotal,1),freqOri)];
        end
        try
            MscohereVmVmTotalTotal=[MscohereVmVmTotalTotal; mean(MscohereVmVmTotal,1)];
        catch
            MscohereVmVmTotalTotal=[MscohereVmVmTotalTotal; interp1(freq,mean(MscohereVmVmTotal,1),freqOri)];
        end
        TmpCounter=TmpCounter+1;
    end
end

close all



h=figure;
hold on
tmpx=find(FmscohereOri>0.5,1);
PlotPatchSEM(FmscohereOri(tmpx:end)',MscohereVrVmTotalTotal(:,tmpx:end),[1 0 0]);
PlotPatchSEM(FmscohereOri(tmpx:end)',MscohereVfVmTotalTotal(:,tmpx:end),[0 0 1]);
title(['MSCoh ' strrep(DirOri2,'_','-')])
set(gca,'box','off')
xlabel('Freq (Hz)')
ylabel('MsCoh (-)')
xlim([0 20])
ylim([0 1])
copy_fig2pptx_opened_blank(1,1000,1000,FigResolution);


h=figure;
hold on
tmpx=find(freqOri>=0.5,1);
tmpy=find(freqOri>=20,1);
PlotPatchSEM(freqOri(tmpx:tmpy),MscohereVrVrTotalTotal(:,tmpx:tmpy),[1 0 0]);
PlotPatchSEM(freqOri(tmpx:tmpy),MscohereVfVfTotalTotal(:,tmpx:tmpy),[0 0 1]);
PlotPatchSEM(freqOri(tmpx:tmpy),MscohereVmVmTotalTotal(:,tmpx:tmpy),[0 1 0]);
title(['Pow. Spec. ' strrep(DirOri2,'_','-')])
set(gca,'box','off')
xlabel('Freq (Hz)')
ylabel('Power/Frequency (dB/Hz)')
xlim([0 20])
copy_fig2pptx_opened_blank(1,1000,1000,FigResolution);


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