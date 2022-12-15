function [peak,lim1]=plotActivityAndPeak(Data,sepdist,trial,cmap,ax1,ax2,varargin)
hold(ax1, 'on');
hold(ax2, 'on');
%PEAK
[peak,c]=find(Data==max(Data));
peak(~diff(c))=[];

%Centroid
% centroid=zeros(size(Data,2),1);
% for columndata=1:size(Data,2)
%     [electrodecentroid, ~]=centroidFspiking(Data(:,columndata),sepdist);
%     centroid(columndata)=electrodecentroid;
% end
%peak=centroid;

if ~isempty(varargin) && strcmp(varargin{1},'laminar')
    [lineOut, ~] = stdshade(Data',0.2,cmap(trial*floor((length(cmap))/5),:),1:size(Data',2),1,ax1,'laminar');
    avgpcurr=mean(peak);
    stdavg=std(peak)./sqrt(length(peak));
    color_current=cmap(trial*floor((length(cmap))/5),:);
    er = errorbar(ax2,trial,avgpcurr,stdavg,stdavg,'vertical');er.Color = [0, 0, 0]; er.LineStyle = 'none';
    scatter(ax2, trial,avgpcurr, [], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', color_current)

    if sepdist==0 && trial==5
        fill(ax1,[ax1.YLim(2) ax1.YLim(1) ax1.YLim(1) ax1.YLim(2)],[14 14 16 16],'r','FaceAlpha', 0.2,'linestyle','none');
        fill(ax1,[ax1.YLim(2) ax1.YLim(1) ax1.YLim(1) ax1.YLim(2)],[22 22 24 24],'r','FaceAlpha', 0.2,'linestyle','none');
        fill(ax2,[ax2.YLim(2) ax2.YLim(1) ax2.YLim(1) ax2.YLim(2)],[14 14 16 16],'r','FaceAlpha', 0.2,'linestyle','none');
        fill(ax2,[ax2.YLim(2) ax2.YLim(1) ax2.YLim(1) ax2.YLim(2)],[22 22 24 24],'r','FaceAlpha', 0.2,'linestyle','none');
        sepdist=5;
    elseif sepdist~=0
        yline(ax1, 16,'Color',cmap(5*floor((length(cmap))/5),:))
        yline(ax1, 16+sepdist+1,'Color',cmap(floor((length(cmap))/5),:))
        yline(ax2, 16,'Color',cmap(5*floor((length(cmap))/5),:))
        yline(ax2, 16+sepdist+1,'Color',cmap(floor((length(cmap))/5),:))
    end

    %ax1
    DATApos1=find(~isnan(lineOut.XData),1,'first');
    DATAposend=find(~isnan(lineOut.XData),1,'last');
    minax=min([(16-DATApos1) (DATAposend-(16+sepdist+1))]);
    if minax>0
        ylim(ax1,[16-minax (16+sepdist+1)+minax])
    else
        ylim(ax1,[DATApos1 DATAposend])
    end
    xt = yticks(ax1);
    xtl=(xt-16)*-50;
    yticklabels(ax1,xtl)
    xlabel(ax1,'Sp/s')
    ylabel(ax1,'Distance from deepest stim elect (\mum)')
    set(ax1,'TickDir','out');
    %set (ax1, 'ydir', 'reverse' )%HORIZONTAL
    %ylim(ax1,[0 0.6])
    hold(ax1,'off')
    %ax2
    if minax>0
        ylim(ax2,[16-minax (16+sepdist+1)+minax])
    else
        ylim(ax2,[DATApos1 DATAposend])
    end
    set(ax2,'ytick',[])
    set(ax2,'yticklabel',[])
    set(ax2,'xtick',[])
    set(ax2,'xticklabel',[])
    if trial==5
        %     text(ax2,1,(16-DATApos1)+9.5,['N=' num2str(size(Data,2)) ' shanks'])
        %     text(ax2,1,(16-DATApos1)+11.5,['N=' num2str(size(Data,2)/4) ' pairs'])
        lim1=(16-DATApos1);
    else
        lim1=0;
    end
else
    [lineOut, ~] = stdshade(Data',0.2,cmap(trial*floor((length(cmap))/5),:),1:size(Data',2),1,ax1);
     avgpcurr=mean(peak);
    stdavg=std(peak)./sqrt(length(peak));
    color_current=cmap(trial*floor((length(cmap))/5),:);
    er = errorbar(ax2,avgpcurr,trial,stdavg,stdavg,'horizontal');er.Color = [0, 0, 0]; er.LineStyle = 'none';
    scatter(ax2, avgpcurr, trial, [], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', color_current)

    if sepdist==0 && trial==5
        fill(ax1,[ax1.YLim(2) ax1.YLim(1) ax1.YLim(1) ax1.YLim(2)],[14 14 16 16],'r','FaceAlpha', 0.2,'linestyle','none');
        fill(ax1,[ax1.YLim(2) ax1.YLim(1) ax1.YLim(1) ax1.YLim(2)],[22 22 24 24],'r','FaceAlpha', 0.2,'linestyle','none');
        fill(ax2,[ax2.YLim(2) ax2.YLim(1) ax2.YLim(1) ax2.YLim(2)],[14 14 16 16],'r','FaceAlpha', 0.2,'linestyle','none');
        fill(ax2,[ax2.YLim(2) ax2.YLim(1) ax2.YLim(1) ax2.YLim(2)],[22 22 24 24],'r','FaceAlpha', 0.2,'linestyle','none');
        sepdist=5;
    elseif sepdist~=0
        xline(ax1, 16,'Color',cmap(5*floor((length(cmap))/5),:))
        xline(ax1, 16+sepdist+1,'Color',cmap(floor((length(cmap))/5),:))
        xline(ax2, 16,'Color',cmap(5*floor((length(cmap))/5),:))
        xline(ax2, 16+sepdist+1,'Color',cmap(floor((length(cmap))/5),:))
    end

    %ax1
    DATApos1=find(~isnan(lineOut.YData),1,'first');
    DATAposend=find(~isnan(lineOut.YData),1,'last');
    minax=min([(16-DATApos1) (DATAposend-(16+sepdist+1))]);
    if minax>0
        xlim(ax1,[16-minax (16+sepdist+1)+minax])
    else
        xlim(ax1,[DATApos1 DATAposend])
    end
    yt = xticks(ax1);
    ytl=(yt-16)*-50;
    xticklabels(ax1,ytl)
    ylabel(ax1,'Sp/s')
    xlabel(ax1,'Distance from deepest stim elect (\mum)')
    set(ax1,'TickDir','out');
    set (ax1, 'xdir', 'reverse' )%HORIZONTAL
    ylim(ax1,[0 1])
    hold(ax1,'off')
    %ax2
    if minax>0
        xlim(ax2,[16-minax (16+sepdist+1)+minax])
    else
        xlim(ax2,[DATApos1 DATAposend])
    end
    set(ax2,'xtick',[])
    set(ax2,'xticklabel',[])
    set(ax2,'ytick',[])
    set(ax2,'yticklabel',[])
    set (ax2, 'xdir', 'reverse' )%HORIZONTAL
    if trial==5
        %     text(ax2,1,(16-DATApos1)+9.5,['N=' num2str(size(Data,2)) ' shanks'])
        %     text(ax2,1,(16-DATApos1)+11.5,['N=' num2str(size(Data,2)/4) ' pairs'])
        lim1=(16-DATApos1);
    else
        lim1=0;
    end
end
%set (ax2, 'ydir', 'reverse' )%HORIZONTAL
hold(ax2,'off')
end