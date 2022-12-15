map.('T2')=nan(9,10);
map.('T3')=nan(9,10);
map.('T4')=nan(9,10);
peakmap.('T2')=nan(9,10);
peakmap.('T3')=nan(9,10);
peakmap.('T4')=nan(9,10);
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=[1 2 3 4 6 8 10] 
        currcheck=['C' num2str(current)];
        for trial=2:4
            trialcheck=['T' num2str(trial)];
            map.(trialcheck)(sepdist, current)=nanmean(saveCplit_all4shanks.(sepcheck).(currcheck).(trialcheck),'all');
            peakmap.(trialcheck)(sepdist, current)=nanmean(peak_all.(sepcheck).(currcheck)(3,:),'all');
        end
    end
end
nonan_map=map.T3;
nonan_map(isnan(nonan_map(:,1)),:)=[];
nonan_map(:,isnan(nonan_map(1,:)))=[];
figure
mesh([1 2 3 4 6 8 10],[300 400 500], nonan_map,'Facecolor','g')
xlabel('Current (\muA)')
ylabel('Sep dist (\mum)')
zlabel('Firing rate (sp/s)')
title('50/50')

hold on
nonan_map=map.T2;
nonan_map(isnan(nonan_map(:,1)),:)=[];
nonan_map(:,isnan(nonan_map(1,:)))=[];
mesh([1 2 3 4 6 8 10],[300 400 500], nonan_map,'Facecolor','r')
nonan_map=map.T4;
nonan_map(isnan(nonan_map(:,1)),:)=[];
nonan_map(:,isnan(nonan_map(1,:)))=[];
mesh([1 2 3 4 6 8 10],[300 400 500], nonan_map,'Facecolor','b')

nonan_map=peakmap.T3;
nonan_map(isnan(nonan_map(:,1)),:)=[];
nonan_map(:,isnan(nonan_map(1,:)))=[];
figure
mesh([1 2 3 4 6 8 10],[300 400 500], nonan_map,'Facecolor','g')
xlabel('Current (\muA)')
ylabel('Sep dist (\mum)')
zlabel('Peak location (\mum)')
title('50/50')

%%
clear peak_all ax1 ax2
vec = [100;80;50;30;15;0];
N = 128;
ploty=0;
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=[1 2 3 4 6 8 10]
        currcheck=['C' num2str(current)];

        for trial=1:5
            trialcheck=['T' num2str(trial)];
            for shank=0:3
                shankcheck=['S' num2str(shank)];
                if ploty==1
                figure(1+sepdist+shank*3)
                end
                if trial==1 && ploty==1
                    axes('Position',[0.13         0.112396822842341                     0.775         0.62])
                    hold on
                    set(gca,'FontSize',12)
                    axcheck=['s' num2str(shank)];
                    ax1.(axcheck)=gca;
                    ylim([0 350])
                    xlim([15 16+9+2])
                    xline(16,'Color',cmap(5*floor((length(cmap))/5),:))
                    xline(16+sepdist+1,'Color',cmap(1*floor((length(cmap))/5),:))
                    title([num2str(shank*200) '\mum'])
                    ylabel('Firing rate (sp/s)')
                    xlabel('Dist from electrode closest to array tip (\mum)')
                    xt = xticks;
                    xtl=(xt-16)*50;
                    xticklabels(xtl)
                    set(gca,'TickDir','out');
                elseif ploty==1
                     axcheck=['s' num2str(shank)];
                    axes(ax1.(axcheck));
                end

                shankcheck=['D' num2str(shank)];
                shanksepdist_dwnsample=saveshanksepdist.(sepcheck).(currcheck).(trialcheck).(shankcheck);%(:, samples_rand.(sepcheck));% downsample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ploty==1
                stdshade(shanksepdist_dwnsample',0.2,cmap(trial*floor((length(cmap))/5),:))
                end

                    
                    dat=shanksepdist_dwnsample;
                    dat(isinf(dat))=nan;
                    dat(sum(~isnan(dat),2)<3,:)=nan;
                    clear pairavgcur erpairavgcur
                    
                    for paircount=1:size(dat,2)
                        datnonan=dat(:,paircount);
                        firstnonan=find(~isnan(datnonan),1,'first');
                        lastnonan=find(~isnan(datnonan),1,'last');
                        datnonan(isnan(datnonan))=[];
                        if length(datnonan)<13
                            continue
                        end
                        FilterLength=5;
                        b=ones(FilterLength,1)./FilterLength;
                        smoothedenvelope=filtfilt(b,1,datnonan);
                        pairavgcur(:,paircount)=[nan(firstnonan-1,1); smoothedenvelope; nan(32-lastnonan,1)];
                    end
                    pairavgcur(:,sum(pairavgcur>0,1)==0)=nan;
                    [peak,c]=find(pairavgcur==max(pairavgcur));

                    peak_all.(sepcheck).(currcheck).(shankcheck)(trial,c)=(peak'-16).*50;
                    avgpcurr=nanmean(peak_all.(sepcheck).(currcheck).(shankcheck)(trial,c),2);
                    stdavg=nanstd(peak_all.(sepcheck).(currcheck).(shankcheck)(trial,c),[],2)./sqrt(sum(~isnan(peak_all.(sepcheck).(currcheck).(shankcheck)(trial,c)),2));

                    if ploty==1
                    figure(1+sepdist+shank*3)
                    if trial==1
                        axes('Position',[ax1.(axcheck).Position(1) .73 ax1.(axcheck).Position(3) .2])
                        ax2.(axcheck)=gca;
                        hold on
                        set(gca,'FontSize',12)
                        yticks([0 1 2 3 4 5 6 7 8])
                        labels=fliplr({'0:100','25:75','50:50','75:25','100:0'});
                        labels=labels([1:5]);
                        yticklabels([{''} labels {''}])
                        ylim([0.5 length([1:5])+0.5])
                        xlim((ax1.(axcheck).XLim-16).*50)
                        set(gca,'TickDir','out');
                        set(gca,'xtick',[])
                        set(gca,'xticklabel',[])
                        set(gca,'ytick',[])
                        set(gca,'yticklabel',[])
                        xline(0,'Color',cmap(5*floor((length(cmap))/5),:))
                        xline((sepdist+1)*50,'Color',cmap(1*floor((length(cmap))/5),:))

                    else
                        axes(ax2.(axcheck));
                    end
                    




                    color_current=cmap(trial*floor((length(cmap))/5),:);


                    er = errorbar(avgpcurr,trial,stdavg,stdavg,'horizontal');er.Color = [0, 0, 0]; er.LineStyle = 'none';
                    scatter(avgpcurr, trial, [], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', color_current)
                    %er = errorbar(avgpcurrc(colorbar),colorbar+5,stdavgc(trials(colorbar)),stdavgc(trials(colorbar)),'horizontal');er.Color = [0, 0, 0]; er.LineStyle = 'none';
                    %scatter(avgpcurrc(colorbar), colorbar+5, [], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', color_current)
                    end
                    if trial==1
                    numpairs.(sepcheck).(currcheck)(shank+1)=size(dat,2);
                    end

            end
        end
    end
end

%%
% is the peak shift significant
AMP=[1 2 3 4 6 8 10];
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for shank=0:3
            shankcheck=['D' num2str(shank)];
            [p,tbl,stats]=anova2(peak_all.(sepcheck).(currcheck).(shankcheck)',1,'off'); % for parametric
            pall.(sepcheck).(shankcheck) (current)=p(1);
            [p,tbl,stats]=friedman(peak_all.(sepcheck).(currcheck).(shankcheck)',1,'off'); %for non-para
            pall2.(sepcheck).(shankcheck) (current)=p;
        end
    end
end