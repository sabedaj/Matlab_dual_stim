function Online_raster(ID,Spike_array)

sp = Spike_array;
trig = loadTrig(0);
TP = loadTrialParams;
%tID = find(cell2mat(TP(:,2)) == ID);
numelect=find(diff(cell2mat(TP(:,1)))~=0,1,'first');
theseTrig = trig(cell2mat(TP(1:numelect:end,2)) == ID)./30;
nT=length(theseTrig);
%% Set up the raster data structure
BIN = [-200 200]; 
SMOOTHING = 1; MAX = 700;
xdata = [];
ydata = [];


for tr = 1:nT
    theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
    for i = 1:length(theseSp)
        xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>.
       %figure(100);hold on; plot((sp(sp(:,1) > theseTrig(tr)+BIN(1) & sp(:,1) < theseTrig(tr)+BIN(2),2:end))')
        ydata = [ydata, tr];
    end
end

    figure; 
        figure;
    p = panel;
    p.pack(2,1)
    p(1,1).select();
    hold on; axM = gca;
    yScale = MAX./nT;
    for i = 1:size(xdata,2)
        line([xdata(i) xdata(i)],[ydata(i)-1/2.3 ydata(i)+1/2.3],'Color','k','LineWidth',2);
    end
   line([200 200],[0 nT],'Color','b');
      xmin=110;
    xmax=290;
        ylim(axM,[0 nT]); xlim(axM,[xmin xmax]);
    % Tidy Up
    xlabel(axM,'Time (ms)');
    ylabel(axM,'Trial');
    xticks(axM,xmin:30:xmax);
    xticklabels(axM,xmin-200:30:xmax-200);

    % Add the convolved spikerate
        p(2,1).select();
         axM2 = gca;
        hold on;
    Z = hist(xdata,0:400); %#ok<HIST>
    window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
    rate = (1000/nT)*conv(Z,window);
    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
    plot(axM2,rate,'LineWidth',2);
     line([200 200],[0 MAX],'Color','b');
    xmin=110;
    xmax=290;
    ylim(axM2,[0 MAX]); xlim(axM2,[xmin xmax]);
    % Tidy Up
    xlabel(axM2,'Time (ms)');
    ylabel(axM2,'Firing Rate (Sp/s)');
    xticks(axM2,xmin:30:xmax);
    xticklabels(axM2,xmin-200:30:xmax-200);
    %beautifyPlot(30,axM);
   set(gca,'TickDir','out');
end