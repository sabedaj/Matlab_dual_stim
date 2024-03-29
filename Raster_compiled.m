function Raster_compiled(ID,chn,ratestruct)
Tnum=['T', num2str(ID)];
spkarray=squeeze(ratestruct.(Tnum)(chn,:,:));
nT=size(spkarray,2);
 MAX = 100; SMOOTHING = 1;
 xdata=[];
 ydata=[];
 for tr=1:nT
     xdata=[xdata;find(spkarray(:,tr)>0)];
     ydata = [ydata; tr*(MAX/nT).*ones(length(find(spkarray(:,tr)>0)),1)];
 end
    figure; hold on; axM = gca;
    yScale = MAX./nT;
    for i = 1:size(xdata,1)
        line([xdata(i) xdata(i)],[ydata(i)-yScale/2.3 ydata(i)+yScale/2.3],'Color','k','LineWidth',2);
    end
   line([91 91],[0 MAX],'Color','b');
    Z = hist(xdata,1:181); %#ok<HIST>
    window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
    rate = (1000/nT)*conv(Z,window);
    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
    plot(axM,rate,'LineWidth',2);


end