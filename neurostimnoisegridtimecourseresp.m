timec=0:10:100;%timecourse
sWin=20;
chns=[115];%channel of interest
arrayresp=cell(length(chns),3);
for chnit=1:length(chns)
for timeit=1:length(timec)%time iterate
%lg=timec(timeit)*0.001;
lg=timec(timeit);
 RF = getRFMap(sTrain_sitmulus, onsetInds_stimulus,onsetInds_baseline, sWin, sBin, lg);
 arrayrasp{chnit,1}(:,:,timeit)=RF;
% arrayrasp{chnit,2}(:,:,timeit)=xax;
% arrayrasp{chnit,3}(:,:,timeit)=yax;
 
%[sta,xax,yax] = d.getSTA('lag',lg,'channels', chns(chnit),'units',1);
% arrayrasp{chnit,1}(:,:,timeit)=sta;
% arrayrasp{chnit,2}(:,:,timeit)=xax;
% arrayrasp{chnit,3}(:,:,timeit)=yax;
end
end

chnchoose=1;
figure;
for timeit=1:length(timec)%time iterate
%imagesc(xax,yax,arrayresp{chnchoose}(:,:,timeit))
imagesc(xs,ys,arrayrasp{chnchoose}(:,:,timeit))
pause(0.5)
end