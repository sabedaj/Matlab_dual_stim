function fwhm(Sp)
Sp=Sp(:,2:end)';
sigmin_hm=(Sp(13,:))./2;

tmp=diff(Sp<sigmin_hm);
flag=false(size(tmp,2),1);
for i=1:size(tmp,2)
    posstart=13-find(flip(tmp(1:13,i)),1,'first')+2;%acctount for diff
    posend=13+find(tmp(14:end,i)==-1,1,'first');%acctount for diff
    if posend-posstart>6
        flag(i)=true;
    end
    
end
end