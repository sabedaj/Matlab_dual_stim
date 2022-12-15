function [centroidpos, mincentroidstart]=centroidFspiking(pairavgcur_input,sepdist)

                centroidpos1=16-find(~isnan(pairavgcur_input),1,'first');
                centroidpos16=abs(16+sepdist+1-find(~isnan(pairavgcur_input),1,'last'));
                mincentroidstart=min(centroidpos1,centroidpos16);
                %rate1=pairavgcur(:,paircount);%rate1=pairavgcur(16:16+sepdist+1,paircount);%
                rate1=pairavgcur_input(16-mincentroidstart:16+sepdist+1+mincentroidstart);
                rate1(rate1<0)=0;
                rate1(isnan(rate1))=[];
                %A = trapz(1:16, rate1);%A = trapz(1:sepdist+2, rate1);% 
                A = trapz(1:sepdist+2+mincentroidstart*2, rate1);
                %B=zeros(16,1);%B=zeros(sepdist+2,1);%
                B=zeros(sepdist+2+mincentroidstart*2,1);%
                for lims = 2:sepdist+2+mincentroidstart*2%2:16 %2:sepdist+2
                    B(lims) =  trapz(1:lims, rate1(1:lims));
                end
                [~,electrodecentroid]=min(abs(B-(A/2)));
                centroidpos=find(~isnan(pairavgcur_input),electrodecentroid,'first')+mincentroidstart;
                centroidpos=centroidpos(end);

end
