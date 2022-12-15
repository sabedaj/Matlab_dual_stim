%% hypoth 1 - additive envelope
clear centroidpos_all peak_all
singlecond=0;%picks the column unless 0 then does average for all
multiply1=0.5;%modify this if you go to half amplitude
AMP=[0 1 2 3 4 6 8 10];%
Csplit_dwnsample=[];
Envelope_add_dwnsmple=[];
seedpoint=65;%20,55
s = RandStream('mlfg6331_64','Seed',seedpoint);
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    figure(sepdist+10)
    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        input1=saveCplit.(checksep).(checkcurr).T1;%./max(saveCplit.(checksep).(checkcurr).T1);
        %input1(isinf(input1))=0;
        input2=saveCplit.(checksep).(checkcurr).T5;%./max(saveCplit.(checksep).(checkcurr).T5);
        %input2(isinf(input2))=0;
        Envelope_add.(checksep).(checkcurr).T3=input1.*multiply1+input2.*multiply1;
        %             for trial=1:5
        %                 check=['T' num2str(trial)];
        dat=saveCplit.(checksep).(checkcurr).T3;
        dat2=Envelope_add.(checksep).(checkcurr).T3;
        dat(isinf(dat))=nan;
        clear pairavgcur erpairavgcur pairavgcursmooth_env pairavgcur_env
        for paircount=1:size(dat,2)/4
            pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*4+1:4* paircount),2);
            pairavgcur_env(:,paircount)=nanmean(dat2(:,( paircount-1)*4+1:4* paircount),2);
            [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur_env(:,paircount),sepdist);
            centroidpos_all.(checksep).(checkcurr)(paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%
            % peak calc
            datnonan=pairavgcur_env(:,paircount);
            firstnonan=find(~isnan(datnonan),1,'first');
            lastnonan=find(~isnan(datnonan),1,'last');
            datnonan(isnan(datnonan))=[];
            FilterLength=5;
            b=ones(FilterLength,1)./FilterLength;
            smoothedenvelope=filtfilt(b,1,datnonan);
            pairavgcursmooth_env(:,paircount)=[nan(firstnonan-1,1); smoothedenvelope; nan(32-lastnonan,1)];

            %centroidpos_all.(checksep).(checkcurr)(paircount)=(centroidpos-16).*50;

        end
        lengthneeded=size(saveCplit.sep7.C1.T1,2)/4;
        [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
        samples_rand.(checksep)=samples_rand1;
        centroidpos_all.(checksep).(checkcurr)=centroidpos_all.(checksep).(checkcurr)(:, samples_rand.(checksep));
        centroidpos_all.(checksep).(checkcurr)(sum(~isnan(centroidpos_all.(checksep).(checkcurr)),2)<10,:)=nan;
        pairavg_dwnsample=pairavgcur(:, samples_rand.(checksep));
        %pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan;
        Csplit_dwnsample.(checksep).(checkcurr).T3=pairavg_dwnsample;

        pairavgcur_env_dwnsmple=pairavgcur_env(:, samples_rand.(checksep));
        %pairavgcur_env_dwnsmple(sum(~isnan(pairavgcur_env_dwnsmple),2)<10,:)=nan;
        Envelope_add_dwnsmple.(checksep).(checkcurr).T3=pairavgcur_env_dwnsmple;

        %peak calc
        pairavgcursmooth_env=pairavgcursmooth_env(:, samples_rand.(checksep));
        pairavgcursmooth_env(:,sum(pairavgcursmooth_env>0,1)==0)=nan;
        [peak,c]=find(pairavgcursmooth_env==max(pairavgcursmooth_env));
        peak_all.(checksep).(checkcurr)=(peak'-16).*50;

        if length(AMP)>1
            subplot(2,4,current)
        end
        if singlecond==0
            %Envelope_add.(checksep).(checkcurr).T3(sum(Envelope_add.(checksep).(checkcurr).T3>0,2)<15,:)=0; %if you want to take out edge conditions
            % Envelope_add.(checksep).(checkcurr).T3(Envelope_add.(checksep).(checkcurr).T3==0)=nan;
            stdshade(Envelope_add_dwnsmple.(checksep).(checkcurr).T3',0.2,'b');
            hold on
            %saveCplit.(checksep).(checkcurr).T3(sum(saveCplit.(checksep).(checkcurr).T3>0,2)<15,:)=0; %if you want to take out edge conditions
            %saveCplit.(checksep).(checkcurr).T3(saveCplit.(checksep).(checkcurr).T3==0)=nan;

           stdshade(Csplit_dwnsample.(checksep).(checkcurr).T3',0.2,'k');


            %plot(mean(saveCplit.(checksep).(checkcurr).T1.*multiply1,2),'Color','k')
            %plot(mean(saveCplit.(checksep).(checkcurr).T5.*multiply1,2),'Color','k')
        else
            plot(Envelope_add_dwnsmple.(checksep).(checkcurr).T3(:,singlecond),'b')
            hold on
            plot(Csplit_dwnsample.(checksep).(checkcurr).T3(:,singlecond),'r')
            plot(Csplit_dwnsample.(checksep).(checkcurr).T1(:,singlecond),'Color','k')
            plot(Csplit_dwnsample.(checksep).(checkcurr).T5(:,singlecond),'Color','k')
        end
        title([num2str(AMP(current)) '\muA'])

        xline(16,'r')
        xline(16+sepdist+1,'r')
        if sepdist==5
            xlim([16-3 16+sepdist+1+3])%([find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'first'), find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'last')])
        elseif sepdist==7
            xlim([16-2 16+sepdist+1+2])
        elseif sepdist==9
            xlim([16-1 16+sepdist+1+1])
        end
        xt = xticks;
        xtl=(xt-16)*50;
        xticklabels(xtl)
        ylim([0 350])
        ylabel('Sp/s')
        xlabel('Distance from deepest stim elect (\mum)')
        set(gca,'TickDir','out');

        %             MSE.(checksep)(current)=(sum((saveCplit.(checksep).(checkcurr).T3(:)-Envelope_add.(checksep).(checkcurr).T3(:)).^2))/length(saveCplit.(checksep).(checkcurr).T3(:));
        %             rsquared.(checksep)(current)=1-(sum((saveCplit.(checksep).(checkcurr).T3(:)-(Envelope_add.(checksep).(checkcurr).T3(:))).^2,'omitnan')./sum((saveCplit.(checksep).(checkcurr).T3(:)-nanmean(saveCplit.(checksep).(checkcurr).T3(:))).^2,'omitnan'));
        %

        y=saveCplit.(checksep).(checkcurr).T3;%realdata
        yinput=saveCplit.(checksep).(checkcurr).T1;%inputdata1
        yfit=Envelope_add.(checksep).(checkcurr).T3;%modeldata
        numparam=0;
        [rsq,AIC,MSE]=modelfitcharacterisation(y,yfit,yinput,numparam);
        rsquared.(checksep)(current)=rsq;
        AIC_all.(checksep)(current)=AIC;
        MSE_all.(checksep)(current)=MSE;



        %             numparam=0;%number of input parameters to the model
        %             NumObservations=sum(~isnan(yinput),'all');
        %             DFE=NumObservations-numparam;%degrees freedom
        %             MSE=(sum((y(:)-yfit(:)).^2,'omitnan'))/length(y);%mean square error
        %             Var = DFE/NumObservations*MSE;
        %             L = -( sum((y-yfit).^2,'omitnan')/Var + NumObservations*log(2*pi) + NumObservations*log(Var))/2;%-ve loglikelihood
        %             AIC=2*(numparam-L);
        %             rsq=1-(sum((y-yfit).^2,'omitnan')./sum((y-nanmean(y)).^2,'omitnan'));


    end
    f=get(gca,'Children');
    %legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
end

%ylim([])
b1centroid_5050.sep5=centroidpos_all.sep5.C6;
b1centroid_5050.sep7=centroidpos_all.sep7.C6;
b1centroid_5050.sep9=centroidpos_all.sep9.C6;

%% hypoth 2/3 - additive envelope current dependent
singlecond=0;%picks the column unless 0 then does average for all
ploty=1;
clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        Envelope_add.(checksep).(checkcurr).T3=saveCplit.(checksep).(checkcurr).T1+saveCplit.(checksep).(checkcurr).T5;
        fun = @(b,X) b(1).*(X(:,1)+X(:,2));
        b0=[1];
        input_eq1=saveCplit.(checksep).(checkcurr).T1;
        %input_eq1(input_eq1==0)=nan;             %input_eq1(sum(input_eq1>0,2)<10,:)=0; if you want to take out edge conditions
        input_eq2=saveCplit.(checksep).(checkcurr).T5;
        %input_eq2(input_eq2==0)=nan;
        input_ans=saveCplit.(checksep).(checkcurr).T3;
        %input_ans(input_ans==0)=nan;
        mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
        mdlcur.(checksep) (current,:)=mdl.Coefficients{1,1};
        if ploty==1
            figure (sepdist)
            subplot(2,4,current)
            if singlecond==0
                plot(mean(Envelope_add.(checksep).(checkcurr).T3,2))
                hold on
                plot(mean(saveCplit.(checksep).(checkcurr).T3,2))
                plot(mean(saveCplit.(checksep).(checkcurr).T1,2),'Color','k')
                plot(mean(saveCplit.(checksep).(checkcurr).T5,2),'Color','k')
            else
                plot(Envelope_add.(checksep).(checkcurr).T3(:,singlecond))
                hold on
                plot(saveCplit.(checksep).(checkcurr).T3(:,singlecond))
                plot(saveCplit.(checksep).(checkcurr).T1(:,singlecond),'Color','k')
                plot(saveCplit.(checksep).(checkcurr).T5(:,singlecond),'Color','k')
            end
            title([num2str(AMP(current)) '\muA'])
            ylim([0, 360])
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')

            figure (sepdist+1)
            subplot(2,4,current)
            if singlecond==0
                plot(mean(Envelope_add.(checksep).(checkcurr).T3.*mdl.Coefficients{1,1},2))
                hold on
                plot(mean(saveCplit.(checksep).(checkcurr).T3,2))
            else
                plot(Envelope_add.(checksep).(checkcurr).T3(:,singlecond).*mdl.Coefficients{1,1})
                hold on
                plot(saveCplit.(checksep).(checkcurr).T3(:,singlecond))
            end
            title([num2str(AMP(current)) '\muA'])
            ylim([0, 360])
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
        end
        MSE.(checksep)(current)=mdl.MSE;
        rsquared.(checksep)(current)=mdl.Rsquared.Ordinary;
        AIC.(checksep)(current)=mdl.ModelCriterion.AICc;
    end
    if ploty==1
        legend('model_d_u_a_l','real_d_u_a_l')
    end
end

%% hypoth 4/5 - additive envelope independent additive variables current dependent
clear mdlcur MSE rsquared AIC_all centroidpos_all peak_all MSE_all
singlecond=0;%2 picks the column unless 0 then does average for all
ploty=0;
AMP=[0 1 2 3 4 6 8 10];
Csplit_dwnsample=[];
Envelope_add_dwnsmple=[];
seedpoint=65;%20,55
s = RandStream('mlfg6331_64','Seed',seedpoint);

for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 1:length(AMP)
         counter=0;
        checkcurr=['C' num2str(AMP(current))];
        randit = datasample(s,1:size(saveCplit.(checksep).(checkcurr).T1,2),size(saveCplit.(checksep).(checkcurr).T1,2),'Replace',false);
        for iterate10=floor(size(saveCplit.(checksep).(checkcurr).T1,2)/100*10):floor(size(saveCplit.(checksep).(checkcurr).T1,2)/100*10):size(saveCplit.(checksep).(checkcurr).T1,2) %cross-validation
            %Envelope_add.(checksep).(checkcurr).T3=saveCplit.(checksep).(checkcurr).T1+saveCplit.(checksep).(checkcurr).T5;
            fun = @(b,X) b(1).*X(:,1)+b(2).*X(:,2);
            b0=[1 1];
            %input_eq1=[nan(13,size(saveCplit.(checksep).(checkcurr).T1,2)); saveCplit.(checksep).(checkcurr).T1(14:16+sepdist+1+2,:); nan(32-19-sepdist,size(saveCplit.(checksep).(checkcurr).T1,2))];
            num_columns=floor(size(saveCplit.(checksep).(checkcurr).T1,2)/100*10);
            input_eq1=saveCplit.(checksep).(checkcurr).T1;
            input_eq1(:,randit(iterate10-num_columns+1:iterate10))=nan;
            %input_eq1(input_eq1==0)=nan;
            %input_eq2=[nan(13,size(saveCplit.(checksep).(checkcurr).T5,2)); saveCplit.(checksep).(checkcurr).T5(14:16+sepdist+1+2,:); nan(32-19-sepdist,size(saveCplit.(checksep).(checkcurr).T5,2))];

            input_eq2=saveCplit.(checksep).(checkcurr).T5;
            input_eq2(:,randit(iterate10-num_columns+1:iterate10))=nan;
            %input_eq2(input_eq2==0)=nan;
            input_ans=saveCplit.(checksep).(checkcurr).T3;
            input_ans(:,randit(iterate10-num_columns+1:iterate10))=nan;
            %input_ans=[nan(13,size(saveCplit.(checksep).(checkcurr).T3,2)); saveCplit.(checksep).(checkcurr).T3(14:16+sepdist+1+2,:); nan(32-19-sepdist,size(saveCplit.(checksep).(checkcurr).T3,2))];
            %input_ans(input_ans==0)=nan;
            %             electnumN=(1.96.*nanstd(input_ans,[],2)./(75)).^2;
            %             electnumperavg=sum(input_ans>=0,2);
            %             input_eq1(electnumperavg<electnumN,:)=nan;
            %             input_eq2(electnumperavg<electnumN,:)=nan;
            %             input_ans(electnumperavg<electnumN,:)=nan;
            counter=counter+1;
            mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
            mdlcur.(checksep).(checkcurr) (counter,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1}];
            Envelope_add.(checksep).(checkcurr).T3(:,randit(iterate10-num_columns+1:iterate10))=mdl.Coefficients{1,1}.*saveCplit.(checksep).(checkcurr).T1(:,randit(iterate10-num_columns+1:iterate10))+mdl.Coefficients{2,1}.*saveCplit.(checksep).(checkcurr).T5(:,randit(iterate10-num_columns+1:iterate10));

            %Envelope_add.(checksep).(checkcurr).T3=mdl.Coefficients{1,1}.*input_eq1+mdl.Coefficients{2,1}.*input_eq2;

            %for trial=1:5
            % check=['T' num2str(trial)];
        end
        Envelope_add.(checksep).(checkcurr).T3(:,randit(iterate10+1:end))=mdl.Coefficients{1,1}.*saveCplit.(checksep).(checkcurr).T1(:,randit(iterate10+1:end))+mdl.Coefficients{2,1}.*saveCplit.(checksep).(checkcurr).T5(:,randit(iterate10+1:end));

        y=saveCplit.(checksep).(checkcurr).T3;%realdata
        yinput=input_eq1;%inputdata1
        yfit=Envelope_add.(checksep).(checkcurr).T3;%modeldata
        numparam=2;
        [rsq,AIC,MSE]=modelfitcharacterisation(y,yfit,yinput,numparam);

            rsquared.(checksep)(current)=rsq;
            AIC_all.(checksep)(current)=AIC;
            MSE_all.(checksep)(current)=MSE;






        dat=saveCplit.(checksep).(checkcurr).T3;
        dat2=Envelope_add.(checksep).(checkcurr).T3;
        dat(isinf(dat))=nan;
        clear pairavgcur erpairavgcur pairavgcursmooth_env




        for paircount=1:size(dat,2)/4
            pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*4+1:4* paircount),2);
            pairavgcur_env(:,paircount)=nanmean(dat2(:,( paircount-1)*4+1:4* paircount),2);
            [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur_env(:,paircount),sepdist);
            centroidpos_all.(checksep).(checkcurr)(paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%


            % peak calc
            datnonan=pairavgcur_env(:,paircount);
            firstnonan=find(~isnan(datnonan),1,'first');
            lastnonan=find(~isnan(datnonan),1,'last');
            datnonan(isnan(datnonan))=[];
            FilterLength=5;
            b=ones(FilterLength,1)./FilterLength;
            smoothedenvelope=filtfilt(b,1,datnonan);
            pairavgcursmooth_env(:,paircount)=[nan(firstnonan-1,1); smoothedenvelope; nan(32-lastnonan,1)];

        end
        lengthneeded=size(saveCplit.sep7.C1.T1,2)/4;
        [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
        samples_rand.(checksep)=samples_rand1;
        centroidpos_all.(checksep).(checkcurr)=centroidpos_all.(checksep).(checkcurr)(:, samples_rand.(checksep));
        centroidpos_all.(checksep).(checkcurr)(sum(~isnan(centroidpos_all.(checksep).(checkcurr)),2)<10,:)=nan;

        pairavg_dwnsample=pairavgcur(:, samples_rand.(checksep));
        pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan;
        Csplit_dwnsample.(checksep).(checkcurr).T3=pairavg_dwnsample;

        pairavgcur_env_dwnsmple=pairavgcur_env(:, samples_rand.(checksep));
        pairavgcur_env_dwnsmple(sum(~isnan(pairavgcur_env_dwnsmple),2)<10,:)=nan;
        Envelope_add_dwnsmple.(checksep).(checkcurr).T3=pairavgcur_env_dwnsmple;

        %peak calc
        pairavgcursmooth_env=pairavgcursmooth_env(:, samples_rand.(checksep));
        pairavgcursmooth_env(:,sum(pairavgcursmooth_env>0,1)==0)=nan;
        [peak,c]=find(pairavgcursmooth_env==max(pairavgcursmooth_env));
        peak_all.(checksep).(checkcurr)=(peak'-16).*50;

        %                         input_eq1(sum(input_eq1>=0,2)<55,:)=nan;
        %             input_eq2(sum(input_eq2>=0,2)<55,:)=nan;
        %             input_ans(sum(input_ans>=0,2)<55,:)=nan;
        if ploty==1
            %                 figure (sepdist)
            %                 if length(AMP)>1
            %                     subplot(2,4,current)
            %                 end
            %                 plot(mean(Envelope_add_dwnsmple.(checksep).(checkcurr).T3,2))
            %                 hold on
            %                 plot(mean(Csplit_dwnsample.(checksep).(checkcurr).T3,2))
            %                 plot(mean(Csplit_dwnsample.(checksep).(checkcurr).T1,2),'k')
            %                 plot(mean(Csplit_dwnsample.(checksep).(checkcurr).T5,2),'k')
            %                 title([num2str(AMP(current)) '\muA'])
            %                 ylim([0, 360])
            %                 xline(16,'r')
            %                 xline(16+sepdist,'r')
            %                 xt = xticks;
            %                 xtl=(xt-16)*50;
            %                 xticklabels(xtl)
            %                 ylabel('Sp/s')
            %                 xlabel('Distance from deepest stim elect (\mum)')
            figure (sepdist+1)
            if length(AMP)>1
                subplot(2,4,current)
            end
            if singlecond==0
                %Envelope_add.(checksep).(checkcurr).T3(sum(Envelope_add.(checksep).(checkcurr).T3>0,2)<15,:)=0; %if you want to take out edge conditions
                %Envelope_add.(checksep).(checkcurr).T3(Envelope_add.(checksep).(checkcurr).T3==0)=nan;
                stdshade(Envelope_add_dwnsmple.(checksep).(checkcurr).T3',0.2,'r');
                hold on
                %saveCplit.(checksep).(checkcurr).T3(sum(saveCplit.(checksep).(checkcurr).T3>0,2)<15,:)=0; %if you want to take out edge conditions
                %saveCplit.(checksep).(checkcurr).T3(saveCplit.(checksep).(checkcurr).T3==0)=nan;
                stdshade(Csplit_dwnsample.(checksep).(checkcurr).T3',0.2,'k');
                ylim([0, 360])

                xlim([find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).T3>=0,2),1,'first')-1, find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).T3>=0,2),1,'last')+1])
            else

                plot(Envelope_add_dwnsmple.(checksep).(checkcurr).T3(:,singlecond),'b')
                hold on
                plot(Csplit_dwnsample.(checksep).(checkcurr).T3(:,singlecond),'r')
                ylim([0, 800])
                xlim([find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).T3(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).T3(:,singlecond)>=0,2),1,'last')+1])
            end
            title([num2str(AMP(current)) '\muA'])
            set(gca,'TickDir','out');
            xline(16,'r')
            xline(16+sepdist+1,'r')
            if sepdist==5
                xlim([16-3 16+sepdist+1+3])%([find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'first'), find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'last')])
            elseif sepdist==7
                xlim([16-2 16+sepdist+1+2])
            elseif sepdist==9
                xlim([16-1 16+sepdist+1+1])
            end
            ylim([0 350])
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            %xlim([13 27])
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
        end
        %             MSE.(checksep)(current)=mdl.MSE;
        %             rsquared.(checksep)(current)=mdl.Rsquared.Ordinary;
        %             AIC.(checksep)(current)=mdl.ModelCriterion.AIC;

        %             figure(14)
        %             subplot(2,4,current)
        %             hold on
        %             plot(electnumN)
        %             plot(electnumperavg)
        %             title([num2str(AMP(current)) '\muA'])
        %             ylabel('Electrode number')
        %             xlabel('Relative position')

    end
    %legend('Number electrodes needed', 'Actual number electrodes')
    if ploty==1
        figure (sepdist+1)
        f=get(gca,'Children');
        legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
    end
    if length(AMP)>1
        figure
        hold on
        meancoef1=zeros(length(AMP),2);
        for itcurr=1:length(AMP)
            checkcurr=['C' num2str(AMP(itcurr))];
            meancoef1(itcurr,1:2)=mean(mdlcur.(checksep).(checkcurr),1);

        end
        plot(AMP,   meancoef1(:,1), 'r')
        plot(AMP,   meancoef1(:,2), 'k')
        ylabel('Scaling Factor')
        xlabel('Current \muA')
        legend({'Factor 1','Factor 2'})
        set(gca,'TickDir','out');
        title(['Separation ' num2str(sepdist*50+50) '\mum'])
        xlim([1 AMP(end)])
        ylim([0 1])
    end
end

b2centroid_5050.sep5=centroidpos_all.sep5.C6;
b2centroid_5050.sep7=centroidpos_all.sep7.C6;
b2centroid_5050.sep9=centroidpos_all.sep9.C6;
%% hypoth 6 - single scaling factor
clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
altogethernumel=0;
condition='T3';
numelarray=zeros(length(AMP),1);
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        altogethernumel=altogethernumel+numel(saveCplit.(checksep).(checkcurr).T1(:));
        numelarray(current)=numel(saveCplit.(checksep).(checkcurr).T1(:));
    end
    T1=zeros(altogethernumel,1);
    T3=zeros(altogethernumel,1);
    T5=zeros(altogethernumel,1);
    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        Envelope_add.(checksep).(checkcurr).(condition)=saveCplit.(checksep).(checkcurr).T1+saveCplit.(checksep).(checkcurr).T5;

        T1(1+sum(numelarray(1:current-1)):numelarray(current)+sum(numelarray(1:current-1)))=saveCplit.(checksep).(checkcurr).T1(:);
        T5(1+sum(numelarray(1:current-1)):numelarray(current)+sum(numelarray(1:current-1)))=saveCplit.(checksep).(checkcurr).T5(:);
        T3(1+sum(numelarray(1:current-1)):numelarray(current)+sum(numelarray(1:current-1)))=saveCplit.(checksep).(checkcurr).(condition)(:);
    end
    fun = @(b,X) b(1).*(X(:,1)+X(:,2));
    b0=[1];
    input_eq1=T1;
    input_eq1(input_eq1==0)=nan;
    input_eq2=T5;
    input_eq2(input_eq2==0)=nan;
    input_ans=T3;
    input_ans(input_ans==0)=nan;
    mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
    mdlcur.(checksep)=mdl.Coefficients{1,1};


    MSE.(checksep)=mdl.MSE;
    rsquared.(checksep)=mdl.Rsquared.Ordinary;
    AIC.(checksep)=mdl.ModelCriterion.AIC;

    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        Envelope_add.(checksep).(checkcurr).T3=mdl.Coefficients{1,1}.*(saveCplit.(checksep).(checkcurr).T1+saveCplit.(checksep).(checkcurr).T5);
        figure (sepdist+1)
        subplot(2,4,current)
        plot(mean(Envelope_add.(checksep).(checkcurr).(condition),2))
        hold on
        plot(mean(saveCplit.(checksep).(checkcurr).(condition),2))
        title([num2str(AMP(current)) '\muA'])
        ylim([0, 360])
        xline(16,'r')
        xline(16+sepdist,'r')
        xt = xticks;
        xtl=(xt-16)*50;
        xticklabels(xtl)
        ylabel('Sp/s')
        xlabel('Distance from deepest stim elect (\mum)')
    end
    legend('model_d_u_a_l','real_d_u_a_l')
end


%% hypoth 7 - complicated additional parameters
clear mdlcur MSE rsquared AIC
singlecond=0;%2 picks the column unless 0 then does average for all
ploty=1;
multiply1=0.5;
multiply2=0.5;
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        input_eq1=saveCplit.(checksep).(checkcurr).T1;
        input_eq2=saveCplit.(checksep).(checkcurr).T5;
        input_ans=saveCplit.(checksep).(checkcurr).T3;
        electnumN=(1.96.*nanstd(input_ans,[],2)./(75)).^2;
        electnumperavg=sum(input_ans>=0,2);

        fun = @(b,X) (b(1).*X(:,1))+(b(2).*X(:,2))+(multiply1.*X(:,1))+(multiply2.*X(:,2));
        b0=[1 1];
        mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
        mdlcur.(checksep) (current,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1} mdl.Coefficients{3,1} mdl.Coefficients{4,1}];
        Envelope_add.(checksep).(checkcurr).T3=mdl.Coefficients{1,1}.*input_eq1+mdl.Coefficients{2,1}.*input_eq2+multiply1.*mdl.Coefficients{3,1}.*input_eq1+multiply2.*mdl.Coefficients{4,1}.*input_eq2;

        if ploty==1
            figure (sepdist)
            subplot(2,4,current)
            plot(mean(Envelope_add.(checksep).(checkcurr).T3,2))
            hold on
            plot(mean(saveCplit.(checksep).(checkcurr).T3,2))
            plot(mean(saveCplit.(checksep).(checkcurr).T1,2),'k')
            plot(mean(saveCplit.(checksep).(checkcurr).T5,2),'k')
            title([num2str(AMP(current)) '\muA'])
            ylim([0, 360])
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
            figure (sepdist+1+20)
            subplot(2,4,current)
            if singlecond==0
                %Envelope_add.(checksep).(checkcurr).T3(sum(Envelope_add.(checksep).(checkcurr).T3>0,2)<15,:)=0; %if you want to take out edge conditions
                %Envelope_add.(checksep).(checkcurr).T3(Envelope_add.(checksep).(checkcurr).T3==0)=nan;
                stdshade(Envelope_add.(checksep).(checkcurr).T3',0.2,'b');
                hold on
                %saveCplit.(checksep).(checkcurr).T3(sum(saveCplit.(checksep).(checkcurr).T3>0,2)<15,:)=0; %if you want to take out edge conditions
                %saveCplit.(checksep).(checkcurr).T3(saveCplit.(checksep).(checkcurr).T3==0)=nan;
                stdshade(input_ans',0.2,'r');
                ylim([0, 360])

                xlim([find(sum(Envelope_add.(checksep).(checkcurr).T3>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).T3>=0,2),1,'last')+1])
            else

                plot(Envelope_add.(checksep).(checkcurr).T3(:,singlecond),'b')
                hold on
                plot(saveCplit.(checksep).(checkcurr).T3(:,singlecond),'r')
                ylim([0, 800])
                xlim([find(sum(Envelope_add.(checksep).(checkcurr).T3(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).T3(:,singlecond)>=0,2),1,'last')+1])
            end
            title([num2str(AMP(current)) '\muA'])

            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
        end
        MSE.(checksep)(current)=mdl.MSE;
        rsquared.(checksep)(current)=mdl.Rsquared.Ordinary;
        AIC.(checksep)(current)=mdl.ModelCriterion.AIC;

        %             figure(14)
        %             subplot(2,4,current)
        %             hold on
        %             plot(electnumN)
        %             plot(electnumperavg)
        %             title([num2str(AMP(current)) '\muA'])
        %             ylabel('Electrode number')
        %             xlabel('Relative position')
    end
    legend('Number electrodes needed', 'Actual number electrodes')
    if ploty==1
        figure (sepdist+1+20)
        f=get(gca,'Children');
        legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
    end

end

%% hypoth 8 - additional parameters
clear mdlcur MSE rsquared AIC
singlecond=0;%2 picks the column unless 0 then does average for all
ploty=0;
multiply1=0.5;
multiply2=0.5;
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        input_eq1=saveCplit.(checksep).(checkcurr).T1;
        input_eq2=saveCplit.(checksep).(checkcurr).T5;
        input_ans=saveCplit.(checksep).(checkcurr).T3;
        electnumN=(1.96.*nanstd(input_ans,[],2)./(75)).^2;
        electnumperavg=sum(input_ans>=0,2);

        fun = @(b,X) (multiply1.*b(1).*X(:,1))+(multiply2.*b(2).*X(:,2));
        b0=[1 1];
        mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
        mdlcur.(checksep) (current,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1}];
        Envelope_add.(checksep).(checkcurr).T3=multiply1.*mdl.Coefficients{1,1}.*input_eq1+multiply2.*mdl.Coefficients{2,1}.*input_eq2;

        for trial=1:5
            check=['T' num2str(trial)];
            dat=saveCplit.(checksep).(checkcurr).T3;
            dat2=Envelope_add.(checksep).(checkcurr).T3;
            dat(isinf(dat))=nan;
            clear pairavgcur erpairavgcur
            for paircount=1:size(dat,2)/4
                pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*4+1:4* paircount),2);
                pairavgcur_env(:,paircount)=nanmean(dat2(:,( paircount-1)*4+1:4* paircount),2);
                %erpairavgcur(:,paircount)=nanstd(dat(:,( paircount-1)*3+1:3* paircount),[],2)./sqrt(sum(~isnan(dat(:,( paircount-1)*3+1:3* paircount)), 2));
            end
            lengthneeded=size(saveCplit.sep7.C1.T1,2)/4;
            if trial==1
                [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
                samples_rand.(checksep)=samples_rand1;
            end
            pairavg_dwnsample=pairavgcur(:, samples_rand.(checksep));
            pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan;
            Csplit_dwnsample.(checksep).(checkcurr).(check)=pairavg_dwnsample;

            pairavgcur_env_dwnsmple=pairavgcur_env(:, samples_rand.(checksep));
            pairavgcur_env_dwnsmple(sum(~isnan(pairavgcur_env_dwnsmple),2)<10,:)=nan;
            Envelope_add_dwnsmple.(checksep).(checkcurr).(check)=pairavgcur_env_dwnsmple;
        end

        if ploty==1
            figure (sepdist+1)
            if length(AMP)>1
                subplot(2,4,current)
            end
            if singlecond==0
                stdshade(Envelope_add_dwnsmple.(checksep).(checkcurr).T3',0.2,'b');
                hold on
                stdshade(Csplit_dwnsample.(checksep).(checkcurr).T3',0.2,'r');
                ylim([0, 360])

                xlim([find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).T3>=0,2),1,'first')-1, find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).T3>=0,2),1,'last')+1])
            else

                plot(Envelope_add_dwnsmple.(checksep).(checkcurr).T3(:,singlecond),'b')
                hold on
                plot(Csplit_dwnsample.(checksep).(checkcurr).T3(:,singlecond),'r')
                ylim([0, 800])
                xlim([find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).T3(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).T3(:,singlecond)>=0,2),1,'last')+1])
            end
            title([num2str(AMP(current)) '\muA'])
            set(gca,'TickDir','out');
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xlim([15 27])
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
        end
        MSE.(checksep)(current)=mdl.MSE;
        rsquared.(checksep)(current)=mdl.Rsquared.Ordinary;
        AIC.(checksep)(current)=mdl.ModelCriterion.AIC;

    end
    if ploty==1
        figure (sepdist+1)
        f=get(gca,'Children');
        legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
    end
    if length(AMP)>1
        figure
        hold on
        plot(AMP,  mdlcur.(checksep)(:,1), 'r')
        plot(AMP,  mdlcur.(checksep)(:,2), 'k')
        ylabel('Scaling Factor')
        xlabel('Current \muA')
        legend({'Factor 1','Factor 2'})
        set(gca,'TickDir','out');
        title(['Separation ' num2str(sepdist*50+50) '\mum'])
        xlim([1 AMP(end)])
        ylim([0 2])
    end

end


%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%75/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% hypoth 1 - pure additive envelope with full currents 75/25
%need to make sure you have the right loopsubdirect ticked!!! Tick
%threequarter_amplitude_single (option 1 below) or
%quarter_amplitude_single (option 2 below)
singlecond=0;
which75_25condition=1;
if which75_25condition==1
    Tcond='T2';
elseif which75_25condition==2
    Tcond='T4';
else
    error('Not valid')
end
AMP=6;%[0 1 2 3 4 6 8 10];
Csplit_dwnsample=[];
Envelope_add_dwnsmple=[];
seedpoint=65;%20,55
s = RandStream('mlfg6331_64','Seed',seedpoint);
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    figure(sepdist)
    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        Envelope_add.(checksep).(checkcurr).(Tcond)=saveCplit.(checksep).(checkcurr).T1+saveCplit.(checksep).(checkcurr).T5;%can just add because the correct amplitudes should have been selected before!!!!!!!!!!!!

        dat=saveCplit.(checksep).(checkcurr).(Tcond);
        dat2=Envelope_add.(checksep).(checkcurr).(Tcond);
        dat(isinf(dat))=nan;
        clear pairavgcur erpairavgcur
        for paircount=1:size(dat,2)/4
            pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*4+1:4* paircount),2);
            pairavgcur_env(:,paircount)=nanmean(dat2(:,( paircount-1)*4+1:4* paircount),2);
            %erpairavgcur(:,paircount)=nanstd(dat(:,( paircount-1)*3+1:3* paircount),[],2)./sqrt(sum(~isnan(dat(:,( paircount-1)*3+1:3* paircount)), 2));
            [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur_env(:,paircount),sepdist);
            centroidpos_all.(checksep).(checkcurr)(paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%

            %centroidpos_all.(checksep).(checkcurr)(paircount)=(centroidpos-16).*50;
        end
        lengthneeded=size(saveCplit.sep7.C1.T1,2)/4;
        [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
        samples_rand.(checksep)=samples_rand1;
        centroidpos_all.(checksep).(checkcurr)=centroidpos_all.(checksep).(checkcurr)(:, samples_rand.(checksep));
        centroidpos_all.(checksep).(checkcurr)(sum(~isnan(centroidpos_all.(checksep).(checkcurr)),2)<10,:)=nan;
        pairavg_dwnsample=pairavgcur(:, samples_rand.(checksep));
        pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan;
        Csplit_dwnsample.(checksep).(checkcurr).(Tcond)=pairavg_dwnsample;

        pairavgcur_env_dwnsmple=pairavgcur_env(:, samples_rand.(checksep));
        pairavgcur_env_dwnsmple(sum(~isnan(pairavgcur_env_dwnsmple),2)<10,:)=nan;
        Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)=pairavgcur_env_dwnsmple;

        if length(AMP)>1
            subplot(2,4,current)
        end
        if singlecond==0
            %Envelope_add.(checksep).(checkcurr).(Tcond)(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>0,2)<15,:)=0; %if you want to take out edge conditions
            % Envelope_add.(checksep).(checkcurr).(Tcond)(Envelope_add.(checksep).(checkcurr).(Tcond)==0)=nan;
            stdshade(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)',0.2,'b');
            hold on
            %saveCplit.(checksep).(checkcurr).(Tcond)(sum(saveCplit.(checksep).(checkcurr).(Tcond)>0,2)<15,:)=0; %if you want to take out edge conditions
            %saveCplit.(checksep).(checkcurr).(Tcond)(saveCplit.(checksep).(checkcurr).(Tcond)==0)=nan;
            stdshade(Csplit_dwnsample.(checksep).(checkcurr).(Tcond)',0.2,'k');
            %plot(mean(saveCplit.(checksep).(checkcurr).T1.*multiply1,2),'Color','k')
            %plot(mean(saveCplit.(checksep).(checkcurr).T5.*multiply1,2),'Color','k')
        else
            plot(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
            hold on
            plot(Csplit_dwnsample.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
            plot(Csplit_dwnsample.(checksep).(checkcurr).T1(:,singlecond),'Color','k')
            plot(Csplit_dwnsample.(checksep).(checkcurr).T5(:,singlecond),'Color','k')
        end
        title([num2str(AMP(current)) '\muA'])
        ylim([0, 375])
        xlim([13 27])
        xline(16,'r')
        xline(16+sepdist+1,'r')
        if sepdist==5
            xlim([16-3 16+sepdist+1+3])%([find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'first'), find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'last')])
        elseif sepdist==7
            xlim([16-2 16+sepdist+1+2])
        elseif sepdist==9
            xlim([16-1 16+sepdist+1+1])
        end
        xt = xticks;
        xtl=(xt-16)*50;
        xticklabels(xtl)
        ylabel('Sp/s')
        xlabel('Distance from deepest stim elect (\mum)')
        set(gca,'TickDir','out');
        ylim([0 350])
        y=saveCplit.(checksep).(checkcurr).(Tcond);%realdata
        yinput=saveCplit.(checksep).(checkcurr).T1;%inputdata1
        yfit=Envelope_add.(checksep).(checkcurr).(Tcond);%modeldata
        numparam=0;
        [rsq,AIC,MSE]=modelfitcharacterisation(y,yfit,yinput,numparam);
        rsquared.(checksep)(current)=rsq;
        AIC_all.(checksep)(current)=AIC;
        MSE_all.(checksep)(current)=MSE;
    end
    f=get(gca,'Children');
    legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
end


%% hypoth 2 - additive envelope with full currents 75/25
%need to make sure you have the right loopsubdirect ticked!!! Do not
%tick any! all 0
clear peak_all
which75_25condition=1;
if which75_25condition==1
    multiply1=3/4;
    multiply2=1/4;
    Tcond='T2';
elseif which75_25condition==2
    multiply1=1/4;
    multiply2=3/4;
    Tcond='T4';
else
    error('Not valid')
end

AMP=6;%[0 1 2 3 4 6 8 10];
Csplit_dwnsample=[];
Envelope_add_dwnsmple=[];
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);

for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    figure(sepdist)
    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        Envelope_add.(checksep).(checkcurr).(Tcond)=saveCplit.(checksep).(checkcurr).T1.*multiply1+saveCplit.(checksep).(checkcurr).T5.*multiply2;
        dat=saveCplit.(checksep).(checkcurr).(Tcond);
        dat2=Envelope_add.(checksep).(checkcurr).(Tcond);
        dat(isinf(dat))=nan;
        clear pairavgcur erpairavgcur pairavgcursmooth_env pairavgcur_env
        for paircount=1:size(dat,2)/4
            pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*4+1:4* paircount),2);
            pairavgcur_env(:,paircount)=nanmean(dat2(:,( paircount-1)*4+1:4* paircount),2);
            %erpairavgcur(:,paircount)=nanstd(dat(:,( paircount-1)*3+1:3* paircount),[],2)./sqrt(sum(~isnan(dat(:,( paircount-1)*3+1:3* paircount)), 2));
            [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur_env(:,paircount),sepdist);
            centroidpos_all.(checksep).(checkcurr)(paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%

            %centroidpos_all.(checksep).(checkcurr)(paircount)=(centroidpos-16).*50;
            % peak calc
            datnonan=pairavgcur_env(:,paircount);
            firstnonan=find(~isnan(datnonan),1,'first');
            lastnonan=find(~isnan(datnonan),1,'last');
            datnonan(isnan(datnonan))=[];
            FilterLength=5;
            b=ones(FilterLength,1)./FilterLength;
            smoothedenvelope=filtfilt(b,1,datnonan);
            pairavgcursmooth_env(:,paircount)=[nan(firstnonan-1,1); smoothedenvelope; nan(32-lastnonan,1)];

        end
        lengthneeded=size(saveCplit.sep7.C1.T1,2)/4;
        [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
        samples_rand.(checksep)=samples_rand1;
        centroidpos_all.(checksep).(checkcurr)=centroidpos_all.(checksep).(checkcurr)(:, samples_rand.(checksep));
        centroidpos_all.(checksep).(checkcurr)(sum(~isnan(centroidpos_all.(checksep).(checkcurr)),2)<10,:)=nan;
        pairavg_dwnsample=pairavgcur(:, samples_rand.(checksep));
        pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan;
        Csplit_dwnsample.(checksep).(checkcurr).(Tcond)=pairavg_dwnsample;

        pairavgcur_env_dwnsmple=pairavgcur_env(:, samples_rand.(checksep));
        pairavgcur_env_dwnsmple(sum(~isnan(pairavgcur_env_dwnsmple),2)<10,:)=nan;
        Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)=pairavgcur_env_dwnsmple;



        %peak calc
        pairavgcursmooth_env=pairavgcursmooth_env(:, samples_rand.(checksep));
        pairavgcursmooth_env(:,sum(pairavgcursmooth_env>0,1)==0)=nan;
        [peak,c]=find(pairavgcursmooth_env==max(pairavgcursmooth_env));
        peak_all.(checksep).(checkcurr)=(peak'-16).*50;

        if length(AMP)>1
            subplot(2,4,current)
        end
        if singlecond==0
            %Envelope_add.(checksep).(checkcurr).(Tcond)(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>0,2)<15,:)=0; %if you want to take out edge conditions
            % Envelope_add.(checksep).(checkcurr).(Tcond)(Envelope_add.(checksep).(checkcurr).(Tcond)==0)=nan;
            stdshade(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)',0.2,'r');
            hold on
            %saveCplit.(checksep).(checkcurr).(Tcond)(sum(saveCplit.(checksep).(checkcurr).(Tcond)>0,2)<15,:)=0; %if you want to take out edge conditions
            %saveCplit.(checksep).(checkcurr).(Tcond)(saveCplit.(checksep).(checkcurr).(Tcond)==0)=nan;

            stdshade(Csplit_dwnsample.(checksep).(checkcurr).(Tcond)',0.2,'k');

            %plot(mean(saveCplit.(checksep).(checkcurr).T1.*multiply1,2),'Color','k')
            %plot(mean(saveCplit.(checksep).(checkcurr).T5.*multiply1,2),'Color','k')
        else
            plot(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
            hold on
            plot(Csplit_dwnsample.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
            plot(Csplit_dwnsample.(checksep).(checkcurr).T1(:,singlecond),'Color','k')
            plot(Csplit_dwnsample.(checksep).(checkcurr).T5(:,singlecond),'Color','k')
        end
        title([num2str(AMP(current)) '\muA'])
        ylim([0, 375])
        xlim([13 27])
        xline(16,'r')
        xline(16+sepdist+1,'r')
        if sepdist==5
            xlim([16-3 16+sepdist+1+3])%([find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'first'), find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'last')])
        elseif sepdist==7
            xlim([16-2 16+sepdist+1+2])
        elseif sepdist==9
            xlim([16-1 16+sepdist+1+1])
        end
        ylim([0 350])
        xt = xticks;
        xtl=(xt-16)*50;
        xticklabels(xtl)
        ylabel('Sp/s')
        xlabel('Distance from deepest stim elect (\mum)')
        set(gca,'TickDir','out');
        y=saveCplit.(checksep).(checkcurr).(Tcond);%realdata
        yinput=saveCplit.(checksep).(checkcurr).T1;%inputdata1
        yfit=Envelope_add.(checksep).(checkcurr).(Tcond);%modeldata
        numparam=0;
        [rsq,AIC,MSE]=modelfitcharacterisation(y,yfit,yinput,numparam);
        rsquared.(checksep)(current)=rsq;
        AIC_all.(checksep)(current)=AIC;
        MSE_all.(checksep)(current)=MSE;

    end
    f=get(gca,'Children');
    legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
end
if which75_25condition==1
    b1centroid_7525.sep5=centroidpos_all.sep5.C6;
    b1centroid_7525.sep7=centroidpos_all.sep7.C6;
    b1centroid_7525.sep9=centroidpos_all.sep9.C6;
else
    b1centroid_2575.sep5=centroidpos_all.sep5.C6;
    b1centroid_2575.sep7=centroidpos_all.sep7.C6;
    b1centroid_2575.sep9=centroidpos_all.sep9.C6;
end
%% hypoth 3 - true additive envelope current dependent scaled
%need to make sure you have the right loopsubdirect ticked!!! Tick
%threequarter_amplitude_single (option 1 below) or
%quarter_amplitude_single (option 2 below)
singlecond=0;
ploty=0;
which75_25condition=1;
if which75_25condition==1
    Tcond='T2';
elseif which75_25condition==2
    Tcond='T4';
else
    error('Not valid')
end
clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        Envelope_add.(checksep).(checkcurr).(Tcond)=saveCplit.(checksep).(checkcurr).T1+saveCplit.(checksep).(checkcurr).T5;
        fun = @(b,X) b(1).*(X(:,1)+X(:,2));
        b0=[1];
        input_eq1=saveCplit.(checksep).(checkcurr).T1;
        %input_eq1(input_eq1==0)=nan;             %input_eq1(sum(input_eq1>0,2)<10,:)=0; if you want to take out edge conditions
        input_eq2=saveCplit.(checksep).(checkcurr).T5;
        %input_eq2(input_eq2==0)=nan;
        input_ans=saveCplit.(checksep).(checkcurr).(Tcond);
        %input_ans(input_ans==0)=nan;
        mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
        mdlcur.(checksep) (current,:)=mdl.Coefficients{1,1};
        if ploty==1
            figure (sepdist+10)
            subplot(2,4,current)
            plot(mean(Envelope_add.(checksep).(checkcurr).(Tcond),2))
            hold on
            plot(mean(saveCplit.(checksep).(checkcurr).(Tcond),2))
            plot(mean(saveCplit.(checksep).(checkcurr).T1,2),'k')
            plot(mean(saveCplit.(checksep).(checkcurr).T5,2),'k')
            title([num2str(AMP(current)) '\muA'])
            ylim([0, 360])
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')

            figure (sepdist+1+10)
            subplot(2,4,current)
            if singlecond==0
                plot(mean(Envelope_add.(checksep).(checkcurr).(Tcond).*mdl.Coefficients{1,1},2))
                hold on
                plot(mean(saveCplit.(checksep).(checkcurr).(Tcond),2))
            else
                plot(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond).*mdl.Coefficients{1,1})
                hold on
                plot(saveCplit.(checksep).(checkcurr).(Tcond)(:,singlecond))
            end
            title([num2str(AMP(current)) '\muA'])
            ylim([0, 650])
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
        end
        MSE.(checksep)(current)=mdl.MSE;
        rsquared.(checksep)(current)=mdl.Rsquared.Ordinary;
        AIC.(checksep)(current)=mdl.ModelCriterion.AIC;
    end
    if ploty==1
        legend('model_d_u_a_l','real_d_u_a_l')
    end
end


%% hypoth 4/5 - additive envelope current dependent scaled equal
%need to make sure you have the right loopsubdirect ticked!!! Do not
%tick any! all 0
ploty=0;
which75_25condition=1;
if which75_25condition==1
    multiply1=3/4;
    multiply2=1/4;
    Tcond='T2';
elseif which75_25condition==2
    multiply1=1/4;
    multiply2=3/4;
    Tcond='T4';
else
    error('Not valid')
end
clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        Envelope_add.(checksep).(checkcurr).(Tcond)=saveCplit.(checksep).(checkcurr).T1.*multiply1+saveCplit.(checksep).(checkcurr).T5.*multiply2;
        fun = @(b,X) b(1).*(multiply1.*X(:,1)+multiply2.*X(:,2));
        b0=[1];
        input_eq1=saveCplit.(checksep).(checkcurr).T1;
        %input_eq1(input_eq1==0)=nan;             %input_eq1(sum(input_eq1>0,2)<10,:)=0; if you want to take out edge conditions
        input_eq2=saveCplit.(checksep).(checkcurr).T5;
        %input_eq2(input_eq2==0)=nan;
        input_ans=saveCplit.(checksep).(checkcurr).(Tcond);
        %input_ans(input_ans==0)=nan;
        mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
        mdlcur.(checksep) (current,:)=mdl.Coefficients{1,1};
        if ploty==1
            figure (sepdist)
            subplot(2,4,current)
            plot(mean(Envelope_add.(checksep).(checkcurr).(Tcond),2))
            hold on
            plot(mean(saveCplit.(checksep).(checkcurr).(Tcond),2))
            plot(mean(saveCplit.(checksep).(checkcurr).T1*multiply1,2),'k')
            plot(mean(saveCplit.(checksep).(checkcurr).T5*multiply2,2),'k')
            title([num2str(AMP(current)) '\muA'])
            ylim([0, 360])
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')

            figure (sepdist+1)
            subplot(2,4,current)
            plot(mean(mdl.Coefficients{1,1}.*Envelope_add.(checksep).(checkcurr).(Tcond),2))
            hold on
            plot(mean(saveCplit.(checksep).(checkcurr).(Tcond),2))
            title([num2str(AMP(current)) '\muA'])
            ylim([0, 360])
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
        end
        MSE.(checksep)(current)=mdl.MSE;
        rsquared.(checksep)(current)=mdl.Rsquared.Ordinary;
        AIC.(checksep)(current)=mdl.ModelCriterion.AIC;
    end
    if ploty==1
        legend('model_d_u_a_l','real_d_u_a_l')
    end
end

%% hypoth 6/7/10 - individual scaling factors additive envelope current dependent
%need to make sure you have the right loopsubdirect ticked!!! need full
%currents
singlecond=0;%3
ploty=0;
which75_25condition=2;
if which75_25condition==1
    Tcond='T2';%7525
elseif which75_25condition==2
    Tcond='T4';%2575
else
    error('Not valid')
end
clear mdlcur MSE rsquared AIC_all peak_all MSE_all
AMP=[1 2 3 4 6 8 10];
Csplit_dwnsample=[];
Envelope_add_dwnsmple=[];
seedpoint=65;%20,55
s = RandStream('mlfg6331_64','Seed',seedpoint);
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        counter=0;
        randit = datasample(s,1:size(saveCplit.(checksep).(checkcurr).T1,2),size(saveCplit.(checksep).(checkcurr).T1,2),'Replace',false);
        for iterate10=floor(size(saveCplit.(checksep).(checkcurr).T1,2)/100*10):floor(size(saveCplit.(checksep).(checkcurr).T1,2)/100*10):size(saveCplit.(checksep).(checkcurr).T1,2) %cross-validation
            counter=counter+1;
            fun = @(b,X) (b(1).*X(:,1))+(b(2).*X(:,2));
            b0=[1 1];
            num_columns=floor(size(saveCplit.(checksep).(checkcurr).T1,2)/100*10);
            input_eq1=saveCplit.(checksep).(checkcurr).T1;
            input_eq1(:,randit(iterate10-num_columns+1:iterate10))=nan;
            input_eq2=saveCplit.(checksep).(checkcurr).T5;
            input_eq2(:,randit(iterate10-num_columns+1:iterate10))=nan;
            input_ans=saveCplit.(checksep).(checkcurr).(Tcond);
            input_ans(:,randit(iterate10-num_columns+1:iterate10))=nan;
            mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
            mdlcur.(checksep).(checkcurr) (counter,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1}];
            Envelope_add.(checksep).(checkcurr).(Tcond)(:,randit(iterate10-num_columns+1:iterate10))=mdl.Coefficients{1,1}.*saveCplit.(checksep).(checkcurr).T1(:,randit(iterate10-num_columns+1:iterate10))+mdl.Coefficients{2,1}.*saveCplit.(checksep).(checkcurr).T5(:,randit(iterate10-num_columns+1:iterate10));


        end
         Envelope_add.(checksep).(checkcurr).(Tcond)(:,randit(iterate10+1:end))=mdl.Coefficients{1,1}.*saveCplit.(checksep).(checkcurr).T1(:,randit(iterate10+1:end))+mdl.Coefficients{2,1}.*saveCplit.(checksep).(checkcurr).T5(:,randit(iterate10+1:end));
            y=saveCplit.(checksep).(checkcurr).(Tcond);%realdata
            yinput=saveCplit.(checksep).(checkcurr).T1;%inputdata1
            yfit=Envelope_add.(checksep).(checkcurr).(Tcond);%modeldata
            
            numparam=2;
            [rsq,AIC,MSE]=modelfitcharacterisation(y,yfit,yinput,numparam);
            rsquared.(checksep)(current)=rsq;
            AIC_all.(checksep)(current)=AIC;
            MSE_all.(checksep)(current)=MSE;
        dat=saveCplit.(checksep).(checkcurr).(Tcond);
        dat2=Envelope_add.(checksep).(checkcurr).(Tcond);
        dat(isinf(dat))=nan;
        clear pairavgcur erpairavgcur pairavgcur_env pairavgcursmooth_env
        for paircount=1:size(dat,2)/4
            pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*4+1:4* paircount),2);
            pairavgcur_env(:,paircount)=nanmean(dat2(:,( paircount-1)*4+1:4* paircount),2);
            %erpairavgcur(:,paircount)=nanstd(dat(:,( paircount-1)*3+1:3* paircount),[],2)./sqrt(sum(~isnan(dat(:,( paircount-1)*3+1:3* paircount)), 2));
            [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur_env(:,paircount),sepdist);
            centroidpos_all.(checksep).(checkcurr)(paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%
            %centroidpos_all.(checksep).(checkcurr)(paircount)=(centroidpos-16).*50;


            % peak calc
            datnonan=pairavgcur_env(:,paircount);
            firstnonan=find(~isnan(datnonan),1,'first');
            lastnonan=find(~isnan(datnonan),1,'last');
            datnonan(isnan(datnonan))=[];
            FilterLength=5;
            b=ones(FilterLength,1)./FilterLength;
            smoothedenvelope=filtfilt(b,1,datnonan);
            pairavgcursmooth_env(:,paircount)=[nan(firstnonan-1,1); smoothedenvelope; nan(32-lastnonan,1)];

        end
        lengthneeded=size(saveCplit.sep7.C1.T1,2)/4;
        [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
        samples_rand.(checksep)=samples_rand1;
        centroidpos_all.(checksep).(checkcurr)=centroidpos_all.(checksep).(checkcurr)(:, samples_rand.(checksep));
        centroidpos_all.(checksep).(checkcurr)(sum(~isnan(centroidpos_all.(checksep).(checkcurr)),2)<10,:)=nan;
        pairavg_dwnsample=pairavgcur(:, samples_rand.(checksep));
        pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan;
        Csplit_dwnsample.(checksep).(checkcurr).(Tcond)=pairavg_dwnsample;

        pairavgcur_env_dwnsmple=pairavgcur_env(:, samples_rand.(checksep));
        pairavgcur_env_dwnsmple(sum(~isnan(pairavgcur_env_dwnsmple),2)<10,:)=nan;
        Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)=pairavgcur_env_dwnsmple;


        %peak calc
        pairavgcursmooth_env=pairavgcursmooth_env(:, samples_rand.(checksep));
        pairavgcursmooth_env(:,sum(pairavgcursmooth_env>0,1)==0)=nan;
        [peak,c]=find(pairavgcursmooth_env==max(pairavgcursmooth_env));
        peak_all.(checksep).(checkcurr)=(peak'-16).*50;
        if ploty==1
            figure (sepdist+1)
            if length(AMP)>1
                subplot(2,4,current)
            end
            if singlecond==0

                stdshade(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)',0.2,'r');
                hold on

                stdshade(Csplit_dwnsample.(checksep).(checkcurr).(Tcond)',0.2,'k');
                ylim([0, 360])

                xlim([find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])
            else

                plot(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                hold on
                plot(Csplit_dwnsample.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                ylim([0, 800])
                xlim([find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
            end
            title([num2str(AMP(current)) '\muA'])
            set(gca,'TickDir','out');
            xline(16,'r')
            xline(16+sepdist+1,'r')
            if sepdist==5
                xlim([16-3 16+sepdist+1+3])%([find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'first'), find(~isnan(pairavgdwnsample_all.(sepcheck).(currcheck)(:,1)), 1, 'last')])
            elseif sepdist==7
                xlim([16-2 16+sepdist+1+2])
            elseif sepdist==9
                xlim([16-1 16+sepdist+1+1])
            end
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            %xlim([13 27])
            ylim([0 350])
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
        end
        %             MSE.(checksep)(current)=mdl.MSE;
        %             rsquared.(checksep)(current)=mdl.Rsquared.Ordinary;
        %             AIC.(checksep)(current)=mdl.ModelCriterion.AIC;

    end
    if ploty==1
        figure (sepdist+1)
        f=get(gca,'Children');
        legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
    end
    if length(AMP)>1
        figure
        hold on
        meancoef1=zeros(length(AMP),2);
        for itcurr=1:length(AMP)
            checkcurr=['C' num2str(AMP(itcurr))];
            meancoef1(itcurr,1:2)=mean(mdlcur.(checksep).(checkcurr),1);

        end
        plot(AMP,   meancoef1(:,1), 'r')
        plot(AMP,   meancoef1(:,2), 'k')
        ylabel('Scaling Factor')
        xlabel('Current \muA')
        legend({'Factor 1','Factor 2'})
        set(gca,'TickDir','out');
        title(['Separation ' num2str(sepdist*50+50) '\mum'])
        xlim([1 AMP(end)])
        ylim([0 1])
    end
end

if which75_25condition==1
    b2centroid_7525.sep5=centroidpos_all.sep5.C6;
    b2centroid_7525.sep7=centroidpos_all.sep7.C6;
    b2centroid_7525.sep9=centroidpos_all.sep9.C6;
else
    b2centroid_2575.sep5=centroidpos_all.sep5.C6;
    b2centroid_2575.sep7=centroidpos_all.sep7.C6;
    b2centroid_2575.sep9=centroidpos_all.sep9.C6;
end

%% hypoth 8 - more complicated dual comparison model
%need to make sure you have the right loopsubdirect ticked!!! need full
%currents
singlecond=0;%3
ploty=1;
if which75_25condition==1
    multiply1=3/4;
    multiply2=1/4;
    Tcond='T2';
elseif which75_25condition==2
    multiply1=1/4;
    multiply2=3/4;
    Tcond='T4';
else
    error('Not valid')
end
clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        Envelope_add.(checksep).(checkcurr).(Tcond)=saveCplit.(checksep).(checkcurr).T1+saveCplit.(checksep).(checkcurr).T5;
        input_eq1=saveCplit.(checksep).(checkcurr).T1;
        input_eq2=saveCplit.(checksep).(checkcurr).T5;
        input_ans=saveCplit.(checksep).(checkcurr).(Tcond);
        fun = @(b,X) (b(1).*X(:,1))+(b(2).*X(:,2))+(multiply1.*b(3).*X(:,1))+(multiply2.*b(4).*X(:,2));
        b0=[1 1 1 1];
        mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
        mdlcur.(checksep) (current,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1} mdl.Coefficients{3,1} mdl.Coefficients{4,1}];
        if ploty==1
            figure (sepdist+10)
            subplot(2,4,current)
            plot(mean(Envelope_add.(checksep).(checkcurr).(Tcond),2))
            hold on
            plot(mean(saveCplit.(checksep).(checkcurr).(Tcond),2))
            plot(mean(saveCplit.(checksep).(checkcurr).T1,2),'k')
            plot(mean(saveCplit.(checksep).(checkcurr).T5,2),'k')
            title([num2str(AMP(current)) '\muA'])
            ylim([0, 360])
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
            %                  input_eq1(sum(input_eq1>0,2)<10,:)=nan;
            %                  input_eq2(sum(input_eq2>0,2)<10,:)=nan;
            %                  input_ans(sum(input_ans>0,2)<10,:)=nan;
            Envelope_add.(checksep).(checkcurr).(Tcond)=mdl.Coefficients{1,1}.*input_eq1+mdl.Coefficients{2,1}.*input_eq2+multiply1.*mdl.Coefficients{3,1}.*input_eq1+multiply2.*mdl.Coefficients{4,1}.*input_eq2;
            figure (sepdist+1)
            subplot(2,4,current)
            if singlecond==0
                %Envelope_add.(checksep).(checkcurr).(Tcond)(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>0,2)<5,:)=0; %if you want to take out edge conditions
                %Envelope_add.(checksep).(checkcurr).(Tcond)(Envelope_add.(checksep).(checkcurr).(Tcond)==0)=nan;
                stdshade(Envelope_add.(checksep).(checkcurr).(Tcond)',0.2,'b');
                hold on
                %saveCplit.(checksep).(checkcurr).(Tcond)(sum(saveCplit.(checksep).(checkcurr).(Tcond)>0,2)<15,:)=0; %if you want to take out edge conditions
                %saveCplit.(checksep).(checkcurr).(Tcond)(saveCplit.(checksep).(checkcurr).(Tcond)==0)=nan;
                %input_ans(sum(input_ans>0,2)<10,:)=nan;
                stdshade(input_ans',0.2,'r');
                ylim([0, 360])
                xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])

            else
                plot(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                hold on
                plot(saveCplit.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
                ylim([0, 800])
            end
            title([num2str(AMP(current)) '\muA'])
            xline(16,'r')
            xline(16+sepdist,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
        end
        MSE.(checksep)(current)=mdl.MSE;
        rsquared.(checksep)(current)=mdl.Rsquared.Adjusted;
        AIC.(checksep)(current)=mdl.ModelCriterion.AIC;
    end
    if ploty==1
        f=get(gca,'Children');
        legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
    end
end


%% hypoth 9 - more parameters
%need to make sure you have the right loopsubdirect ticked!!! need full
%currents
singlecond=0;%3
ploty=1;
if which75_25condition==1
    multiply1=3/4;
    multiply2=1/4;
    Tcond='T2';
elseif which75_25condition==2
    multiply1=1/4;
    multiply2=3/4;
    Tcond='T4';
else
    error('Not valid')
end
clear mdlcur MSE rsquared AIC
AMP=6;%[0 1 2 3 4 6 8 10];
Csplit_dwnsample=[];
Envelope_add_dwnsmple=[];
seedpoint=65;%20,55
s = RandStream('mlfg6331_64','Seed',seedpoint);
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        Envelope_add.(checksep).(checkcurr).(Tcond)=saveCplit.(checksep).(checkcurr).T1+saveCplit.(checksep).(checkcurr).T5;
        input_eq1=saveCplit.(checksep).(checkcurr).T1;
        input_eq2=saveCplit.(checksep).(checkcurr).T5;
        input_ans=saveCplit.(checksep).(checkcurr).(Tcond);
        fun = @(b,X) (multiply1.*b(1).*X(:,1))+(multiply2.*b(2).*X(:,2));
        b0=[1 1];
        mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
        mdlcur.(checksep) (current,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1}];
        Envelope_add.(checksep).(checkcurr).(Tcond)=multiply1.*mdl.Coefficients{1,1}.*input_eq1+multiply2.*mdl.Coefficients{2,1}.*input_eq2;
        for trial=1:5
            check=['T' num2str(trial)];
            dat=saveCplit.(checksep).(checkcurr).(Tcond);
            dat2=Envelope_add.(checksep).(checkcurr).(Tcond);
            dat(isinf(dat))=nan;
            clear pairavgcur erpairavgcur
            for paircount=1:size(dat,2)/4
                pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*4+1:4* paircount),2);
                pairavgcur_env(:,paircount)=nanmean(dat2(:,( paircount-1)*4+1:4* paircount),2);
            end
            lengthneeded=size(saveCplit.sep7.C1.T1,2)/4;
            if trial==1
                [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
                samples_rand.(checksep)=samples_rand1;
            end
            pairavg_dwnsample=pairavgcur(:, samples_rand.(checksep));
            pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan;
            Csplit_dwnsample.(checksep).(checkcurr).(check)=pairavg_dwnsample;

            pairavgcur_env_dwnsmple=pairavgcur_env(:, samples_rand.(checksep));
            pairavgcur_env_dwnsmple(sum(~isnan(pairavgcur_env_dwnsmple),2)<10,:)=nan;
            Envelope_add_dwnsmple.(checksep).(checkcurr).(check)=pairavgcur_env_dwnsmple;
        end


        if ploty==1
            figure (sepdist+1)
            if length(AMP)>1
                subplot(2,4,current)
            end
            if singlecond==0

                stdshade(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)',0.2,'b');
                hold on

                stdshade(Csplit_dwnsample.(checksep).(checkcurr).(Tcond)',0.2,'r');
                ylim([0, 360])

                xlim([find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])
            else

                plot(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                hold on
                plot(Csplit_dwnsample.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                ylim([0, 800])
                xlim([find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add_dwnsmple.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
            end
            title([num2str(AMP(current)) '\muA'])
            set(gca,'TickDir','out');
            xline(16,'r')
            xline(16+sepdist+1,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            xlim([13 27])
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
        end
        MSE.(checksep)(current)=mdl.MSE;
        rsquared.(checksep)(current)=mdl.Rsquared.Ordinary;
        AIC.(checksep)(current)=mdl.ModelCriterion.AIC;

    end
    if ploty==1
        figure (sepdist+1)
        f=get(gca,'Children');
        legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
    end
    if length(AMP)>1
        figure
        hold on
        plot(AMP,  mdlcur.(checksep)(:,1), 'r')
        plot(AMP,  mdlcur.(checksep)(:,2), 'k')
        ylabel('Scaling Factor')
        xlabel('Current \muA')
        legend({'Factor 1','Factor 2'})
        set(gca,'TickDir','out');
        title(['Separation ' num2str(sepdist*50+50) '\mum'])
        xlim([1 AMP(end)])
        ylim([0 2])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% altogether model #1
input_ans=[];
input_eq1=[];
input_eq2=[];
weights1=[];
weights2=[];
singlecond=0;%3
ploty=0;

clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        for trial=2:4
            checktrial=['T' num2str(trial)];
            if trial==2
                multiply1=3/4;
                multiply2=1/4;
            elseif trial==3
                multiply1=1/2;
                multiply2=1/2;
            elseif trial==4
                multiply1=1/4;
                multiply2=3/4;
            end
            input_ans=[input_ans saveCplit.(checksep).(checkcurr).(checktrial)];
            input_eq1=[input_eq1 saveCplit.(checksep).(checkcurr).T1];
            input_eq2=[input_eq2 saveCplit.(checksep).(checkcurr).T5];
            weights1=[weights1 ones(size(saveCplit.(checksep).(checkcurr).T1)).*multiply1];
            weights2=[weights2 ones(size(saveCplit.(checksep).(checkcurr).T5)).*multiply2];
        end

    end

end
fun = @(b,X) (X(:,3).*b(1).*X(:,1))+(X(:,4).*b(2).*X(:,2));
b0=[1 1];
mdl = fitnlm([input_eq1(:) input_eq2(:) weights1(:) weights2(:)], input_ans(:),fun,b0);
mdlcur.(checksep) (current,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1}];

MSE=mdl.MSE;
rsquared=mdl.Rsquared.Adjusted;
AIC=mdl.ModelCriterion.AIC;
%%%%%%%% PLOTTING
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for trial=2:4
        Tcond=['T' num2str(trial)];
        if trial==2
            multiply1=3/4;
            multiply2=1/4;
        elseif trial==3
            multiply1=1/2;
            multiply2=1/2;
        elseif trial==4
            multiply1=1/4;
            multiply2=3/4;
        end
        for current = 1:length(AMP)
            checkcurr=['C' num2str(AMP(current))];
            Envelope_add.(checksep).(checkcurr).(Tcond)=multiply1.*mdl.Coefficients{1,1}.*saveCplit.(checksep).(checkcurr).T1+multiply2.*mdl.Coefficients{2,1}.*saveCplit.(checksep).(checkcurr).T5;

            if ploty==1

                figure (sepdist+trial*5)
                subplot(2,4,current)
                if singlecond==0
                    stdshade(Envelope_add.(checksep).(checkcurr).(Tcond)',0.2,'b');
                    hold on
                    stdshade(saveCplit.(checksep).(checkcurr).(Tcond)',0.2,'r');
                    ylim([0, 360])
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])

                else
                    plot(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                    hold on
                    plot(saveCplit.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
                    ylim([0, 800])
                end
                title([num2str(AMP(current)) '\muA'])
                xline(16,'r')
                xline(16+sepdist,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
            end
        end
        if ploty==1
            f=get(gca,'Children');
            legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
        end
    end
end
%% altogether model #2
input_ans=[];
input_eq1=[];
input_eq2=[];
index=[];
singlecond=0;%3
ploty=1;
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
clear mdlcur MSE rsquared AIC_all MSE_all mdlcoef
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        for trial=2:4
            checktrial=['T' num2str(trial)];
            input_ans=[input_ans saveCplit.(checksep).(checkcurr).(checktrial)];
            input_eq1=[input_eq1 saveCplit.(checksep).(checkcurr).T1];
            input_eq2=[input_eq2 saveCplit.(checksep).(checkcurr).T5];
            temp=strings(size(saveCplit.(checksep).(checkcurr).(checktrial)));
            temp(:)=[checksep ' ' checkcurr ' ' checktrial ' '];
            tempnum=1:size(temp,2);
            temp(1,:)=strcat(temp(1,:), string(tempnum));
            index=[index temp];
             Envelope_add.(checksep).(checkcurr).(checktrial)=[];
        end

    end

end


            counter=0;
            randit = datasample(s,1:size(input_eq1,2),size(input_eq1,2),'Replace',false);
            for iterate10=floor(size(input_eq1,2)/100*10):floor(size(input_eq1,2)/100*10):size(input_eq1,2) %cross-validation
                counter=counter+1;
                fun = @(b,X) (b(1).*X(:,1))+(b(2).*X(:,2));
                b0=[1 1];
                num_columns=floor(size(input_eq1,2)/100*10);
                in1=input_eq1;
                in1(:,randit(iterate10-num_columns+1:iterate10))=nan;
                in2=input_eq2;
                in2(:,randit(iterate10-num_columns+1:iterate10))=nan;
                inans=input_ans;
                inans(:,randit(iterate10-num_columns+1:iterate10))=nan;
                mdl = fitnlm([in1(:) in2(:)], inans(:),fun,b0);
                mdlcoef(counter,1:2)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1}];
                for itall=iterate10-num_columns+1:iterate10
                    temp=index(1,randit(itall));
                    temp=split(temp);
                    checksep=temp{1};
                    checkcurr=temp{2};
                    Tcond=temp{3};
                    col=temp{4};
                    Envelope_add.(checksep).(checkcurr).(Tcond)(:,str2double(col))=mdl.Coefficients{1,1}.*input_eq1(:,randit(itall))+mdl.Coefficients{2,1}.*input_eq2(:,randit(itall));
                end
                
            end

            for itendcols=iterate10+1:size(input_eq1,2)
                    temp=index(1,randit(itendcols));
                    temp=split(temp);
                    checksep=temp{1};
                    checkcurr=temp{2};
                    Tcond=temp{3};
                    col=temp{4};
                    Envelope_add.(checksep).(checkcurr).(Tcond)(:,str2double(col))=mdl.Coefficients{1,1}.*input_eq1(:,randit(itall))+mdl.Coefficients{2,1}.*input_eq2(:,randit(itall));
            end
%      MSE=mdl.MSE;
%      rsquared=mdl.Rsquared.Adjusted;
%      AIC=mdl.ModelCriterion.AIC;

%             y=input_ans;%realdata
%             yinput=input_eq1;%inputdata1
%             yfit=mdl.Coefficients{1,1}.*input_eq1+mdl.Coefficients{2,1}.*input_eq2;%modeldata
%             numparam=2;
%             [rsq,AIC]=modelfitcharacterisation(y,yfit,yinput,numparam);


%%%%%%%% PLOTTING
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for trial=2:4
        Tcond=['T' num2str(trial)];
        for current = 2:length(AMP)
            checkcurr=['C' num2str(AMP(current))];
            y=saveCplit.(checksep).(checkcurr).(Tcond);%realdata
            yinput=saveCplit.(checksep).(checkcurr).T1;%inputdata1
            yfit=Envelope_add.(checksep).(checkcurr).(Tcond);%modeldata
            numparam=2;
            [rsq,AIC,MSE]=modelfitcharacterisation(y,yfit,yinput,numparam);
            rsquared.(Tcond).(checksep)(current)=rsq;
            AIC_all.(Tcond).(checksep)(current)=AIC;
            MSE_all.(Tcond).(checksep)(current)=MSE;
            if ploty==1

                figure (sepdist+trial*5)
                subplot(2,4,current)
                if singlecond==0
                    stdshade(Envelope_add.(checksep).(checkcurr).(Tcond)',0.2,'b');
                    hold on
                    stdshade(saveCplit.(checksep).(checkcurr).(Tcond)',0.2,'r');
                    ylim([0, 360])
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])

                else
                    plot(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                    hold on
                    plot(saveCplit.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
                    ylim([0, 800])
                end
                title([num2str(AMP(current)) '\muA'])
                xline(16,'r')
                xline(16+sepdist,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
            end
        end
        if ploty==1
            f=get(gca,'Children');
            legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
        end
    end
end

%% altogether model #3
input_ans=[];
input_eq1=[];
input_eq2=[];
weights1=[];
weights2=[];
singlecond=0;%3
ploty=0;

clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        for trial=2:4
            checktrial=['T' num2str(trial)];
            if trial==2
                multiply1=3/4;
                multiply2=1/4;
            elseif trial==3
                multiply1=1/2;
                multiply2=1/2;
            elseif trial==4
                multiply1=1/4;
                multiply2=3/4;
            end
            input_ans=[input_ans saveCplit.(checksep).(checkcurr).(checktrial)];
            input_eq1=[input_eq1 saveCplit.(checksep).(checkcurr).T1];
            input_eq2=[input_eq2 saveCplit.(checksep).(checkcurr).T5];
            weights1=[weights1 ones(size(saveCplit.(checksep).(checkcurr).T1)).*multiply1];
            weights2=[weights2 ones(size(saveCplit.(checksep).(checkcurr).T5)).*multiply2];
        end

    end

end
fun = @(b,X) (b(1).*X(:,1))+(b(2).*X(:,2))+(X(:,3).*b(3).*X(:,1))+(X(:,4).*b(4).*X(:,2));
b0=[1 1 1 1];
mdl = fitnlm([input_eq1(:) input_eq2(:) weights1(:) weights2(:)], input_ans(:),fun,b0);
mdlcur.(checksep) (current,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1} mdl.Coefficients{3,1} mdl.Coefficients{4,1}];

MSE=mdl.MSE;
rsquared=mdl.Rsquared.Adjusted;
AIC=mdl.ModelCriterion.AIC;
%%%%%%%% PLOTTING
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for trial=2:4
        Tcond=['T' num2str(trial)];
        if trial==2
            multiply1=3/4;
            multiply2=1/4;
        elseif trial==3
            multiply1=1/2;
            multiply2=1/2;
        elseif trial==4
            multiply1=1/4;
            multiply2=3/4;
        end
        for current = 1:length(AMP)
            checkcurr=['C' num2str(AMP(current))];
            Envelope_add.(checksep).(checkcurr).(Tcond)=mdl.Coefficients{1,1}.*saveCplit.(checksep).(checkcurr).T1+mdl.Coefficients{2,1}.*saveCplit.(checksep).(checkcurr).T5+multiply1.*mdl.Coefficients{3,1}.*saveCplit.(checksep).(checkcurr).T1+multiply2.*mdl.Coefficients{4,1}.*saveCplit.(checksep).(checkcurr).T5;

            if ploty==1

                figure (sepdist+trial*5)
                subplot(2,4,current)
                if singlecond==0
                    stdshade(Envelope_add.(checksep).(checkcurr).(Tcond)',0.2,'b');
                    hold on
                    stdshade(saveCplit.(checksep).(checkcurr).(Tcond)',0.2,'r');
                    ylim([0, 360])
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])

                else
                    plot(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                    hold on
                    plot(saveCplit.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
                    ylim([0, 800])
                end
                title([num2str(AMP(current)) '\muA'])
                xline(16,'r')
                xline(16+sepdist,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
            end
        end
        if ploty==1
            f=get(gca,'Children');
            legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
        end
    end
end

%% altogether model #4
input_ans=[];
input_eq1=[];
input_eq2=[];
weights1=[];
weights2=[];
singlecond=0;%3
ploty=0;

clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        for trial=2:4
            checktrial=['T' num2str(trial)];
            if trial==2
                multiply1=3/4;
                multiply2=1/4;
            elseif trial==3
                multiply1=1/2;
                multiply2=1/2;
            elseif trial==4
                multiply1=1/4;
                multiply2=3/4;
            end
            input_ans=[input_ans saveCplit.(checksep).(checkcurr).(checktrial)];
            input_eq1=[input_eq1 saveCplit.(checksep).(checkcurr).T1];
            input_eq2=[input_eq2 saveCplit.(checksep).(checkcurr).T5];
            weights1=[weights1 ones(size(saveCplit.(checksep).(checkcurr).T1)).*multiply1];
            weights2=[weights2 ones(size(saveCplit.(checksep).(checkcurr).T5)).*multiply2];
        end

    end

end
fun = @(b,X) (b(1).*X(:,1))+(b(2).*X(:,2))+(X(:,3).*b(1).*X(:,1))+(X(:,4).*b(2).*X(:,2));
b0=[1 1];
mdl = fitnlm([input_eq1(:) input_eq2(:) weights1(:) weights2(:)], input_ans(:),fun,b0);
mdlcur.(checksep) (current,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1}];

MSE=mdl.MSE;
rsquared=mdl.Rsquared.Adjusted;
AIC=mdl.ModelCriterion.AIC;
%%%%%%%% PLOTTING
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for trial=2:4
        Tcond=['T' num2str(trial)];
        if trial==2
            multiply1=3/4;
            multiply2=1/4;
        elseif trial==3
            multiply1=1/2;
            multiply2=1/2;
        elseif trial==4
            multiply1=1/4;
            multiply2=3/4;
        end
        for current = 1:length(AMP)
            checkcurr=['C' num2str(AMP(current))];
            Envelope_add.(checksep).(checkcurr).(Tcond)=mdl.Coefficients{1,1}.*saveCplit.(checksep).(checkcurr).T1+mdl.Coefficients{2,1}.*saveCplit.(checksep).(checkcurr).T5+multiply1.*mdl.Coefficients{1,1}.*saveCplit.(checksep).(checkcurr).T1+multiply2.*mdl.Coefficients{2,1}.*saveCplit.(checksep).(checkcurr).T5;

            if ploty==1

                figure (sepdist+trial*5)
                subplot(2,4,current)
                if singlecond==0
                    stdshade(Envelope_add.(checksep).(checkcurr).(Tcond)',0.2,'b');
                    hold on
                    stdshade(saveCplit.(checksep).(checkcurr).(Tcond)',0.2,'r');
                    ylim([0, 360])
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])

                else
                    plot(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                    hold on
                    plot(saveCplit.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
                    ylim([0, 800])
                end
                title([num2str(AMP(current)) '\muA'])
                xline(16,'r')
                xline(16+sepdist,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
            end
        end
        if ploty==1
            f=get(gca,'Children');
            legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
        end
    end
end

%% altogether model #5
input_ans=[];
input_eq1=[];
input_eq2=[];
weights1=[];
weights2=[];
singlecond=0;%3
ploty=0;

clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        for trial=2:4
            checktrial=['T' num2str(trial)];
            if trial==2
                multiply1=3/4;
                multiply2=1/4;
            elseif trial==3
                multiply1=1/2;
                multiply2=1/2;
            elseif trial==4
                multiply1=1/4;
                multiply2=3/4;
            end
            input_ans=[input_ans saveCplit.(checksep).(checkcurr).(checktrial)];
            input_eq1=[input_eq1 saveCplit.(checksep).(checkcurr).T1];
            input_eq2=[input_eq2 saveCplit.(checksep).(checkcurr).T5];
            weights1=[weights1 ones(size(saveCplit.(checksep).(checkcurr).T1)).*multiply1];
            weights2=[weights2 ones(size(saveCplit.(checksep).(checkcurr).T5)).*multiply2];
        end

    end

end
fun = @(b,X) (X(:,3).*b(1).*X(:,1))+(X(:,4).*b(2).*X(:,2))+b(3);
b0=[1 1 1];
mdl = fitnlm([input_eq1(:) input_eq2(:) weights1(:) weights2(:)], input_ans(:),fun,b0);
mdlcur.(checksep) (current,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1}];

MSE=mdl.MSE;
rsquared=mdl.Rsquared.Adjusted;
AIC=mdl.ModelCriterion.AIC;
%%%%%%%% PLOTTING
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for trial=2:4
        Tcond=['T' num2str(trial)];
        if trial==2
            multiply1=3/4;
            multiply2=1/4;
        elseif trial==3
            multiply1=1/2;
            multiply2=1/2;
        elseif trial==4
            multiply1=1/4;
            multiply2=3/4;
        end
        for current = 1:length(AMP)
            checkcurr=['C' num2str(AMP(current))];
            Envelope_add.(checksep).(checkcurr).(Tcond)=multiply1.*mdl.Coefficients{1,1}.*saveCplit.(checksep).(checkcurr).T1+multiply2.*mdl.Coefficients{2,1}.*saveCplit.(checksep).(checkcurr).T5;

            if ploty==1

                figure (sepdist+trial*5)
                subplot(2,4,current)
                if singlecond==0
                    stdshade(Envelope_add.(checksep).(checkcurr).(Tcond)',0.2,'b');
                    hold on
                    stdshade(saveCplit.(checksep).(checkcurr).(Tcond)',0.2,'r');
                    ylim([0, 360])
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])

                else
                    plot(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                    hold on
                    plot(saveCplit.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
                    ylim([0, 800])
                end
                title([num2str(AMP(current)) '\muA'])
                xline(16,'r')
                xline(16+sepdist,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
            end
        end
        if ploty==1
            f=get(gca,'Children');
            legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
        end
    end
end
%% altogether model #6
input_ans=[];
input_eq1=[];
input_eq2=[];
weights1=[];
weights2=[];
singlecond=0;%3
ploty=1;

clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        for trial=2:4
            checktrial=['T' num2str(trial)];
            if trial==2
                multiply1=3/4;
                multiply2=1/4;
            elseif trial==3
                multiply1=1/2;
                multiply2=1/2;
            elseif trial==4
                multiply1=1/4;
                multiply2=3/4;
            end
            input_ans=[input_ans saveCplit.(checksep).(checkcurr).(checktrial)];
            input_eq1=[input_eq1 saveCplit.(checksep).(checkcurr).T1];
            input_eq2=[input_eq2 saveCplit.(checksep).(checkcurr).T5];
            weights1=[weights1 ones(size(saveCplit.(checksep).(checkcurr).T1)).*multiply1];
            weights2=[weights2 ones(size(saveCplit.(checksep).(checkcurr).T5)).*multiply2];
        end

    end

end

%%%%%%%% PLOTTING
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for trial=2:4
        Tcond=['T' num2str(trial)];
        if trial==2
            multiply1=3/4;
            multiply2=1/4;
        elseif trial==3
            multiply1=1/2;
            multiply2=1/2;
        elseif trial==4
            multiply1=1/4;
            multiply2=3/4;
        end
        for current = 1:length(AMP)
            checkcurr=['C' num2str(AMP(current))];
            Envelope_add.(checksep).(checkcurr).(Tcond)=multiply1.*saveCplit.(checksep).(checkcurr).T1+multiply2.*saveCplit.(checksep).(checkcurr).T5;

            if ploty==1

                figure (sepdist+trial*5)
                subplot(2,4,current)
                if singlecond==0
                    stdshade(Envelope_add.(checksep).(checkcurr).(Tcond)',0.2,'b');
                    hold on
                    stdshade(saveCplit.(checksep).(checkcurr).(Tcond)',0.2,'r');
                    ylim([0, 360])
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])

                else
                    plot(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                    hold on
                    plot(saveCplit.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
                    ylim([0, 800])
                end
                title([num2str(AMP(current)) '\muA'])
                xline(16,'r')
                xline(16+sepdist,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
            end
        end
        if ploty==1
            f=get(gca,'Children');
            legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
        end
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check scaling factor bias
input_ans=[];
input_eq1=[];
input_eq2=[];
preferredelectresp=0;
nonpreferredelectresp=0;
singlecond=0;%3
preferredelectrespNUM=0;
nonpreferredelectrespNUM=0;

clear mdlcur MSE rsquared AIC
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];

    for current = 2:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        sumelectpreferred1.(checksep).(checkcurr)=0;
        sumelectpreferred2.(checksep).(checkcurr)=0;
        sumNUMelectpreferred1.(checksep).(checkcurr)=0;
        sumNUMelectpreferred2.(checksep).(checkcurr)=0;
        for Pair=1:size(saveCplit.(checksep).(checkcurr).T1,2)
            checkpair=['P' num2str(Pair)];

            rate1=saveCplit.(checksep).(checkcurr).T1(:,Pair)-saveCplit.(checksep).(checkcurr).T5(:,Pair);
            electpreferred=(sum(saveCplit.(checksep).(checkcurr).T1(:,Pair)-saveCplit.(checksep).(checkcurr).T5(:,Pair)>0))/sum(~isnan(saveCplit.(checksep).(checkcurr).T1(:,Pair)-saveCplit.(checksep).(checkcurr).T5(:,Pair)));
            sumNUMelectpreferred1.(checksep).(checkcurr)=sumNUMelectpreferred1.(checksep).(checkcurr)+ double(electpreferred>0.5);
            sumNUMelectpreferred2.(checksep).(checkcurr)=sumNUMelectpreferred2.(checksep).(checkcurr)+double(electpreferred<0.5);
            sumelectpreferred1.(checksep).(checkcurr)=sumelectpreferred1.(checksep).(checkcurr)+double(nansum(rate1)>0);
            sumelectpreferred2.(checksep).(checkcurr)=sumelectpreferred2.(checksep).(checkcurr)+double(nansum(rate1)<0);
            %                   firstnonnan=find(~isnan(saveCplit.(checksep).(checkcurr).(checktrial)(:,1)),1,'first');
            %                   lastnonnan=find(~isnan(saveCplit.(checksep).(checkcurr).(checktrial)(:,1)),1,'last');
            %                 rate1(rate1<0)=0;
            %                 A = trapz(firstnonnan:lastnonnan, rate1(firstnonnan:lastnonnan));
            %
            %                 B=nan(32,1);
            %                 B(1)=0;
            %                 for lims = firstnonnan+1:lastnonnan
            %                     B(lims) =  trapz(firstnonnan:lims, rate1(firstnonnan:lims));
            %                 end
            %
            %                 [~,electrodecentroid]=min(abs(B-(A/2)));
            %
            %                 electpos=[16 16+sepdist+1];
            %                 [~,preferredElect]=min(abs(electpos-electrodecentroid));

            input_ans= saveCplit.(checksep).(checkcurr).(Tcond)(:,Pair);
            input_eq1=saveCplit.(checksep).(checkcurr).T1(:,Pair);
            input_eq2=saveCplit.(checksep).(checkcurr).T5(:,Pair);

            fun = @(b,X) (b(1).*X(:,1))+(b(2).*X(:,2));
            b0=[1 1];
            mdl = fitnlm([input_eq1(:) input_eq2(:)], input_ans(:),fun,b0);
            mdlcur.(checksep) (current,:)=[mdl.Coefficients{1,1} mdl.Coefficients{2,1}];
            if ((mdl.Coefficients{1,1}>mdl.Coefficients{2,1}) && nansum(rate1)>0) || ((mdl.Coefficients{1,1}<mdl.Coefficients{2,1}) && nansum(rate1)<0)
                preferredelectresp=preferredelectresp+1;

            else
                nonpreferredelectresp=nonpreferredelectresp+1;
            end

            if ((mdl.Coefficients{1,1}>mdl.Coefficients{2,1}) &&  electpreferred>0.5) || ((mdl.Coefficients{1,1}<mdl.Coefficients{2,1}) &&  electpreferred<0.5)
                preferredelectrespNUM=preferredelectrespNUM+1;

            else
                nonpreferredelectrespNUM=nonpreferredelectrespNUM+1;
            end


        end

    end

end



%% mean of scaling factors

algethermeanmdl=[mean(mdlcur_2575(2:end,:));mean(mdlcur_5050(2:end,:));mean(mdlcur_7525(2:end,:))];
currratio=[0.25 0.25; 0.5 0.5; 0.75 0.75];
altogetherstdmdl=[nanstd(mdlcur_2575(2:end,:),[],1)./sqrt(sum(~isnan(mdlcur_2575(2:end,:)),1));nanstd(mdlcur_5050(2:end,:),[],1)./sqrt(sum(~isnan(mdlcur_5050(2:end,:)),1));nanstd(mdlcur_7525(2:end,:),[],1)./sqrt(sum(~isnan(mdlcur_7525(2:end,:)),1))];
figure;  hold on; scatter(currratio(:,1),algethermeanmdl(:,1),'filled','b'); scatter(currratio(:,2),algethermeanmdl(:,2),'filled','r');
errorbar(currratio(:,1),algethermeanmdl(:,1), altogetherstdmdl(:,1), 'LineStyle','none','Color','b');
errorbar(currratio(:,2),algethermeanmdl(:,2), altogetherstdmdl(:,2), 'LineStyle','none','Color','r');
xticks([0.25 0.5 0.75])
xlim([0 1])
xlabel('Current ratio')
ylabel('Scaling factor')
ploty=0;
c1 = polyfit(currratio(:,1),algethermeanmdl(:,1),1);
c2 = polyfit(currratio(:,2),algethermeanmdl(:,2),1);
y_est1 = polyval(c1,currratio(:,1));
y_est2 = polyval(c2,currratio(:,2));
plot(currratio(:,1),y_est1,'b')
plot(currratio(:,2),y_est2,'r')
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for trial=2:4
        Tcond=['T' num2str(trial)];
        if trial==2
            multiply1=algethermeanmdl(3,1);
            multiply2=algethermeanmdl(3,2);
        elseif trial==3
            multiply1=algethermeanmdl(2,1);
            multiply2=algethermeanmdl(2,2);
        elseif trial==4
            multiply1=algethermeanmdl(1,1);
            multiply2=algethermeanmdl(1,2);
        end
        for current = 1:length(AMP)
            checkcurr=['C' num2str(AMP(current))];
            Envelope_add.(checksep).(checkcurr).(Tcond)=multiply1.*saveCplit.(checksep).(checkcurr).T1+multiply2.*saveCplit.(checksep).(checkcurr).T5;

            if ploty==1

                figure (sepdist+trial*5)
                subplot(2,4,current)
                if singlecond==0
                    stdshade(Envelope_add.(checksep).(checkcurr).(Tcond)',0.2,'b');
                    hold on
                    stdshade(saveCplit.(checksep).(checkcurr).(Tcond)',0.2,'r');
                    ylim([0, 360])
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)>=0,2),1,'last')+1])

                else
                    plot(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond),'b')
                    hold on
                    plot(saveCplit.(checksep).(checkcurr).(Tcond)(:,singlecond),'r')
                    xlim([find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'first')-1, find(sum(Envelope_add.(checksep).(checkcurr).(Tcond)(:,singlecond)>=0,2),1,'last')+1])
                    ylim([0, 800])
                end
                title([num2str(AMP(current)) '\muA'])
                xline(16,'r')
                xline(16+sepdist,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
            end
        end
        if ploty==1
            f=get(gca,'Children');
            legend([f(3),f(5)],'real_d_u_a_l','model_d_u_a_l')
        end
    end
end
%% USING SIG CURVES TO RECREATE MODEL
%% Make curves
%spike rate curves
sigcurve_pool=[];
AMP=[0 1 2 3 4 6 8 10];
dat=[];
figure (1)
hold on
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for  elect=1:2:5
        for current=1:length(AMP)
            checkcurr=['C' num2str(AMP(current))];
            dat=saveCplit.(checksep).(checkcurr).(['T' num2str(elect)]);
            dat(dat==0)=nan;
            sigcurve_pool.(checksep)(elect,current)=nanmean(dat,'all'); %curves
        end
        plot(AMP, sigcurve_pool.(checksep)(elect,:))
    end
end

%% Invert and scale
figure
hold on
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    weights_theoretical.(checksep)=(sigcurve_pool.(checksep)./max(sigcurve_pool.(checksep),[],2))+1;
    plot(AMP,1./weights_theoretical.(checksep) (1,:));
    plot(AMP,1./weights_theoretical.(checksep) (5,:));
end












%%
clear mdlcur
for current = 1:length(AMP)
    checkcurr=['C' num2str(AMP(current))];
    fun = @(b,X) b(1).*(X(:,1)+X(:,2))+b(2);
    b0=[1 1];
    mdl = lsqcurvefit(fun,b0,[[saveCplit.sep5.(checkcurr).T1(:);saveCplit.sep7.(checkcurr).T1(:);saveCplit.sep9.(checkcurr).T1(:)], [saveCplit.sep5.(checkcurr).T5(:);saveCplit.sep7.(checkcurr).T5(:);saveCplit.sep9.(checkcurr).T5(:)]], [saveCplit.sep5.(checkcurr).(Tcond)(:);saveCplit.sep7.(checkcurr).(Tcond)(:);saveCplit.sep9.(checkcurr).(Tcond)(:)]);
    mdlcur(current,:)=mdl;
end

%% hypoth 2 - additive scaling
%individual spike rate curves
clear sigcurve

AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        for row=1:32
            checkrow=['row' num2str(row)];
            for col=1:size(saveCplit.(checksep).(checkcurr).T1,2)
                checkcol=['col' num2str(col)];
                sigcurve.(checksep).(checkrow).(checkcol)(1,current)=saveCplit.(checksep).(checkcurr).T1(row,col);
                sigcurve.(checksep).(checkrow).(checkcol)(2,current)=AMP(current);
                if current==length(AMP)
                    sigcurve.(checksep).(checkrow).(checkcol)(1,:)=sigcurve.(checksep).(checkrow).(checkcol)(1,:).*2./max(sigcurve.(checksep).(checkrow).(checkcol)(1,:));
                end
            end
        end

    end
end

%% pool based on row
clear sigcurve  sigcurve_avgepos sigcurve2 sigcurve_avgepos2
AMP=[0 1 2 3 4 6 8 10];
figure
hold on
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];
        for row=1:32
            checkrow=['row' num2str(row)];
            for col=1:size(saveCplit.(checksep).(checkcurr).T1,2)
                checkcol=['col' num2str(col)];
                sigcurve.(checksep).(checkrow)(col,current)=saveCplit.(checksep).(checkcurr).T1(row,col);
            end

            for col=1:size(saveCplit.(checksep).(checkcurr).T5,2)
                checkcol=['col' num2str(col)];
                sigcurve2.(checksep).(checkrow)(col,current)=saveCplit.(checksep).(checkcurr).T5(row,col);
            end
            sigcurve.(checksep).(checkrow)(sigcurve.(checksep).(checkrow)==0)=nan;
            sigcurve2.(checksep).(checkrow)(sigcurve2.(checksep).(checkrow)==0)=nan;
            sigcurve_avgepos.(checksep)(row,current)=nanmean(sigcurve.(checksep).(checkrow)(:,current));
            sigcurve_avgepos2.(checksep)(row,current)=nanmean(sigcurve2.(checksep).(checkrow)(:,current));
        end
    end
    plot(AMP,nanmean(sigcurve_avgepos.(checksep)))
    plot(AMP,nanmean(sigcurve_avgepos2.(checksep)))
end


%% %% fit each row exponential
clear mdl_params
for sepdist=5:2:9
    for row=1:32
        checkrow=['row' num2str(row)];
        checksep=['sep' num2str(sepdist)];
        for col=1:size(saveCplit.(checksep).(checkcurr).T1,2)
            checkcol=['col' num2str(col)];
            fun = @(b,X) (b(1)-(b(1)-b(2)).*exp(-b(3).*X));
            b0=[200 1 1];
            if sum(isnan(sigcurve.(checksep).(checkrow)(col,:)))==0
                mdl = lsqcurvefit(fun,b0,AMP, sigcurve.(checksep).(checkrow)(col,:));
                if mdl(1)<700 && isreal(mdl(1))
                    mdl_params.(checksep)(row,:)=mdl;
                end
            end
            if sum(isnan(sigcurve2.(checksep).(checkrow)(col,:)))==0
                mdl = lsqcurvefit(fun,b0,AMP, sigcurve2.(checksep).(checkrow)(col,:));
                if mdl(1)<700 && isreal(mdl(1))
                    mdl_params.(checksep)(row+50,:)=mdl;
                end
            end
        end
    end
    mdl_params.(checksep)(mdl_params.(checksep)==0)=nan;
end

%% check params exponential
clear curve_norm
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    YM=nanmean(mdl_params.(checksep)(:,1));
    Y0=nanmean(mdl_params.(checksep)(:,2));
    k=nanmean(mdl_params.(checksep)(:,3));
    x=[AMP];
    curve=YM-(YM-Y0).*exp(-k.*x);
    %     figure
    %     plot(x,curve)
    curve_norm.(checksep)=curve.*2./max(curve);
    figure
    plot(x,curve_norm.(checksep))
end

%%
clear Envelope_add Envelope_normalise MSE rsquared
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    figure
    for current = 1:length(AMP)
        checkcurr=['C' num2str(AMP(current))];

        Envelope_add.(checksep).(checkcurr).T3=(saveCplit.(checksep).(checkcurr).T1+saveCplit.(checksep).(checkcurr).T5)./ curve_norm.(checksep)(current);%./(AMP(current)^(1/sigcurve_normalise(current)));
        subplot(2,4,current)
        plot(mean(Envelope_add.(checksep).(checkcurr).T3,2))
        hold on
        plot(mean(saveCplit.(checksep).(checkcurr).T3,2))
        title([num2str(AMP(current)) '\muA'])
        ylim([0, 300])
        xline(16,'r')
        xline(16+sepdist,'r')
        xt = xticks;
        xtl=(xt-16+(-sepdist/2))*50;
        xticklabels(xtl)
        ylabel('Sp/s')
        xlabel('Distance from deepest stim elect (\mum)')
        RealDat=(mean(saveCplit.(checksep).(checkcurr).T3,2));
        ModelledDat=mean(Envelope_add.(checksep).(checkcurr).T3,2);
        RealDatNO0=(RealDat((RealDat==0)+(ModelledDat==0)~=2));
        ModelledDatNO0=(ModelledDat((RealDat==0)+(ModelledDat==0)~=2));
        MSE.(checksep)(current)=(sum((RealDatNO0-ModelledDatNO0).^2))/length(ModelledDat);
        rsquared.(checksep)(current)=1-(sum((RealDatNO0-ModelledDatNO0).^2)./sum((RealDatNO0-mean(RealDatNO0)).^2));
    end
end


%% hypoth 3 - normalisation?

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% everything I tried

%% this works
alpha=199.5;%max(sigcurvepool_avgall);
gamma=1.98;%%AMP((abs(sigcurvepool_avgall-(max(sigcurvepool_avgall)./2))==min(abs(sigcurvepool_avgall-(max(sigcurvepool_avgall)./2)))));
delta=30.5;%sigcurvepool_avgall(1);
x=AMP;
beta=0.8;%4
curve=(alpha.*x.^beta./((x.^beta) + (gamma.^beta)))+ delta;
curve=curve.*2./max(curve);%*1./max(curve)+1;
figure;plot(AMP, curve)
fun = @(b,X) (b(1).*X.^b(2)./((X.^b(2)) + (b(3).^b(2))))+ b(4);
b0=[200 1 3 8];
mdl = lsqcurvefit(fun,b0,AMP,nanmean([sigcurve_avgepos.sep5]));

%% exponential plateau
YM=199.5;
Y0=30.5;
k=2;
x=AMP;
curve=YM-(YM-Y0).*exp(-k*x);
curve=curve.*2./max(curve);
figure;plot(AMP, curve)


%%
alpha=mdl(1);
beta=mdl(2);
gamma=mdl(3);
delta=mdl(4);
x=[AMP];
curve=(alpha.*x.^beta./((x.^beta) + (gamma.^beta)))+ delta;
figure
plot(x,curve)
curve_norm=curve.*2./max(curve);
figure
plot(x,curve_norm)
%%
%spike rate curves
sigcurve_pool=[];
AMP=[0 1 2 3 4 6 8 10];
dat=[];
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    for  elect=1:4:5
        for current=1:length(AMP)
            checkcurr=['C' num2str(AMP(current))];
            dat=saveCplit.(checksep).(checkcurr).(['T' num2str(elect)]);
            dat(dat==0)=nan;
            sigcurve_pool.(checksep)(elect,current)=nanmean(dat,'all');
        end
    end
end
%sigcurve_normalise=sigcurve_pool.*2./max(sigcurve_pool);
%% %% fit each row sigmoid
clear mdl_params
for sepdist=5:2:9
    for row=1:32
        checkrow=['row' num2str(row)];
        checksep=['sep' num2str(sepdist)];
        for col=1:size(saveCplit.(checksep).(checkcurr).T1,2)
            checkcol=['col' num2str(col)];
            fun = @(b,X) (b(1).*X.^b(2)./((X.^b(2)) + (b(3).^b(2))))+ b(4);
            b0=[200 1 3 8];
            if sum(isnan(sigcurve.(checksep).(checkrow)(col,:)))==0
                mdl = lsqcurvefit(fun,b0,AMP, sigcurve.(checksep).(checkrow)(col,:));
                if mdl(1)<700 && isreal(mdl(1)) && (sum(mdl<0)==0) && mdl(2)<=10
                    mdl_params.(checksep)(row,:)=mdl;
                end
            end
            if sum(isnan(sigcurve2.(checksep).(checkrow)(col,:)))==0
                mdl = lsqcurvefit(fun,b0,AMP, sigcurve2.(checksep).(checkrow)(col,:));
                if mdl(1)<700 && isreal(mdl(1)) && (sum(mdl<0)==0) && mdl(2)<=10
                    mdl_params.(checksep)(row+50,:)=mdl;
                end
            end
        end
    end
    mdl_params.(checksep)(mdl_params.(checksep)==0)=nan;
end
%% check params
clear curve_norm
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    alpha=nanmean(mdl_params.(checksep)(:,1));
    beta=nanmean(mdl_params.(checksep)(:,2));
    gamma=nanmean(mdl_params.(checksep)(:,3));
    delta=nanmean(mdl_params.(checksep)(:,4));
    x=[AMP];
    curve=(alpha.*x.^beta./((x.^beta) + (gamma.^beta)))+ delta;
    figure
    plot(x,curve)
    curve_norm.(checksep)=curve.*2./max(curve);
    % figure
    % plot(x,curve_norm)
end


%% %% fit each row Gompertz growth
clear mdl_params
for sepdist=5:2:9
    for row=1:32
        checkrow=['row' num2str(row)];
        checksep=['sep' num2str(sepdist)];
        for col=1:size(saveCplit.(checksep).(checkcurr).T1,2)
            checkcol=['col' num2str(col)];
            fun = @(b,X) b(1).*((b(2)/b(1)).^exp(-b(3).*X));
            b0=[200 1 1];
            if sum(isnan(sigcurve.(checksep).(checkrow)(col,:)))==0
                mdl = lsqcurvefit(fun,b0,AMP, sigcurve.(checksep).(checkrow)(col,:));
                if mdl(1)<700 && isreal(mdl(1))
                    mdl_params.(checksep)(row,:)=mdl;
                end
            end
            if sum(isnan(sigcurve2.(checksep).(checkrow)(col,:)))==0
                mdl = lsqcurvefit(fun,b0,AMP, sigcurve2.(checksep).(checkrow)(col,:));
                if mdl(1)<700 && isreal(mdl(1))
                    mdl_params.(checksep)(row+50,:)=mdl;
                end
            end
        end
    end
    mdl_params.(checksep)(mdl_params.(checksep)==0)=nan;
end

%% check params exponential
clear curve_norm
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
    YM=nanmean(mdl_params.(checksep)(:,1));
    Y0=nanmean(mdl_params.(checksep)(:,2));
    k=nanmean(mdl_params.(checksep)(:,3));
    x=[AMP];
    curve=YM.*((Y0/YM).^exp(-k.*x));
    %     figure
    %     plot(x,curve)
    curve_norm.(checksep)=curve.*2./max(curve);
    figure
    plot(x,curve_norm.(checksep))
end

%%


sigcurvepool_avgall=nanmean([sigcurve_avgepos.sep5; sigcurve_avgepos.sep7; sigcurve_avgepos.sep9; sigcurve_avgepos2.sep5; sigcurve_avgepos2.sep7; sigcurve_avgepos2.sep9;]);
% alpha=199.5;%max(sigcurvepool_avgall);
% gamma=2.45;%%AMP((abs(sigcurvepool_avgall-(max(sigcurvepool_avgall)./2))==min(abs(sigcurvepool_avgall-(max(sigcurvepool_avgall)./2)))));
% delta=30.5;%sigcurvepool_avgall(1);
% x=AMP;
% beta=0.8;

