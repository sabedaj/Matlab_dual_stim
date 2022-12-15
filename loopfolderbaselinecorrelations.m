clear
folders=dir;
matsave=[];
%% correlations
for loopfold=3:size(folders,1)
    try
        cd([folders(loopfold).folder filesep folders(loopfold).name])
        load('ElectLayerClass.mat','ElectLayerClass')
    catch
        continue
    end
    
    
    folders_rat=dir;
    timesmall=3;
    for loopfoldrat=4:size(folders_rat,1)
        if ((str2double(folders_rat(loopfoldrat).name(end-5:end))<str2double(folders_rat(timesmall).name(end-5:end))) || isnan(str2double(folders_rat(timesmall).name(end-5:end)))) && ~strcmpi(folders_rat(loopfoldrat).name(1:5),'flash')
        timesmall=loopfoldrat;
        end
    end
    cd([folders_rat(timesmall).folder filesep folders_rat(timesmall).name])
    maxcor=noisecorrelations(ElectLayerClass);
    %%
    shankpoint=nan(4,136);
    
 for shank=1:4
     index=1;
    r=1+(shank*16-16);
        for c=1+(shank*16-16):(shank*16)
            shankpoint(shank,index:index+length((c-r+1):16)-1)=maxcor(c,(c-r+1)+(shank*16-16):16+(shank*16-16));
            index=index+length((c-r+1):16);
        end
 end
 %% var per shank per electrode excluding En corr En
 S = spdiags(nan(64,64),0,maxcor);
 maxcor_noselfcor=full(S);
 var_matrix=nan(64,4);
 for elect=1:64
  for shank=1:4
       var_matrix(elect,shank)=std(maxcor_noselfcor(elect, (shank-1)*16+1:shank*16),'omitnan');
  end
 end
 cd([folders(loopfold).folder filesep folders(loopfold).name])
 save('VarianceCorrelation.mat','var_matrix')
 
 %%
 
shankpoint(:, [1 1+16 1+16+15 1+16+15+14 1+16+15+14+13 1+16+15+14+13+12 1+16+15+14+13+12+11 1+16+15+14+13+12+11+10 1+16+15+14+13+12+11+10+9 1+16+15+14+13+12+11+10+9+8 1+16+15+14+13+12+11+10+9+8+7 1+16+15+14+13+12+11+10+9+8+7+6 1+16+15+14+13+12+11+10+9+8+7+6+5 1+16+15+14+13+12+11+10+9+8+7+6+5+4 1+16+15+14+13+12+11+10+9+8+7+6+5+4+3 1+16+15+14+13+12+11+10+9+8+7+6+5+4+3+2])=[];
saveshankVAR.(folders(loopfold).name)=shankpoint;
variance_animal(loopfold-2)=var(shankpoint,[],'all','omitnan');
halfcormatrix=tril(maxcor',-1);
halfcormatrix=halfcormatrix(halfcormatrix~=0);
matsave=[matsave halfcormatrix];
variance_animal_halfcor(loopfold-2)=var(halfcormatrix,[],'all','omitnan');
end

