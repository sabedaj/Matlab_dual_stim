function allExtractBR(dName,nChn,FS,T,par)
%% This function takes a blanked Blackrock datafile (ns6_dn)
% and constructs mu, lfp, and sp datafiles
%% Variables
threshfac = -3.5;
sp = cell(1, nChn);
thresh = cell(1, nChn);
SCALEFACTOR = 10;
chk = 1; N = 0; time = 0;
NSp = zeros(1,nChn);
Ninput = 1e6;
ARTMAX = 1.5e3;
if isempty(dir('*muBR.dat'))
    justMu = 1;
else
    justMu = 0;
end
%% Load in the amplifier waveform data
try
    fileinfo = dir([dName '*_dnBR.dat']);
    nSam = fileinfo.bytes/(nChn * 2);
    v_fid = fopen([dName '_dnBR.dat'], 'r');
catch
    fprintf('No amplifier data found. Aborting . . .\n');
    return; 
end
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);
%% Generate filters
[Mufilt,Lfpfilt] = generate_Filters;
MuNf = length(Mufilt);
LfpNf = length(Lfpfilt);
%% Calculate thresholds
dispstat('','init');
dispstat(sprintf('Processing thresholds . . .'),'keepthis','n');
for iChn = 1:nChn
    sp{iChn} = zeros(ntimes * double(T) * 20, FS * 1.6 / 1e3 + 1 + 1);
end
artchk = zeros(1,nChn);
munoise = cell(1,nChn);
nRun = 0;
tRun = ceil(nSam / FS / 10) + 1; % Number of 10 second chunks of time in the data
while sum(artchk) < nChn
    v = fread(v_fid, [nChn, FS*20], 'int16') * 0.195;    
    nRun = nRun + 1;
    dispstat(sprintf('Progress %03.2f%%',(100*(nRun/tRun))),'timestamp');
    if ~isempty(v)
        % Checks for artifact
        for iChn = 1:nChn
            if isempty(thresh{iChn})
                munoise{iChn} = [];
                mu = conv(v(iChn,:),Mufilt);
                if max(abs(mu(MuNf+1:end-MuNf))) < ARTMAX
                    artchk(iChn) = 1;
                    sd = median(abs(mu(MuNf+1:end-MuNf)))./0.6745;
                    thresh{iChn} = threshfac*sd;
                else
                    munoise{iChn} = [munoise{iChn} mu(MuNf+1:end-MuNf)];
                end
            end
        end
    else
        for iChn = 1:nChn
            if isempty(thresh{iChn}) % If threshold is still 0 - noisy channel
                artchk(iChn) = 1;
                sd = median(abs(munoise{iChn}))./0.6745;
                thresh{iChn} = threshfac*sd;
            end
        end
    end
end
disp(['Total recording time: ' num2str(nSam / FS) ' seconds.']);
disp(['Time analysed per loop: ' num2str(T) ' seconds.']);
fseek(v_fid,0,'bof'); % Returns the data pointer to the beginning of the file
%% Loop through the data
if (justMu)
    mu_fid = fopen([dName '.muBR.dat'],'W');
    dispstat('','init');
    dispstat(sprintf('Processing MU . . .'),'keepthis','n');
    while (chk && N < Ninput)
        N = N + 1;        
        dispstat(sprintf('Progress %03.2f%%',100*((N-1)/ntimes)),'timestamp');        
        data = fread(v_fid, [nChn, FS*T], 'int16') * 0.195;
        if (size(data,2))
            Ndata = size(data,2);
            MuOut = cell(1,nChn);
            if (par)
                parfor iChn = 1:nChn
                    flip_data = fliplr(data(iChn,:));
                    tmp = conv(flip_data,Mufilt);
                    mu = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
                    MuOut{iChn} = mu;
                end
            else
                for iChn = 1:nChn
                    flip_data = fliplr(data(iChn,:));
                    tmp = conv(flip_data,Mufilt);
                    mu = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
                    MuOut{iChn} = mu;
                end
            end
            mu2 = zeros(nChn,Ndata);
            for iChn = 1:nChn
                MuOut{iChn} = MuOut{iChn} - mean(MuOut{iChn});
                mu2(iChn,:) = MuOut{iChn};
            end
            fwrite(mu_fid,SCALEFACTOR*mu2,'short');
            time = time + T*1e3;
        end
        if (size(data,2) < FS * T)
            chk = 0;
        end
        dispstat(sprintf('Progress %03.2f%%',100*((N)/ntimes)),'timestamp');
    end
    fclose(mu_fid);
    clear data mu mu2 tmp flip_data
end
%% Calculate LFP and SP
m_fid = fopen([dName '.muBR.dat'],'r');
lfp_fid = fopen([dName '.lfpBR.dat'],'W');
dispstat('','init');
dispstat(sprintf('Processing SPLFP . . .'),'keepthis','n');
chk = 1; N = 0; time = 0;
fseek(v_fid,0,'bof');
while (chk && N < Ninput)
    N = N + 1;
    dispstat(sprintf('Progress %03.2f%%',100*((N-1)/ntimes)),'timestamp');
    mu = fread(m_fid,[nChn, FS*T],'short') ./ 10;        
    datalfp = fread(v_fid, [nChn, FS*T], 'int16') * 0.195;    
    if (size(datalfp,2))
        Ndata = size(datalfp,2);
        LfpOut = cell(1,nChn);
        if (par)
            parfor iChn = 1:nChn
                tmp = conv(datalfp(iChn,:),Lfpfilt);
                lfp = tmp(1,LfpNf/2:FS/1e3:Ndata+LfpNf/2-1);
                [Sp_tmp, Spktimes_tmp] = spikeextract(mu(iChn,:), thresh{iChn}, FS);                
                NSp_tmp = length(Spktimes_tmp);
                if NSp_tmp > 1
                    SpMat = [Spktimes_tmp+double(time) Sp_tmp];
                    NSp_tmp = size(SpMat,1);
                    sp{iChn}(NSp(iChn)+1:NSp(iChn)+NSp_tmp,:) = SpMat;
                    NSp(iChn) = NSp(iChn) + NSp_tmp;
                end
                LfpOut{iChn} = lfp;
            end
        else
            for iChn = 1:nChn
                tmp = conv(datalfp(iChn,:),Lfpfilt);
                lfp = tmp(1,LfpNf/2:FS/1e3:Ndata+LfpNf/2-1);
                [Sp_tmp, Spktimes_tmp] = spikeextract(mu(iChn,:), thresh{iChn}, FS);
                NSp_tmp = length(Spktimes_tmp);
                if NSp_tmp > 1
                    SpMat = [Spktimes_tmp+double(time) Sp_tmp];
                    NSp_tmp = size(SpMat,1);
                    sp{iChn}(NSp(iChn)+1:NSp(iChn)+NSp_tmp,:) = SpMat;
                    NSp(iChn) = NSp(iChn) + NSp_tmp;
                end
                LfpOut{iChn} = lfp;
            end
        end
        lfp2 = zeros(nChn,size(LfpOut{1},2));
        for iChn = 1:nChn
            lfp2(iChn,:) = LfpOut{iChn};
        end
        fwrite(lfp_fid,lfp2,'float');
        time = time + T*1e3;
    end
    if (size(datalfp,2) < FS * T)
        chk = 0;
    end
    dispstat(sprintf('Progress %03.2f%%',100*((N)/ntimes)),'timestamp');
end
fclose(lfp_fid);
fclose(v_fid);
%% Calculate a zero-condition spike template and r2t
if ~isempty(dir('*_exp_datafile_*'))
    trig = loadTrigBR(0); loadNREP;
    TrialParams = loadTrialParams;
    TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == 1));
    TrialParams = TrialParams(1:n_REP);
    trig = trig(TrialParams);
    nTrig = length(trig);
    spWN = [-300 300]; template = cell(1,nChn); r2t = cell(1,nChn); spt = cell(1,nChn);
    for c = 1:nChn
        for n = 1:nTrig
            spt{c} = [spt{c}; sp{c}(sp{c}(:,1) > trig(n)/30+spWN(1) & sp{c}(:,1) < trig(n)/30+spWN(2),2:end)];
        end
        template{c} = mean(cell2mat(spt(c)),1);
        r2t{c} = mean((cell2mat(spt(c)) - template{c}).^2);
    end
    %% Save the datafiles
    disp('Saving spikes');
    for iChn = 1:nChn
        sp{iChn} = sp{iChn}(1:NSp(iChn),:);
    end
    save([dName '.spBR.mat'],'sp','thresh','threshfac','template','r2t','-v7.3');
else
    disp('Saving spikes');
    for iChn = 1:nChn
        sp{iChn} = sp{iChn}(1:NSp(iChn),:);
    end
    save([dName '.spBR.mat'],'sp','thresh','threshfac','-v7.3');
end