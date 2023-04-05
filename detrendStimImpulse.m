function detrendStimImpulse(filepath,dName,T,par)
%% Initial Setup
name = filepath;
name = strsplit(name, filesep);
name = name{end};
name = name(1:end-14);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;
if nChn>32
    E_Mapnumber=1;
else
    E_Mapnumber=0;
end
%% Blank stimulation artefact
trig = loadTrig(0);
theseTrig = trig;
if isempty(dir([name '_DT.mu.dat']))
    if ~isempty(theseTrig)
        info = dir([filepath filesep dName '.dat']);
        info = info.bytes/2;
        nL = (ceil(info / (nChn*FS*double(T)))+1);
        vFID = fopen([filepath filesep name '.mu_sab.dat'],'r');
        vdnFID = fopen([filepath filesep name '_DT.mu.dat'],'W');
        N = 1;
        CHK = split(filepath,'_');
        blank = 1;
        for i = 1:size(CHK,1)
            if strcmp(CHK{i},'CSD') || strcmp(dName,'analogin') || isempty(dir('*exp_datafile*.mat'))
                blank = 0;
            end
        end
        dispstat('','init');
        if (blank)
            dispstat(sprintf('Detrending artefact . . .'),'keepthis','n');
        else
            dispstat(sprintf('Skipping detrending . . .'),'keepthis','n');
        end
        BREAK = 1;
        while (N) && (BREAK)
            
            v = fread(vFID,[nChn, (FS * T)],'short')./10;
            
            if ~size(v,2)
                BREAK = 0;
            end
            dispstat(sprintf('Progress %03.2f%%',100*((N-1)/nL)),'timestamp');
            % Append saved data to the start of the next loop
            if (N ~= 1)
                v = [hv v]; %#ok<*AGROW>
            end
            if (blank)
                if (par)
                    parfor iChn = 1:nChn
                        v(iChn,:) = simpleDetrend(v(iChn,:),N,T,theseTrig);
                    end
                else
                    d = Depth(E_Mapnumber);
                    for iChn = 1:nChn
                        if d(iChn) == 13
                            pause(0.01);
                        end
                        v(iChn,:) = simpleDetrend(v(iChn,:),N,T,theseTrig);
                    end
                end
            end
            % Save a little data at the start of each loop
            hv = v(:,end-FS+1:end);
            % Don't save the last second, we'll deal with that next time
            if (BREAK)
                v = v(:,1:end-FS);
            end
            fwrite(vdnFID,v.*10,'short');
            dispstat(sprintf('Progress %03.2f%%',100*((N)/nL)),'timestamp');
            N = N + 1;
        end
    end
    %% Clean up
    fclose(vFID);
    fclose(vdnFID);
end


end