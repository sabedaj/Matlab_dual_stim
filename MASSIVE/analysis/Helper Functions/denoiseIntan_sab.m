function denoiseIntan_sab(filepath,dName,T,par,varargin)
dbstop if error
%% This function denoises an Intan datafile
%% Initial Setup
if size(varargin,2)==0
    Startpoint_analyse=0;
    Overall_time_to_analyse=0;
else
    Startpoint_analyse=cell2mat(varargin(1));
    Overall_time_to_analyse=cell2mat(varargin(2));
end
if strcmp(dName,'analogin')
    nChn = 1;
elseif strcmp(dName,'amplifier')
    fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
    fileinfo = dir([filepath, filesep, 'info.rhs']);
    if (datetime(fileinfo.date) < fourShank_cutoff)
        nChn=32;
        E_Mapnumber=0;
    else
        E_Mapnumber=loadMapNum;
        if E_Mapnumber>0
            nChn=64;
        else
            nChn=32;
        end
    end
    E_MAP = Depth;
end
shortbytes=2;
FS = 30000;
%% Deal with the digital lines
cleanTrig_sabquick;
%% Blank stimulation artefact
trig = loadTrig(0);
theseTrig = trig;
if isempty(dir([dName '_dn_sab.dat']))
    if ~isempty(theseTrig)
        info = dir([filepath filesep dName '.dat']);
        info = info.bytes/2;
        nL = (ceil(info / (nChn*FS*double(T)))+1);
        %lvFID = fopen([filepath '\' dName '.dat'],'r');
        vFID = fopen([filepath filesep dName '.dat'],'r');
        vdnFID = fopen([filepath filesep dName '_dn_sab.dat'],'W');
        %lbFID = fopen([filepath '\' dName '_lfpraw.dat'],'w');
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
            dispstat(sprintf('Blanking artefact . . .'),'keepthis','n');
        else
            dispstat(sprintf('Skipping blanking . . .'),'keepthis','n');
        end
        BREAK = 1;
        while (N) && (BREAK)
            if strcmp(dName,'analogin')
                v = fread(vFID,[nChn, (FS * T)],'uint16');
                v = (v - 32768) .* 0.0003125;
            elseif strcmp(dName,'amplifier')
                if Startpoint_analyse~=0
                    offset=Startpoint_analyse*FS*nChn*shortbytes;%offset from beginning of file
                    ftell(vFID)
                    fseek(vFID,offset,'bof');
                    ftell(vFID)
                end
                v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195;
            end
            if ~size(v,2)
                BREAK = 0;
            end
            %dispstat(sprintf('Progress %03.2f%%',100*((N-1)/nL)),'timestamp');
            % Append saved data to the start of the next loop
            if (N ~= 1)
                v = [hv v]; %#ok<*AGROW>
            end
            if (blank)
                if (par)
                    parfor iChn = 1:nChn
                        v(iChn,:) = simpleBlank(v(iChn,:),N,T,theseTrig,1);
                    end
                else
                    d = Depth(E_Mapnumber);
                    for iChn = 1:nChn
                        if d(iChn) == 13
                            pause(0.01);
                        end
                        v(iChn,:) = simpleBlank(v(iChn,:),N,T,theseTrig,1);
                    end
                end
            end
            % Save a little data at the start of each loop
            hv = v(:,end-FS+1:end);
            % Don't save the last second, we'll deal with that next time
            if (BREAK)
                v = v(:,1:end-FS);
            end
            if Overall_time_to_analyse~=0
                BREAK=0;
            end
            fwrite(vdnFID,v./0.195,'int16');
            dispstat(sprintf('Progress %03.2f%%',100*((N)/nL)),'timestamp');
            N = N + 1;
        end
    end
    %% Clean up
    fclose(vFID);
    fclose(vdnFID);
end
end