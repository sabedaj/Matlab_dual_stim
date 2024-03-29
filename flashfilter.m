
%%for analysing data on massive and saving denoised files

% Parameters to alter
Startpoint_analyse=0; %set to 0 for no input
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input

filepath = pwd;
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;
fs=FS;
dName='amplifier';
mem_check=dir('amplifier.dat');
T = mem_check.bytes ./ (2 * nChn * fs);
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * fs);
if t_len < (Overall_time_to_analyse+Startpoint_analyse)
    error('Time to analyse exceeds the data recording period')
end
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256;
elseif Overall_time_to_analyse~=0
    T=Overall_time_to_analyse;
end


%% 3. LFP
% Loop through the data
chk = 1; N = 0; Ninput = 1e6; time = 0; SCALEFACTOR = 10;

E_Mapnumber=1;
par=0;
E_MAP = Depth(E_Mapnumber);

filepath = pwd;
[~,name,~] = fileparts(filepath);
name = name(1:end-14);
[~,Lfpfilt] = generate_Filters;
LfpNf = length(Lfpfilt);
dName='amplifier';
fileinfo = dir([filepath filesep dName '.dat']);
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);
if isempty(dir('*.lfp_sab.dat'))
    lfp_fid = fopen([name '.lfp_sab.dat'],'W');
     v_fid = fopen([dName '.dat'], 'r');
    dispstat('','init');
    dispstat(sprintf('Processing LFP . . .'),'keepthis','n');
    while (chk && N < Ninput)
        N = N + 1;        
        dispstat(sprintf('Progress %03.2f%%',100*((N-1)/ntimes)),'timestamp');
        if strcmp(dName,'analogin')
            data = fread(v_fid, [nChn, FS*T], 'uint16');
            data = (data - 32768) .* 0.0003125;
        elseif strcmp(dName,'amplifier')
            data = fread(v_fid, [nChn, FS*T], 'int16')* 0.195;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if (size(data,2))
            Ndata = size(data,2);
            LFPOut = cell(1,nChn);
            if (par)
                parfor iChn = 1:nChn
                    flip_data = fliplr(data(iChn,:));
                    tmp = conv(flip_data,Lfpfilt);
                    lfp = fliplr(tmp(1,LfpNf/2:Ndata+LfpNf/2-1));
                    LFPOut{iChn} = lfp;
                    %LFPOut{iChn} = data(iChn,:);
                end
            else
                for iChn = 1:nChn
                    flip_data = fliplr(data(iChn,:));
                    tmp = conv(flip_data,Lfpfilt);
                    lfp = fliplr(tmp(1,LfpNf/2:Ndata+LfpNf/2-1));
                    LFPOut{iChn} = lfp;
                    %LFPOut{iChn} = data(iChn,:);
                end
            end
            lfp2 = zeros(nChn,Ndata);
            for iChn = 1:nChn
                lfp2(iChn,:) = LFPOut{iChn};
            end
            fwrite(lfp_fid,SCALEFACTOR*lfp2,'short');
            time = time + T*1e3;
            dispstat(sprintf('Progress %03.2f%%',100*((N)/ntimes)),'timestamp');
        end
        if (size(data,2) < FS * T)
            chk = 0;
        end
    end
    fclose('all');
    clear data lfp lfp2 tmp flip_data
end



%% 4. meanLFP

if isempty(dir('meanlfp.mat'))
fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);
stimDig = flip(find(digital_in == 0)); % Fix for finding the falling line instead of the rising line

visDig = flip(find(digital_in == 2)); % Fix for finding the falling line instead of the rising line
dt = diff(stimDig);
kill = dt == -1;
stimDig(kill) = [];
dt = diff(visDig);
kill = dt == -1;
visDig(kill) = [];
nStamps = max([length(stimDig); length(visDig)]);
time_stamps = nan(1,nStamps);
time_stamps(1,1:length(stimDig)) = flip(stimDig);
time_stamps(1)=[];
trig=time_stamps;

%load([name '.lfp_sab.mat'])
for chsp=1:nChn
    %check(chsp,:)=['T', num2str(chsp)];
    
LFPstruct{chsp}=[];
end
n_REP=nStamps-1;

        nTrig = n_REP;
        fprintf('Percentage trig processed: \n')
        for indT=1:nTrig
                fileID=fopen([name '.lfp_sab.dat'],'r');
                shortbytes=2;
                offset=trig(indT)*nChn*shortbytes-0.25*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
                %ftell(fileID)
                fseek(fileID,offset,'bof');
                %ftell(fileID)
                vlfp = fread(fileID,[nChn, (0.25*FS+0.5*FS)],'short')./10; %plots 500ms from trigger and 250ms brefore
                fclose(fileID);
                %vlfp = vlfp(:,1:30:end);
                fprintf('%0.8f%% \n', (indT*100/nTrig))
           for chsp=1:nChn
                %check=['T', num2str(chsp)];
                LFPstruct{chsp}(indT,:) =  vlfp(chsp,:);
            end
        end
        meanlfpstruct=zeros(nChn,size(vlfp,2));
        for chsp=1:1:nChn
            check=['T', num2str(chsp)];
            meanlfpstruct(chsp,:)=mean(LFPstruct{chsp},1);
        end
        
        
save('meanlfp.mat','meanlfpstruct')
%%
end


% %%
% if isempty(dir('meanresp.mat'))
% fileinfo = dir('digitalin.dat');
% nSam = fileinfo.bytes/2;
% digin_fid = fopen('digitalin.dat','r');
% digital_in = fread(digin_fid, nSam, 'uint16');
% fclose(digin_fid);
% stimDig = flip(find(digital_in == 0)); % Fix for finding the falling line instead of the rising line
% 
% visDig = flip(find(digital_in == 2)); % Fix for finding the falling line instead of the rising line
% dt = diff(stimDig);
% kill = dt == -1;
% stimDig(kill) = [];
% dt = diff(visDig);
% kill = dt == -1;
% visDig(kill) = [];
% nStamps = max([length(stimDig); length(visDig)]);
% time_stamps = nan(1,nStamps);
% time_stamps(1,1:length(stimDig)) = flip(stimDig);
% time_stamps(1)=[];
% trig=time_stamps;
% 
% 
% %load([name '.lfp_sab.mat'])
% LFPstruct=[];
% n_REP=nStamps-1;
% 
%         nTrig = n_REP;
%         fprintf('Percentage trig processed: \n')
%         for indT=1:nTrig
%                 fileID=fopen('amplifier.dat','r');
%                 shortbytes=2;
%                 offset=trig(indT)*nChn*shortbytes-0.25*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
%                 %ftell(fileID)
%                 fseek(fileID,offset,'bof');
%                 %ftell(fileID)
%                 vlfp = fread(fileID,[nChn, (0.25*FS+1.5*FS)],'short')./10; %plots 500ms from trigger and 250ms brefore
%                 fclose(fileID);
%                 %vlfp = vlfp(:,1:30:end);
%                 fprintf('%0.8f%% \n', (indT*100/nTrig))
%            for chsp=1:1:nChn
%                 check=['T', num2str(chsp)];
%                 if isfield(LFPstruct,check)
%                     LFPstruct.(check) = [LFPstruct.(check); vlfp(chsp,:)];
%                 else
%                     LFPstruct.(check)= vlfp(chsp,:);
%                 end
%             end
%         end
%         meanrespstruct=zeros(nChn,size(vlfp,2));
%         for chsp=1:1:nChn
%             check=['T', num2str(chsp)];
%             meanrespstruct(chsp,:)=mean(LFPstruct.(check),1);
%         end
% save('meanresp.mat','meanrespstruct')
% end
% 
% 
% cd(folder)
% fprintf(['End of Analysis for: ' SubDir_Path newline])


