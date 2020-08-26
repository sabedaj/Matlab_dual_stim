function data = simpleBlank(data,N,T,trig,mode)
TrialParams = loadTrialParams;
%%
TrialParams = cell2mat(TrialParams(:,2))';
nChn = size(data,1); FS = 30000;
if mode == 1 % for HPF
    BIN_b = -30; % Amount to blank in samples
    missed_trials = TrialParams(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
    missed_trig = trig(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
    missed_trials((missed_trig==-500))=[];
    missed_trig((missed_trig==-500))=[];
    
    trials = TrialParams(trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000);
    trig = trig(trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000);
    if ~isempty(missed_trig)
        trials = [missed_trials trials];
        trig = [missed_trig trig];
    end
    if (N == 1)
        trig = trig - ((N-1)*T*FS);
    else
        trig = trig - ((N-1)*T*FS) + FS;
    end
    for t = 1:length(trig) %used to go from 1
        % Don't blank zero trials
%         if trials(t) == 1
%             continue;
%         end


        if trials(t) == 5
            pause(0.01);
        end
        thisTrig = trig(t);
        for c = 1:nChn
            % Find the appropriate range
            ra = [1 46];
            range_max = 150;
            while range(data(c,thisTrig+ra(1):thisTrig+ra(2))) > range_max
                ra = ra + 1;
                if ra(1) > 180
                    ra(1) = 180;
                    break;
                end
            end
            ra(1) = ra(1) + 10;
            if ra(1) > 180
                ra(1) = 180;
            end
            % Shift the entire data array to minimize artefact. Convenient

            %original code
             % place to reintroduce artefact is +100 msec
            data(c,thisTrig+BIN_b:thisTrig+ra(1)) = interpolate(data(c,thisTrig+BIN_b:thisTrig+ra(1)),1);
            shift = diff([data(c,thisTrig+BIN_b),data(c,thisTrig+ra(1))]);
            data(c,thisTrig+ra(1):thisTrig+3000) = data(c,thisTrig+ra(1):thisTrig+3000) - shift;
            data(c,thisTrig+BIN_b:thisTrig+ra(1)+1) = interpolate(data(c,thisTrig+BIN_b:thisTrig+ra(1)+1),1);
            
%data(c,thisTrig+BIN_b) = detrend(data(c,thisTrig+BIN_b),1);
        

        end
    end
elseif mode == 2 % for LPF
    bBIN = [5, 250*30]; % Amount to blank in samples
    % Check for missed trigger lines
    missed_trig = trig(trig > ((N-1)*T*1000)-550 & trig < ((N-1)*T*1000)+550);
    trig = trig(trig >= ((N-1)*T*1000)+550 & trig <= (N*T*1000)-550);
    if ~isempty(missed_trig)
        trig = [missed_trig trig];
    end
    if (N == 1)
        trig = trig - ((N-1)*T*1000);
    else
        trig = trig - ((N-1)*T*1000) + 1000;
    end
    trig = cast(trig,'int64');
    for t = 1:length(trig)
        for c = 1:size(data,1)
            %distant_channel = 13 + c; % Get as far away as possible
            %if distant_channel > 26
            %    distant_channel = distant_channel - 26;
            %end
            %if distant_channel == 2 || distant_channel == 9
            %    distant_channel = distant_channel + 1; % Deal with dead channels
            %end
            trend1 = interpolate(data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)),1);
            wave = (data(c,trig(t)*30 + bBIN(2):trig(t)*30 + bBIN(2) + diff(bBIN)));
            trend2 = interpolate(wave,1);
            wave = wave - trend2 + trend1;
            data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)) = wave;
            %data(c,trig(t)*30 - 5:trig(t)*30 + 5) = interpolate(data(c,trig(t)*30 - 5:trig(t)*30 + 5),1);
        end
    end
end
end