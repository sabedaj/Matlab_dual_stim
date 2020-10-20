%% single Channel Processing


D = dir; % D is a struct ... first elements are '.' and '..' used for navigation.
folder = fileparts(which('MASSIVE_ANALYSIS_LOOP')); % Determines filepath to folder containing your .m file even if your is not stored in same dir as long as path is added
addpath(genpath(folder)); % adds folder and all subfolders to path for all functions
fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'
cd([fileptochange filesep 'Data' filesep 'Rat_006' filesep 'Pen2_laminar_006_200925_143840'])%change working directory to where data is stored - currently manually input
D_data=dir;
fprintf([D_data(1).folder '\n'])
MASSIVE_ANALYSIS_LOOP(D_data(1).folder); %Perform data analysis and save files in appropriate folders
