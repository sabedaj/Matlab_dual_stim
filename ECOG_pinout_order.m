% ECOG electrode order working from the bottom RH corner up each column
% e.g. columns 1-8 then 16-23 in the diagram

ECOG_order=[1 15 2 14 3 13 4 12 5 11 6 10 7 9 8,   16 29 17 28 18 27 19 26 20 25 21 24 22 23,   30 43 31 42 32 41 33 40 34 39 35 38 36 37,   44 57 45 56 46 55 47 54 48 53 49 52 50 51,   58 71 59 70 60 69 61 68 62 67 63 66 64 65,  ...
    72 85 73 84 74 83 75 82 76 81 77 80 78 79,    86 99 87 98 88 97 89 96 90 95 91 94 92 93,    100 113 101 112 102 111 103 110 104 109 105 108 106 107,    114 128 115 127 116 126 117 125 118 124 119 123 120 122 121]; % order of each column
nnconnector=[63 50 47 34 31 18 15 2 49 64 33 48 17 32 1 16 14 3 30 19 46 35 62 51 4 13 20 29 36 45 52 61,    8 9 24 25 40 41 56 57 10 7 26 23 42 39 58 55 53 60 37 44 21 28 5 12 59 54 43 38 27 22 11 6,...
    128 113 112 97 96 81 80 65 114 127 98 111 82 95 66 79 77 68 93 84 109 100 125 116 67 78 83 94 99 110 115 126,   69 76 85 92 101 108 117 124 75 70 91 86 107 102 123 118 120 121 104 105 88 89 72 73 122 119 106 103 90 87 74 71]; %order on headstage connector according to intan map
order=zeros(128,1);
for i=1:128
    order(i)=find(ECOG_order(i)==nnconnector); %indicates which row of the .dat file belongs in each position e.g. data in row 15 of the intan data belongs to electrode 1 which is column 1 row 1 on the ECOG diagram
    % another example: row 28 of intan data corresponds to electrode 17 which is column 2 row 2 on the ECOG diagram
end

%% stimulating probe order from shank tip to base
NN_flexible=[9,8,10,7,11,6,13,4,14,3,15,2,24,25,23,26,22,27,20,29,19,30,18,31,17,32,21,28,16,1,12,5]; %already unscrambled
% Examples:
% data in row 9 corresponds with electrode 1 at shank tip
%data in row 5 corresponds to electrode 32 at shank base
