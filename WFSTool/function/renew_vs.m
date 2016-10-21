function [namelist,vsx,vsy,isActive] = renew_vs(handles)
% Function description: read music file name into namelist, virtual sound
% source position into vsx and vsy. If no_file = 0, means that there will
% be at least on file in the table; if no_file = -1, means that there is no
% file in the table at all; if no_file = 1, means that there are files in
% the table but not active. 
% Example: 
% namelist = {'test.wav','car.wav','dog.wav'}
% vsx = [0,0,0]
% vsy = [0,2,4]

global isTable;
isTable = 0;

% Cell array list stores data read from the uitable
list = get(handles.table,'data');
row = size(list,1);

% Data initialisation
namelist = {};
vsx = [];
vsy = [];
isActive = [];

% Empty table
if row == 0
    isTable = -1;
    return;
end

% Read and fill the data
for i = 1:row
    namelist(end+1) = list(i,1);
    vsx(end+1) = cell2mat(list(i,2));
    vsy(end+1) = cell2mat(list(i,3));
    isActive(end+1) = cell2mat(list(i,4));
    if cell2mat(list(i,4)) == 1
        isTable = isTable + 1;
    end
end
end
