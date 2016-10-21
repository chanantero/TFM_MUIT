function [endPoint,music] = read_music_file(namelist)
global signal_list pageSize;
% Function description: read data related to every music file in the
% namelist. endPoint will be the sample data size of each file. fileChan 
% records how many channels each file has, fileCount stores the number of 
% how many file in the namelist. 
% Example: 
% endPoint = [416500,184300,279300]
% fileChan = [2,2,2]
% fileCount = 3

% Initialization of data
fileCount = size(namelist,2);
endPoint = zeros(1,fileCount);
music = cell(1,fileCount);

% Read every file in the namelist
for i = 1:fileCount
    if isempty(strfind(cell2mat(namelist(i)),'wav')) == 1
        c = strsplit(cell2mat(namelist(i)),'signal');
        m = cell2mat(signal_list(str2num(cell2mat(c(2)))));
    else
        m = audioread(cell2mat(namelist(i)));
    end
    fileChan = size(m,2);
    if  fileChan > 1
        m = sum(m,2) / fileChan;
    end
    m = [zeros(pageSize,1);m];
    endPoint(i) = size(m,1);
    music(i) = mat2cell(m,size(m,1),size(m,2));
end

end

