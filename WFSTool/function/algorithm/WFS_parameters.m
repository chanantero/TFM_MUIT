function [parray_act,an,tn] = WFS_parameters(playDeviceID,vsx,vsy)

% Function description: According to different playDeviceID, the function
% will call different algorithm to calculate three parameters related to
% the WFS implementation. parray_act represents active loudspeakers related
% to different virtual sound source
% Examples: 
% parray_act = {[1,2,3,4,5,6,7,8],[],[]}
% an = 3x96 matrix 
%      each row shows different an coefficients for each file
%      each column shows an coefficients for each loudspeakers
% tn = 3x96 matrix
%      similar to an

% Read number of resources
source_no = size(vsx,2);

% Set number of loudspeakers according to device ID
if playDeviceID == 3
    speakers_no = 2;
else
    speakers_no = 96;
end

% Initialization of data
parray_act = cell(1,source_no);
an = zeros(source_no,speakers_no);
tn = zeros(source_no,speakers_no);

% For each sound source, select algorithm to calculate three return
% parameters based on device ID. 
for i = 1:source_no
    if playDeviceID == 3
        [array,an(i,:),tn(i,:)] = Headphone_test(vsx(i),vsy(i));
    else
        [array,an(i,:),tn(i,:)] = WFS_DrivingSignals(vsx(i),vsy(i));
    end
    parray_act(i) = mat2cell(array,size(array,1),size(array,2));
end

end

