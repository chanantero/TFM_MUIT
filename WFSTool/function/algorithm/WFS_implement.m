function data = WFS_implement(y,parray_act,an,tn,pageSize,buffer)

% Function description: used to calculate the data in each page based on
% different attenuation and time delay. To calculate delay influence, a
% buffer is used to store last page info. If first page, buffer is filled
% with all zeros

nact = size(parray_act,2);
data = zeros(pageSize,nact);

% May changed depend on different value unit of tn. Finally rn represents
% how many samples are delayed
rn = round(tn);

% Calculate operated data for each active loudspeakers
for i = 1:nact
    % If no delay
    if rn(i) == 0
        data(:,i) = an(i)*y;
    else
        data(:,i) = [an(parray_act(i))*buffer(pageSize-rn(parray_act(i))+1:pageSize,1);...
            an(parray_act(i)) * y(1:pageSize-rn(parray_act(i)),1)];
    end
end

end

