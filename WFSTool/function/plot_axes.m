function [ ] = plot_axes(parray_act,handles,vsx,vsy)
global isTable isActive;
% Function description: plot the data about active loudspeakers array and
% sound source position on the axes

% Get the coordinates of loudspeakers
alt = generate_array();

% Clean the axes
cla(handles.axes1);

% Plot the loudspeakers first with blue circle
plot(handles.axes1,alt(1,:),alt(2,:),'ob');
xlim(handles.axes1,[-2,6]);
ylim(handles.axes1,[-2,8]);
grid(handles.axes1,'on');

% For each virtual source, plot it on the axes using green circle, while
% its active loudspeakers array in red circle
if isTable > 0
    for i = 1:size(parray_act,2)
        array = cell2mat(parray_act(i));
        column = size(array,2);
        ac_speakers = alt(:,array(1:column));
        hold(handles.axes1,'on');
        if isActive(i) == 1
            plot(handles.axes1,ac_speakers(1,:),ac_speakers(2,:),'or',vsx(i),vsy(i),'og');
        end
    end
end

end

