function varargout = wave_input(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_input_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_input_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function wave_input_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
uiwait(handles.wave_input_figure);

function varargout = wave_input_OutputFcn(hObject, eventdata, handles) 

function input_fs_Callback(hObject, eventdata, handles)

function input_fs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function input_a_Callback(hObject, eventdata, handles)

function input_a_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in menu_cs.
function menu_cs_Callback(hObject, eventdata, handles)

function menu_cs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function input_we_Callback(hObject, eventdata, handles)

function input_we_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function but_ok_Callback(hObject, eventdata, handles)
global signal_list sg_count;

signal = [];
if get(handles.select_sc,'value') == 1
    wave_fs = str2double(get(handles.input_fs,'string'));
    a = str2double(get(handles.input_a,'string'));
    % error information
    if isnan(a) || isnan(wave_fs)
        fprintf('Invalid Input of FS or Amplitude! \n');
        h = dialog('name','Error! ','position',[640,480,200,70]);  
        uicontrol('parent',h,'style','text','string','Invalid Input of FS or Amplitude! ','position',[15 10 180 50],'fontsize',10);  
        uicontrol('parent',h,'style','pushbutton','string','OK','position',[75 5 50 30],'callback','delete(gcbf)');
        return;  
    end
    wave_select = get(handles.menu_cs,'value');
    if wave_select == 1
        signal = a * sin((1:wave_fs)*2*pi*1000/wave_fs)';
    else
        signal = a * cos((1:wave_fs)*2*pi*1000/wave_fs)';
    end
end

if get(handles.select_custom,'value') == 1
    str = get(handles.input_we,'string');
    signal = evalin('base',str);
end

signal_list(end+1) = mat2cell(signal,size(signal,1),size(signal,2));
sg_count = sg_count + 1;
uiresume(handles.wave_input_figure);
delete(handles.wave_input_figure);

function but_plot_Callback(hObject, eventdata, handles)

if get(handles.select_sc,'value') == 1
    wave_fs = str2double(get(handles.input_fs,'string'));
    a = str2double(get(handles.input_a,'string'));
    wave_select = get(handles.menu_cs,'value');
    if wave_select == 1
        signal = a * sin((1:wave_fs)*2*pi*1000/wave_fs)';
    else
        signal = a * cos((1:wave_fs)*2*pi*1000/wave_fs)';
    end
end

if get(handles.select_custom,'value') == 1
    str = get(handles.input_we,'string');
    signal = evalin('base',str);
end

x = 1:size(signal,1);
plot(handles.wave_axes,x,signal);
xlim(handles.wave_axes,[0,size(signal,1)]);
ylim(handles.wave_axes,[-size(signal,2),size(signal,2)]);
grid(handles.wave_axes,'on');
