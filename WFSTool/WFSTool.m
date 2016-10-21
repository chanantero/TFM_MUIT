function varargout = WFSTool(varargin)
%% Begin initialization code - DO NOT EDIT
% Input arguments:
% - gui Callback. Optional.

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WFSTool_OpeningFcn, ...
                   'gui_OutputFcn',  @WFSTool_OutputFcn, ...
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

% Initialize the file route
addpath('function');
addpath('function/algorithm');
addpath('function/check_variable');
addpath('examples');
addpath('playrec');
addpath('gui_element');
addpath('gui_element/Timers');
addpath('gui_element/Menu');
addpath('gui_element/Buttons');
% End initialization code - DO NOT EDIT

% --- Executes just before WFSTool is made visible.

% --- Outputs from this function are returned to the command line.
function varargout = WFSTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in menu_driver.
function menu_driver_Callback(hObject, eventdata, handles)
% hObject    handle to menu_driver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function menu_driver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_driver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in menu_buffersize.
function menu_buffersize_Callback(hObject, eventdata, handles)
% hObject    handle to menu_buffersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function menu_buffersize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_buffersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in menu_device.
function menu_device_Callback(hObject, eventdata, handles)
% hObject    handle to menu_device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function menu_device_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function file_Callback(hObject, eventdata, handles)

function Exit_Callback(h, eventdata, handles, varargin)

% --- Executes on mouse press over axes1 background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function wave_Callback(hObject, eventdata, handles)
% hObject    handle to wave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ed_pattern_Callback(hObject, eventdata, handles)
% hObject    handle to ed_pattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_pattern as text
%        str2double(get(hObject,'String')) returns contents of ed_pattern as a double


% --- Executes during object creation, after setting all properties.
function ed_pattern_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_pattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_angle_Callback(hObject, eventdata, handles)
% hObject    handle to ed_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_angle as text
%        str2double(get(hObject,'String')) returns contents of ed_angle as a double


% --- Executes during object creation, after setting all properties.
function ed_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_z_Callback(hObject, eventdata, handles)
% hObject    handle to ed_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_z as text
%        str2double(get(hObject,'String')) returns contents of ed_z as a double


% --- Executes during object creation, after setting all properties.
function ed_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_patternz_Callback(hObject, eventdata, handles)
% hObject    handle to ed_patternz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_patternz as text
%        str2double(get(hObject,'String')) returns contents of ed_patternz as a double


% --- Executes during object creation, after setting all properties.
function ed_patternz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_patternz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
