function varargout = findprobesHD_GUI(varargin)
% findprobesHD_GUI MATLAB code for findprobesHD_GUI.fig
%      findprobesHD_GUI, by itself, creates a new findprobesHD_GUI or raises the existing
%      singleton*.
%
%      H = findprobesHD_GUI returns the handle to a new findprobesHD_GUI or the handle to
%      the existing singleton*.
%
%      findprobesHD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in findprobesHD_GUI.M with the given input arguments.
%
%      findprobesHD_GUI('Property','Value',...) creates a new findprobesHD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before findprobesHD_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to findprobesHD_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help findprobesHD_GUI

% Last Modified by GUIDE v2.5 27-Jun-2012 11:09:09
%
% AUTHOR:
% Bobby Kent 2012

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @findprobesHD_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @findprobesHD_GUI_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before findprobesHD_GUI is made visible.
function findprobesHD_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to findprobesHD_GUI (see VARARGIN)

% Browse for a .fasta file
[fastafile, filepath] = uigetfile({'*.fasta';'*.fa';'*.txt'},...
                                     'Choose a .fasta file', 'D:\0Raj Lab');
cd(filepath)
handles.filename = fastafile;
set(handles.outfilename_box,'String',fastafile(1:end-6));

% Choose default command line output for findprobesHD_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes findprobesHD_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = findprobesHD_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in repeatmask_box.
function repeatmask_box_Callback(hObject, eventdata, handles)
% hObject    handle to repeatmask_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of repeatmask_box


% --- Executes on button press in pseudogenemask_box.
function pseudogenemask_box_Callback(hObject, eventdata, handles)
% hObject    handle to pseudogenemask_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pseudogenemask_box


% --- Executes on button press in blastmask_box.
function blastmask_box_Callback(hObject, eventdata, handles)
% hObject    handle to blastmask_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blastmask_box


% --- Executes on button press in miRNAmask_box.
function miRNAmask_box_Callback(hObject, eventdata, handles)
% hObject    handle to miRNAmask_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of miRNAmask_box


% --- Executes on button press in humanfiltermask_box.
function humanfiltermask_box_Callback(hObject, eventdata, handles)
% hObject    handle to humanfiltermask_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of humanfiltermask_box


% --- Executes on button press in genomemask_box.
function genomemask_box_Callback(hObject, eventdata, handles)
% hObject    handle to genomemask_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of genomemask_box


% --- Executes on button press in GCrunmask_box.
function GCrunmask_box_Callback(hObject, eventdata, handles)
% hObject    handle to GCrunmask_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GCrunmask_box


% --- Executes on button press in GCmask_box.
function GCmask_box_Callback(hObject, eventdata, handles)
% hObject    handle to GCmask_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GCmask_box


% --- Executes on button press in gobutton.
function gobutton_Callback(hObject, eventdata, handles)
% hObject    handle to gobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Checkboxes
repeatmask = logical(get(handles.repeatmask_box,'Value'));
pseudogenemask = logical(get(handles.pseudogenemask_box,'Value'));
blastmask = logical(get(handles.blastmask_box,'Value'));
miRNAmask = logical(get(handles.miRNAmask_box,'Value'));
humanfiltermask = logical(get(handles.humanfiltermask_box,'Value'));
genomemask = logical(get(handles.genomemask_box,'Value'));
GCrunmask = logical(get(handles.GCrunmask_box,'Value'));
GCmask = logical(get(handles.GCmask_box,'Value'));

% Parameters
noligos = str2num(get(handles.noligos_box,'String'));
targetTM = str2num(get(handles.targetTM_box,'String'));
spacerlength = str2num(get(handles.spacerlength_box,'String'));
oligolength = str2num(get(handles.oligolength_box,'String'));
outfilename = get(handles.outfilename_box,'String');

n = get(handles.species_box,'Value');
speclist = get(handles.species_box,'String');
species = speclist{n};

close findprobesHD_GUI;

% Run findprobesHD with the input arguments
findprobesHD(handles.filename,noligos,... 
             'targetTM',targetTM,... 
             'spacerlength',spacerlength,... 
             'oligolength',oligolength,... 
             'outfilename',outfilename,...
             'species',species,...
             'repeatmask',repeatmask,...
             'pseudogenemask',pseudogenemask,...
             'blastmask',blastmask,...
             'miRNAmask',miRNAmask,...
             'humanfiltermask',humanfiltermask,...
             'genomemask',genomemask,...
             'GCrunmask',GCrunmask,...
             'GCmask',GCmask);
return;

function noligos_box_Callback(hObject, eventdata, handles)
% hObject    handle to noligos_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noligos_box as text
%        str2double(get(hObject,'String')) returns contents of noligos_box as a double


% --- Executes during object creation, after setting all properties.
function noligos_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noligos_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function targetTM_box_Callback(hObject, eventdata, handles)
% hObject    handle to targetTM_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetTM_box as text
%        str2double(get(hObject,'String')) returns contents of targetTM_box as a double


% --- Executes during object creation, after setting all properties.
function targetTM_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetTM_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spacerlength_box_Callback(hObject, eventdata, handles)
% hObject    handle to spacerlength_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spacerlength_box as text
%        str2double(get(hObject,'String')) returns contents of spacerlength_box as a double


% --- Executes during object creation, after setting all properties.
function spacerlength_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spacerlength_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function oligolength_box_Callback(hObject, eventdata, handles)
% hObject    handle to oligolength_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oligolength_box as text
%        str2double(get(hObject,'String')) returns contents of oligolength_box as a double


% --- Executes during object creation, after setting all properties.
function oligolength_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oligolength_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outfilename_box_Callback(hObject, eventdata, handles)
% hObject    handle to outfilename_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outfilename_box as text
%        str2double(get(hObject,'String')) returns contents of outfilename_box as a double


% --- Executes during object creation, after setting all properties.
function outfilename_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outfilename_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in species_box.
function species_box_Callback(hObject, eventdata, handles)
% hObject    handle to species_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns species_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from species_box


% --- Executes during object creation, after setting all properties.
function species_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to species_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
