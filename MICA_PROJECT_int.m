function varargout = MICA_PROJECT_int(varargin)
% MICA_PROJECT_INT MATLAB code for MICA_PROJECT_int.fig
%      MICA_PROJECT_INT, by itself, creates a new MICA_PROJECT_INT or raises the existing
%      singleton*.
%
%      H = MICA_PROJECT_INT returns the handle to a new MICA_PROJECT_INT or the handle to
%      the existing singleton*.
%
%      MICA_PROJECT_INT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MICA_PROJECT_INT.M with the given input arguments.
%
%      MICA_PROJECT_INT('Property','Value',...) creates a new MICA_PROJECT_INT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MICA_PROJECT_int_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MICA_PROJECT_int_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MICA_PROJECT_int

% Last Modified by GUIDE v2.5 18-May-2018 12:20:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MICA_PROJECT_int_OpeningFcn, ...
                   'gui_OutputFcn',  @MICA_PROJECT_int_OutputFcn, ...
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


% --- Executes just before MICA_PROJECT_int is made visible.
function MICA_PROJECT_int_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MICA_PROJECT_int (see VARARGIN)

% Choose default command line output for MICA_PROJECT_int
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MICA_PROJECT_int wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MICA_PROJECT_int_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pqrst_box.
function pqrst_box_Callback(hObject, eventdata, handles)
% hObject    handle to pqrst_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pqrst_box
