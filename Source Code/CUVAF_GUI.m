function varargout = CUVAF_GUI(varargin)
% CUVAF_GUI MATLAB code for CUVAF_GUI.fig
%      CUVAF_GUI, by itself, creates a new CUVAF_GUI or raises the existing
%      singleton*.
%
%      H = CUVAF_GUI returns the handle to a new CUVAF_GUI or the handle to
%      the existing singleton*.
%
%      CUVAF_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CUVAF_GUI.M with the given input arguments.
%
%      CUVAF_GUI('Property','Value',...) creates a new CUVAF_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CUVAF_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CUVAF_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CUVAF_GUI

% Last Modified by GUIDE v2.5 01-Dec-2015 13:36:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CUVAF_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CUVAF_GUI_OutputFcn, ...
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


% --- Executes just before CUVAF_GUI is made visible.
function CUVAF_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CUVAF_GUI (see VARARGIN)

% Choose default command line output for CUVAF_GUI
handles.output = hObject;

%% ###### START UP CODE #####
clc;
clear global;

%% Initialize Variables
handles.CUVAF_version = 'v2';
set(handles.UVAF_area_text, 'String','0');
set(handles.ROI_checkbox, 'Value', 0);
set(handles.comments_edit, 'String', '');
set(handles.axes1, 'Visible', 'off');
set(handles.axes2, 'Visible', 'off');
set(handles.Image_display,'SelectedObject',[]);
set(handles.ChooseAxes_uipanel,'SelectedObject',[]);
handles.selected_points = [];                           % coordinates of points to remove regions

%% Calibration factor calculation
im_pix_length = 3008;       % pixels
im_mm_length = 24;         % mm
resize_factor = 0.5;

handles.calibration = [num2str(im_mm_length) 'mm2 = ' num2str(im_pix_length) 'pix'];

calibfactor = im_mm_length/im_pix_length;   % mm/pixel
calibfactor = calibfactor/resize_factor;

handles.calibfactor = calibfactor;
handles.resize_factor = resize_factor;

%% Set Parameters for analysis
win_ratio = 0.15;
handles.win_ratio = win_ratio;
area_filt = 0.03;       % regions < 0.03 mm^2 will be excluded (deemed insignificant)
handles.area_filt = area_filt;

%% Initialisation of POI Libs
% Add Java POI Libs to matlab javapath
javaaddpath('poi_library/poi-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('poi_library/xmlbeans-2.3.0.jar');
javaaddpath('poi_library/dom4j-1.6.1.jar');
javaaddpath('poi_library/stax-api-1.0.1.jar');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CUVAF_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CUVAF_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Maximize GUI figure
jf = get(handle(gcf), 'JavaFrame');
jf.setMaximized(true);

% --- Executes on button press in browse_pushbutton.
function browse_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to browse_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.save_pushbutton,'Enable','on');

%% reset variables - allow for multiple runs without closing the GUI
set(handles.UVAF_area_text, 'String', '0');
set(handles.axes1_text, 'String', '');
set(handles.axes2_text, 'String', '');
set(handles.threshold_edit_text, 'String', 0);
set(handles.Image_display,'SelectedObject',[]);
set(handles.ChooseAxes_uipanel,'SelectedObject',[]);
set(handles.folder_listbox,'string','');
set(handles.ROI_checkbox, 'Value', 0);
set(handles.subjective_checkbox, 'Value', 0);
set(handles.weak_checkbox, 'Value', 0);
cla(handles.axes1);
cla(handles.axes2);
set(handles.measure_pushbutton,'Enable','on');

%% Get pathname of folder/images
directory = get(handles.directory_edit, 'String');
folder_name = uigetdir(directory);
if folder_name == 0  %user cancelled/closed
    return;
end
handles.ImageFolder = folder_name;

% find out if it is a max or windows pc
handles.delimiter = filesep;

directory = [];
path_parts = regexp(folder_name, handles.delimiter, 'split');
num_parts = length(path_parts);
for i = 1: num_parts-1
	directory = [directory path_parts{i} handles.delimiter];
end
set(handles.directory_edit, 'String', directory);

set(handles.foldername_edit_text, 'String', path_parts{num_parts});
handles.excel_filename = [path_parts{num_parts} '-CUVAF_measurements.xls'];      % name of excel file

% Update handles structure
guidata(hObject, handles);

% display image name in folder in listbox
set(handles.folder_listbox,'Value',1);
handles = LoadImageList(handles);

%% Create/open excel spreadsheet
filename = handles.excel_filename;
path = [folder_name handles.delimiter filename];
    
% Sheet1 headings
headings1 = {'ID' 'Resize Factor' 'Overlap Thresh' 'Area (mm2)', 'Notes' 'Date Analyzed', 'ROI Date', 'Subjective?', 'Initials', 'Comments'};  
sub1 = {'OD nasal'; 'OD temporal'; 'OS nasal'; 'OS temporal'};
sub2 = {'PARAMETERS'; 'CUVAF Version'; 'Calibration'; 'Window Width Ratio'; 'Area Filter Size'};
sub3 = {'IMAGE STATS', 'Mean', 'Standard Dev', 'Min', 'Max', 'Mean G', 'SD G', 'Mean B', 'SD B'};
xlwrite(path, headings1, 'Results', 'A1');
xlwrite(path, sub1, 'Results', 'A2');   
xlwrite(path, sub2, 'Results', 'A10');
xlwrite(path, sub3, 'Results', 'M1');

sub2_values = {handles.CUVAF_version; handles.calibration};
xlwrite(path, sub2_values, 'Results', 'B11');
sub2_values2 = [handles.win_ratio; handles.area_filt];
xlwrite(path, sub2_values2, 'Results', 'B13');

% Sheet2 headings
xlwrite(path, {'OD nasal'}, 'ROI Coordinates', 'A1');
xlwrite(path, {'OD temporal'}, 'ROI Coordinates', 'D1');
xlwrite(path, {'OS nasal'}, 'ROI Coordinates', 'G1');
xlwrite(path, {'OS temporal'}, 'ROI Coordinates', 'J1');

% Sheet3 headings
headings3 = {'OD nasal' '' 'OD temporal' ''  'OS nasal' '' 'OS temporal'};
xlwrite(path, headings3, 'Removed Regions', 'A1');    

% Update handles structure
guidata(hObject, handles);
    
%=====================================================================
% http://au.mathworks.com/matlabcentral/answers/171611-load-folder-s-images-names-to-a-matlab-listbox
% --- Load up the listbox with tif files in folder handles.handles.ImageFolder
function handles=LoadImageList(handles)        
    ListOfImageNames = {};
    if ~isempty(handles.ImageFolder) 
		ImageFiles = dir([handles.ImageFolder handles.delimiter '*.*']);
        for Index = 1:length(ImageFiles)
            baseFileName = ImageFiles(Index).name;
            [folder, name, extension] = fileparts(baseFileName);
            extension = upper(extension);
            switch lower(extension)
            case {'.png', '.bmp', '.jpg', '.tif', '.avi', '.xls'}
                ListOfImageNames = [ListOfImageNames baseFileName];
            otherwise
            end
        end
        set(handles.folder_listbox,'string',ListOfImageNames);
        return
	else
		msgbox('No folder specified as input for function LoadImageList.');
		return;
    end   
%=====================================================================

% --- Executes on selection change in folder_listbox.
function folder_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to folder_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns folder_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from folder_listbox

folder_name = handles.ImageFolder;
imagename = get(handles.folder_listbox, 'String');
index = get(handles.folder_listbox,'Value');
pathname = [folder_name handles.delimiter imagename{index}];

im = imread(pathname);
imshow(im, 'Parent', handles.axes1);

set(handles.axes1_text, 'String', imagename{index})
handles.imagename = get(handles.axes1_text, 'String');

handles.pathname = pathname;
handles.im = im;

set(handles.ROI_checkbox, 'Value', 0);
set(handles.subjective_checkbox, 'Value', 0);
set(handles.weak_checkbox, 'Value', 0);
set(handles.comments_edit, 'String', '');
set(handles.UVAF_area_text, 'String', '0');
set(handles.measure_pushbutton,'Enable','on');

[m,n,o] = size(im);
imG = im(:,:,2);
imB = im(:,:,3);
meanG = mean2(double(imG));

gamma = (1/50)*meanG + 0.6;              % min gamma = y-intercept
if gamma > 1.5
    gamma = 1.5;
end

tol = stretchlim(imB, [0.02 0.99]);
im2 = imadjust(im, tol, [0;1], gamma);
imG2 = im2(:,:,2);
imB2 = im2(:,:,3);

disp_im = cat(3, zeros(m,n), imG2, imB2);
imshow(disp_im, 'Parent', handles.axes2);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function folder_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to folder_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ROI_checkbox.
function ROI_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ROI_checkbox
checkbox = get(handles.ROI_checkbox, 'Value');

if checkbox == 1
    set(handles.ROI_checkbox, 'Value', 1);
   
    % check if ROI coordinates exists
    foldername = [handles.ImageFolder handles.delimiter];                    %ie. C:\Users\emily huynh\Desktop\160814\10331\
    pathname = handles.pathname;
    path_parts = regexp(pathname, handles.delimiter, 'split');
    num_parts = length(path_parts);

    imagename = path_parts{num_parts};
    path_parts = regexp(imagename, '\.', 'split');
    imagename = path_parts{1};          %ie. OS temporal UV_AF
    
    filename = handles.excel_filename;
    path = [foldername filename];
    if exist(path, 'file')
        if isempty(regexpi(imagename, 'OD NASAL')) == 0
            ROI_coords = xlsread(path , 'ROI Coordinates', 'A2:B60');                    
            resize_factor = xlsread(path, 'Results', 'B2');
        elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
            ROI_coords = xlsread(path , 'ROI Coordinates', 'D2:E60');
            resize_factor = xlsread(path, 'Results', 'B3');
        elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
            ROI_coords = xlsread(path , 'ROI Coordinates', 'G2:H60');
            resize_factor = xlsread(path, 'Results', 'B4');
        elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
            ROI_coords = xlsread(path , 'ROI Coordinates', 'J2:K60');
            resize_factor = xlsread(path, 'Results', 'B5');
        else
            msgbox('No ROI coordinates stored for this image')
            return
        end
    else
        msgbox('No data for this image')
        return
    end  
    if isempty(ROI_coords) == 1 
        msgbox('No ROI coordinates stored for this image')
        return
    end
    im = handles.im;
    im = imresize(im, resize_factor);
    roi_border = roipoly(im, ROI_coords(:,1), ROI_coords(:,2));
    roi_border = bwperim(roi_border);
    im(roi_border) = 255;
    imshow(im, 'Parent', handles.axes1);

elseif checkbox == 0
    set(handles.ROI_checkbox, 'Value', 0);
    im = handles.im;
    imshow(im, 'Parent', handles.axes1);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in measure_pushbutton.
function measure_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to measure_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.measure_pushbutton,'Enable','off');
pathname = handles.pathname;
set(handles.save_pushbutton,'Enable','on');

[UV_area_mm, im, disp_im, image, UV_roi_mask, smoothed_flag edge_UV, overlap_thresh, UV_region] = Analyze_CUVAF(handles); %Analyze_CUVAF(pathname, checkbox, excel, excel_results);

if isnumeric(UV_area_mm) == 1
    set(handles.UVAF_area_text, 'String', UV_area_mm)
    set(handles.threshold_edit_text, 'String', overlap_thresh);
    set(handles.threshold_slider, 'Value', overlap_thresh);

    handles.im_resized = im;                    % resized original image
    handles.red_free = disp_im;                 % red-free display image
    handles.disp_im = disp_im;
    handles.image = image;                      % grayscale image
    handles.UV_roi_mask = UV_roi_mask;          % user-defined ROI
    handles.smoothed_flag = smoothed_flag;      % double-type image of local thresholding - (masked, no threshold applied, Gaussian smoothed)
    handles.edge_UV = edge_UV;                  % delineated UV region
    handles.overlap_thresh = overlap_thresh;    % overlap color map
    handles.UV_region = UV_region;              % UV region (binary) - (after masking, and overlap threshold applied)
    handles.current_axes = 2;                   % specify which axes to change display options

    set(handles.ChooseAxes_uipanel,'SelectedObject',handles.axes2_radiobutton);
    set(handles.axes2_text, 'String', handles.imagename);

    %% display image
    disp_im(handles.edge_UV) = 255;
    imshow(disp_im, 'Parent', handles.axes2);

    set(handles.Image_display,'SelectedObject',handles.redfree_im_radiobutton);

    %% set slider properties
    max_val = ceil(max(smoothed_flag(:)));
    slider_step = 1/max_val;

    set(handles.threshold_slider, 'Min', 0);
    set(handles.threshold_slider, 'Max', max_val);
    set(handles.threshold_slider, 'SliderStep', [slider_step , slider_step+2]);

    set(handles.measure_pushbutton,'Enable','on');
    
    % Update handles structure
    guidata(hObject, handles);
end


function UVAF_area_text_Callback(hObject, eventdata, handles)
% hObject    handle to UVAF_area_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UVAF_area_text as text
%        str2double(get(hObject,'String')) returns contents of UVAF_area_text as a double


% --- Executes during object creation, after setting all properties.
function UVAF_area_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UVAF_area_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function foldername_edit_text_Callback(hObject, eventdata, handles)
% hObject    handle to foldername_edit_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of foldername_edit_text as text
%        str2double(get(hObject,'String')) returns contents of foldername_edit_text as a double

% --- Executes during object creation, after setting all properties.
function foldername_edit_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to foldername_edit_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in ChooseAxes_uipanel.
function ChooseAxes_uipanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ChooseAxes_uipanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'String')
    case 'Axes1'
        handles.current_axes = 1;
        set(handles.ChooseAxes_uipanel,'SelectedObject',handles.axes1_radiobutton);
        set(handles.Image_display,'SelectedObject', []);
    case 'Axes2'
        handles.current_axes = 2;
        set(handles.ChooseAxes_uipanel,'SelectedObject',handles.axes2_radiobutton);
        set(handles.Image_display,'SelectedObject', []);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in Image_display.
function Image_display_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Image_display 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'String') 
    case 'Original Image'
        original_im_radiobutton_Callback(hObject, eventdata, handles);
    case 'Red-free Image'
        redfree_im_radiobutton_Callback(hObject, eventdata, handles);
    case 'Grayscale Image'
        grayscale_im_radiobutton_Callback(hObject, eventdata, handles);
    case 'Colour Photo'
        colourphoto_radiobutton_Callback(hObject, eventdata, handles);
end

% --- Executes on button press in original_im_radiobutton.
function original_im_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to original_im_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of original_im_radiobutton

disp_im= handles.im_resized;
handles.disp_im = disp_im;
    
if handles.current_axes == 1 
    imshow(disp_im, 'Parent', handles.axes1);
elseif handles.current_axes == 2 
    disp_im(handles.edge_UV) = 255;
    imshow(disp_im, 'Parent', handles.axes2);
else
    msgbox('ERROR');
end

set(handles.Image_display,'SelectedObject',handles.original_im_radiobutton);
set(handles.axes1_text, 'String', handles.imagename);
set(handles.axes2_text, 'String', handles.imagename);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in redfree_im_radiobutton.
function redfree_im_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to redfree_im_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of redfree_im_radiobutton

%% Display image
disp_im = handles.red_free;
handles.disp_im = disp_im;

if handles.current_axes == 1 
    imshow(disp_im, 'Parent', handles.axes1);
elseif handles.current_axes == 2 
    disp_im(handles.edge_UV) = 255;
    imshow(disp_im, 'Parent', handles.axes2);
else
    msgbox('ERROR');
end

set(handles.Image_display,'SelectedObject',handles.redfree_im_radiobutton);
set(handles.axes1_text, 'String', handles.imagename);
set(handles.axes2_text, 'String', handles.imagename);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in grayscale_im_radiobutton.
function grayscale_im_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to grayscale_im_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of grayscale_im_radiobutton

disp_im = handles.image;
handles.disp_im = disp_im;

if handles.current_axes == 1 
    imshow(disp_im, 'Parent', handles.axes1);
elseif handles.current_axes == 2 
    disp_im(handles.edge_UV) = 255;
    imshow(disp_im, 'Parent', handles.axes2);
else
    msgbox('ERROR');
end

set(handles.Image_display,'SelectedObject',handles.grayscale_im_radiobutton);
set(handles.axes1_text, 'String', handles.imagename);
set(handles.axes2_text, 'String', handles.imagename);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in colourphoto_radiobutton.
function colourphoto_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to original_im_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of original_im_radiobutton

imagename = get(handles.axes1_text, 'String');    %ie. OS temporal UV_AF
parts = regexp(imagename, ' ', 'split');
LR = parts{1};
NT = parts{2};

ListImageNames = get(handles.folder_listbox,'string');
num_im = size(ListImageNames,1);
for i = 1:num_im
    name = ListImageNames{i};
    parts = regexp(name, ' ', 'split');
    num_parts = length(parts);
    if num_parts >=3
        if isempty(regexpi(parts{3}, 'colour')) == 0 || isempty(regexpi(parts{3}, 'color')) == 0
            if strcmpi(LR, parts{1}) == 1 && strcmpi(NT, parts{2}) == 1
                colourphoto_index = i;
                break;
            end
        end
    end
end

folder_name = handles.ImageFolder;
colour_pathname = [folder_name handles.delimiter ListImageNames{colourphoto_index}];

disp_im = imread(colour_pathname);
% disp_im = imresize(disp_im, handles.resize_factor);
% handles.disp_im = disp_im;
if handles.current_axes == 1 
    imshow(disp_im, 'Parent', handles.axes1);
    set(handles.axes1_text, 'String', ListImageNames{colourphoto_index});
elseif handles.current_axes == 2 
    imshow(disp_im, 'Parent', handles.axes2);
    set(handles.axes2_text, 'String', ListImageNames{colourphoto_index});
else
    msgbox('ERROR');
end

set(handles.Image_display,'SelectedObject',handles.colourphoto_radiobutton);

% Update handles structure
guidata(hObject, handles);

% --- Executes on slider movement.
function threshold_slider_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

overlap_thresh = get(handles.threshold_slider, 'Value');
overlap_thresh = round(overlap_thresh);
handles.overlap_thresh = overlap_thresh;
set(handles.threshold_edit_text, 'String', overlap_thresh);

handles.selected_points = [];

[UV_area_mm, edge_UV, UV_region] = Refine_UVAF(handles); %Refine_UVAF(uv_flag, UV_roi_mask, overlap_thresh, resize_factor);

% update handles and display in gui
handles.UV_region = UV_region;
handles.edge_UV = edge_UV;

set(handles.UVAF_area_text, 'String', UV_area_mm)
handles.edge_UV = edge_UV;

%% display image
disp_im = handles.disp_im;
disp_im(handles.edge_UV) = 255;
imshow(disp_im, 'Parent', handles.axes2);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function threshold_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function threshold_edit_text_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_edit_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_edit_text as text
%        str2double(get(hObject,'String')) returns contents of threshold_edit_text as a double

overlap_thresh = get(handles.threshold_edit_text, 'String');
overlap_thresh = round(str2num(overlap_thresh));
handles.overlap_thresh = overlap_thresh;
set(handles.threshold_slider, 'Value', overlap_thresh);

handles.selected_points = [];

[UV_area_mm, edge_UV, UV_region] = Refine_UVAF(handles);

% update handles and display in gui
handles.UV_region = UV_region;
handles.edge_UV = edge_UV;

set(handles.UVAF_area_text, 'String', UV_area_mm)
handles.edge_UV = edge_UV;

%% display image
disp_im = handles.disp_im;
disp_im(handles.edge_UV) = 255;
imshow(disp_im, 'Parent', handles.axes2);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function threshold_edit_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_edit_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reset_pushbutton.
function reset_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to reset_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

overlap_thresh = get(handles.threshold_edit_text, 'String');
overlap_thresh = round(str2num(overlap_thresh));
handles.overlap_thresh = overlap_thresh;
set(handles.threshold_slider, 'Value', overlap_thresh);

handles.selected_points = [];

[UV_area_mm, edge_UV, UV_region] = Refine_UVAF(handles);

% update handles and display in gui
handles.UV_region = UV_region;
handles.edge_UV = edge_UV;

set(handles.UVAF_area_text, 'String', UV_area_mm)
handles.edge_UV = edge_UV;

%% display image
disp_im = handles.disp_im;
disp_im(handles.edge_UV) = 255;
imshow(disp_im, 'Parent', handles.axes2);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in remove_region_pushbutton.
function remove_region_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to remove_region_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UV_region = handles.UV_region;
selected_points = handles.selected_points;
% resize_factor = handles.resize_factor;

disp_im = handles.disp_im;
disp_im(handles.edge_UV) = 255;
figure; imshow(disp_im);
title('Remove regions')
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
[x,y] = ginput;
close;

for i = 1:size(x,1)
    if UV_region(round(y(i)),round(x(i))) > 0
        selected_points = [selected_points; y(i) x(i)];
    else
        msgbox('Not a valid region to remove.', 'Remove Region', 'warn')
        return
    end
end

handles.selected_points = selected_points;

[UV_area_mm, edge_UV, new_UV_region] = RemoveRegions(handles);  %RemoveRegions(UV_region, resize_factor, selected_points);

handles.UV_region = new_UV_region;
handles.edge_UV = edge_UV;
set(handles.UVAF_area_text, 'String', UV_area_mm)

%% display image
disp_im = handles.disp_im;
disp_im(handles.edge_UV) = 255;
imshow(disp_im, 'Parent', handles.axes2);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in reproduce_pushbutton.
function reproduce_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to reproduce_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pathname = handles.pathname;

[UV_area_mm, edge_UV, im, overlap_thresh, subjective, user_initials, comments] = Reproduce_Results(handles);

set(handles.UVAF_area_text, 'String', UV_area_mm)
set(handles.threshold_edit_text, 'String', overlap_thresh);
set(handles.threshold_slider, 'Value', overlap_thresh);
set(handles.UserInitials_edit, 'String', user_initials);
set(handles.comments_edit, 'String', comments);

if strcmp(subjective, '***') == 1
    set(handles.subjective_checkbox, 'Value', 1);
else
    set(handles.subjective_checkbox, 'Value', 0);
end
    
% NOTE: does not overwrite UV_region, UV_flag... just for display purposes.

%% image title          % WHAT TO PUT WHEN THERE IS AN ERROR... ##################
path_parts = regexp(pathname, handles.delimiter, 'split');
num_parts = length(path_parts);
imagename = path_parts{num_parts};

set(handles.axes2_text, 'String', imagename);

%% display image
[m,n,o] = size(im);
imG = im(:,:,2);
imB = im(:,:,3);
meanG = mean2(double(imG));

gamma = (1/50)*meanG + 0.6;              % min gamma = y-intercept
if gamma > 1.5
    gamma = 1.5;
end

tol = stretchlim(imB, [0.02 0.99]);
im2 = imadjust(im, tol, [0;1], gamma);
imG2 = im2(:,:,2);
imB2 = im2(:,:,3);

disp_im = cat(3, zeros(m,n), imG2, imB2);
handles.disp_im = disp_im;

disp_im(edge_UV) = 255;
imshow(disp_im, 'Parent', handles.axes2);

figure; imshow(disp_im);

set(handles.Image_display,'SelectedObject',handles.redfree_im_radiobutton);
set(handles.ChooseAxes_uipanel,'SelectedObject',handles.axes2_radiobutton);

% --- Executes on button press in save_pushbutton.
function save_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

foldername = handles.ImageFolder;
filename = handles.excel_filename;
path = [foldername handles.delimiter filename];
imagename = get(handles.axes1_text, 'String');

%% get user initials
initials = get(handles.UserInitials_edit, 'String');

if isempty(initials) == 1 || strcmp(initials, {''}) == 1
    initials = inputdlg('Please enter initials:');
    initials = initials{1};
    set(handles.UserInitials_edit, 'String', initials);
end

%% save comments
comments = get(handles.comments_edit, 'String');
if isempty(comments)
    comments = repmat({' '},1,1);
end
% notes = repmat({' '},1,1);

% Variables and Area  (sheet 1)
resize_factor = handles.resize_factor;
overlap_thresh = str2num(get(handles.threshold_edit_text, 'String'));
UV_area_mm = str2num(get(handles.UVAF_area_text, 'String'));

if isempty(regexpi(imagename, 'OD NASAL')) == 0
    xlwrite(path, resize_factor, 'Results', 'B2');
    xlwrite(path, overlap_thresh, 'Results', 'C2');
    xlwrite(path, UV_area_mm, 'Results', 'D2');
%     xlwrite(path, {notes}, 'Results', 'E2');
    xlwrite(path, {datestr(now)}, 'Results', 'F2');
    [num, ROI_date, data] = xlsread(path, 'ROI Coordinates', 'C1');
    xlwrite(path, ROI_date, 'Results', 'G2');   
    xlwrite(path, {initials}, 'Results', 'I2');
    xlwrite(path, {comments}, 'Results', 'J2');
    
elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
    xlwrite(path, resize_factor, 'Results', 'B3');
    xlwrite(path, overlap_thresh, 'Results', 'C3');
    xlwrite(path, UV_area_mm, 'Results', 'D3');
%     xlwrite(path, {notes}, 'Results', 'E3');
    xlwrite(path, {datestr(now)}, 'Results', 'F3')
    [num, ROI_date, data] = xlsread(path, 'ROI Coordinates', 'F1');
    xlwrite(path, ROI_date, 'Results', 'G3');   
    xlwrite(path, {initials}, 'Results', 'I3');
    xlwrite(path, {comments}, 'Results', 'J3');

elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
    xlwrite(path, resize_factor, 'Results', 'B4');
    xlwrite(path, overlap_thresh, 'Results', 'C4');
    xlwrite(path, UV_area_mm, 'Results', 'D4');
%     xlwrite(path, {notes}, 'Results', 'E4');
    xlwrite(path, {datestr(now)}, 'Results', 'F4');
    [num, ROI_date, data] = xlsread(path, 'ROI Coordinates', 'I1');
    xlwrite(path, ROI_date, 'Results', 'G4');    
    xlwrite(path, {initials}, 'Results', 'I4');
    xlwrite(path, {comments}, 'Results', 'J4');
    
elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
    xlwrite(path, resize_factor, 'Results', 'B5');
    xlwrite(path, overlap_thresh, 'Results', 'C5');
    xlwrite(path, UV_area_mm, 'Results', 'D5');
%     xlwrite(path, {notes}, 'Results', 'E5');
    xlwrite(path, {datestr(now)}, 'Results', 'F5');
    [num, ROI_date, data] = xlsread(path, 'ROI Coordinates', 'L1');
    xlwrite(path, ROI_date, 'Results', 'G5');
    xlwrite(path, {initials}, 'Results', 'I5');
    xlwrite(path, {comments}, 'Results', 'J5');
    
else
    msgbox('Invalid image name.');
    return;
end

% Calculate total UVAF
area_array = xlsread(path, 'Results', 'D2:D5');
total_area = 0;
num_entries = size(area_array,1);
for i = 1:num_entries
    if isnan(area_array(i)) == 0
        total_area = total_area+area_array(i);
    end
end
xlwrite(path, {'Total Area'}, 'Results', 'C7');
xlwrite(path, total_area, 'Results', 'D7');

%% Points to remove (sheet3)
selected_points = handles.selected_points;
clear_cells = nan(10,2);
if isempty(selected_points) == 0
    if isempty(regexpi(imagename, 'OD NASAL')) == 0
        xlwrite(path, clear_cells, 'Removed Regions', 'A2');
        xlwrite(path, selected_points, 'Removed Regions', 'A2');

    elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
        xlwrite(path, clear_cells, 'Removed Regions', 'C2');
        xlwrite(path, selected_points, 'Removed Regions', 'C2');   

    elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
        xlwrite(path, clear_cells, 'Removed Regions', 'E2');
        xlwrite(path, selected_points, 'Removed Regions', 'E2');

    elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0 
        xlwrite(path, clear_cells, 'Removed Regions', 'G2');
        xlwrite(path, selected_points, 'Removed Regions', 'G2');
    end
end

%% Subjective image
subjective_checkbox = get(handles.subjective_checkbox, 'Value');

if subjective_checkbox == 1
    if isempty(regexpi(imagename, 'OD NASAL')) == 0
        xlwrite(path, {'***'}, 'Results', 'H2');
    elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
        xlwrite(path, {'***'}, 'Results', 'H3');  
    elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
        xlwrite(path, {'***'}, 'Results', 'H4');
    elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
       xlwrite(path, {'***'}, 'Results', 'H5');
    else
        msgbox('Invalid image name.');
        return;
    end

elseif subjective_checkbox == 0 
    clear_cell = repmat({' '},1,1);
    if isempty(regexpi(imagename, 'OD NASAL')) == 0
        xlwrite(path, clear_cell, 'Results', 'H2');
    elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
        xlwrite(path, clear_cell, 'Results', 'H3');
    elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
        xlwrite(path, clear_cell, 'Results', 'H4');
    elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
       xlwrite(path, clear_cell, 'Results', 'H5');
    else
        msgbox('Invalid image name.');
        return;
    end
end

%% Weak image
weak_checkbox = get(handles.weak_checkbox, 'Value');

if weak_checkbox == 1
    if isempty(regexpi(imagename, 'OD NASAL')) == 0
        xlwrite(path, {'Weak CUVAF'}, 'Results', 'E2');
    elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
        xlwrite(path, {'Weak CUVAF'}, 'Results', 'E3');  
    elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
        xlwrite(path, {'Weak CUVAF'}, 'Results', 'E4');
    elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
       xlwrite(path, {'Weak CUVAF'}, 'Results', 'E5');
    else
        msgbox('Invalid image name.');
        return;
    end

elseif weak_checkbox == 0 
    clear_cell = repmat({' '},1,1);
    if isempty(regexpi(imagename, 'OD NASAL')) == 0
        xlwrite(path, clear_cell, 'Results', 'E2');
    elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
        xlwrite(path, clear_cell, 'Results', 'E3');
    elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
        xlwrite(path, clear_cell, 'Results', 'E4');
    elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
       xlwrite(path, clear_cell, 'Results', 'E5');
    else
        msgbox('Invalid image name.');
        return;
    end
end

%% display image
disp_im = handles.disp_im;
edge_UV = handles.edge_UV;
edge_UV = bwmorph(edge_UV, 'dilate', 1);
disp_im(edge_UV) = 255;

figure; imshow(disp_im)

save_path = [foldername handles.delimiter 'marked ' handles.CUVAF_version '-' imagename];
imwrite(disp_im, save_path)
handles = LoadImageList(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in NoUVAF_pushbutton.
function NoUVAF_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to NoUVAF_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

foldername = handles.ImageFolder;
filename = handles.excel_filename;
path = [foldername handles.delimiter filename];

imagename = get(handles.axes1_text, 'String');

%% Read Image
pathname = handles.pathname;
im=imread(pathname);

im = imresize(im, handles.resize_factor);
[m,n,o] = size(im);
imG = im(:,:,2);
imB = im(:,:,3);

meanG = mean(double(imG(:)));
sdG = std(double(imG(:)));
meanB = mean(double(imB(:)));
sdB = std(double(imB(:)));

gamma = (1/50)*meanG + 0.6;              % min gamma = y-intercept
if gamma > 1.5
    gamma = 1.5;
end

tol = stretchlim(imB, [0.02 0.99]);
im2 = imadjust(im, tol, [0;1], gamma);
imG2 = im2(:,:,2);
imB2 = im2(:,:,3);

disp_im = cat(3, zeros(m,n), imG2, imB2);

save_path = [foldername handles.delimiter 'marked ' handles.CUVAF_version '-' imagename];
imwrite(disp_im, save_path)

image = rgb2gray(disp_im);

%% image stats - tophat
max_int_pp = max(image(:));
min_int_pp = min(image(:));

im_mean_pp = mean(double(image(:)));
im_sd_pp = std(double(image(:)));

image_stats = [im_mean_pp im_sd_pp min_int_pp max_int_pp meanG sdG meanB sdB];

%% get user initials
initials = get(handles.UserInitials_edit, 'String');

if isempty(initials) == 1 || strcmp(initials, {''}) == 1
    initials = inputdlg('Please enter initials:');
    initials = initials{1};
    set(handles.UserInitials_edit, 'String', initials);
end

%% save comments
comments = get(handles.comments_edit, 'String');

%% Variables and Area  (sheet 1)
resize_factor = {'N/A'};
overlap_thresh = {'N/A'};
UV_area_mm = 0;
% notes = {'No UVAF'};

if isempty(regexpi(imagename, 'OD NASAL')) == 0
    % write sheet1 
    clear_cells = nan(1,3);
    xlwrite(path, clear_cells, 'Results', 'B2');
    xlwrite(path, resize_factor, 'Results', 'B2');
    xlwrite(path, overlap_thresh, 'Results', 'C2');
    xlwrite(path, UV_area_mm, 'Results', 'D2');
%     xlwrite(path, {notes}, 'Results', 'E2');
    xlwrite(path, {datestr(now)}, 'Results', 'F2');
    xlwrite(path, {'N/A'}, 'Results', 'G2');
    xlwrite(path, {initials}, 'Results', 'I2');
    xlwrite(path, {comments}, 'Results', 'J2');
    xlwrite(path, image_stats, 'Results', 'N2');
    
    %write sheet2
    clear_cells = nan(60,2);
    xlwrite(path, clear_cells, 'ROI Coordinates', 'A2');
    clear_cells = nan(1,2);
    xlwrite(path, clear_cells, 'ROI Coordinates', 'B1');
    
    % write sheet3
    clear_cells = nan(10,2);
    xlwrite(path, clear_cells, 'Removed Regions', 'A2');

elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
    % write sheet1 
    clear_cells = nan(1,3);
    xlwrite(path, clear_cells, 'Results', 'B3');
    xlwrite(path, resize_factor, 'Results', 'B3');
    xlwrite(path, overlap_thresh, 'Results', 'C3');
    xlwrite(path, UV_area_mm, 'Results', 'D3');
%     xlwrite(path, {notes}, 'Results', 'E3');
    xlwrite(path, {datestr(now)}, 'Results', 'F3');
    xlwrite(path, {'N/A'}, 'Results', 'G3');
    xlwrite(path, {initials}, 'Results', 'I3'); 
    xlwrite(path, {comments}, 'Results', 'J3');
    xlwrite(path, image_stats, 'Results', 'N3');
    
    % write sheet2
    clear_cells = nan(60,2);
    xlwrite(path, clear_cells, 'ROI Coordinates', 'D2');
    clear_cells = nan(1,2);
    xlwrite(path, clear_cells, 'ROI Coordinates', 'E1');
    
    % write sheet3
    clear_cells = nan(10,2);
    xlwrite(path, clear_cells, 'Removed Regions', 'C2');
  
elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
    % write sheet1
    clear_cells = nan(1,3);
    xlwrite(path, clear_cells, 'Results', 'B4');
    xlwrite(path, resize_factor, 'Results', 'B4');
    xlwrite(path, overlap_thresh, 'Results', 'C4');
    xlwrite(path, UV_area_mm, 'Results', 'D4');
%     xlwrite(path, {notes}, 'Results', 'E4');
    xlwrite(path, {datestr(now)}, 'Results', 'F4');
    xlwrite(path, {'N/A'}, 'Results', 'G4');
    xlwrite(path, {initials}, 'Results', 'I4'); 
    xlwrite(path, {comments}, 'Results', 'J4');
    xlwrite(path, image_stats, 'Results', 'N4');
    
    % write sheet2
    clear_cells = nan(60,2);
    xlwrite(path, clear_cells, 'ROI Coordinates', 'G2');
    clear_cells = nan(1,2);
    xlwrite(path, clear_cells, 'ROI Coordinates', 'H1');
    
    % write sheet3
    clear_cells = nan(10,2);
    xlwrite(path, clear_cells, 'Removed Regions', 'E2');

elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
    % write sheet1
    clear_cells = nan(1,3);
    xlwrite(path, clear_cells, 'Results', 'B5');
    xlwrite(path, resize_factor, 'Results', 'B5');
    xlwrite(path, overlap_thresh, 'Results', 'C5');
    xlwrite(path, UV_area_mm, 'Results', 'D5');
%     xlwrite(path, {notes}, 'Results', 'E5');
    xlwrite(path, {datestr(now)}, 'Results', 'F5');
    xlwrite(path, {'N/A'}, 'Results', 'G5');
    xlwrite(path, {initials}, 'Results', 'I5'); 
    xlwrite(path, {comments}, 'Results', 'J5');
    xlwrite(path, image_stats, 'Results', 'N5');
    
    % write sheet2
    clear_cells = nan(60,2);
    xlwrite(path, clear_cells, 'ROI Coordinates', 'J2');
    clear_cells = nan(1,2);
    xlwrite(path, clear_cells, 'ROI Coordinates', 'K1');
    
    % write sheet3 
    clear_cells = nan(10,2);
    xlwrite(path, clear_cells, 'Removed Regions', 'G2');
    
else
    msgbox('Invalid image name.');
    return;
end

% Calculate total UVAF
area_array = xlsread(path, 'Results', 'D2:D5');
total_area = 0;
num_entries = size(area_array,1);
for i = 1:num_entries
    if isnan(area_array(i)) == 0
        total_area = total_area+area_array(i);
    end
end
xlwrite(path, {'Total Area'}, 'Results', 'C7');
xlwrite(path, total_area, 'Results', 'D7');

%% Subjective image
subjective_checkbox = get(handles.subjective_checkbox, 'Value');

if subjective_checkbox == 1
    if isempty(regexpi(imagename, 'OD NASAL')) == 0
        xlwrite(path, {'***'}, 'Results', 'H2');
    elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
        xlwrite(path, {'***'}, 'Results', 'H3');  
    elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
        xlwrite(path, {'***'}, 'Results', 'H4');
    elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
       xlwrite(path, {'***'}, 'Results', 'H5');
    else
        msgbox('Invalid image name.');
        return;
    end

elseif subjective_checkbox == 0    
    clear_cell = repmat({' '},1,1);
    if isempty(regexpi(imagename, 'OD NASAL')) == 0
        xlwrite(path, clear_cell, 'Results', 'H2');
    elseif isempty(regexpi(imagename, 'OD TEMPORAL')) == 0
        xlwrite(path, clear_cell, 'Results', 'H3');
    elseif isempty(regexpi(imagename, 'OS NASAL')) == 0
        xlwrite(path, clear_cell, 'Results', 'H4');
    elseif isempty(regexpi(imagename, 'OS TEMPORAL')) == 0
       xlwrite(path, clear_cell, 'Results', 'H5');
    else
        msgbox('Invalid image name.');
        return;
    end
end

set(handles.save_pushbutton,'Enable','off');
handles = LoadImageList(handles);

% Update handles structure
guidata(hObject, handles);

msgbox('Data for "No CUVAF" saved.');


function directory_edit_Callback(hObject, eventdata, handles)
% hObject    handle to directory_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of directory_edit as text
%        str2double(get(hObject,'String')) returns contents of directory_edit as a double


% --- Executes during object creation, after setting all properties.
function directory_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to directory_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in subjective_checkbox.
function subjective_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to subjective_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subjective_checkbox
subjective_checkbox = get(handles.subjective_checkbox, 'Value');

if subjective_checkbox == 1
    set(handles.subjective_checkbox, 'Value', 1);
elseif subjective_checkbox == 0
    set(handles.subjective_checkbox, 'Value', 0);
end

% --- Executes on button press in weak_checkbox.
function weak_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to weak_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of weak_checkbox
weak_checkbox = get(handles.weak_checkbox, 'Value');

if weak_checkbox == 1
    set(handles.weak_checkbox, 'Value', 1);
elseif weak_checkbox == 0
    set(handles.weak_checkbox, 'Value', 0);
end

function UserInitials_edit_Callback(hObject, eventdata, handles)
% hObject    handle to UserInitials_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UserInitials_edit as text
%        str2double(get(hObject,'String')) returns contents of UserInitials_edit as a double

% --- Executes during object creation, after setting all properties.
function UserInitials_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UserInitials_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function comments_edit_Callback(hObject, eventdata, handles)
% hObject    handle to comments_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of comments_edit as text
%        str2double(get(hObject,'String')) returns contents of comments_edit as a double


% --- Executes during object creation, after setting all properties.
function comments_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to comments_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function RenameFile_Callback(hObject, eventdata, handles)
% hObject    handle to RenameFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

folder_path = handles.ImageFolder;
source = handles.pathname;          %path of current image
old_name = handles.imagename;
new_name = inputdlg('Please enter new file name:', 'Rename File', 1, {old_name});    
new_name = new_name{1};
destination = fullfile(folder_path, new_name);
movefile(source, destination);

handles = LoadImageList(handles);

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function Listbox_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Listbox_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

