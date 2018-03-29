function varargout = showMe(varargin)
% A general purpose viewer for multi-dimensional arrays
% Usage: showMe(M, colorMap)
% Where M is either a 2D, 3D, 4D, or 5D matrix & colorMap is optional.
% The dimensions are interpreted as: X, Y, Z, Time, RGB
% A special case occurs when the last dimension has length 3 and is interpreted as RGB channels.
% e.g. A size of: 100x100x2 is interpreted as 2 slices of 100x100 images.
% 		   		  100x100x3 is interpreted as a single slice of a 100x100 RGB image.
% For multi-dimensional data, you can click on the image, then use arrow keys.
% Up/Down keys will scroll through slices. Left/Right go through time.

% Note: Matlab might stop automatically creating new figures after running
% this GUI. In other words, if no figure is open then the command plot(1:10)
% will not make a new figure. To return this behaviour, try the command: close all 

% Last Modified by GUIDE v2.5 04-Apr-2015 20:18:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0; % Setting to 0 allows multiple instances
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showMe_OpeningFcn, ...
                   'gui_OutputFcn',  @showMe_OutputFcn, ...
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


% --- Executes just before showMe is made visible.
function showMe_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for showMe
handles.output = hObject;
% Decide input type and display appropriate image
handles.imageData = varargin{1};
if nargin > 4
	handles.cMap = varargin{2};
else
    handles.cMap = colormap('jet');
end
S = size(handles.imageData);
handles.imMode = length(S);
handles.isRGB = 0;
handles.haveSlices = 0;
handles.haveTime = 0;
switch handles.imMode
    case 2 					% 2D data
    	[handles.sX, handles.sY] = size(handles.imageData);
        imshow(handles.imageData, []);
        set(handles.sSliceIndex,'visible','off')
        set(handles.sTimeIndex,'visible','off')
    case 3 		
    	[handles.sX, handles.sY, handles.sZ] = size(handles.imageData);
    	if handles.sZ == 3 	% 2D RGB data
        	imshow(handles.imageData(:,:,:), []);
        	handles.sZ = 0;
        	set(handles.sSliceIndex,'visible','off')
        	set(handles.sTimeIndex,'visible','off')
        	handles.isRGB = 1;
        else 				% 3D data
        	imshow(handles.imageData(:,:,1),[]);
        	set(handles.sSliceIndex,'Min',1);
	        set(handles.sSliceIndex,'Max',handles.sZ);
	        set(handles.sSliceIndex,'SliderStep',[1/handles.sZ,5/handles.sZ]);
	        handles.sliceIndex = 1;
	        set(handles.sSliceIndex,'Value',handles.sliceIndex);
	        set(handles.sTimeIndex,'visible','off')
	        handles.haveSlices = 1;
        end
    case 4
    	[handles.sX, handles.sY, handles.sZ, handles.sT] = size(handles.imageData);
        if handles.sT == 3 	% 3D RGB data
    		imshow(squeeze(handles.imageData(:,:,1,:)),[]);
	        set(handles.sSliceIndex,'Min',1);
	        set(handles.sSliceIndex,'Max',handles.sZ);
	        set(handles.sSliceIndex,'SliderStep',[1/handles.sZ,5/handles.sZ]);
	        handles.sliceIndex = 1;
	        set(handles.sSliceIndex,'Value',handles.sliceIndex);
        	set(handles.sTimeIndex,'visible','off')
        	handles.isRGB = 1;
        	handles.haveSlices = 1;
        else 				% 4D data
	        imshow(squeeze(handles.imageData(:,:,1,1)),[]);
	        set(handles.sSliceIndex,'Min',1);
	        set(handles.sSliceIndex,'Max',handles.sZ);
	        set(handles.sSliceIndex,'SliderStep',[1/handles.sZ,5/handles.sZ]);
	        handles.sliceIndex = 1;
	        set(handles.sSliceIndex,'Value',handles.sliceIndex);
	        set(handles.sTimeIndex,'Min',1);
	        set(handles.sTimeIndex,'Max',handles.sT);
	        set(handles.sTimeIndex,'SliderStep',[1/handles.sT,5/handles.sT]);
	        handles.timeIndex = 1;
	        set(handles.sTimeIndex,'Value',handles.timeIndex);
	        handles.haveSlices = 1;
	        handles.haveTime = 1;
    	end
    case 5 					% 4D RGB data
    	[handles.sX, handles.sY, handles.sZ, handles.sT, grbg] = size(handles.imageData);
    	if grbg~= 3
    		error('Can not handle this kind of 5D input.')
    	end
    	imshow(squeeze(handles.imageData(:,:,1,1,:)),[]);
        set(handles.sSliceIndex,'Min',1);
        set(handles.sSliceIndex,'Max',handles.sZ);
        set(handles.sSliceIndex,'SliderStep',[1/handles.sZ,5/handles.sZ]);
        handles.sliceIndex = 1;
        set(handles.sSliceIndex,'Value',handles.sliceIndex);
        set(handles.sTimeIndex,'Min',1);
        set(handles.sTimeIndex,'Max',handles.sT);
        set(handles.sTimeIndex,'SliderStep',[1/handles.sT,5/handles.sT]);
        handles.timeIndex = 1;
        set(handles.sTimeIndex,'Value',handles.timeIndex);
        handles.haveSlices = 1;
        handles.haveTime = 1;
    otherwise
        error('Can not handle input. Only 2D, 3D or 4D data allowed.');
end
colormap(gca,handles.cMap)
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes showMe wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = showMe_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function updateDisplay
	handles=guidata(gcbf);  % Load all the handles
	switch handles.imMode
		case 2
			curData = handles.imageData;
		case 3
			if handles.isRGB == 1
				curData = handles.imageData(:,:,1:3);
			else
				curData = handles.imageData(:,:,handles.sliceIndex);
			end
		case 4
			if handles.isRGB
				curData = squeeze(handles.imageData(:,:,handles.sliceIndex,1:3));
			else
				curData = handles.imageData(:,:,handles.sliceIndex, handles.timeIndex);
			end
		case 5
			curData = squeeze(handles.imageData(:,:,handles.sliceIndex, handles.timeIndex,1:3));
	end
	% The x/ylim business is there so that matlab doesn't reset the user's zoom
	XL=xlim;
	YL=ylim;
	imshow(curData,[])
	colormap(gca,handles.cMap)
	xlim(XL);
	ylim(YL);
    % This next section is for custom colourmap limits
    if get(handles.cAutoCMap,'Value')==1
        caxis auto;
        curVals = caxis;
        set(handles.tCMin,'String',num2str(curVals(1,1)));
        set(handles.tCMax,'String',num2str(curVals(1,2)));
    else
        cMin = get(handles.tCMin,'String');
        cMax = get(handles.tCMax,'String');
        caxis([str2num(cMin) str2num(cMax)]);
    end

% --- Executes during object creation, after setting all properties.
% I think this gives an error, not sure why, but everything seems to work
function mainAxes_CreateFcn(hObject, eventdata, handles)
    set(gcf,'toolbar','figure');  % Shows toolbar (for displaying colorbar or zooming)
    set(gcf, 'KeyPressFcn', @(h, eventdata) GetKey(eventdata) );

% This allows user to scroll through data using arrow keys
function GetKey (eventdata)
    handles=guidata(gcbf);
    if (strcmp(eventdata.Key, 'rightarrow')) && handles.haveTime
        handles.timeIndex = handles.timeIndex + 1;
        if handles.timeIndex > handles.sT
            handles.timeIndex = handles.sT;
        end
        guidata(gcbf, handles);
        updateTimeIndex;
    elseif (strcmp(eventdata.Key, 'leftarrow')) && handles.haveTime
        handles.timeIndex = handles.timeIndex - 1;
        if handles.timeIndex < 1
            handles.timeIndex = 1;
        end
        guidata(gcbf, handles);
        updateTimeIndex;
    elseif (strcmp(eventdata.Key, 'uparrow')) && handles.haveSlices
        handles.sliceIndex = handles.sliceIndex + 1;
        if handles.sliceIndex > handles.sZ
            handles.sliceIndex = handles.sZ;
        end
        guidata(gcbf, handles);
        updateSliceIndex;
    elseif (strcmp(eventdata.Key, 'downarrow')) && handles.haveSlices
        handles.sliceIndex = handles.sliceIndex - 1;
        if handles.sliceIndex < 1
            handles.sliceIndex = 1;
        end
        guidata(gcbf, handles);
        updateSliceIndex;
    end

function updateSliceIndex
    handles=guidata(gcbf);
    set(handles.sSliceIndex,'Value',handles.sliceIndex);
    updateDisplay;

function updateTimeIndex
    handles=guidata(gcbf);  % Load all the handles
    set(handles.sTimeIndex,'Value',handles.timeIndex);
    updateDisplay; 

% --- Executes on slider movement.
function sSliceIndex_Callback(hObject, eventdata, handles)
	handles.sliceIndex=round(get(hObject,'Value'));  % Get new value
	guidata(hObject, handles);
	updateDisplay

% --- Executes on slider movement.
function sTimeIndex_Callback(hObject, eventdata, handles)
	handles.timeIndex=round(get(hObject,'Value'));  % Get new value
	guidata(hObject, handles);
	updateDisplay
    
function bColourMap_Callback(hObject, eventdata, handles)
    colormapeditor;     % Launch matlab's colourmap editor

function cAutoCMap_Callback(hObject, eventdata, handles)
    updateDisplay

function tCMin_Callback(hObject, eventdata, handles)
    set(handles.cAutoCMap,'Value',0);
    updateDisplay;
    
function tCMax_Callback(hObject, eventdata, handles)
    set(handles.cAutoCMap,'Value',0);
    updateDisplay;


% --- Executes during object creation, after setting all properties.
function sSliceIndex_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sTimeIndex_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function tCMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tCMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function tCMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tCMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
