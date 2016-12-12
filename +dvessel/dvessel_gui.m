function [varargout] = dvessel_gui(varargin)
%3D vessel
%Andreas S�derlund

if nargin==0
  varargin = {'init'};
end

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init(force) %#ok<DEFNU,INUSD>
%------------
% Initialize GUI
%
global DATA

if nargin == 0
  force = false;
end

% Launch GUI in a mygui object.
DATA.GUI.DVESSELGUI = mygui('+dvessel/dvessel_gui.fig');

initializegui;

%---------------------
function initializegui
%---------------------
%Initialize the gui

global DATA SET
gui = DATA.GUI.DVESSELGUI;

gui.loadedimage = true; % only one instance of it. Should it be removed?
gui.selectedImgs = 1;
gui.smoothObjectsCell = {'All data', 'Good data', 'Along splines'};

% Plot the first image selected automatically
selectbox_Callback;

% Set the names for the possible image series to choose from.
selBoxCell = cell(max(size(SET)), 1);
for i = 1:size(selBoxCell, 1)
    selBoxCell(i) = {sprintf('%d. %s', i, SET(i).ImageViewPlane)};
end;
set(gui.handles.selectbox, 'string', selBoxCell);

% Set the possible methods to choose from in making the spline
gui.tMethodsCell = {'uniform', 'centripetal', 'chordlength', 'affineInvariant', 'random'};
set(gui.handles.tMethod, 'string', gui.tMethodsCell);

% Set the possible methods to choose from in making the spline
% gui.tMethodsCell = {'uniform', 'centripetal', 'chordlength', 'affineInvariant', 'random'};
% set(gui.handles.tMethod, 'string', gui.tMethodsCell);
% 
% % Set the possible export data types the user can export to creo
% gui.exportCreoCell = {'splines', 'projected'};
% set(gui.handles.exportCreo, 'string', gui.exportCreoCell);

if ispc
    gui.homeFolder = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    gui.homeFolder = getenv('HOME');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nAverages defines how may times repeated average should be run on the
% centroids. Has to be a non-negative integer (0,1,2,...) 0 is to not have
% repeated average.
%
% aaMethod defines the digital filter method to be used. (aa stands for
% anti-aliasing)
%
% convexTolerance defines at which threshold the curves be ignored. The
% smaller it is the more convex the curves has to be. Choose it to be
% positive
%
% tMethod is the parametrization method for the B-spline creation.
%
% sliceLength is the minimum length of consecutive curves satisfying the
% convexTolerance to use in finding the centerline.
%
% If the code doesn't work ('singular matrix') then change nAverages,
% aaMethod or tMethod
gui.nAverages = 1;
gui.aaMethod = {[1 2 1]}; % Choose between [1 2 1], [1 1 1], [1 1 1 1 1] or any combination thereof
gui.convexTolerance = 1.2e-2;
gui.circleTolerance = 1e-1;
gui.cosineTolerance = 4e-1;
gui.distanceTolerance = 1.5;
gui.tMethod = 'chordlength'; %choose between 'uniform', 'centripetal', 'chordlength', 'affineInvariant', 'random'
gui.sliceLength = 4;
gui.exportPath = fullfile(gui.homeFolder, getfield(dir(fullfile(gui.homeFolder,'*ocument*')), 'name'), 'creofile.ibl');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------
function selectbox_Callback
%------------------------
%Listback callback 

global DATA
gui = DATA.GUI.DVESSELGUI;
% Initialize everything
gui.goodCell = {};
gui.goodHandleCell = {};
gui.slicesCell = {};
gui.centroidHandleCell = {};
gui.splinesCell = {};
gui.selectedImgs = get(gui.handles.selectbox,'value');
gui.sizeSelected = size(gui.selectedImgs);

set(gui.handles.smoothbox, 'value', 1);
set(gui.handles.smoothbox, 'string', gui.smoothObjectsCell{1});

% Store the data in the gui
gui.data = savedata(gui.selectedImgs);

clearplot;
colors = {':r', ':g', ':m', ':c'};
hold on;
for i = 1:max(gui.sizeSelected)
    tmpData = gui.data{i}.rlapfh();
    for j = 1:size(tmpData, 4)
        plotdata(tmpData(:, :, :, j), colors{j});
    end;
end;
hold off;
if 1 == 0
    dat = gui.data{1}.rlapfh();
    save ../Documents/tst.mat dat
end;
%Turn on surface option
set(gui.handles.surf, 'Enable', 'on');

%------------------------------
function optionsbutton_Callback %#ok<DEFNU>
%------------------------------
% Choose all the options for the method
global DATA
gui = DATA.GUI.DVESSELGUI;

promptOpts = {'Convex error:', 'Circular error:', 'Cosine error:', ...
    '# slices allowed:', 'Closeness of slices allowed:', 'Smoothing method:', ...
    }; %'Spline method used'};
title = 'Options';
defaultOpts = {num2str(gui.convexTolerance), num2str(gui.circleTolerance), ...
    num2str(gui.cosineTolerance), num2str(gui.sliceLength), num2str(gui.distanceTolerance), ...
    sprintf('[%s]',num2str(gui.aaMethod{1}))};
answer = inputdlg(promptOpts, title, 1, defaultOpts);
if ~isempty(answer)
    gui.convexTolerance = str2num(answer{1});
    gui.circleTolerance = str2num(answer{2});
    gui.cosineTolerance = str2num(answer{3});
    gui.sliceLength = str2num(answer{4});
    gui.distanceTolerance = str2num(answer{5});
    gui.aaMethod{1} = str2num(answer{6}); %#ok<*ST2NM>
end;

%----------------------------
function exportpath_Callback%#ok<DEFNU>
%----------------------------
% Not used
global DATA
gui = DATA.GUI.DVESSELGUI;
gui.exportPath = get(gui.handles.exportPath,'string');

%--------------------------
function smoothbox_Callback
%--------------------------
% Smooth the data
global DATA
gui = DATA.GUI.DVESSELGUI;
smoothselect = get(gui.handles.smoothbox,'value');
switch smoothselect
    case 1
        % Smooth all data
        for i = 1:max(size(gui.data))
            %%%%%%%%%%%%%%%%%
            % IMPORTANT! See if this is good
            gui.data{i}.curveaveraging(gui.aaMethod, gui.nAverages);
            % It might throw away important data?
            %%%%%%%%%%%%%%%%%
        end;
        clearplot;
        colors = {':r', ':g', ':m', ':c'};
        hold on;
        for i = 1:max(gui.sizeSelected)
            tmpData = gui.data{i}.rlapfh();
            for j = 1:size(tmpData, 4)
                plotdata(tmpData(:, :, :, j), colors{j});
            end;
        end;
        hold off;
    case 2
        % Smooth only good data
        for i = 1:max(size(gui.goodCell))
            %%%%%%%%%%%%%%%%%
            % IMPORTANT! See if this is good
            gui.goodCell{i}.curveaveraging(gui.aaMethod, gui.nAverages);
            % It might throw away important data?
            %%%%%%%%%%%%%%%%%
        end;
        try
            gui.goodHandleCell;
        catch
            gui.goodHandleCell = cell(size(gui.goodCell));
        end;
        hold on;
        for i = 1:max(size(gui.data))
            try
                delete(gui.goodHandleCell{i});
            catch
            end;
            gui.goodHandleCell{i} = plotdata(gui.goodCell{i}.rlapfh(), 'k');
        end;
        hold off;
    case 3
        % Smooth along splines
        for i = 1:max(size(gui.slicesCell))
            % Create smoothed curves from the good curves
            %%%%%%%%%%%%%%%%%
            % OBS se om detta �r bra!
            gui.slicesCell{i}.vert_averaging(gui.aaMethod, gui.nAverages);
            % Sl�nger bort viktig data kanske?
            %%%%%%%%%%%%%%%%%
        end;
        tmethod_Callback;
end;

%------------------------
function tmethod_Callback %#ok<DEFNU>
%------------------------
global DATA
gui = DATA.GUI.DVESSELGUI;
gui.tMethod = gui.tMethodsCell{get(gui.handles.tMethod,'Value')};

if isempty(gui.slicesCell)
    gui.slicesCell = {};
    for i = 1:max(gui.sizeSelected)
        tmpCell = gui.goodCell{i}.groupclosest('rlapfh', gui.distanceTolerance, gui.sliceLength);
        emptyCells = cellfun('isempty', tmpCell); 
        tmpCell(emptyCells) = [];
        for j = 1:max(size(tmpCell))
            gui.slicesCell{end + 1} = dvessel.curvesclass(gui.goodCell{i}.dataNumber, tmpCell{j});
        end;
    end;
end;

gui.splinesCell = cell(size(gui.slicesCell));
for i = 1:max(size(gui.slicesCell))
    gui.splinesCell{i} = dvessel.splineclass(gui.slicesCell{i}.centroids('rlapfh'), 3, [0 1], gui.tMethod);
end;
try
    gui.splinesHandleCell;
catch
    gui.splinesHandleCell = cell(size(gui.splinesCell));
end;
hold on;
for i = 1:max(size(gui.splinesCell))
    try
        delete(gui.splinesHandleCell{i});
    catch
    end;
    try
        delete(gui.slicesHandleCell{i});
    catch
    end;
    try
        delete(gui.goodHandleCell{i});
    catch
    end;
    try
        delete(gui.projectedHandleCell{i});
    catch
    end;
    % Create the points to be plotted
    points = fnplt(gui.splinesCell{i}.spline());
    % Plot the points and save the handles.
    gui.slicesHandleCell{i} = plotdata(gui.slicesCell{i}.rlapfh(), 'k');
    gui.splinesHandleCell{i} = plot3(points(1, :), points(2, :), points(3, :), 'r');
end;
% Allow for smoothing along splines
set(gui.handles.smoothbox, 'value', 3);
set(gui.handles.smoothbox, 'string', gui.smoothObjectsCell(1:3));

%-------------------------
function findgood_Callback %#ok<DEFNU>
%-------------------------
global DATA
gui = DATA.GUI.DVESSELGUI;

% Check if the datapoints have been imported or not. Should not trigger as
% the data is imported in the initialize phase.
try
    gui.data;
catch
    warning('Images not imported. Images imported automatically');
    selectbox_Callback;
end;

% create a cell containing the curves that both satisfies the convex error
% and the circle error.
gui.goodCell = cell(size(gui.data));
for i = 1:max(size(gui.data))
    gui.goodCell{i} = dvessel.dataclass(gui.data{i}.dataNumber, gui.data{i}.getcurvesbymatr('xyz', (gui.data{i}.compc2conv('mm') > 1 - gui.convexTolerance).*(gui.data{i}.compc2circ('mm') > 1 - gui.circleTolerance)));
end;

% Plot the good slices in black and store the handles in slicesHandleCell
try
    gui.goodHandleCell;
catch
    gui.goodHandleCell = cell(size(gui.goodCell));
end;
hold on;
for i = 1:max(size(gui.data))
    try
        delete(gui.goodHandleCell{i});
    catch
    end;
    gui.goodHandleCell{i} = plotdata(gui.goodCell{i}.rlapfh(), 'k');
end;
hold off;
%Allow for smoothing on all data and good data only.
set(gui.handles.smoothbox, 'value', 2);
set(gui.handles.smoothbox, 'string', gui.smoothObjectsCell(1:2));

%---------------------------
function exportcreo_Callback%#ok<DEFNU>
%---------------------------
global DATA
gui = DATA.GUI.DVESSELGUI;
switch get(gui.handles.exportCreo,'value')
    case 1
        gui.exportData = gui.splinesCell;
        tmpDat = cell(size(gui.exportData));
        for i = 1:max(size(tmpDat))
            tmpDat{i} = gui.exportData{i}.interppts;
        end;
        gui.exportData = tmpDat;
    case 2
        gui.exportData = gui.projectedCell;
    otherwise
        error('Not implemented yet');
end;
save2creo(gui.exportPath, gui.exportData);

%-----------------------
function rotate_Callback%#ok<DEFNU>
%-----------------------
% Button to rotate.
global DATA
gui = DATA.GUI.DVESSELGUI;
rotate3d(gui.handles.axes1);

%-----------------------
function zoomon_Callback%#ok<DEFNU>
%-----------------------
% Is a bit difficult to use
global DATA
gui = DATA.GUI.DVESSELGUI;
h = zoom(gui.handles.axes1);
h.Enable = 'on';

%-----------------------
function zoomoff_Callback%#ok<DEFNU>
%-----------------------
% A bit difficult to use
global DATA
gui = DATA.GUI.DVESSELGUI;
h = zoom(gui.handles.axes1);
h.Enable = 'off';

%---------------------
function save_Callback%#ok<DEFNU>
%---------------------
global DATA
gui = DATA.GUI.DVESSELGUI;
% Opens up a save window where the user can specify which data so save and
% to where. Is saved as a file readable by Creo.

title = 'Export data as Creo file';
datatypes = {'*.ibl', 'Center splines (*.ibl)'; '*.ibl', 'Projected outlines (*.ibl)'; '*.ibl', 'Non projected data (*.ibl)'};
[fileName, path, index] = uiputfile(datatypes, title, gui.exportPath);
% If the save is canceled, exit without doing anything.
if index == 0
    return;
end;
gui.exportPath = fullfile(path, fileName);
switch index
    case 1
        gui.exportData = gui.splinesCell;
        tmpDat = cell(size(gui.exportData));
        for i = 1:max(size(tmpDat))
            tmpDat{i} = gui.exportData{i}.interppts;
        end;
        gui.exportData = tmpDat;
    case 2
        gui.exportData = gui.projectedCell;
    case 3
        gui.exportData = cell(1,1);
        size(gui.exportData)
        for i = 1:size(gui.data, 2)
            gui.exportData{i} = gui.data{i}.rlapfh();
        end;
        tmpDat = cell(1, 4*size(gui.exportData, 2));
        for i = 1:size(gui.exportData, 2)
            for j = 1:4
                tmpDat{4*(i - 1) + j} = gui.exportData{i}(:, :, :, j);
            end;
        end;
        gui.exportData = tmpDat;
    otherwise
        error('Not implemented yet');
end;
save2creo(gui.exportPath, gui.exportData);

%------------------------
function connect_Callback%#ok<DEFNU>
%------------------------
% Connect two splines together to form one long spline.
global DATA
gui = DATA.GUI.DVESSELGUI;
newCrv = {};
newSpl = {};
averArea = [];
oldCrvIndex = [];
for i = 1:max(size(gui.slicesCell))
    for j = i + 1:max(size(gui.slicesCell))
        
        [cosComp, connectedCrvs, dataNumber] = curvesconnect(gui.slicesCell{i}, gui.slicesCell{j});
        
        if all(cosComp >= 1 - gui.cosineTolerance)
            newCrv{end + 1} = dvessel.curvesclass(dataNumber, connectedCrvs);
            newSpl{end + 1} = dvessel.splineclass(newCrv{end}.centroids('rlapfh'), 3, [0 1], gui.tMethod);
            averArea(end + 1) = newCrv{end}.averarea('mm');
            oldCrvIndex(end + 1, :) = [i j];
        end;
    end;
end;
[~, ind] = sort(averArea, 'descend');
if isempty(ind)
    error('The program cannot find two center lines to connect to a main vessel.');
end;
gui.splinesCell(oldCrvIndex(ind(1), :)) = [];
gui.slicesCell(oldCrvIndex(ind(1), :)) = [];

gui.splinesCell{end + 1} = newSpl{ind(1)};
gui.slicesCell{end + 1} = newCrv{ind(1)};

tmethod_Callback;
% Surface plot does not work with connected vessels
set(gui.handles.surf, 'Enable', 'off');
% Smoothing along spline does not work with connected vessels
set(gui.handles.smoothbox, 'value', 2);
set(gui.handles.smoothbox, 'string', gui.smoothObjectsCell(1:2));

%------------------------------
function projectcurves_Callback%#ok<DEFNU>
%------------------------------
global DATA
gui = DATA.GUI.DVESSELGUI;
gui.projectedCell = cell(size(gui.splinesCell));
for i = 1:max(size(gui.splinesCell))
    sliceMatr = gui.slicesCell{i}.rlapfh();
    derivMatr = gui.splinesCell{i}.nderivs(1, 'tau');
    for k = 1:size(gui.slicesCell{i}.xyz(), 3)
        gui.projectedCell{i}(:, :, k) = pointsprojection(sliceMatr(:, :, k), derivMatr(k, :, 2), derivMatr(k, :, 1));
    end;
end;
try
    gui.projectedHandleCell;
catch
    gui.projectedHandleCell = cell(size(gui.projectedCell));
end;
hold on;
for i = 1:max(size(gui.projectedCell))
    try
        delete(gui.projectedHandleCell{i});
    catch
    end;
    try
        delete(gui.slicesHandleCell{i});
    catch
    end;

    gui.projectedHandleCell{i} = plotdata(gui.projectedCell{i}, 'k');
end;
hold off;

%----------------------------
function plotsurface_Callback %#ok<DEFNU>
%----------------------------
global DATA
gui = DATA.GUI.DVESSELGUI;
try
    gui.projectedCell;
catch
    projectedcurves_Callback;
end;

try
    gui.surfHandleCell;
catch
    gui.surfHandleCell = cell(size(gui.projectedCell));
end;
hold on;
for i = 1:max(size(gui.projectedCell))
    try
        delete(gui.surfHandleCell{i});
    catch
    end;
    % plot the surface along the smoothed out projected curves
    XQ = permute(gui.projectedCell{i}(:,1,:),[3 1 2]);
    YQ = permute(gui.projectedCell{i}(:,2,:),[3 1 2]);
    ZQ = permute(gui.projectedCell{i}(:,3,:),[3 1 2]);
    gui.surfHandleCell{i} = surf(XQ, YQ, ZQ, ones(size(ZQ)), 'faceAlpha', 0.5, 'faceColor', 'b');
end;
hold off;

%----------------------------------%
%----------------------------------%
%%%%%%    Helper Functions    %%%%%%
%----------------------------------%
%----------------------------------%

%-----------------
function clearplot
%-----------------
hold off;
plot3(0,0,0);

%---------------------------------------------
function plothandle = plotdata(data, varargin)
%---------------------------------------------
if iscell(data)
    cLen = size(data, 2);
    for i = 1:cLen
        plotdata(data{i}, varargin{:});
    end;
elseif isfloat(data)
    XQ = permute(data(:,1,:),[3 1 2]);
    YQ = permute(data(:,2,:),[3 1 2]);
    ZQ = permute(data(:,3,:),[3 1 2]);
    plothandle = plot3(XQ', YQ', ZQ', varargin{:});
else
    error('The data is neither a matrix or a cell structure. Plotting doesnt work for other structures.');
end;
axis equal;

%-------------------------------------
function dataCell = savedata(dataList)
%-------------------------------------
global SET
dataCell = cell(size(dataList));
for i = 1:length(dataCell)
    % Create lvEndo
    if max(size(SET(dataList(i)).EndoX)) > 0
        % Add x and y coordinates
        lvEndo = cat(2, SET(dataList(i)).EndoX(:,1,:), SET(dataList(i)).EndoY(:,1,:));
        for j = 1:SET(dataList(i)).ZSize
            % Add z coordinates
            lvEndo(:, 3, j) = j*ones(size(SET(dataList(i)).EndoX,1),1);
        end;
    else
        lvEndo = NaN(80,3,SET(dataList(i)).ZSize);
    end;
    % Create lvEpi
    if max(size(SET(dataList(i)).EpiX)) > 0
        % Add x and y coordinates
        lvEpi = cat(2, SET(dataList(i)).EpiX(:,1,:), SET(dataList(i)).EpiY(:,1,:));
        for j = 1:SET(dataList(i)).ZSize
            % Add z coordinates
            lvEpi(:, 3, j) = j*ones(size(SET(dataList(i)).EpiX,1),1);
        end;
    else
        lvEpi = NaN(80,3,SET(dataList(i)).ZSize);
    end;
    % Create rvEndo
    if max(size(SET(dataList(i)).RVEndoX)) > 0
        % Add x and y coordinates
        rvEndo = cat(2, SET(dataList(i)).RVEndoX(:,1,:), SET(dataList(i)).RVEndoY(:,1,:));
        for j = 1:SET(dataList(i)).ZSize
            % Add z coordinates
            rvEndo(:, 3, j) = j*ones(size(SET(dataList(i)).RVEndoX,1),1);
        end;
    else
        rvEndo = NaN(80,3,SET(dataList(i)).ZSize);
    end;
    % Create rvEpi
    if max(size(SET(dataList(i)).RVEpiX)) > 0
        % Add x and y coordinates
        rvEpi = cat(2, SET(dataList(i)).RVEpiX(:,1,:), SET(dataList(i)).RVEpiY(:,1,:));
        for j = 1:SET(dataList(i)).ZSize
            % Add z coordinates
            rvEpi(:, 3, j) = j*ones(size(SET(dataList(i)).RVEpiX,1),1);
        end;
    else
        rvEpi = NaN(80,3,SET(dataList(i)).ZSize);
    end;
    dataMatr = cat(4, lvEndo, lvEpi, rvEndo, rvEpi);
    dataCell{i} = dvessel.dataclass(dataList(i), dataMatr);
end;

%-------------------------------------------
function save2creo(fileName, data, varargin)
%-------------------------------------------

curveNo = 0;
fileID = fopen(fileName,'w');
fprintf(fileID,'%s\r\n','open');
fprintf(fileID,'%s\r\n','arclength');

for j = 1:size(data, 2)
    for i = 1:size(data{j}, 3)
        if ~any(any(isnan(data{j}(:, :, i))))
            curveNo = curveNo + 1;
            fprintf(fileID,'%s %i\r\n','begin section !', j);
            fprintf(fileID,'%s %i\r\n','begin curve !',curveNo);
            fprintf(fileID,'%f %f %f \r\n', data{j}(:, :, i)');
            fprintf(fileID,'\r\n','');
        end;
    end;
end;
fclose(fileID);

%-------------------------------------------------------
function [cosComp, dataCrvs, dataNr] = curvesconnect(crv1, crv2)
%-------------------------------------------------------
% Try to connect crv1 and crv2
global DATA
gui = DATA.GUI.DVESSELGUI;

pts1 = crv1.centroids('rlapfh');
pts2 = crv2.centroids('rlapfh');
xyz1 = crv1.xyz();
xyz2 = crv2.xyz();

% Find the endpoints of the crvs that are closest to eachother
ptsToCompare = [pts1([1 end], :); pts2([1 end], :)];
dstBetweenPts = pdist(ptsToCompare);
triangDst = tril(squareform(dstBetweenPts));
smllst = min(dstBetweenPts(2:end - 1));
[dstIndx(:, 1), dstIndx(:, 2)] = find(ismember(triangDst, smllst));
dstIndx = sort(dstIndx);

% Four possible pairings for closest endpoints
if dstIndx(1,1) == 1
    if dstIndx(1,2) == 3
        % Combine the endpoints and each of their closest neighbour to do 
        % cosine comparison.
        cosineCompPts = [pts1([2 1], :); pts2([1 2], :)];
        % Combine curves and dataNr in correct order
        dataCrvs = cat(3, xyz1(:, :, end:-1:1), xyz2);
        dataNr = [crv1.dataNumber crv2.dataNumber];
    else
        % Combine the endpoints and each of their closest neighbour to do 
        % cosine comparison.
        cosineCompPts = [pts2([end-1 end], :); pts1([1 2], :)];
        % Combine curves and dataNr in correct order
        dataCrvs = cat(3, xyz2, xyz1);
        dataNr = [crv2.dataNumber crv1.dataNumber];
    end;
else
    if dstIndx(1,2) == 3
        % Combine the endpoints and each of their closest neighbour to do 
        % cosine comparison.
        cosineCompPts = [pts1([end-1 end], :); pts2([1 2], :)];
        % Combine curves and dataNr in correct order
        dataCrvs = cat(3, xyz1, xyz2);
        dataNr = [crv1.dataNumber crv2.dataNumber];
    else
        % Combine the endpoints and each of their closest neighbour to do 
        % cosine comparison.
        cosineCompPts = [pts1([end-1 end], :); pts2([end end-1], :)];
        % Combine curves and dataNr in correct order
        dataCrvs = cat(3, xyz1, xyz2(:, :, end:-1:1));
        dataNr = [crv1.dataNumber crv2.dataNumber];
    end;
end;
% Do cosine comparison.
cosComp = gui.data{1}.pointsimilarity(cosineCompPts);

%-------------------------------------------------------------------
function projectedpts = pointsprojection(points, vector, planepoint)
%-------------------------------------------------------------------
% Function to project the points in a plane defined by a point (planepoint)
% and vector normal to the plane
nvector = vector/norm(vector); %Normalize the vector.
projectedpts = zeros(size(points));
for i= 1:length(points)
    projectedpts(i,:) = points(i,:) - dot(points(i,:) - planepoint, nvector)*nvector;
end;

function indx = findclosest(iCell, jCell)
shortest = NaN(length(iCell), length(jCell));
for i = 1:length(iCell)
    ipts = fnval(iCell{i}.spline(), 0:1e-3:1);
    for j = 1:length(jCell)
        jpts = fnval(jCell{j}.spline(), 0:1e-3:1);
        [~, dst] = knnsearch(ipts, jpts);
        shortest(i, j) = min(dst);
    end;
end;
[~, in] = min(shortest(:));
[ind, jnd] = ind2sub(size(shortest),in);
indx = [ind, jnd];
