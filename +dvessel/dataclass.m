classdef dataclass < dvessel.dataparentclass
    %curvesclass Main class for the data.
    %   Detailed explanation goes here
    %Andreas S�derlund
    
    properties
        %For which data series it comes
        dataNumber;
        
        %The n x 3 x m x 4 array containing the data
        dataMatrix;
    end
    
    methods
        %-------------------------------------------
        function obj = dataclass(dataNr, dataMatr)
        %-------------------------------------------
            global SET
            
            % Check if dataNr is valid
            if (dataNr > max(size(SET))) | (dataNr < 1)
                error('dataNumber = %d is an invalid data series. must be between 1 and %d.', dataNr, max(size(SET)));
            end;
            
            if (size(dataMatr, 2) ~= 3) | (size(dataMatr, 4) ~= 4)
                error('dataMatr must be on the form n x 3 x m x 4, the input on the form [%d x %d x %d x %d] is not valid', size(dataMatr));
            end;
            
            %Assign the data
            obj.dataNumber = dataNr;
            obj.dataMatrix = dataMatr;
        end;
        
        %---------------------------------------
        function obj = set.dataMatrix(obj, data)
        %---------------------------------------
            obj.dataMatrix = data;
        end;
        
        %-----------------------------
        function res = resolution(obj)
        %-----------------------------
            % Returns the resolution of the slices in the form
            % [xr;
            %  yr;
            %  zr]
            global SET
            res = [SET(obj.dataNumber).ResolutionX; SET(obj.dataNumber).ResolutionY; ([SET(obj.dataNumber).SliceThickness] + [SET(obj.dataNumber).SliceGap])];
        end;
        
        %----------------------
        function pts = xyz(obj)
        %----------------------
            % returns the data in the original segment (xyz) format. The
            % data is assumed to already be in this format when stored.
            pts = obj.dataMatrix;
        end;
        
        %---------------------
        function pts = mm(obj)
        %---------------------
            % returns the data in the mm format.
            global SET
            pts = obj.xyz();
            res = obj.resolution();
            % For every slice
            for e = 1:size(pts, 4)
                for j = 1:size(pts, 3)
                    % For every (x,y,z) point in the slice
                    for k = 1:size(pts, 1)
                        % Do pointwise multiplication with the resolution
                        % vector.
                        pts(k, :, j, e) = pts(k, :, j, e) .* res';
                    end;
                end;
            end;
        end;
        
        %-------------------------
        function pts = rlapfh(obj)
        %-------------------------
            % returns the data in the rlapfh format.
            global SET %#ok<*NUSED>
            pts = obj.xyz();
            
            for i = 1:size(pts, 4)
                % For every slice
                for j = 1:size(pts, 3)
                    % Transform slice points using built in function.
                    pts(:, :, j, i) = calcfunctions('xyz2rlapfh', obj.dataNumber, pts(:, 1, j, i), pts(:, 2, j, i), pts(:, 3, j, i));
                end;
            end;
        end;
        
        %-----------------------------------
        function xyz = mm2xyz(obj, dataMatr)
        %-----------------------------------
            xyz = ones(size(dataMatr));
            res = obj.resolution();
            for i = 1:size(dataMatr, 3)
                for j = 1:length(res)
                    xyz(:, j, i) = dataMatr(:, j, i)/res(j);
                end;
            end;
        end;
        
        %------------------------------------------------------
        function pts = centroids(obj, representation, varargin)
        %------------------------------------------------------
            switch representation
                case 'xyz'
                    data = obj.xyz();
                case 'mm'
                    data = obj.mm();
                case 'rlapfh'
                    data = obj.rlapfh();
                otherwise
                    types = 'xyz, mm or rlapfh';
                    error('The representation: %s is not supported. choose between %s', dataType, types);
            end;
            if max(size(varargin)) < 1
                pts = ones(size(obj.dataMatrix, 3), size(obj.dataMatrix, 2), 4);
                for k = 1:size(pts, 3)
                    for i = 1:size(obj.dataMatrix, 3)
                        len = length(data(:, 1, i, k));
                        pts(i, :, k) = [sum(data(:, 1, i, k)) sum(data(:, 2, i, k)) sum(data(:, 3, i, k))]/len;
                    end;
                end;
            elseif max(size(varargin)) < 2
                pts = ones(size(obj.dataMatrix, 3), size(obj.dataMatrix, 2), max(size(varargin{1})));
                for k = length(varargin{1})
                    for i = 1:size(obj.dataMatrix, 3)
                        len = length(data(:, 1, i, varargin{1}(k)));
                        pts(i, :, k) = [sum(data(:, 1, i, varargin{1}(k))) sum(data(:, 2, i, varargin{1}(k))) sum(data(:, 3, i, varargin{1}(k)))]/len;
                    end;
                end;
            elseif max(size(varargin)) < 3
                pts = ones(max(size(varargin{2})), size(obj.dataMatrix, 2), max(size(varargin{1})));
                for k = length(varargin{1})
                    for i = length(varargin{2})
                        len = length(data(:, 1, varargin{2}(i), varargin{1}(k)));
                        pts(i, :, k) = [sum(data(:, 1, varargin{2}(i), varargin{1}(k))) sum(data(:, 2, varargin{2}(i), varargin{1}(k))) sum(data(:, 3, varargin{2}(i), varargin{1}(k)))]/len;
                    end;
                end;
            else
                error('error');
            end;
        end;
        
        %--------------------------------------------------------
        function areas = curvearea(obj, representation, varargin)
        %--------------------------------------------------------
            switch representation
                case 'xyz'
                    data = obj.xyz();
                case 'mm'
                    data = obj.mm();
                otherwise
                    types = 'xyz or mm';
                    error('The representation: %s is not supported. choose between %s', dataType, types);
            end;
            if max(size(varargin)) < 1
                areas = ones(size(obj.dataMatrix, 3), 4);
                for k = 1:size(areas, 2)
                    for i = 1:size(obj.dataMatrix, 3)
                        areas(i, k) = polyarea(data(:, 1, i, k), data(:, 2, i, k));
                    end;
                end;
            elseif max(size(varargin)) < 2
                areas = ones(size(obj.dataMatrix, 3), max(size(varargin{1})));
                for k = length(varargin{1})
                    for i = 1:size(obj.dataMatrix, 3)
                        areas(i, k) = polyarea(data(:, 1, i, varargin{1}(k)), data(:, 2, i, varargin{1}(k)));
                    end;
                end;
            elseif max(size(varargin)) < 3
                areas = ones(max(size(varargin{2})), max(size(varargin{1})));
                for k = length(varargin{1})
                    for i = length(varargin{2})
                        areas(i, k) = polyarea(data(:, 1, varargin{2}(i), varargin{1}(k)), data(:, 2, varargin{2}(i), varargin{1}(k)));
                    end;
                end;
            else
                error('error');
            end;
        end;
        
        %----------------------------------------------------------------
        function circ = curvecircumference(obj, representation, varargin)
        %----------------------------------------------------------------
            switch representation
                case 'xyz'
                    data = obj.xyz();
                case 'mm'
                    data = obj.mm();
                case 'rlapfh'
                    data = obj.rlapfh();
                otherwise
                    types = 'xyz, mm or rlapfh';
                    error('The representation: %s is not supported. choose between %s', dataType, types);
            end;
            if max(size(varargin)) < 1
                circ = ones(size(obj.dataMatrix, 3), 4);
                for k = 1:size(circ, 2)
                    for i = 1:size(obj.dataMatrix, 3)
                        circ(i, k) = obj.curvelength(data(:, :, i, k));
                    end;
                end;
            elseif max(size(varargin)) < 2
                circ = ones(size(obj.dataMatrix, 3), max(size(varargin{1})));
                for k = length(varargin{1})
                    for i = 1:size(obj.dataMatrix, 3)
                        circ(i, k) = obj.curvelength(data(:, :, i, varargin{1}(k)));
                    end;
                end;
            elseif max(size(varargin)) < 3
                circ = ones(max(size(varargin{2})), max(size(varargin{1})));
                for k = length(varargin{1})
                    for i = length(varargin{2})
                        circ(i, k) = obj.curvelength(data(:, :, varargin{2}(i), varargin{1}(k)));
                    end;
                end;
            else
                error('error');
            end;
        end;
        
        %----------------------------------------------------
        function boolvar = vert_averaging(obj, methodCell, n)
        %----------------------------------------------------
            data = obj.mm();
            for j = 1:size(data, 1)
                obj.dataMatrix(j, :, :) = obj.mm2xyz(permute(obj.repeatedaaaverage(permute(data(j, :, :), [3 2 1]), methodCell, n), [3 2 1]));
            end;
            boolvar = true;
        end;
        
        %----------------------------------------------------
        function boolvar = curveaveraging(obj, methodCell, n)
        %----------------------------------------------------
            data = obj.mm();
            % Go through each outline color endolv, endorv, epilv, epirv
            for i = 1:size(data, 4)
                % Go through each curve.
                for j = 1:size(data, 3)
                    obj.dataMatrix(:, :, j, i) = obj.mm2xyz(obj.repeatedaaaverage(data(:, :, j, i), methodCell, n));
                end;
            end;
            boolvar = true;
        end;
        
        %--------------------------------------------------------
        function pts = getcurvesbymatr(obj, representation, matr)
        %--------------------------------------------------------
            switch representation
                case 'xyz'
                    data = obj.xyz();
                case 'mm'
                    data = obj.mm();
                case 'rlapfh'
                    data = obj.rlapfh();
                otherwise
                    types = 'xyz, mm or rlapfh';
                    error('The representation: %s is not supported. choose between %s', dataType, types);
            end;
            if size(data, 3) ~= size(matr, 1) | size(data, 4) ~= size(matr, 2)
                error('matr has to have the correct size. Has to be %d x %d', size(data, 3), size(data, 4));
            end;
            matr = double(matr); matr(matr == 0) = NaN;
            matr = reshape(kron(matr, ones(size(data, 1)*size(data, 2),1)), size(data));
            pts = matr.*data;
        end;
        
        %--------------------------------------------------------
        function newPoints = aafilter(~, oldPoints, methodVector)
        %--------------------------------------------------------
            % A filter that acts as a way of averaging the points in oldPoints.
            % methodVector is the method, can be ]1 2 1[ or ]1 1 1[ for example
            methodLen = length(methodVector);
            ptsLen = size(oldPoints, 1);

            % If the methodVector is of even length, make a convolution to get a method
            % that cancels the same frequencies but is of odd length.
            if mod(methodLen, 2) == 0
                methodVector = conv(methodVector, fliplr(methodVector));
                methodLen = length(methodVector);
            end;

            methodMat = repmat(methodVector', [1 size(oldPoints, 2)]);
            methodSiz = (methodLen - 1)/2;

            newPoints = oldPoints;
            for i = 1 + methodSiz:ptsLen - methodSiz
                newPoints(i, :) = sum(methodMat.*oldPoints(i - methodSiz:i + methodSiz, :),1)/sum(methodVector);
            end;
        end;

        %---------------------------------------------------------------
        function newPoints = repeatedaaaverage(obj, oldPoints, methodCell, n)
        %---------------------------------------------------------------
            newPoints = oldPoints;
            % Repeat the methods in methodCell n times
            for i = 1:n
                % Go through each method specified in the cell.
                for j = methodCell
                    newPoints = obj.aafilter(newPoints, j{1});
                end;
            end;
        end;
        
        %-------------------------------------------------------
        function cosTheta = pointsimilaritynew(obj, representation, xVals, yVals)
        %-------------------------------------------------------
            centds = obj.centroids(representation, xVals, yVals);
            if any(isnan(centds))
                error('Some curves you wish to access contain at least one NaN value.');
            end;
            cosTheta = obj.pointsimilarity(centds);
            % cosTheta = pointsimilarity@dvessel.dataparentclass(centds);
        end;
        
        %---------------------------------------------------
        function ratioMatr = compc2conv(obj, representation)
        %---------------------------------------------------
            ratioMatr = NaN(size(obj.dataMatrix, 3), size(obj.dataMatrix, 4));
            switch representation
                case 'xyz'
                    data = obj.xyz();
                case 'mm'
                    data = obj.mm();
                otherwise
                    types = 'xyz or mm';
                    error('The representation: %s is not supported. choose between %s', dataType, types);
            end;
            for i = 1:size(ratioMatr, 2)
                for j = 1:size(ratioMatr, 1)
                    if any(any(isnan(data(:, :, j, i))))
                        ratioMatr(j, i) = NaN;
                    else
                        ratioMatr(j, i) = obj.compchullper(data(:, :, j, i));
                    end;
                end;
            end;
        end;
        
        %---------------------------------------------------
        function ratioMatr = compc2circ(obj, representation)
        %---------------------------------------------------
            circC = 2*sqrt(pi)*sqrt(obj.curvearea(representation));
            ratioMatr = circC./obj.curvecircumference(representation);
        end;
        
        %--------------------------------------------------------
        function centr = centroidsbyxy(obj, representation, X, Y)
        %--------------------------------------------------------
            tmpMatr = obj.centroids(representation);
            centr = ones(max(size(X)), 3);
            for i = 1:max(size(X))
                centr(i, :) = tmpMatr(Y(i), :, X(i));
            end;
        end;
        
        %----------------------------------------------------------------------
        function closestCell = groupclosest(obj, representation, lenTol, nSize)
        %----------------------------------------------------------------------
            % Creates a cell with a number of entries where the entries
            % contains at least nSize curves from obj.xyz() where each
            % curve center (centroid) is at least
            % lenTol*obj.resolution(3,1) or shorter away from at least one
            % other centroid in the cell entry.
            centroids = obj.centroids(representation);
            centroids = cat(1, centroids(:, :, 1), centroids(:, :, 2), centroids(:, :, 3), centroids(:, :, 4));
            % create a list of every piecewise distance between the points
            % in centroids
            B = pdist(centroids);
            % turn it into lower triangular matrix
            D = tril(squareform(B));
            % res(3, 1) is the minimum distance between points
            res = obj.resolution()*lenTol;
            % remove distances that are too large
            D = D.*(D < 1.5*res(3, 1));
            % Find the indexes for the distances. Here centroids(Y(k, 1), :) contains the
            % starting point and centroids(Y(k, 2), :) contain the ending
            % point for a distance that is short enough.
            [Y(:, 1), Y(:, 2)] = find(ismember(D, B));
            % group them into sets where the points are close enough
            % together
            lstCell = {Y(1, :)};
            for i = 2:max(size(Y))
                % if the point is not close enough, a new cell entry is
                % needed.
                newCellNeeded = true;
                % for each pair of points
                for j = 1:max(size(lstCell))
                    % if starting or ending entry is in it then the other
                    % point is close enough. Add both, duplicates will be
                    % removed later.
                    if any(Y(i, 1) == lstCell{j}) | any(Y(i, 2) == lstCell{j})
                        newCellNeeded = false;
                        lstCell{j} = [lstCell{j} Y(i, :)];
                    end;
                end;
                % Create a new cell entry if needed
                if newCellNeeded
                    lstCell{end + 1} = Y(i, :);
                end;
            end;
            % Cell containing the indexes for the groups that are large enough
            nSizeCell = {};
            for i = 1:max(size(lstCell))
                % Remove duplicate entries
                lstCell{i} = unique(lstCell{i});
                % If large enough, add to new cell
                if max(size(lstCell{i})) >= nSize
                    nSizeCell{end + 1} = lstCell{i};
                end;
            end;
            % retrieve the curves
            data = obj.xyz();
            data = cat(3, data(:, :, :, 1), data(:, :, :, 2), data(:, :, :, 3), data(:, :, :, 4));
            %cell containing the curves in xyz format
            closestCell = cell(size(nSizeCell));
            for i = 1:max(size(nSizeCell))
                closestCell{i} = data(:, :, nSizeCell{i});
                % sort the data on height
                [~, I] = sort(closestCell{i}(1, 3, :));
                closestCell{i} = closestCell{i}(:, :, I);
            end;
        end;
        
        %-------------------------------
        function pos = removeNaN(~, values)
        %-------------------------------
            pos = cell(1, size(values, 4));
            for i = 1:size(values, 4)
                pos{i} = values(:, :, :, i);
            end;
            
            for i = 1:size(pos, 2)
                pos{i}(:,:,any(isnan(pos{i}(1,1,:)),2)) = [];
            end;
        end;
        
        %----------------------------------------------
        function curvesCell = saveslices(obj, interval)
        %----------------------------------------------
            % Interval is a cell of the same length as the number of colors
            % (lvepi, lvendo, rvepi, rvendo) selected. It contains a 
            % matrix at each cell that 
            % The data with removed NaN entries
            noNaN = obj.removeNaN(obj.xyz());
            
            curvesCell = cell(size(noNaN));
            for i = 1:max(size(noNaN))
                % How many intervals to save
                cLen = size(interval{i}, 1);
                
                % Divide the positions into cell elements on the intervals
                % defined in interval{i}
                positionCell = cell(1, cLen);
                for j = 1:cLen
                    % Save the data as a curvesclass object.
                    positionCell{j} = dvessel.curvesclass(obj.dataNumber, noNaN{i}(:,:,interval{i}(j, 1):interval{i}(j, 2)));
                end;
                % Store the cell in curvesCell
                curvesCell{i} = positionCell;
            end;
        end;
            
    end
    
end

