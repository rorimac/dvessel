classdef curvesclass < dvessel.dataparentclass
    %curvesclass Class for the individual curves
    %   Detailed explanation goes here
    %Andreas S�derlund
    
    properties
        %For which data series it comes
        dataNumber;
        
        %The n x 3 x m array containing the data
        dataMatrix;
    end
    
    methods
        %-------------------------------------------
        function obj = curvesclass(dataNr, dataMatr)
        %-------------------------------------------
            global SET DATA
            
            % Check if dataNr is valid
            if ( max(dataNr) > max(size(SET))) | (min(dataNr) < 1)
                error('dataNumber = %d is an invalid data series. must be between 1 and %d.', dataNr, max(size(SET)));
            end;
            
            if (size(dataMatr, 2) ~= 3)
                error('dataMatr must be on the form n x 3 x m, the input on the form [%d x %d x %d] is not valid', size(dataMatr));
            end;
            
            %Assign the data
            if numel(dataNr) == 1
                obj.dataNumber = dataNr*ones(1, size(dataMatr, 3));
            else
                obj.dataNumber = dataNr;
            end;
            obj.dataMatrix = dataMatr;
        end;
        
        function obj = set.dataMatrix(obj, data)
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
            res = ones(3, 1, max(size(obj.dataNumber)));
            for i = 1:max(size(obj.dataNumber))
                res(:, :, i) = [SET(obj.dataNumber(i)).ResolutionX; SET(obj.dataNumber(i)).ResolutionY; ([SET(obj.dataNumber(i)).SliceThickness] + [SET(obj.dataNumber(i)).SliceGap])];
            end;
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
            for j = 1:size(pts, 3)
                % For every (x,y,z) point in the slice
                for k = 1:size(pts, 1)
                    % Do pointwise multiplication with the resolution
                    % vector.
                    pts(k, :, j) = pts(k, :, j) .* res(:, :, j)';
                end;
            end;
        end;
        
        %-------------------------
        function pts = rlapfh(obj)
        %-------------------------
            % returns the data in the rlapfh format.
            global SET %#ok<*NUSED>
            pts = obj.xyz();
            
            % For every slice
            for j = 1:size(pts, 3)
                % Transform slice points using built in function.
                pts(:,:,j) = calcfunctions('xyz2rlapfh', obj.dataNumber(j), pts(:, 1, j), pts(:, 2, j), pts(:, 3, j));
            end;
        end;
        
        %------------------------------------
        function deleteslices(obj, indexVect)
        %------------------------------------
            global SET
            % Checking if sorted
            tmpVect = sort(indexVect);
            if ~all(indexVect == tmpVect)
                warning('The vector is not sorted. Sorting performed automatically.');
                indexVect = tmpVect;
            end;
            
            if (indexVect(1) < 1) | (indexVect(end) > max(size(SET)))
                error('The index values are out of bounds. Has to be between 1 and %d.', max(size(SET)));
            end;
            obj.dataMatrix(:,:,indexVect) = [];
        end;
        
        function addslices(obj, newSlices, dataType, varargin)
            switch dataType
                case 'xyz'
                    pts = newSlices;
                case 'mm'
                    pts = mm2xyz(newSlices);
                case 'rlapfh'
                    pts = rlapfh2xyz(newSlices);
                otherwise
                    types = 'xyz, mm or rlapfh';
                    error('The dataType: %s is not supported. choose between %s', dataType, types);
            end
            if max(size(varargin)) == 0
                obj.dataMatrix = cat(3, obj.dataMatrix, pts);
            elseif varargin{1} == 1
                obj.dataMatrix = cat(3, pts, obj.dataMatrix);
            else
                obj.dataMatrix = cat(3, cat(3, obj.dataMatrix(1:varargin{1}-1), pts), obj.dataMatrix(varargin{1}:end));
            end;
            
        end;
        
        %-----------------------------------
        function xyz = mm2xyz(obj, dataMatr)
        %-----------------------------------
            xyz = ones(size(dataMatr));
            res = obj.resolution();
            for i = 1:size(dataMatr, 3)
                for j = 1:size(res, 1)
                    xyz(:, j, i) = dataMatr(:, j, i)/res(j, 1, i);
                end;
            end;
        end;
        
        %---------------------------------------
        function xyz = rlapfh2xyz(obj, dataMatr)
        %---------------------------------------
            xyz = calcfunctions('rlapfh2xyz', obj.dataNumber, dataMatr(:, 1), dataMatr(:, 2), dataMatr(:, 3))';
        end;
        
        %--------------------------------------------
        function pts = centroids(obj, representation)
        %--------------------------------------------
            % function for computing the centroids for the curves in the
            % representation defined.
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
            pts = ones(size(obj.dataMatrix, 3), size(obj.dataMatrix, 2));
            % calculate the centroids as the average of the points defining
            % each curve.
            for i = 1:size(obj.dataMatrix, 3)
                len = length(data(:,1,i));
                pts(i,:) = [sum(data(:,1,i)) sum(data(:,2,i)) sum(data(:,3,i))]/len;
            end;
        end;
        
        %----------------------------------------------------
        function boolvar = vert_averaging(obj, methodCell, n)
        %----------------------------------------------------
            % average the points defining the curves in the plane that the
            % points lie in.
            data = obj.mm();
            for j = 1:size(data, 1)
                obj.dataMatrix(j, :, :) = obj.mm2xyz(permute(obj.repeatedaaaverage(permute(data(j, :, :), [3 2 1]), methodCell, n), [3 2 1]));
            end;
            boolvar = true;
        end;
            
        
        %-----------------------------------------------------
        function newPoints = aafilter(~, oldPoints, methodVector)
        %-----------------------------------------------------
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
        function cosTheta = pointsimilarity(obj, representation)
        %-------------------------------------------------------
            switch representation
                case 'xyz'
                    cosTheta = pointsimilarity@plottest.dataparentclass(obj.xyz());
                case 'mm'
                    cosTheta = pointsimilarity@plottest.dataparentclass(obj.mm());
                case 'rlapfh'
                    cosTheta = pointsimilarity@plottest.dataparentclass(obj.rlapfh());
                otherwise
                    types = 'xyz, mm or rlapfh';
                    error('The representation: %s is not supported. choose between %s', dataType, types);
            end;
        end;
        
        %----------------------------------------------
        function areas = curvearea(obj, representation)
        %----------------------------------------------
            switch representation
                case 'xyz'
                    data = obj.xyz();
                case 'mm'
                    data = obj.mm();
                otherwise
                    types = 'xyz or mm';
                    error('The representation: %s is not supported. choose between %s', dataType, types);
            end;
            areas = ones(size(data, 3), 1);
            for i = 1:size(areas, 2)
                areas(i, 1) = polyarea(data(:, 1, i), data(:, 2, i));
            end;
        end;
        
        %--------------------------------------------
        function aver = averarea(obj, representation)
        %--------------------------------------------
            aver = mean(obj.curvearea(representation));
        end;
        
    end
    
end

