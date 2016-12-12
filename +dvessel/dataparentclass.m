classdef (Abstract) dataparentclass < handle
    %dataparentclass Parent class that the other classes builds upon. Can
    %not be initialized.
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % Functions that have to be written in child classes.
        %-----------------------------
        function res = resolution(obj)
        %-----------------------------
        end
        
        %----------------------
        function pts = xyz(obj)
        %----------------------
        end
        
        %---------------------
        function pts = mm(obj)
        %---------------------
        end
        
        %-------------------------
        function pts = rlapfh(obj)
        %-------------------------
        end
        
        %--------------------------------------------
        function pts = centroids(obj)
        %--------------------------------------------
        end
        
        %------------------------------------
        function ratio = compchullper(obj, points)
        %------------------------------------
            % If the points are 2-dimensional
            try
                cHull = convhull(points);
                ratio = obj.curvelength(points(cHull))/obj.curvelength(points);
            % If the points are 3-dimensional, drop the z coordinates
            catch
                cHull = convhull(points(:,1:2));
                % For debugging, remove later
                if true == 0
                    figure
                    plot(points(:,1),points(:,2), 'b', points(cHull,1), points(cHull, 2), 'r')
                    obj.curvelength(points(:,1:2))
                    obj.curvelength(points(cHull,1:2))
                end;
                ratio = obj.curvelength(points(cHull,1:2))/obj.curvelength(points(:,1:2));
            end;
        end;
        
        %------------------------------------
        function length = curvelength(~, points)
        %------------------------------------
            length = 0;
            for i = 2:size(points, 1)
                length = length + norm(points(i, :) - points(i - 1, :),2);
            end;
        end;
        
        %-------------------------------
        function pos = removeNaN(~, values)
        %-------------------------------
        % Remove NaN values from data
            pos = values;
            for i = 1:max(size(pos))
                pos{i}(:,:,any(isnan(pos{i}(1,1,:)),2)) = [];
            end;
        end;
        
        %----------------------------------------------------
        function cosTheta = cosinesimilarity(~, vect1, vect2)
        %----------------------------------------------------
        % Compute cosine similarity
            cosTheta = dot(vect1, vect2)/(norm(vect1)*norm(vect2));
        end;
        
        %-----------------------------------------------
        function cosTheta = pointsimilarity(obj, points)
        %-----------------------------------------------
            if ~(size(points, 2) == 3)
                error('points has to be on the format [n 3 m] to work');
            end;
            cosTheta = ones(size(points, 1) - 2, 1, size(points, 3));
            for i = 1:size(cosTheta, 3)
                for j = 1:size(cosTheta, 1)
                    vect1 = points(j + 1, :, i) - points(j, :, i);
                    vect2 = points(j + 2, :, i) - points(j + 1, :, i);
                    cosTheta(j, 1, i) = obj.cosinesimilarity(vect1, vect2);
                end;
            end;
        end;
        
        function boolvar = save2creo(obj, fileName, representation, dataObject, varargin)
            if isa(obj,'dvessel.splineclass')
                error('%s not supported.', class(obj));
            end;
            switch dataObject
                case 'centroids'
                    data = obj.centroids(representation);
                case 'curves'
                    switch representation
                        case 'xyz'
                            data = obj.xyz();
                        case 'mm'
                            data = obj.mm();
                        case 'rlapfh'
                            data = obj.rlapfh();
                        otherwise
                            types = 'xyz, mm or rlapfh';
                            error('The representation: %s is not supported. choose between %s', representation, types);
                    end;
                case 'projection'
                    data = obj.vectorprojection(representation, varargin{1}, varargin{2});
                otherwise
                    types = 'centroids, curves or projection';
                    error('The dataObject: %s is not supported. choose between %s', dataObject, types);
            end;
            curveNo = 0;
            fileID = fopen(fileName,'w');
            fprintf(fileID,'%s\r\n','open');
            fprintf(fileID,'%s\r\n','arclength');
            
            for i = 1:size(data, 4)
                for j = 1:size(data, 3)
                    switch dataObject
                        case 'centroids'
                            if any(any(~isnan(data(:, :, j, i))))
                                curveNo = curveNo + 1;
                                fprintf(fileID,'%s %i\r\n','begin section !',obj.dataNumber);
                                fprintf(fileID,'%s %i\r\n','begin curve !',curveNo);
                                for k = 1:size(data, 1)
                                    if ~any(isnan(data(k, :, j, i)))
                                        fprintf(fileID,'%f %f %f \r\n', data(k, :, j, i)');
                                    end;
                                end;
                            end;
                        case {'curves', 'projection'}
                            if ~any(any(isnan(data(:, :, j, i))))
                                curveNo = curveNo + 1;
                                fprintf(fileID,'%s %i\r\n','begin section !',obj.dataNumber);
                                fprintf(fileID,'%s %i\r\n','begin curve !',curveNo);
                                fprintf(fileID,'%f %f %f \r\n', data(:, :, j, i)');
        %                         for k = 1:size(data, 1)
        %                             fprintf(fileID,'%f %f %f \r\n', data(k, 1, j, i), data(k, 2, j, i), data(k, 3, j, i) );
        %                         end;
                                fprintf(fileID,'\r\n','');
                            end;
                    end;
                end;
            end;
            fclose(fileID);
            boolvar = true;
        end;
        
        %-------------------------------------------------------------------
        function projectedpts = pointsprojection(~, points, vector, planepoint)
        %-------------------------------------------------------------------
            % Function to project the points in a plane defined by a point (planepoint)
            % and vector normal to the plane
            nvector = vector/norm(vector); %Normalize the vector.
            projectedpts = zeros(size(points));
            for i= 1:length(points)
                projectedpts(i,:) = points(i,:) - dot(points(i,:) - planepoint, nvector)*nvector;
            end;
        end;
        
        %------------------------------------------------------------------------
        function data = vectorprojection(obj, representation, vectMatr, varargin)
        %------------------------------------------------------------------------
            if isa(obj,'dvessel.splineclass')
                error('%s not supported.', class(obj));
            end;
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
            if size(data, 2) ~= size(vectMatr, 2) | size(data, 3) ~= size(vectMatr, 1) | size(data, 4) ~= size(vectMatr, 3)
                error('You must supply exactly one vector for each curve in your data');
            end;
            if length(varargin) < 1
                ptsMatr = obj.centroids(representation);
            else
                if ~any(size(vectMatr) == size(varargin{1}))
                    error('You need as many points as vectors and they need to be the same dimension');
                end;
                ptsMatr = varargin{1};
            end;
            for i = 1:size(data, 4)
                for j = 1:size(data, 3)
                    data(:, :, j, i) = obj.pointsprojection(data(:, :, j, i), vectMatr(j, :, i), ptsMatr(j, :, i));
                end;
            end;
            
        end;
            
        
        
        
    end
end

