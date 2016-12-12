classdef splineclass < dvessel.dataparentclass
    %splineclass This class is used to create the center lines of the
    %vessels.
    %   Detailed explanation goes here
    
    properties
        % The points that the curve will pass through
        interppts;
        % The control points of the spline
        ctrlpts;
        % The interpolation vector, the parameter values at which the
        % spline will pass through interppts.
        tau;
        % The knot vector of the spline
        knots;
        % The degree of the spline
        degree;
        % Domain of the spline. Common is [0, 1].
        domain;
    end
    
    methods
        function obj = splineclass(interppoints, degree, domain, taumethod, varargin)
            if size(interppoints, 3) > 1 | size(interppoints, 2) ~= 3
                error('The dimensions of interppoints provided are wrong. You provided %d, but it should be [a 3 1]', size(interppoints));
            end;
            if degree < 1
                error('You provided a value of %d for degree. It has to be 1 or larger.', degree);
            end;
            if mod(degree, 1)
                error('degree has to take an integer value');
            end;
            if size(interppoints, 1) < degree - 1
                error('The number of interpolation points has to be larger than the degree');
            end;
            obj.degree = degree;
            obj.interppts = interppoints;
            
            if ~all(size(domain) == [1 2])
                error('domain is of wrong dimensions. Has to be [1 2]');
            end;
            if domain(2) <= domain(1)
                error('domain(2) has to be larger than domain(1)');
            end;
            obj.domain = domain;
            
            obj.tau = obj.spacing(taumethod, varargin{:});
            obj.knots = obj.knotaveraging();
            obj.ctrlpts = obj.interpolate();
        end;
        
        % Andreas function
        function knotvect = knotaveraging(obj)
            % From A. S�derlund, Parameter selection and derivative conditions for 
            % B-splines applied to gas turbine blade modeling

            % Averaging the knot vector according to de Boor from the parameter vector.
            % Will create a clamped knot vector of degree 'degree'. Quadratic splines
            % are assumed to be of degree 3 in matlab.

            tlen = length(obj.tau);
            if tlen <= obj.degree
                error('Not enough parameter values . Has to be strictly more than the degree')
            end;

            starting = ones(1, obj.degree) * obj.domain(1);
            ending   = ones(1, obj.degree) * obj.domain(2);
            center   = ones(1, tlen - obj.degree);

            len = obj.domain(2) - obj.domain(1);
            for i = 1:length(center)
                s = i + 1; e = i + obj.degree - 1;
                center(i) = len/(obj.degree - 1)*sum(obj.tau(s:e)); 
            end;
            knotvect = [starting center ending];
        end;

        % Andreas function
        function tau = spacing(obj, varargin)
            % Function producing the parameter vector. Supply it with a string and
            % optional extra arguments depending on which method used. Implementer
            % methods are:
            %
            % 'uniform'. No extra input.
            % 'centripetal'. Takes a real valued number. Default is 0.5
            % 'chordlength'. No extra input.

            switch varargin{1}
                case 'uniform'
                    hvector = obj.uniform_h();
                case 'centripetal'
                    if length(varargin) < 2
                        % If alpha has not been provided.
                        alpha = 0.5;
                    else
                        alpha = varargin{2};
                    end
                    hvector = obj.centripetal_h(alpha);
                case {'chordlength', 'chordLength'}
                    % Chordlength is a special case of centripetal.
                    hvector = obj.centripetal_h(1);
                case {'affineInvariant', 'affineinvariant'}
                    hvector = obj.affineinvariant_h();
                case 'random'
                    % If random is chosen, return tau as a randomly chosen
                    % vector and function is done.
                    tau = sort(rand(1, length(obj.interppts)));
                    return
                otherwise
                    error('some error');
            end;
            tau = ones(1, length(obj.interppts));
            
            % First and last already known.
            tau(1) = obj.domain(1); 
            tau(end) = obj.domain(2);
            %Go through all except for first and last
            for i = 1:length(tau) - 2
                tau(i + 1) = tau(i) + hvector(i);
            end;
        end;

        %Andreas function
        function hvector = uniform_h(obj)
            % From A. S�derlund, Parameter selection and derivative conditions for 
            % B-splines applied to gas turbine blade modeling

            % Creating a uniformly distributed h vector.
            hvector = ones(1, length(obj.interppts) - 1) * (obj.domain(2) - obj.domain(1))/(length(obj.interppts) - 1);
        end;

        %Andreas function
        function hvector = centripetal_h(obj, alpha)
            % From A. S�derlund, Parameter selection and derivative conditions for 
            % B-splines applied to gas turbine blade modeling

            % Creating a centripetally distributed h vector.
            hvector = ones(1, length(obj.interppts) - 1);
            for i = 1:length(hvector)
                hvector(i) = norm(obj.interppts(i + 1,:) - obj.interppts(i,:))^alpha;
            end;
            hvector = hvector/sum(hvector)*(obj.domain(2) - obj.domain(1));
        end;

        % Andreas function
        function A = affineinvariantnorm(obj)
            % Following (Wendelin L.E Degen, Volker Milbrandt, 1997).
            % Returns the Nielsen Affine Invariant Norm matrix A of the associated
            % points.

            % centre of gravity, the centre of the points. Is a n dim vector for n dim points.
            c = sum(obj.interppts)/length(obj.interppts);
            % V matrix to create the affine norm matrix A
            V = ones(size(obj.interppts));
            for i = 1:length(obj.interppts)
                V(i,:) = obj.interppts(i,:) - c;
            end;
            % matrix A is the norm of the metric
            A = inv(V.'*V);
        end;

        function hvector = affineinvariant_h(obj)
            % Following (T. Foley, G. Nielson, 1989).
            
            % Creating a affine invariant distrinuted h vector
            A = obj.affineinvariantnorm();
            hvector = ones(1, length(obj.interppts) - 1);

            for i = 1:length(hvector)
                hvector(i) = sqrt( (obj.interppts(i,:) - obj.interppts(i + 1,:))*A*(obj.interppts(i,:) - obj.interppts(i + 1,:)).' );
            end;

            hvector = hvector/sum(hvector)*(obj.domain(2) - obj.domain(1));
        end;
        
        %----------------------------------
        function ctrlpts = interpolate(obj)
        %----------------------------------
            colmat = spcol(obj.knots, obj.degree, obj.tau);
            ctrlpts = (colmat\obj.interppts).';
        end;
        
        function vals = pts(obj, tvect)
            % Evaluate spline at tvect.
            vals = fnval(obj.spline(), tvect)';
        end;
        
        %--------------------------
        function spln = spline(obj)
        %--------------------------
            spln = spmak(obj.knots, obj.ctrlpts);
        end;
        
        function dermatr = nderivs(obj, n, vector)
            % Check if vector is specified as tau (or begins with t)
            if any(vector == 'tau')
                vector = obj.tau;
            end;
            dermatr(:, :, 1) = obj.pts(vector);
            deriv = fnder(obj.spline());
            for i = 1:n
                dermatr(:, :, end + 1) = fnval(deriv, vector)';
                deriv = fnder(deriv);
            end;
                
        end;
        
        %---------------------------------------
        function cosTheta = pointsimilarity(obj)
        %---------------------------------------
            cosTheta = pointsimilarity@dvessel.dataparentclass(obj, obj.interppts);
        end;
        
        %-------------------------------
        function vects = endvectors(obj)
        %-------------------------------
            % Returns the endvectors of the spline. Will point out of the
            % spline, not along it.
            vects(1, :) = -(obj.ctrlpts(:, 2)'   - obj.ctrlpts(:, 1)');
            vects(2, :) = obj.ctrlpts(:, end)' - obj.ctrlpts(:, end - 1)';
        end;
        
        %----------------------------
        function pts = endpoints(obj)
        %----------------------------
            % Return endpts of spline.
            pts(1, :) = obj.ctrlpts(:, 1)';
            pts(2, :) = obj.ctrlpts(:, end)';
        end;
            

    end
    
end

