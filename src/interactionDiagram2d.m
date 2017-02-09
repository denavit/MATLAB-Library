classdef interactionDiagram2d
    %interactionDiagram2d Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        idX
        idY 
        q
    end
    
    methods
        function obj = interactionDiagram2d(idX,idY)
            assert(isvector(idX),'idX must be a vector')
            assert(isvector(idY),'idY must be a vector')
            if isrow(idX) 
                idX = idX';
            end
            if isrow(idY) 
                idY = idY';
            end                 

            % Compute Angles
            [q,~]   = cart2pol(idX,idY);
            q       = mod(q,2*pi);
            
            % Sort
            [q,IX]  = sort(q);
            idX     = idX(IX);
            idY     = idY(IX);
            
            % Repeat end data
            idX = vertcat(idX(end),idX,idX(1));
            idY = vertcat(idY(end),idY,idY(1));
            q   = vertcat(q(end)-2*pi,q,q(1)+2*pi);
            
            % Store
            obj.idX = idX;
            obj.idY = idY;
            obj.q   = q;
        end
        
        function d = radial_distance(obj,angles)
            angles = mod(angles,2*pi);
            d = nan(size(angles));
                       
            for i = 1:length(angles)
                ind = find( obj.q >= angles(i) ,1,'first');
                if isempty(ind)
                    d(i) = nan; 
                elseif ind==1
                    if (obj.q(1) == angles(i))
                        d(i) = hypot(obj.idX(1),obj.idY(1));
                    else                    
                        d(i) = 0;
                    end
                else
                    Ax = obj.idX(ind);
                    Ay = obj.idY(ind);
                    Bx = obj.idX(ind-1);
                    By = obj.idY(ind-1);
                    [Ix,Iy] = find_intersection_between_two_lines(Ax,Ay,Bx,By,...
                        0,0,cos(angles(i)),sin(angles(i)));
                    if isempty(Ix)
                        error('Bad interaction diagram');
                    else
                        d(i) = hypot(Ix,Iy);
                    end
                end
            end
        end

        % Outside/Inside Functions        
        function errors = compareTwo(obj,test_id,angles)
            d_base = obj.radial_distance(angles);
            d_test = test_id.radial_distance(angles);
            errors = (d_base-d_test)./d_base;
        end
        function max_ratio = checkPoints(obj,pointsX,pointsY)
            [q_pts,d_pts] = cart2pol(pointsX,pointsY);
            d = obj.radial_distance(q_pts);
            max_ratio = max(d_pts./d);
        end
        function [X,Y,ind,x] = findIntersection(obj,pathX,pathY)
            [q_path,d_path] = cart2pol(pathX,pathY);
            d = obj.radial_distance(q_path);
            [ind,x] = find_limit_point_in_vector(d_path-d,0);
            if isempty(ind)
                X = nan;
                Y = nan;
            else
                X = interpolate_vector(pathX,ind,x);
                Y = interpolate_vector(pathY,ind,x);
            end
            if nargout < 4
                clear ind x
            end
        end
        function X = findXgivenY(obj,Y,signX)
            switch lower(signX)
                case {'pos','positive'}
                    peakX = 1.1*max(obj.idX);
                case {'neg','negative'}
                    peakX = 1.1*min(obj.idX);
                otherwise
                    error('Unknown signX');
            end
            npts = 1000;
            pathX = linspace(0,peakX,npts);
            pathY = Y*ones(1,npts);
            [X,~] = obj.findIntersection(pathX,pathY);
        end
        function Y = findYgivenX(obj,X,signY)
            switch lower(signY)
                case {'pos','positive'}
                    peakY = 1.1*max(obj.idY);
                case {'neg','negative'}
                    peakY = 1.1*min(obj.idY);
                otherwise
                    error('Unknown signY');
            end
            npts = 1000;
            pathX = X*ones(1,npts);
            pathY = linspace(0,peakY,npts);
            [~,Y] = obj.findIntersection(pathX,pathY);
        end
        
        % Plotting Functions
        function plot(obj,varargin)
            plot(obj.idX,obj.idY,varargin{:})
        end                
    end
end
