classdef piecewiseLinear
    
    properties
        xp
        yp
    end
    
    methods
        function obj = piecewiseLinear(xp,yp)
            % Constructor
            obj.xp = xp;
            obj.yp = yp; 
        end
        function res = ismonotonic(obj)
            diff = obj.xp(2:end)-obj.xp(1:end-1);
            if ( min(diff) >= 0 )
                % Incresing Monotonic
                res = true;
            elseif ( max(diff) <= 0 )
                % Decreasing Monotonic
                res = true;
            else
                res = false;
            end
        end
        function y = val(obj,x)
            % Evaluate the piecewise linear function at a value of x
            ind = find(obj.xp >= x,1,'first');
            if (ind == 1)
                if (x == obj.xp(1))
                    y = obj.yp(1);
                else                
                    error('x out of range');
                end
            elseif isempty(ind)
                error('x out of range');
            else
                y = obj.yp(ind-1) + (obj.yp(ind)-obj.yp(ind-1))*...
                    ((x-obj.xp(ind-1)))/((obj.xp(ind)-obj.xp(ind-1)));
            end
        end 
        function x = inverse(obj,y)
            % Evaluate the inverse piecewise linear function at a value of y
            ind = find(obj.yp >= y,1,'first');
            if (ind == 1)
                if (y == obj.yp(1))
                    x = obj.xp(1);
                else
                    error('y out of range')
                end
            elseif isempty(ind) 
                error('y out of range')
            else
                x = obj.xp(ind-1) + (y-obj.yp(ind-1))*...
                    ((obj.xp(ind)-obj.xp(ind-1))/(obj.yp(ind)-obj.yp(ind-1)));
            end
        end         
        function A = area(obj)
            % Find the area under the piecewise linear function
            A = trapz(obj.xp,obj.yp);
        end
        function y = peak(obj,opt)
            % Find the peak value of the piecewise linear function
            if nargin < 2
                opt = 'standard';
            end  
            switch lower(opt)
                case 'standard'
                    y = max(obj.yp);
                case 'absolute'
                    y = max(abs(obj.yp));
                case 'signedabsolute'
                    ymax = max(obj.yp);
                    ymin = min(obj.yp);
                    if abs(xmax) > abs(xmin)
                        y = ymax;
                    else
                        y = ymin;
                    end
                otherwise
                    error('Unknown peak option');
            end
        end
        function x = xAtPeak(obj,opt)
            % Find the x-coordinate (absicissa) of the peak value
            if nargin < 2
                opt = 'standard';
            end  
            switch lower(opt)
                case 'standard'
                    [~,ind] = max(obj.yp);
                    x = obj.xp(ind);
                case 'absolute'
                    [~,ind] = max(abs(obj.yp));
                    x = obj.xp(ind);
                otherwise
                    error('Unknown peak option');
            end            
        end
        function e = initialTangent(obj,x)
            % Find the initial tangent
            if nargin < 2
                % True inital tangent
                e = (obj.yp(2)-obj.yp(1))/(obj.xp(2)-obj.xp(1));
            else
                % Secant inital tangent
                yval = x*obj.peak;
                xval = obj.inverse(yval);
                e = yval/xval;
            end
        end
        function [xOffset,yOffset] = offsetStress(obj,ratio,offset)
            slope = obj.initialTangent(ratio);
            vals = obj.yp - slope*(obj.xp-offset);
            ind = find(vals < 0,1,'first');
            if (ind == 1)
                error('out of range')
            elseif isempty(ind) 
                error('out of range')
            else
                [xOffset,yOffset] = find_intersection_between_two_lines(...
                    obj.xp(ind-1),obj.yp(ind-1),obj.xp(ind),obj.yp(ind),...
                    offset,0,1+offset,slope);
            end
        end
        function metrics = loadDeformationMetrics(obj)
            metrics.peakLoad              = obj.peak;
            metrics.deformationAtPeakLoad = obj.xAtPeak;
            metrics.initialTangent        = obj.initialTangent(1/3);
            metrics.areaUnderCurve        = obj.area;
        end
    end
end
