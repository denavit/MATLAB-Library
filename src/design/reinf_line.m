classdef reinf_line < reinf
    
    properties
        xi
        yi
        xj
        yj
        nspaces
        bar_at_start = true;
        bar_at_end = true;
    end
    
    methods
        function obj = reinf_line(xi,yi,xj,yj,nspaces,Ab)
            obj.xi = xi;
            obj.yi = yi;
            obj.xj = xj;
            obj.yj = yj;
            obj.nspaces = nspaces;
            obj.Ab = Ab;
        end
        function [x,y] = coordinates(obj)
            x = linspace(obj.xi,obj.xj,obj.nb);
            y = linspace(obj.yi,obj.yj,obj.nb);
            if ~obj.bar_at_start
                x = x(2:end);
                y = y(2:end);
            end
            if ~obj.bar_at_end
                x = x(1:(end-1));
                y = y(1:(end-1));
            end
        end
        function n = num_bars(obj)
            n = obj.nspaces + 1;
            if ~obj.bar_at_start
                n = n-1;
            end
            if ~obj.bar_at_end
                n = n-1;
            end
        end
    end
    
end