classdef reinf

    properties (Hidden = true)
        color_steelFill     = [0.75 0.75 0.75];
        color_concreteFill  = [0.90 0.90 0.90];
    end
    
    properties
        Fy
        Es
        Ab
        db = nan;
    end
    
    methods
        function i = I(obj,axis)
            [x,y] = obj.coordinates();
            switch lower(axis)
                case 'x'
                    i = sum(obj.Ab*y.^2);
                case 'y'
                    i = sum(obj.Ab*x.^2);
                otherwise
                    error('Unknown axis: %s',axis);
            end
        end
        function plotSection(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            angles = linspace(0,2*pi,25);
            if isnan(obj.db)
                db_plt = sqrt(4/pi*obj.Ab);
            else
                db_plt = obj.db;
            end
            circX = (db_plt/2*cos(angles));
            circY = (db_plt/2*sin(angles));
            [xb,yb] = obj.coordinates();
            for i = 1:length(yb)
                x =  xb(i) + circX;
                y =  yb(i) + circY;
                fill(x,y,obj.color_steelFill,'LineStyle','none')
                plot(x,y,'k-','LineWidth',lineWidth);
            end
        end        
    end

end