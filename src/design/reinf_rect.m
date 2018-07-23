classdef reinf_rect < reinf
    
    properties
        Bx      % center-to-center distance between outside bars in x-direction
        By      % center-to-center distance between outside bars in y-direction
        xc      % x coordiante of the center of bar pattern
        yc      % y coordiante of the center of bar pattern
        nbx     % number of bars in the x-direction
        nby     % number of bars in the y-direction
    end
    
    methods
        function obj = reinf_rect(Bx,By,xc,yc,nbx,nby,Ab)
            obj.Bx  = Bx;
            obj.By  = By;
            obj.xc  = xc;
            obj.yc  = yc;
            obj.nbx = nbx;
            obj.nby = nby;
            obj.Ab  = Ab;
        end
        function [x,y] = coordinates(obj)
            x1 = linspace(-0.5*obj.Bx+obj.xc, 0.5*obj.Bx+obj.xc,obj.nbx);
            x2 = linspace( 0.5*obj.Bx+obj.xc, 0.5*obj.Bx+obj.xc,obj.nby);
            x3 = linspace( 0.5*obj.Bx+obj.xc,-0.5*obj.Bx+obj.xc,obj.nbx);
            x4 = linspace(-0.5*obj.Bx+obj.xc,-0.5*obj.Bx+obj.xc,obj.nby);
            y1 = linspace( 0.5*obj.By+obj.yc, 0.5*obj.By+obj.yc,obj.nbx);
            y2 = linspace( 0.5*obj.By+obj.yc,-0.5*obj.By+obj.yc,obj.nby);
            y3 = linspace(-0.5*obj.By+obj.yc,-0.5*obj.By+obj.yc,obj.nbx);
            y4 = linspace(-0.5*obj.By+obj.yc, 0.5*obj.By+obj.yc,obj.nby);
            x = [x1(1:end-1) x2(1:end-1) x3(1:end-1) x4(1:end-1)];
            y = [y1(1:end-1) y2(1:end-1) y3(1:end-1) y4(1:end-1)];
        end
        function n = nb(obj)
            n = 2*(obj.nbx+obj.nby)-4;
        end
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
            db = sqrt(4/pi*obj.Ab);
            circX = (db/2*cos(angles));
            circY = (db/2*sin(angles));
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