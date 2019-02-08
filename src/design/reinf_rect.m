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
        function n = num_bars(obj)
            n = 2*(obj.nbx+obj.nby)-4;
        end
    end
end