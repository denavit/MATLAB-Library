classdef Circle_Shape < geometric_shape
        
    properties
        D
        
        plot_fill_color = [0.90 0.90 0.90];
        plot_num_angles = 80;
    end
    
    methods
        function obj = Circle_Shape(D)
            obj.D = D;
        end 
        function tf = is_section_valid(obj)
            tf = obj.D > 0;    % D should be positive
        end
        function d = depth(obj,axis)
            d = obj.D;            
        end        
        function a = A(obj)
            a = (pi/4)*obj.D^2;
        end
        function i = I(obj,axis)
            i = (pi/64)*obj.D^4;
        end
        function [x,y,r] = boundary_points(obj)
            x = 0;
            y = 0;
            r = obj.D/2;
        end        
        function plotSection(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            hold all
            angles = linspace(0,2*pi,obj.plot_num_angles);
            x = (obj.D/2)*cos(angles);
            y = (obj.D/2)*sin(angles);
            fill(x,y,obj.plot_fill_color,'LineStyle','none')
            plot(x,y,'k-','LineWidth',lineWidth);
            axis equal
        end        
        function add_to_fiber_section(obj,fiber_section,matID)
            fiber_section.addPatch('circ',matID,...
                0.0,0.0,0.0,obj.D/2);
        end        
    end
end

