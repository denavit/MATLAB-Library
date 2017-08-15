classdef Circle_Shape < geometric_shape
        
    properties
        D
    end
    
    methods
        function obj = Circle_Shape(D)
            obj.D = D;
        end 
        function tf = is_section_valid(obj)
            tf = obj.D > 0;    % D should be positive
        end
        function a = A(obj)
            a = (pi/4)*obj.D^2;
        end
        function i = I(obj,axis)
            i = (pi/64)*obj.D^4;
        end        
    end
end

