classdef Circular_Tube_Shape
        
    properties
        D
        t
    end
    
    methods
        function obj = Circular_Tube_Shape(D,t)
            obj.D = D;
            obj.t = t;
        end 
        function tf = is_section_valid(obj)
            tf1 = obj.D > 0;    % D should be positive
            tf2 = obj.t > 0;    % t should be positive
            tf3 = obj.t < obj.D/2;
            tf = all([tf1 tf2 tf3]);
        end
        function a = A(obj)
            circo = Circle_Shape(obj.D);
            circi = Circle_Shape(obj.D-2*obj.t);
            a = circo.A - circi.A;
        end
        function i = I(obj,axis)
            circo = Circle_Shape(obj.D);
            circi = Circle_Shape(obj.D-2*obj.t);
            i = circo.I(axis) - circi.I(axis);
        end       
    end
    
end

