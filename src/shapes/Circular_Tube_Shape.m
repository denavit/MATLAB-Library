classdef Circular_Tube_Shape < geometric_shape
        
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
            ro = obj.D/2;
            ri = ro-obj.t;
            a = pi*(ro^2-ri^2);
        end
        function i = I(obj,axis)
            ro = obj.D/2;
            ri = ro-obj.t;
            i = pi/4*(ro^4-ri^4);
        end
        function s = S(obj,axis)
            ro = obj.D/2;
            ri = ro-obj.t;
            s = pi/4*(ro^4-ri^4)/ro;
        end
        function j = J(obj)
            ro = obj.D/2;
            ri = ro-obj.t;
            j = pi/2*(ro^4-ri^4);
        end
        function z = Z(obj,axis)
            do = obj.D;
            di = do-2*obj.t;
            z = (do^3-di^3)/6;
        end
    end 
end
