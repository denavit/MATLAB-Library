classdef Rectangular_Tube_Shape < geometric_shape
        
    properties
        H
        B
        t
        ro = 0;
    end
    
    methods
        function obj = Rectangular_Tube_Shape(H,B,t,ro)
            obj.H = H;
            obj.B = B;
            obj.t = t;
            if nargin > 2
                obj.ro = ro;
            end
        end 
        function tf = is_section_valid(obj)
            tf1 = obj.B > 0;    % B should be positive
            tf2 = obj.H > 0;    % H should be positive
            tf3 = obj.t > 0;    % t should be positive
            tf4 = obj.t < min([obj.B obj.H])/2;
            tf5 = obj.ro >= 0;  % ro should be positive or zero
            tf6 = obj.ro < min([obj.B obj.H])/2;
            tf = all([tf1 tf2 tf3 tf4 tf5 tf6]);
        end
        function ri = ri(obj)
            ri = max(0,obj.ro-obj.t);
        end
        function a = A(obj)
            recto = Rectangle_Shape(obj.H,obj.B,obj.ro);
            recti = Rectangle_Shape(obj.H-2*obj.t,obj.B-2*obj.t,obj.ri);
            a = recto.A - recti.A;
        end
        function i = I(obj,axis)
            recto = Rectangle_Shape(obj.H,obj.B,obj.ro);
            recti = Rectangle_Shape(obj.H-2*obj.t,obj.B-2*obj.t,obj.ri);
            i = recto.I(axis) - recti.I(axis);
        end
        function s = S(obj,axis)
            I = obj.I(axis);
            switch lower(axis)
                case {'x','z','major','strong'}
                    s = I/(obj.H/2);
                case {'y','minor','weak'}
                    s = I/(obj.B/2);
                otherwise
                    error('Bad axis');
            end 
        end
        function z = Z(obj,axis)
            recto = Rectangle_Shape(obj.H,obj.B,obj.ro);
            recti = Rectangle_Shape(obj.H-2*obj.t,obj.B-2*obj.t,obj.ri);
            z = recto.Z(axis) - recti.Z(axis);
        end
        function j = J(obj)
            % This equation is from "Torsional Section Properties of Steel Shapes"
            % published by the Canadian Institute of Steel Construction, 2002 
            rc = mean([obj.ro,obj.ri]);
            p = 2 *((obj.H - obj.t) + (obj.B - obj.t)) - 2 * rc * (4 - pi);
            Ap = (obj.H - obj.t)*(obj.B - obj.t) - rc * (4 - pi);
            j  = 4* Ap^2*obj.t/p;
        end
    end
    
end

