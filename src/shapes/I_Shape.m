classdef I_Shape < geometric_shape
        
    properties
        d
        tw
        bf
        tf
        k = 0
    end
    
    methods
        function obj = I_Shape(d,tw,bf,tf,k)
            obj.d  = d;
            obj.tw = tw;
            obj.bf = bf;
            obj.tf = tf;
            if nargin > 4
                obj.k = k;
            end
        end 
        function tf = is_section_valid(obj)
            tf1 = obj.d  > 0;   
            tf2 = obj.tw > 0;   
            tf3 = obj.bf > 0;   
            tf4 = obj.tf > 0;   
            tf5 = obj.tw < obj.bf;
            tf6 = obj.tf < obj.d/2;
            if obj.k == 0
                tf = all([tf1 tf2 tf3 tf4 tf5 tf6]);
            else
                tf7 = obj.k > obj.tf;
                tf8 = obj.k < obj.d/2;
                tf = all([tf1 tf2 tf3 tf4 tf5 tf6 tf7 tf8]);
            end
        end
        function dw = dw(obj)
            dw = obj.d - 2*obj.tf;
        end
        function r = r_fillet(obj)
            if obj.k == 0
                r = 0;
            else
                r = obj.k - obj.tf;
            end
        end
        function a = A_fillet(obj)
            a = (1-pi/4) * obj.r_fillet^2;
        end
        function i = I_fillet(obj)
            i = (176-84*pi+9*pi^2)/144/(4-pi) * obj.r_fillet^4;
        end
        function y_bar = Y_Bar_fillet(obj) 
            y_bar = 2/(12-3*pi) * obj.r_fillet;
        end
        
        function a = A(obj)
            a = 2*obj.bf*obj.tf + obj.dw*obj.tw + 4*obj.A_fillet;
        end
        function i = I(obj,axis)
            switch lower(axis)
                case {'z','x','major','strong'}
                    i = (1/12)*(obj.bf*obj.d^3 - (obj.bf-obj.tw)*obj.dw^3) ...
                        + 4*(obj.I_fillet+obj.A_fillet*(obj.dw/2-obj.r_fillet+obj.Y_Bar_fillet)^2);
                case {'y','minor','weak'}
                    i = (1/12)*(2*obj.tf*obj.bf^3 + obj.dw*obj.tw^3) ...
                        + 4*(obj.I_fillet+obj.A_fillet*(obj.tw/2+obj.r_fillet-obj.Y_Bar_fillet)^2);
                otherwise
                    error('Unknown axis: %s',axis);
            end 
        end
        function s = S(obj,axis)
            I = obj.I(axis);
            switch lower(axis)
                case {'z','x','major','strong'}
                    s = I/(obj.d/2);
                case {'y','minor','weak'}
                    s = I/(obj.bf/2);
                otherwise
                    error('Unknown axis: %s',axis);
            end 
        end
        function z = Z(obj,axis)
            switch lower(axis)
                case {'z','x','major','strong'}
                    z = 2*(...
                        (obj.tf*obj.bf)*(obj.d/2-obj.tf/2) + ...
                        (obj.dw/2*obj.tw)*(obj.dw/4) + ...
                        2*(obj.A_fillet)*(obj.dw/2-obj.r_fillet+obj.Y_Bar_fillet));
                case {'y','minor','weak'}
                    z = 2*(...
                        2*(obj.tf*obj.bf/2)*(obj.bf/4) + ...
                        (obj.dw*obj.tw/2)*(obj.tw/4) + ...
                        2*(obj.A_fillet)*(obj.tw/2+obj.r_fillet-obj.Y_Bar_fillet));
                otherwise
                    error('Unknown axis: %s',axis);
            end 
        end
        function j = J(obj)
            r = obj.r_fillet;
            alpha1 = -0.0420 + 0.220*(obj.tw/obj.tf) + 0.136*(r/obj.tf) - 0.0865*(obj.tw*r/obj.tf^2) - 0.0725*(obj.tw/obj.tf)^2;
            D1 = ((obj.tf+r)^2+obj.tw*(r+obj.tw/4)) / (2*r+obj.tf);
            j = (2/3)*obj.bf*obj.tf^3 + obj.tw^3*(obj.d-2*obj.tf)/3 + 2*alpha1*D1^4 - 0.420*obj.tf^4;
        end
    end
end

