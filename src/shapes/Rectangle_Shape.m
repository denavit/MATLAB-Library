classdef Rectangle_Shape < geometric_shape
        
    properties
        H
        B
        r = 0;
    end
    
    methods
        function obj = Rectangle_Shape(H,B,r)
            obj.H = H;
            obj.B = B;
            if nargin > 2
                obj.r = r;
            end
        end 
        function tf = is_section_valid(obj)
            tf1 = obj.B > 0;    % B should be positive
            tf2 = obj.H > 0;    % H should be positive
            tf3 = obj.r >= 0;   % r should be positive or zero
            tf4 = obj.r < min([obj.B obj.H])/2;
            tf = all([tf1 tf2 tf3 tf4]);
        end
        function a = A(obj)
            a = obj.H*obj.B - (4-pi)*obj.r^2;
        end
        function i = I(obj,axis)
            switch lower(axis)
                case {'major','strong'}
                    if obj.r == 0
                        i = (1/12)*obj.B*obj.H^3;
                    else
                        i = (1/12)*obj.B*obj.H^3 ...
                            - 4 * ((1/12)*obj.r^4 + obj.r^2*(obj.H/2-obj.r/2)^2) ...
                            + 4 * ((pi/16-4/(9*pi))*obj.r^4 ...
                            + (pi/4)*obj.r^2*(obj.H/2-(obj.r-(4*obj.r)/(3*pi)))^2);
                    end
                case {'minor','weak'}
                    if obj.r == 0
                        i = (1/12)*obj.H*obj.B^3;
                    else
                        i = (1/12)*obj.H*obj.B^3 ...
                            - 4 * ((1/12)*obj.r^4 + obj.r^2*(obj.B/2-obj.r/2)^2) ...
                            + 4 * ((pi/16-4/(9*pi))*obj.r^4 ...
                            + (pi/4)*obj.r^2*(obj.B/2-(obj.r-(4*obj.r)/(3*pi)))^2);
                    end
                otherwise
                    error('Bad axis');
            end
        end
        function s = S(obj,axis)
            I = obj.I(axis);
            switch lower(axis)
                case {'major','strong'}
                    s = I/(obj.H/2);
                case {'minor','weak'}
                    s = I/(obj.B/2);
                otherwise
                    error('Bad axis');
            end 
        end
        function z = Z(obj,axis)
            switch lower(axis)
                case {'major','strong'}
                    if obj.r == 0
                        z = obj.B*obj.H^2/4;
                    else
                        z = obj.B*obj.H^2/4 ...
                            - 4*((1-pi/4)*obj.r^2)*(obj.H/2-((10-3*pi)/(12-3*pi))*obj.r);
                    end
                case {'minor','weak'}
                    if obj.r == 0
                        z = obj.H*obj.B^2/4;
                    else
                        z = obj.H*obj.B^2/4 ...
                            - 4*((1-pi/4)*obj.r^2)*(obj.B/2-((10-3*pi)/(12-3*pi))*obj.r);
                    end
                otherwise
                    error('Bad axis');
            end
        end        
    end
    
end

