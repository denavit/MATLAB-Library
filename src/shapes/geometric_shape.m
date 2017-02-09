classdef geometric_shape
         
    methods   
        function r = r(obj,axis)
            r = sqrt(obj.I(axis)/obj.A);
        end
        function z = Z(obj,axis)
            error('Z not implemented for class: %s',class(obj));
        end      
    end
end
