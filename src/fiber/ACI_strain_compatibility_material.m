classdef ACI_strain_compatibility_material
    properties 
        id;
    end
    
    methods (Abstract)
        stress = getStress(obj,strain)
    end    
    methods 
        function obj = set.id(obj,id)
            assert(isscalar(id) && isnumeric(id),'id must be a numeric scalar')
            obj.id = id;
        end
    end    
end