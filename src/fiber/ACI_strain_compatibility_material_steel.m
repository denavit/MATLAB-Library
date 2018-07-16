classdef ACI_strain_compatibility_material_steel < ACI_strain_compatibility_material

    properties (SetAccess = immutable)
        Fy      % Steel yield stress
        Es      % Modulus of elasticity
    end

    methods
        function obj = ACI_strain_compatibility_material_steel(id,Fy,Es)
            % Constructor
            obj.id  = id;
            obj.Fy  = Fy;
            obj.Es  = Es;
        end
        function ey = ey(obj)
            ey = obj.Fy/obj.Es;
        end
        function stress = getStress(obj,strain)
            stress = nan(size(strain));
            ind = strain >= obj.ey;
            stress(ind) = obj.Fy;
            ind = strain <= -obj.ey;
            stress(ind) = -obj.Fy;
            ind = strain < obj.ey & strain > -obj.ey;
            stress(ind) = strain(ind)*obj.Es;
        end
    end
end
