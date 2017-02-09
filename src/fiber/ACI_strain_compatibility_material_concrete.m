classdef ACI_strain_compatibility_material_concrete < ACI_strain_compatibility_material

    properties (SetAccess = immutable)
        fc      % Steel yield stress
        units   % Units
        
        strainAtExtremeConcreteFiber = -0.003;
    end

    methods
        function obj = ACI_strain_compatibility_material_concrete(id,fc,units)
            % Constructor
            obj.id      = id;
            obj.fc      = fc;
            obj.units   = units;
        end
        function b1 = beta1(obj)
            switch obj.units
                case 'US'
                    if obj.fc <= 4.0
                        b1 = 0.85;
                    elseif obj.fc <= 8.0
                        b1 = 1.05 - 0.05*obj.fc;
                    else
                        b1 = 0.65;
                    end
                otherwise
                    error('Unknown units: %s',obj.units);
            end
        end
        function stress = getStress(obj,strain)
            stress = nan(size(strain));
            
            ecr = obj.strainAtExtremeConcreteFiber*(1-obj.beta1);
            
            ind = strain  > ecr;
            stress(ind) = 0;
            ind = strain <= ecr;
            stress(ind) = -0.85*obj.fc;
        end
    end
end