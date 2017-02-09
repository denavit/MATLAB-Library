classdef structural_shape
    
    properties
        Lx = 0;         % Unbraced length of the member for major axis buckling
        Ly = 0;         % Unbraced length of the member for minor axis buckling
        Kx = 1          % Effective length factor in the strong axis
        Ky = 1          % Effective length factor in the weak axis
        units = ' ';    % unit system
                        %  - US = United States customary units (i.e., kips, inches, ksi)
                        %  - SI = International System of Units (i.e., N, mm, MPa)
    end
    
    properties (Hidden = true)
        color_steelFill     = [0.75 0.75 0.75];
        color_concreteFill  = [0.90 0.90 0.90];
    end
       
    methods       
        function l = L(obj,axis)
            switch lower(axis)
                case {'major','strong'}
                    l = obj.Lx;
                case {'minor','weak'}
                    l = obj.Ly;
                otherwise
                    error('Bad axis');
            end
        end
        function k = K(obj,axis)
            switch lower(axis)
                case {'major','strong'}
                    k = obj.Kx;
                case {'minor','weak'}
                    k = obj.Ky;
                otherwise
                    error('Bad axis');
            end
        end
        function ei = EI(obj,axis,type)
            [E,~,I] = obj.sectionPropertiesForElasticAnalysis2d(axis,type);
            ei = E*I;
        end
        function loe = lambdaoe(obj,axis)            
            loe = obj.K(axis)*obj.L(axis)*sqrt(obj.Pnco/obj.EIeff(axis))/pi;
        end
        function loe1 = lambdaoe1(obj,axis)
            loe1 = obj.L(axis)*sqrt(obj.Pnco/obj.EIeff(axis))/pi;
        end
        function lp = Lp(obj,axis,Li)
            error('Lp not implemented for class: %s',class(obj));
        end
        function scACI = strainCompatibilityAciObject(obj)
            error('strainCompatibilityAciObject not implemented for class: %s',class(obj));
        end         
    end
end
