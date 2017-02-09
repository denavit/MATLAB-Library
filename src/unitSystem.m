classdef unitSystem
    
    properties
        lengthUnits
        curvatureUnits
        timeUnits
        forceUnits
        momentUnits
        pressureUnits
        densityUnits
        angleUnits
        volumeUnits
    end
    
    methods
        function obj = unitSystem(type)
            switch type
                case 'US'
                    obj.lengthUnits     = 'in';
                    obj.curvatureUnits  = '1/in';
                    obj.timeUnits       = 's';
                    obj.forceUnits      = 'kip';
                    obj.momentUnits     = 'k-in';
                    obj.pressureUnits   = 'ksi';
                    obj.densityUnits    = '';
                    obj.angleUnits      = 'rad';
                    obj.volumeUnits     = 'cbin';
                case 'SI'
                    obj.lengthUnits     = 'mm';
                    obj.curvatureUnits  = '1/mm';
                    obj.timeUnits       = 's';
                    obj.forceUnits      = 'N';
                    obj.momentUnits     = 'N-mm';
                    obj.pressureUnits   = 'MPa';
                    obj.densityUnits    = '';  
                    obj.angleUnits      = 'rad';                    
                    obj.volumeUnits     = 'cbmm';
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function disp(obj)
            fprintf('Unit System \n');
            fprintf('  Length Units     : %s\n',obj.lengthUnits);
            fprintf('  Curvature Units  : %s\n',obj.curvatureUnits);
            fprintf('  Time Units       : %s\n',obj.timeUnits);
            fprintf('  Force Units      : %s\n',obj.forceUnits);
            fprintf('  Moment Units     : %s\n',obj.momentUnits);
            fprintf('  Pressure Units   : %s\n',obj.pressureUnits);
            fprintf('  Density Units    : %s\n',obj.densityUnits);
            fprintf('  Angle Units      : %s\n',obj.angleUnits);
            fprintf('  Volumn Units     : %s\n',obj.volumeUnits);
        end
    end
end
