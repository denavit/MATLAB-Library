function newValue = unitConvert(type,varargin)
% conversionFactor = unitConvert(type,oldUnit,newUnit)
% newValue         = unitConvert(type,oldValue,oldUnit,newUnit)
%

%% Input
% Parse input data
switch nargin
    case 3
        oldValue    = 1;
        oldUnit     = varargin{1};
        newUnit     = varargin{2};
    case 4
        oldValue    = varargin{1};
        oldUnit     = varargin{2};
        newUnit     = varargin{3};
    otherwise
        error('invalid input')
end

% Check for quick return
if isempty(oldValue)
    newValue = oldValue;
    return
end
if isscalar(oldValue) && isnan(oldValue)
    newValue = oldValue;
    return;
end

% Check input data 
assert(ischar(type),...
    'type must be a charcter string')
assert(ischar(oldUnit)||isa(oldUnit,'unitSystem'),...
    'oldUnit must be a charcter string or unitSystem object')
assert(ischar(newUnit)||isa(newUnit,'unitSystem'),...
    'newUnit must be a charcter string or unitSystem object')
assert(isnumeric(oldValue),...
    'oldValue must be a numeric')

%% Conversions
% Length (Base Unit: inches, in)
length_base2in = 1;
length_base2ft = 1/12;
length_base2m  = 0.0254;
length_base2cm = 2.54;
length_base2mm = 25.4;

% Time (Base Unit: seconds, s)
time_base2s = 1;
time_base2m = 60;
time_base2h = 60*60;

% Force (Base Unit: kip)
force_base2kip = 1;
force_base2lbf = 1000;
force_base2metricton = 0.453592;
force_base2longton = 1/2.24;
force_base2shortton = 0.5;
force_base2kN = 4.448222;
force_base2MN = 4.448222/1000;
force_base2N = 4448.222;

% Moment (Base Unit: kin)
moment_base2kin = 1;
moment_base2kft = 1/12;
moment_base2kNm = 4.448222*0.0254;
moment_base2MNm = 4.448222*0.0254/1000;
moment_base2tfm = 0.453592*0.0254;
moment_base2Nmm = 4448.222*25.4;

% Pressure (Base Unit: ksi)
pressure_base2psi = 1000;
pressure_base2psf = 144000;
pressure_base2ksi = 1;
pressure_base2kPa = 6.89476e3;
pressure_base2MPa = 6.89476;
pressure_base2GPa = 6.89476e-3;
pressure_base2kgscm = 70.306986;
pressure_base2tscm = 0.070306986;
pressure_base2ltsi = 1/2.24;
pressure_base2stsi = 0.5;

% Density (Base Unit: pcf)
density_base2pcf = 1;
density_base2kgcm = 16.018463;

% Angle (Base Unit: rad)
angle_base2rad = 1;
angle_base2deg = 180/pi;

% Volume (Base Unit: cubic inch, cbin)
volume_base2cbin = 1;
volume_base2cbft = (1/12)^3;
volume_base2cbyd = (1/36)^3;
volume_base2cbm  = 0.0254^3;
volume_base2cbmm = 25.4^3;

%% Convert Units
switch lower(type)
    
    case 'length'
        if isa(oldUnit,'unitSystem')
            oldUnit = oldUnit.lengthUnits;
        end
        if isa(newUnit,'unitSystem')
            newUnit = newUnit.lengthUnits;
        end
        switch lower(oldUnit)
            case 'in';      toBase = 1/length_base2in;
            case 'ft';      toBase = 1/length_base2ft;
            case 'm';       toBase = 1/length_base2m;
            case 'cm';      toBase = 1/length_base2cm;
            case 'mm';      toBase = 1/length_base2mm;
            otherwise;      error('old %s unit not recgonized: %s',type,oldUnit);
        end
        switch lower(newUnit)
            case 'in';      fromBase = length_base2in;
            case 'ft';      fromBase = length_base2ft;
            case 'm';       fromBase = length_base2m;
            case 'cm';      fromBase = length_base2cm;
            case 'mm';      fromBase = length_base2mm;
            otherwise;      error('new %s unit not recgonized: %s',type,newUnit);
        end
        
    case 'curvature'
        if isa(oldUnit,'unitSystem')
            oldUnit = oldUnit.curvatureUnits;
        end
        if isa(newUnit,'unitSystem')
            newUnit = newUnit.curvatureUnits;
        end
        switch lower(oldUnit)
            case '1/in';    toBase = length_base2in;
            case '1/mm';    toBase = length_base2mm;
            otherwise;      error('old %s unit not recgonized: %s',type,oldUnit);
        end
        switch lower(newUnit)
            case '1/in';    fromBase = 1/length_base2in;
            case '1/mm';    fromBase = 1/length_base2mm;
            otherwise;      error('new %s unit not recgonized: %s',type,newUnit);
        end
        
    case 'time'
        if isa(oldUnit,'unitSystem')
            oldUnit = oldUnit.timeUnits;
        end
        if isa(newUnit,'unitSystem')
            newUnit = newUnit.timeUnits;
        end
        switch lower(oldUnit)
            case {'s','sec'}; toBase = 1/time_base2s;
            case {'m','min'}; toBase = 1/time_base2m;
            case {'h','hr'};  toBase = 1/time_base2h;
            otherwise;        error('old %s unit not recgonized: %s',type,oldUnit);
        end
        switch lower(newUnit)
            case {'s','sec'}; fromBase = time_base2s;
            case {'m','min'}; fromBase = time_base2m;
            case {'h','hr'};  fromBase = time_base2h;
            otherwise;        error('new %s unit not recgonized: %s',type,newUnit);
        end
        
    case 'force'
        if isa(oldUnit,'unitSystem')
            oldUnit = oldUnit.forceUnits;
        end
        if isa(newUnit,'unitSystem')
            newUnit = newUnit.forceUnits;
        end
        switch lower(oldUnit)
            case {'kip','kips'};        toBase = 1/force_base2kip;
            case 'lbf';                 toBase = 1/force_base2lbf;
            case 'kn';                  toBase = 1/force_base2kN;
            case 'mn';                  toBase = 1/force_base2MN;
            case 'n';                   toBase = 1/force_base2N;
            case {'metricton','tonne'}; toBase = 1/force_base2metricton;
            case 'longton';             toBase = 1/force_base2longton;
            case 'shortton';            toBase = 1/force_base2shortton;
            otherwise;                  error('old %s unit not recgonized: %s',type,oldUnit);
        end
        switch lower(newUnit)
            case {'kip','kips'};        fromBase = force_base2kip;
            case 'lbf';                 fromBase = force_base2lbf;
            case 'kn';                  fromBase = force_base2kN;
            case 'mn';                  fromBase = force_base2MN;
            case 'n';                   fromBase = force_base2N;
            case {'metricton','tonne'}; fromBase = force_base2metricton;
            case 'longton';             fromBase = force_base2longton;
            case 'shortton';            fromBase = force_base2shortton;
            otherwise;                  error('new %s unit not recgonized: %s',type,newUnit);
        end
        
    case {'moment','energy'}
        if isa(oldUnit,'unitSystem')
            oldUnit = oldUnit.momentUnits;
        end
        if isa(newUnit,'unitSystem')
            newUnit = newUnit.momentUnits;
        end
        switch lower(oldUnit)
            case {'kin','k-in'};     
                toBase = 1/moment_base2kin;
            case {'kft','k-ft'};     
                toBase = 1/moment_base2kft;
            case {'nmm','n-mm'};     
                toBase = 1/moment_base2Nmm;
            case {'knm','kn-m'};     
                toBase = 1/moment_base2kNm;
            case {'mnm','mn-m'};     
                toBase = 1/moment_base2MNm;
            case {'tfm','tf-m'};     
                toBase = 1/moment_base2tfm;
            otherwise;                  
                error('old %s unit not recgonized: %s',type,oldUnit);
        end
        switch lower(newUnit)
            case {'kin','k-in'};     
                fromBase = moment_base2kin;
            case {'kft','k-ft'};     
                fromBase = moment_base2kft;
            case {'nmm','n-mm'};     
                fromBase = moment_base2Nmm;
            case {'knm','kn-m'};     
                fromBase = moment_base2kNm;
            case {'mnm','mn-m'};     
                fromBase = moment_base2MNm;                
            case {'tfm','tf-m'};     
                fromBase = moment_base2tfm;
            otherwise;              
                error('new %s unit not recgonized: %s',type,newUnit);
        end
        
    case {'pressure','stress'}
        if isa(oldUnit,'unitSystem')
            oldUnit = oldUnit.pressureUnits;
        end
        if isa(newUnit,'unitSystem')
            newUnit = newUnit.pressureUnits;
        end
        switch lower(oldUnit)
            case 'psi';     toBase = 1/pressure_base2psi;
            case 'psf';     toBase = 1/pressure_base2psf;
            case 'ksi';     toBase = 1/pressure_base2ksi;
            case 'kpa';     toBase = 1/pressure_base2kPa;
            case 'mpa';     toBase = 1/pressure_base2MPa;
            case 'gpa';     toBase = 1/pressure_base2GPa;
            case 'kgscm';   toBase = 1/pressure_base2kgscm;
            case {'tscm','metricton/cm^2'}
                            toBase = 1/pressure_base2tscm;
            case 'longton/in^2'; 
                            toBase = 1/pressure_base2ltsi;
            case 'shortton/in^2'; 
                            toBase = 1/pressure_base2stsi;                         
            otherwise;      error('old %s unit not recgonized: %s',type,oldUnit);
        end
        switch lower(newUnit)
            case 'psi';     fromBase = pressure_base2psi;
            case 'psf';     fromBase = pressure_base2psf;
            case 'ksi';     fromBase = pressure_base2ksi;
            case 'kpa';     fromBase = pressure_base2kPa;
            case 'mpa';     fromBase = pressure_base2MPa;
            case 'gpa';     fromBase = pressure_base2GPa;
            case 'kgscm';   fromBase = pressure_base2kgscm;
            case {'tscm','metricton/cm^2'};
                            fromBase = pressure_base2tscm;
            case 'longton/in^2';
                            fromBase = pressure_base2ltsi;
            case 'shortton/in^2';
                            fromBase = pressure_base2stsi;                            
            otherwise;      error('new %s unit not recgonized: %s',type,newUnit);
        end
        
    case 'density'
        if isa(oldUnit,'unitSystem')
            oldUnit = oldUnit.densityUnits;
        end
        if isa(newUnit,'unitSystem')
            newUnit = newUnit.densityUnits;
        end
        switch lower(oldUnit)
            case 'pcf';     toBase = 1/density_base2pcf;
            case 'kgcm';    toBase = 1/density_base2kgcm;
            otherwise;      error('old %s unit not recgonized: %s',type,oldUnit);
        end
        switch lower(newUnit)
            case 'pcf';     fromBase = density_base2pcf;
            case 'kgcm';    fromBase = density_base2kgcm;
            otherwise;      error('new %s unit not recgonized: %s',type,newUnit);
        end
        
    case 'angle'
        if isa(oldUnit,'unitSystem')
            oldUnit = oldUnit.angleUnits;
        end
        if isa(newUnit,'unitSystem')
            newUnit = newUnit.angleUnits;
        end
        switch lower(oldUnit)
            case 'rad';     toBase = 1/angle_base2rad;
            case 'deg';     toBase = 1/angle_base2deg;
            otherwise;      error('old %s unit not recgonized: %s',type,oldUnit);
        end
        switch lower(newUnit)
            case 'rad';     fromBase = angle_base2rad;
            case 'deg';     fromBase = angle_base2deg;
            otherwise;      error('new %s unit not recgonized: %s',type,newUnit);
        end
        
    case 'volume'
        if isa(oldUnit,'unitSystem')
            oldUnit = oldUnit.volumeUnits;
        end
        if isa(newUnit,'unitSystem')
            newUnit = newUnit.volumeUnits;
        end
        switch lower(oldUnit)
            case 'cbin';    toBase = 1/volume_base2cbin;
            case 'cbft';    toBase = 1/volume_base2cbft;
            case 'cbyd';    toBase = 1/volume_base2cbyd;
            case 'cbm';     toBase = 1/volume_base2cbm;
            case 'cbmm';    toBase = 1/volume_base2cbmm;
            otherwise;      error('old %s unit not recgonized: %s',type,oldUnit);
        end
        switch lower(newUnit)
            case 'cbin';    fromBase = volume_base2cbin;
            case 'cbft';    fromBase = volume_base2cbft;
            case 'cbyd';    fromBase = volume_base2cbyd;
            case 'cbm';     fromBase = volume_base2cbm;
            case 'cbmm';    fromBase = volume_base2cbmm;
            otherwise;      error('new %s unit not recgonized: %s',type,newUnit);
        end
        
    otherwise
        error('Unknown unit type: %s',type);
end

newValue = oldValue*toBase*fromBase;

end

