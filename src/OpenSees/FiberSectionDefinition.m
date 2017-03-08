function [definition,numMat] = FiberSectionDefinition(sectionData,bendingType,sectionID,startMatID,options)
% FiberSectionDefinition is a function to parse MATLAB data
% and returns a string that can be used to define a fiber section in
% OpenSees using the OpenSeesComposite tcl package

% Determine section type
if isa(sectionData,'structural_shape')
    sectionType = sectionData.memberType;
elseif isfield(sectionData,'sectionType')
    sectionType = sectionData.sectionType;
else
    error('Cannot determine section type');
end

% Massage Input Data
if isnumeric(sectionID)
    sectionID = sprintf('%i',sectionID);
end

if strcmpi(bendingType,'strong')
    bendingType = '2dStrong';
end

if strcmpi(bendingType,'weak')
    bendingType = '2dWeak';
end

% Create the section definition
switch lower(sectionType)    
    case 'ccft'
        [definition,numMat] = FiberSectionDefinition_ccft(...
            sectionData,bendingType,sectionID,startMatID,options);
    case 'rcft'
        [definition,numMat] = FiberSectionDefinition_rcft(...
            sectionData,bendingType,sectionID,startMatID,options);
    case 'src'
        [definition,numMat] = FiberSectionDefinition_src(...
            sectionData,bendingType,sectionID,startMatID,options);
    case 'wf'
        [definition,numMat] = FiberSectionDefinition_wf(...
            sectionData,bendingType,sectionID,startMatID,options);
    case 'recthss'
        [definition,numMat] = FiberSectionDefinition_recthss(...
            sectionData,bendingType,sectionID,startMatID,options);
    case 'roundhss'
        [definition,numMat] = FiberSectionDefinition_roundhss(...
            sectionData,bendingType,sectionID,startMatID,options);
    case 'elastic'
        [definition,numMat] = FiberSectionDefinition_elastic(...
            sectionData,bendingType,sectionID,startMatID,options);
    otherwise
        error('Unknown sectionType')
end

% Include Package Definition (if requested)
if isfield(options,'includePackageDefinition')
    includePackageDefinition = options.includePackageDefinition;
else
    includePackageDefinition = false;
end

if includePackageDefinition
    packageDefinition = {...
        'package require OpenSeesComposite',...
        'namespace import OpenSeesComposite::*'}';
    if ischar(definition)
        definition = {definition};
    end
    definition = vertcat(packageDefinition,definition);
end

if nargout < 2
    clear numMat
end

end
    

function [definition,numMat] = FiberSectionDefinition_ccft(sectionData,bendingType,sectionID,startMatID,options)

[nf1,nf2] = parseFiberDensity(bendingType,options);

if isa(sectionData,'CCFT')
    units = sectionData.units;
    D  = sprintf('%g',sectionData.D);
    t  = sprintf('%g',sectionData.t);
    Fy = sprintf('%g',sectionData.Fy);
    if isempty(sectionData.Fu)
        Fu = 'calc';
    else
        Fu = sprintf('%g',sectionData.Fu);
    end
    Es = sprintf('%g',sectionData.Es);
    fc = sprintf('%g',sectionData.fc);
elseif isstruct(sectionData)
    units = sectionData.units;
    D  = parseFromStruct(sectionData,'D','%g');
    t  = parseFromStruct(sectionData,'t','%g');
    Fy = parseFromStruct(sectionData,'Fy','%g');
    Fu = parseFromStruct(sectionData,'Fu','%g','default','calc','replace',-1,'calc');
    Es = parseFromStruct(sectionData,'Es','%g','default','calc','replace',-1,'calc');
    fc = parseFromStruct(sectionData,'fc','%g');
else
    error('Unknown type for sectionData: %s',class(sectionData))
end

definition = sprintf(...
    'ccftSection %s %i %s %s %s %s %s %s %s %s %s',...
    sectionID,startMatID,nf1,nf2,units,D,t,Fy,Fu,Es,fc);

definition = parseSteelMaterialType(definition,options);
definition = parseConcreteMaterialType(definition,options);

% Added Elastic Stiffness
definition = parseAddedElasticStiffness(definition,bendingType,options);

% Number of materials
numMat = 2;

end

function [definition,numMat] = FiberSectionDefinition_rcft(sectionData,bendingType,sectionID,startMatID,options)

[nf1,nf2] = parseFiberDensity(bendingType,options);

if isa(sectionData,'RCFT')
    units = sectionData.units;
    H  = sprintf('%g',sectionData.H);
    B  = sprintf('%g',sectionData.B);
    t  = sprintf('%g',sectionData.t);
    Fy = sprintf('%g',sectionData.Fy);
    if isempty(sectionData.Fu)
        Fu = 'calc';
    else
        Fu = sprintf('%g',sectionData.Fu);
    end
    Es = sprintf('%g',sectionData.Es);
    fc = sprintf('%g',sectionData.fc);
elseif isstruct(sectionData)
    units = sectionData.units;
    H  = parseFromStruct(sectionData,'H','%g');
    B  = parseFromStruct(sectionData,'B','%g');
    t  = parseFromStruct(sectionData,'t','%g');
    Fy = parseFromStruct(sectionData,'Fy','%g');
    Fu = parseFromStruct(sectionData,'Fu','%g','replace',-1,'calc');
    Es = parseFromStruct(sectionData,'Es','%g','replace',-1,'calc');
    fc = parseFromStruct(sectionData,'fc','%g');
else
    error('Unknown type for sectionData: %s',class(sectionData))
end

definition = sprintf(...
    'rcftSection %s %i %s %s %s %s %s %s %s %s %s %s',...
    sectionID,startMatID,nf1,nf2,units,H,B,t,Fy,Fu,Es,fc);

definition = parseSteelMaterialType(definition,options);
definition = parseConcreteMaterialType(definition,options);

% Added Elastic Stiffness
definition = parseAddedElasticStiffness(definition,bendingType,options);

% Number of materials
numMat = 3;

end

function [definition,numMat] = FiberSectionDefinition_src(sectionData,bendingType,sectionID,startMatID,options)

[nf1,nf2] = parseFiberDensity(bendingType,options);

if isa(sectionData,'SRC')
    units = sectionData.units;
    B  = sprintf('%g',sectionData.B);
    H  = sprintf('%g',sectionData.H);
    fc = sprintf('%g',sectionData.fc);
    d  = sprintf('%g',sectionData.d);
    tw = sprintf('%g',sectionData.tw);
    bf = sprintf('%g',sectionData.bf);
    tf = sprintf('%g',sectionData.tf);
    Fy = sprintf('%g',sectionData.Fy);
    if isempty(sectionData.Fu)
        Fu = 'calc';
    else
        Fu = sprintf('%g',sectionData.Fu);
    end
    Es = sprintf('%g',sectionData.Es);
    
    % Reinforcement
    reinf = sprintf('"%s %g %g %s %s %g %g %g %g"',...
        sectionData.rebarConfig,...
        sectionData.db,...
        sectionData.Fylr,...
        'calc',...
        'calc',...
        sectionData.dbTies,...
        sectionData.s,...
        sectionData.Fytr,...
        sectionData.cover);
    
elseif isstruct(sectionData)
    units = sectionData.units;
    B  = parseFromStruct(sectionData,'B','%g');
    H  = parseFromStruct(sectionData,'H','%g');
    fc = parseFromStruct(sectionData,'fc','%g');
    
    d  = parseFromStruct(sectionData,'d','%g');
    tw = parseFromStruct(sectionData,'tw','%g');
    bf = parseFromStruct(sectionData,'bf','%g');
    tf = parseFromStruct(sectionData,'tf','%g');
    
    Fy = parseFromStruct(sectionData,'Fy','%g');
    Fu = parseFromStruct(sectionData,'Fu','%g','replace',-1,'calc');
    Es = parseFromStruct(sectionData,'Es','%g','replace',-1,'calc');
    
    % Reinforcement
    if strcmpi(sectionData.rebarConfig,'none')
        reinf = 'none';
    else
        Fylr = parseFromStruct(sectionData,'Fylr','%g');
        Fulr = parseFromStruct(sectionData,'Fulr','%g','replace',-1,'calc');
        Eslr = parseFromStruct(sectionData,'Eslr','%g','replace',-1,'calc');
        
        reinf = sprintf('"%s %g %s %s %s %g %g %g %g"',...
            sectionData.rebarConfig,...
            sectionData.db,...
            Fylr,Fulr,Eslr,...
            sectionData.dbTies,...
            sectionData.s,...
            sectionData.Fytr,...
            sectionData.cover);
    end
    
else
    error('Unknown type for sectionData: %s',class(sectionData))
end

definition = sprintf(...
    'srcSection %s %i %s %s %s %s %s %s %s %s %s %s %s %s %s %s',...
    sectionID,startMatID,nf1,nf2,units,B,H,fc,d,tw,bf,tf,Fy,Fu,Es,reinf);

definition = parseSteelMaterialType(definition,options);
definition = parseConcreteMaterialType(definition,options);
definition = parseReinforcementMaterialType(definition,options);

% Added Elastic Stiffness
definition = parseAddedElasticStiffness(definition,bendingType,options);

% Number of materials
numMat = 14;

end


function [definition,numMat] = FiberSectionDefinition_pec(sectionData,bendingType,sectionID,startMatID,options)

[nf1,nf2] = parseFiberDensity(bendingType,options);

if isstruct(sectionData)
    units = sectionData.units;
    d  = parseFromStruct(sectionData,'d','%g');
    tw = parseFromStruct(sectionData,'tw','%g');
    bf = parseFromStruct(sectionData,'bf','%g');
    tf = parseFromStruct(sectionData,'tf','%g');
    
    Fy = parseFromStruct(sectionData,'Fy','%g');
    Fu = parseFromStruct(sectionData,'Fu','%g','replace',-1,'calc');
    Es = parseFromStruct(sectionData,'Es','%g','replace',-1,'calc');
    
    fc = parseFromStruct(sectionData,'fc','%g');
    
    % Reinforcement
    if strcmpi(rebarConfig,'none')
        reinf = 'none';
    else
        error('Unknown reinf configuration');
    end
    
else
    error('Unknown type for sectionData: %s',class(sectionData))
end

definition = sprintf(...
    'pecSection %s %i %s %s %s %s %s %s %s %s %s %s %s %s',...
    sectionID,startMatID,nf1,nf2,units,d,tw,bf,tf,Fy,Fu,Es,fc,reinf);

definition = parseSteelMaterialType(definition,options);
definition = parseConcreteMaterialType(definition,options);
definition = parseReinforcementMaterialType(definition,options);

% Added Elastic Stiffness
definition = parseAddedElasticStiffness(definition,bendingType,options);

% Number of materials
numMat = 14;

end

function [definition,numMat] = FiberSectionDefinition_recthss(sectionData,bendingType,sectionID,startMatID,options)

[nf1,nf2] = parseFiberDensity(bendingType,options);

if isa(sectionData,'RectangularHSS')
    units = sectionData.units;
    H  = sprintf('%g',sectionData.H);
    B  = sprintf('%g',sectionData.B);
    t  = sprintf('%g',sectionData.t);
    Fy = sprintf('%g',sectionData.Fy);
    if isempty(sectionData.Fu)
        Fu = 'calc';
    else
        Fu = sprintf('%g',sectionData.Fu);
    end
    Es = sprintf('%g',sectionData.Es);
elseif isstruct(sectionData)
    units = sectionData.units;
    H  = parseFromStruct(sectionData,'H','%g');
    B  = parseFromStruct(sectionData,'B','%g');
    t  = parseFromStruct(sectionData,'t','%g');
    Fy = parseFromStruct(sectionData,'Fy','%g');
    Fu = parseFromStruct(sectionData,'Fu','%g','replace',-1,'calc');
    Es = parseFromStruct(sectionData,'Es','%g','replace',-1,'calc');
else
    error('Unknown type for sectionData: %s',class(sectionData))
end

definition = sprintf(...
    'recthssSection %s %i %s %s %s %s %s %s %s %s %s',...
    sectionID,startMatID,nf1,nf2,units,H,B,t,Fy,Fu,Es);

definition = parseSteelMaterialType(definition,options);

% Added Elastic Stiffness
definition = parseAddedElasticStiffness(definition,bendingType,options);

% Number of materials
numMat = 2;

end

function [definition,numMat] = FiberSectionDefinition_roundhss(sectionData,bendingType,sectionID,startMatID,options)

[nf1,nf2] = parseFiberDensity(bendingType,options);

if isa(sectionData,'RoundHSS')
    units = sectionData.units;
    D  = sprintf('%g',sectionData.D);
    t  = sprintf('%g',sectionData.t);
    Fy = sprintf('%g',sectionData.Fy);
    if isempty(sectionData.Fu)
        Fu = 'calc';
    else
        Fu = sprintf('%g',sectionData.Fu);
    end
    Es = sprintf('%g',sectionData.Es);
elseif isstruct(sectionData)
    units = sectionData.units;
    D  = parseFromStruct(sectionData,'D','%g');
    t  = parseFromStruct(sectionData,'t','%g');
    Fy = parseFromStruct(sectionData,'Fy','%g');
    Fu = parseFromStruct(sectionData,'Fu','%g','replace',-1,'calc');
    Es = parseFromStruct(sectionData,'Es','%g','replace',-1,'calc');
else
    error('Unknown type for sectionData: %s',class(sectionData))
end

definition = sprintf(...
    'roundhssSection %s %i %s %s %s %s %s %s %s %s',...
    sectionID,startMatID,nf1,nf2,units,D,t,Fy,Fu,Es);

definition = parseSteelMaterialType(definition,options);

% Added Elastic Stiffness
definition = parseAddedElasticStiffness(definition,bendingType,options);

% Number of materials
numMat = 1;

end

function [definition,numMat] = FiberSectionDefinition_wf(sectionData,bendingType,sectionID,startMatID,options)

[nf1,nf2] = parseFiberDensity(bendingType,options);

if isa(sectionData,'WF')
    units = sectionData.units;
    d  = sprintf('%g',sectionData.d);
    tw = sprintf('%g',sectionData.tw);
    bf = sprintf('%g',sectionData.bf);
    tf = sprintf('%g',sectionData.tf);
    k  = sprintf('%g',sectionData.k);
    Fy = sprintf('%g',sectionData.Fy);
    Fu = sprintf('%g',sectionData.Fu);
    eu = sprintf('%g',sectionData.eu);
    Es = sprintf('%g',sectionData.Es);
elseif isstruct(sectionData)
    units = sectionData.units;
    d  = parseFromStruct(sectionData,'d','%g');
    bf = parseFromStruct(sectionData,'bf','%g');
    tw = parseFromStruct(sectionData,'tw','%g');
    tf = parseFromStruct(sectionData,'tf','%g');
    Fy = parseFromStruct(sectionData,'Fy','%g');
    Fu = parseFromStruct(sectionData,'Fu','%g');
    eu = parseFromStruct(sectionData,'eu','%g');
    Es = parseFromStruct(sectionData,'Es','%g');
else
    error('Unknown type for sectionData: %s',class(sectionData))
end

definition = sprintf(...
    'wfSection %s %s %s %s %s %s %s',...
    sectionID,nf1,nf2,d,tw,bf,tf);

% Fillet
if isfield(options,'includeFillet')
    includeFillet = options.includeFillet;
else
    includeFillet = false;
end
if includeFillet
    definition = sprintf('%s -Fillet %s',definition,k);
end

% Residual Stress
nResidualStressSectors = ...
    parseFromStruct(sectionData,'nResidualStressSectors','%i',...
    'default','10');
if isfield(options,'residualStressRatio')
    residualStressRatio = options.residualStressRatio;
else
    residualStressRatio = 0.3;
end
frc = parseFromStruct(sectionData,'frc','%g',...
    'default',sprintf('%g',-residualStressRatio*sscanf(Fy,'%g')));
definition = sprintf('%s -Lehigh %s %s',...
    definition,frc,nResidualStressSectors);

% Material Models
if isfield(options,'SteelMaterialType')
    materialType = options.SteelMaterialType;
else
    materialType = 'Shen';
end
switch materialType
    case 'Shen'
        definition = sprintf('%s -ShenSteel %i %s %s %s %s %s',...
            definition,startMatID,Es,Fy,Fu,eu,units);
    case 'ShenDegrade'
        definition = sprintf('%s -ShenSteelDegrade %i %s %s %s %s %s',...
            definition,startMatID,Es,Fy,Fu,eu,units);
    case 'ElasticPP'
        definition = sprintf('%s -ElasticPP %i %s %s',...
            definition,startMatID,Es,Fy);
    case 'Steel02'
        if isfield(options,'StrainHardeningRatio')
            b = options.StrainHardeningRatio;
        else
            b = 0.003;
        end
        definition = sprintf('%s -Steel02 %i %s %s %g',...
            definition,startMatID,Es,Fy,b);
    case 'Elastic'
        definition = sprintf('%s -Elastic %i %s',...
            definition,startMatID,Es);    
    case 'ElasticSmallStiffness'
        definition = sprintf('%s -ElasticSmallStiffness %i %s %s',...
            definition,startMatID,Es,Fy);
    otherwise
        error('Unknown material type')
end

% Added Elastic Stiffness
definition = parseAddedElasticStiffness(definition,bendingType,options);

% Number of materials
numMat = nResidualStressSectors;

end


function [definition,numMat] = FiberSectionDefinition_elastic(sectionData,bendingType,sectionID,startMatID,options)

switch lower(bendingType)
    case '3d'
        if isa(sectionData,'general_elastic_section')
            E  = sprintf('%g',sectionData.E);
            A  = sprintf('%g',sectionData.A);
            Iz = sprintf('%g',sectionData.Iz);
            Iy = sprintf('%g',sectionData.Iy);
            GJ = sprintf('%g',sectionData.GJ);
        else
            error('Unknown type for sectionData: %s',class(sectionData))
        end
        
        definition{1} = sprintf(...
            'uniaxialMaterial Elastic %i %s',...
            startMatID,E);
        definition{2} = sprintf(...
            'fourFiberSectionGJ %s %i %g %g %g %g',...
            sectionID,startMatID,A,Iy,Iz,GJ);
        
    case '2dstrong'
        if isa(sectionData,'general_elastic_section')
            E  = sprintf('%g',sectionData.E);
            A  = sprintf('%g',sectionData.A);
            I  = sprintf('%g',sectionData.Iz);
        else
            error('Unknown type for sectionData: %s',class(sectionData))
        end
        
        definition{1} = sprintf(...
            'uniaxialMaterial Elastic %i %s',...
            startMatID,E);
        definition{2} = sprintf(...
            'twoFiberSection %s %i %s %s',...
            sectionID,startMatID,A,I);
        
    case '2dweak'
        if isa(sectionData,'general_elastic_section')
            E  = sprintf('%g',sectionData.E);
            A  = sprintf('%g',sectionData.A);
            I  = sprintf('%g',sectionData.Iy);
        else
            error('Unknown type for sectionData: %s',class(sectionData))
        end
        
        definition{1} = sprintf(...
            'uniaxialMaterial Elastic %i %s',...
            startMatID,E);
        definition{2} = sprintf(...
            'twoFiberSectionGJ %s %i %g %g',...
            sectionID,startMatID,A,I);
        
    otherwise
        error('Unknown bending type')
end

% Number of materials
numMat = 1;

end


function definition = parseSteelMaterialType(definition,options)
if isfield(options,'SteelMaterialType')
    type = options.SteelMaterialType;
    switch type
        case 'ModifiedAbdelRahman'
            definition = sprintf('%s -SteelMaterialType %s %g %g',definition,type,...
                options.AbdelRahmanResidualStressParameter,...
                options.AbdelRahmanHardeningRatio);
        otherwise
            definition = sprintf('%s -SteelMaterialType %s',definition,type);
    end
end
end

function definition = parseConcreteMaterialType(definition,options)
if isfield(options,'ConcreteMaterialType')
    type = options.ConcreteMaterialType;
    definition = sprintf('%s -ConcreteMaterialType %s',definition,type);
end
end

function definition = parseReinforcementMaterialType(definition,options)
if isfield(options,'ReinforcementMaterialType')
    type = options.ReinforcementMaterialType;
    definition = sprintf('%s -ReinforcementMaterialType %s',definition,type);
end
end

function str = parseFromStruct(obj,field,format,varargin)

defaultValue = '';
checkForReplace = false;
replaceValue = [];
replaceWith = '';

i = 1;
while i < length(varargin)
    switch lower(varargin{i})
        case 'default'
            defaultValue = varargin{i+1};
            i = i+1;
        case 'replace'
            checkForReplace = true;
            replaceValue = varargin{i+1};
            replaceWith  = varargin{i+2};
            i = i+2;
        otherwise
            error('Unknown Option');
    end
    i = i+1;
end

if isfield(obj,field);
    value = obj.(field);
    if checkForReplace && isequal(value,replaceValue)
        str = replaceWith;
    else
        str = sprintf(format,value);
    end
else
    if ~isempty(defaultValue)
        str = defaultValue;
    else
        error('Field "%s" does not exist in structure',field);
    end
end


end

function [nf1,nf2] = parseFiberDensity(bendingType,options)
default_nf1 = 30;
default_nf2 = 30;

% nf1
if isfield(options,'nf1')
    nf1 = sprintf('%i',options.nf1);
else
    nf1 = sprintf('%i',default_nf1);
end
% nf2
switch lower(bendingType)
    case '2dstrong'
        nf2 = 'strong';
    case '2dweak'
        nf2 = 'weak';
    case '3d'
        if isfield(options,'nf2')
            nf2 = sprintf('%i',options.nf2);
        else
            nf2 = sprintf('%i',default_nf2);
        end
    otherwise
        error('Unknown bendingType');
end
end


function definition = parseAddedElasticStiffness(definition,bendingType,options)

includeAddedElasticStiffness = false;

if isfield(options,'includeAddedElasticStiffness')
    includeAddedElasticStiffness = options.includeAddedElasticStiffness;
elseif isfield(options,'addedElasticStiffness_include')
    includeAddedElasticStiffness = options.addedElasticStiffness_include;
end

if includeAddedElasticStiffness
    switch lower(bendingType)
        case {'2dstrong','2dweak'}
            EA = parseFromStruct(options,'addedElasticStiffness_EA','%g');
            EI = parseFromStruct(options,'addedElasticStiffness_EI','%g');
            definition = sprintf('%s -AddedElastic %s %s',...
                definition,EA,EI);
        case '3d'
            EA  = parseFromStruct(options,'addedElasticStiffness_EA','%g');
            EIz = parseFromStruct(options,'addedElasticStiffness_EIz','%g');
            EIy = parseFromStruct(options,'addedElasticStiffness_EIy','%g');
            GJ  = parseFromStruct(options,'addedElasticStiffness_GJ','%g');
            definition = sprintf('%s -AddedElastic %s %s %s %s',...
                definition,EA,EIz,EIy,GJ);
        otherwise
            error('Unknown bendingType');
    end
end

end
