function [definition,numMat] = FiberSectionDefinition(sectionData,bendingType,sectionID,startMatID,options)
% FiberSectionDefinition is a function to parse MATLAB data
% and returns a string that can be used to define a fiber section in
% OpenSees using the OpenSeesComposite tcl package

% Determine section type
if isa(sectionData,'structural_shape')
    sectionType = sectionData.memberType;
elseif isfield(sectionData,'sectionType')
    sectionType = sectionData.sectionType;
elseif isfield(sectionData,'section_type')
    sectionType = sectionData.section_type;
else
    error('Cannot determine section type');
end

% Massage Input Data
if isnumeric(sectionID)
    sectionID = sprintf('%i',sectionID);
end

if any(strcmpi(bendingType,{'strong','x','z'}))
    bendingType = '2dx';
end

if any(strcmpi(bendingType,{'weak','y'}))
    bendingType = '2dy';
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
    case 'rc'
        [definition,numMat] = FiberSectionDefinition_RC(...
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
    if isfield(options,'StrengthReduction')
        Fy = sprintf('%g',options.StrengthReduction*sectionData.Fy);
        if isempty(sectionData.Fu)
            Fu = 'calc';
        else
            Fu = sprintf('%g',options.StrengthReduction*sectionData.Fu);
        end
    else
        Fy = sprintf('%g',sectionData.Fy);
        if isempty(sectionData.Fu)
            Fu = 'calc';
        else
            Fu = sprintf('%g',sectionData.Fu);
        end
    end
    if isfield(options,'StiffnessReduction')
        E = sprintf('%g',options.StiffnessReduction*sectionData.E);
    else
        E = sprintf('%g',sectionData.E);
    end    
elseif isstruct(sectionData)
    units = sectionData.units;
    H  = parseFromStruct(sectionData,'H','%g');
    B  = parseFromStruct(sectionData,'B','%g');
    t  = parseFromStruct(sectionData,'t','%g');
    if isfield(options,'StrengthReduction')
        Fy = parseFromStruct(sectionData,'Fy','%g','reduce',options.StrengthReduction);
        Fu = parseFromStruct(sectionData,'Fu','%g','reduce',options.StrengthReduction,'replace',-1,'calc');
    else
        Fy = parseFromStruct(sectionData,'Fy','%g');
        Fu = parseFromStruct(sectionData,'Fu','%g','replace',-1,'calc');
    end
    if isfield(options,'StiffnessReduction')
        E = parseFromStruct(sectionData,'E','%g','reduce',options.StiffnessReduction,'replace',-1,'calc');
    else
        E = parseFromStruct(sectionData,'E','%g','replace',-1,'calc');
    end
else
    error('Unknown type for sectionData: %s',class(sectionData))
end

definition = sprintf(...
    'recthssSection %s %i %s %s %s %s %s %s %s %s %s',...
    sectionID,startMatID,nf1,nf2,units,H,B,t,Fy,Fu,E);

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
    if isfield(options,'StrengthReduction')
        Fy = sprintf('%g',options.StrengthReduction*sectionData.Fy);
        if isempty(sectionData.Fu)
            Fu = 'calc';
        else
            Fu = sprintf('%g',options.StrengthReduction*sectionData.Fu);
        end
    else
        Fy = sprintf('%g',sectionData.Fy);
        if isempty(sectionData.Fu)
            Fu = 'calc';
        else
            Fu = sprintf('%g',sectionData.Fu);
        end
    end
    if isfield(options,'StiffnessReduction')
        E = sprintf('%g',options.StiffnessReduction*sectionData.E);
    else
        E = sprintf('%g',sectionData.E);
    end 
elseif isstruct(sectionData)
    units = sectionData.units;
    D  = parseFromStruct(sectionData,'D','%g');
    t  = parseFromStruct(sectionData,'t','%g');
    if isfield(options,'StrengthReduction')
        Fy = parseFromStruct(sectionData,'Fy','%g','reduce',options.StrengthReduction);
        Fu = parseFromStruct(sectionData,'Fu','%g','reduce',options.StrengthReduction,'replace',-1,'calc');
    else
        Fy = parseFromStruct(sectionData,'Fy','%g');
        Fu = parseFromStruct(sectionData,'Fu','%g','replace',-1,'calc');
    end
    if isfield(options,'StiffnessReduction')
        E = parseFromStruct(sectionData,'E','%g','reduce',options.StiffnessReduction,'replace',-1,'calc');
    else
        E = parseFromStruct(sectionData,'E','%g','replace',-1,'calc');
    end
else
    error('Unknown type for sectionData: %s',class(sectionData))
end

definition = sprintf(...
    'roundhssSection %s %i %s %s %s %s %s %s %s %s',...
    sectionID,startMatID,nf1,nf2,units,D,t,Fy,Fu,E);

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
    if isfield(options,'StrengthReduction')
        Fy = sprintf('%g',options.StrengthReduction*sectionData.Fy);
        Fu = sprintf('%g',options.StrengthReduction*sectionData.Fu);
    else
        Fy = sprintf('%g',sectionData.Fy);
        Fu = sprintf('%g',sectionData.Fu);
    end
    eu = sprintf('%g',sectionData.eu);
    if isfield(options,'StiffnessReduction')
        E = sprintf('%g',options.StiffnessReduction*sectionData.E);
    else
        E = sprintf('%g',sectionData.E);
    end
elseif isstruct(sectionData)
    units = sectionData.units;
    d  = parseFromStruct(sectionData,'d','%g');
    bf = parseFromStruct(sectionData,'bf','%g');
    tw = parseFromStruct(sectionData,'tw','%g');
    tf = parseFromStruct(sectionData,'tf','%g');
    if isfield(options,'StrengthReduction')
        Fy = parseFromStruct(sectionData,'Fy','%g','reduce',options.StrengthReduction);
        Fu = parseFromStruct(sectionData,'Fu','%g','reduce',options.StrengthReduction);
    else
        Fy = parseFromStruct(sectionData,'Fy','%g');
        Fu = parseFromStruct(sectionData,'Fu','%g');
    end
    eu = parseFromStruct(sectionData,'eu','%g');
    if isfield(options,'StiffnessReduction')
        E = parseFromStruct(sectionData,'E','%g','reduce',options.StiffnessReduction);
    else
        E = parseFromStruct(sectionData,'E','%g');
    end
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
if isfield(options,'include_residual_stress')
    include_residual_stress = options.include_residual_stress;
else
    include_residual_stress = true;
end
if include_residual_stress
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
end

% Material Models
if isfield(options,'SteelMaterialType')
    materialType = options.SteelMaterialType;
else
    materialType = 'Shen';
end
switch materialType
    case 'Shen'
        definition = sprintf('%s -ShenSteel %i %s %s %s %s %s',...
            definition,startMatID,E,Fy,Fu,eu,units);
    case 'ShenDegrade'
        definition = sprintf('%s -ShenSteelDegrade %i %s %s %s %s %s',...
            definition,startMatID,E,Fy,Fu,eu,units);
    case 'ElasticPP'
        definition = sprintf('%s -ElasticPP %i %s %s',...
            definition,startMatID,E,Fy);
    case 'Steel02'
        if isfield(options,'StrainHardeningRatio')
            b = options.StrainHardeningRatio;
        else
            b = 0.003;
        end
        definition = sprintf('%s -Steel02 %i %s %s %g',...
            definition,startMatID,E,Fy,b);
    case 'Elastic'
        definition = sprintf('%s -Elastic %i %s',...
            definition,startMatID,E);    
    case 'ElasticSmallStiffness'
        definition = sprintf('%s -ElasticSmallStiffness %i %s %s',...
            definition,startMatID,E,Fy);
    otherwise
        error('Unknown material type')
end

% Added Elastic Stiffness
definition = parseAddedElasticStiffness(definition,bendingType,options);

% Number of materials
if include_residual_stress
    numMat = nResidualStressSectors;
else
    numMat = 1;
end

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
        elseif isstruct(sectionData)
            E  = parseFromStruct(sectionData,'E','%g');
            A  = parseFromStruct(sectionData,'A','%g');
            Iz = parseFromStruct(sectionData,'Iz','%g');
            Iy = parseFromStruct(sectionData,'Iy','%g');
            GJ = parseFromStruct(sectionData,'GJ','%g');
        else
            error('Unknown type for sectionData: %s',class(sectionData))
        end
        
        definition{1} = sprintf(...
            'uniaxialMaterial Elastic %i %s',...
            startMatID,E);
        definition{2} = sprintf(...
            'fourFiberSectionGJ %s %i %s %s %s %s',...
            sectionID,startMatID,A,Iy,Iz,GJ);
        
    case '2dx'
        if isa(sectionData,'general_elastic_section')
            E  = sprintf('%g',sectionData.E);
            A  = sprintf('%g',sectionData.A);
            I  = sprintf('%g',sectionData.Iz);
        elseif isstruct(sectionData)
            E  = parseFromStruct(sectionData,'E','%g');
            A  = parseFromStruct(sectionData,'A','%g');
            I  = parseFromStruct(sectionData,'I','%g');
        else
            error('Unknown type for sectionData: %s',class(sectionData))
        end
        
        definition{1} = sprintf(...
            'uniaxialMaterial Elastic %i %s',...
            startMatID,E);
        definition{2} = sprintf(...
            'twoFiberSection %s %i %s %s',...
            sectionID,startMatID,A,I);
        
    case '2dy'
        if isa(sectionData,'general_elastic_section')
            E  = sprintf('%g',sectionData.E);
            A  = sprintf('%g',sectionData.A);
            I  = sprintf('%g',sectionData.Iy);
        elseif isstruct(sectionData)
            E  = parseFromStruct(sectionData,'E','%g');
            A  = parseFromStruct(sectionData,'A','%g');
            I  = parseFromStruct(sectionData,'I','%g');
        else
            error('Unknown type for sectionData: %s',class(sectionData))
        end
        
        definition{1} = sprintf(...
            'uniaxialMaterial Elastic %i %s',...
            startMatID,E);
        definition{2} = sprintf(...
            'twoFiberSectionGJ %s %i %s %s',...
            sectionID,startMatID,A,I);
        
    otherwise
        error('Unknown bending type')
end

% Number of materials
numMat = 1;

end


function [definition,numMat] = FiberSectionDefinition_RC(sectionData,bendingType,sectionID,startMatID,options)

[nf1,nf2] = parseFiberDensity(bendingType,options);

if isa(sectionData,'RC')

    switch class(sectionData.conc_cross_section)
        case 'Rectangle_Shape'
            assert(length(sectionData.reinforcement)==1 && isa(sectionData.reinforcement{1},'reinf_rect'),...
                'reinforcement should consist of only a single reinf_rect object');            
            assert(sectionData.reinforcement{1}.xc == 0 && sectionData.reinforcement{1}.yc == 0,...
                'reinforcement should be centered');
            assert(sectionData.conc_cross_section.B - sectionData.reinforcement{1}.Bx == sectionData.conc_cross_section.H - sectionData.reinforcement{1}.By,...
                'cover should be equal on both sides');
                        
            error('This section still needs work')
            definition = sprintf(...
                'RectangularRC %s %i %s %s %s %g %g %g %g %g %g NA %g %g %g %i %i %g %g %g %i %i %i %s %s',...
                sectionID,startMatID,nf1,nf2,sectionData.units,...
                sectionData.conc_cross_section.B,sectionData.conc_cross_section.H,...
                sectionData.fc,sectionData.Ec,dp,...
                sectionData.fy,sectionData.Es,...
                sectionData.reinforcement{1}.db,sectionData.reinforcement{1}.Ab,...
                sectionData.reinforcement{1}.nbx,sectionData.reinforcement{1}.nby,...
                sectionData.fyt,sectionData.dbt,sectionData.Abt,...
                sectionData.nLegX,sectionData.nLegY,...
                sectionData.s,...
                options.conc_material,options.steel_material);

        case 'Circle_Shape'
            error('This section still needs work')
            definition = sprintf(...
                'CircularRC %s %i %s %s %s %g %g %g %g %g %g %g %g %g %g %i %i %g %g %g %i %i %i %s %s',...
                sectionID,startMatID,nf1,nf2,sectionData.units,...
                sectionData.conc_cross_section.D,...
                sectionData.fc,sectionData.Ec,sectionData.cover,...
                sectionData.fy,sectionData.fu,sectionData.Es,...
                sectionData.db,sectionData.Ab,...
                sectionData.nBar,...
                sectionData.fyt,sectionData.dbt,sectionData.Abt,...
                sectionData.s,sectionData.transverse_reinf_type,...
                options.conc_material,options.steel_material);

        otherwise
            error('Unknown type for sectionData.conc_cross_section: %s',class(sectionData.conc_cross_section));
    end
    
elseif isstruct(sectionData)

    switch lower(sectionData.shape)
        case {'square','rectangle'}
            definition = sprintf(...
                'RectangularRC %s %i %s %s %s %g %g %g %g %g %g %g %g %g %g %i %i %g %g %g %i %i %i %s %s',...
                sectionID,startMatID,nf1,nf2,sectionData.units,...
                sectionData.B,sectionData.H,...
                sectionData.fc,sectionData.Ec,sectionData.cover,...
                sectionData.fy,sectionData.fu,sectionData.Es,...
                sectionData.db,sectionData.Ab,...
                sectionData.nBarX,sectionData.nBarY,...
                sectionData.fyt,sectionData.dbt,sectionData.Abt,...
                sectionData.nLegX,sectionData.nLegY,...
                sectionData.s,...
                options.conc_material,options.steel_material);

        case 'circle'
            definition = sprintf(...
                'CircularRC %s %i %s %s %s %g %g %g %g %g %g %g %g %g %i %g %g %g %g %s %s %s',...
                sectionID,startMatID,nf1,nf2,sectionData.units,...
                sectionData.D,...
                sectionData.fc,sectionData.Ec,sectionData.cover,...
                sectionData.fy,sectionData.fu,sectionData.Es,...
                sectionData.db,sectionData.Ab,...
                sectionData.num_bars,...
                sectionData.fyt,sectionData.dbt,sectionData.Abt,...
                sectionData.s,sectionData.transverse_reinf_type,...
                options.conc_material,options.steel_material);

        otherwise
            error('Unknown shape: %s',sectionData.shape);
    end    
else
    error('Unknown type for sectionData: %s',class(sectionData));
end

% Added Elastic Stiffness
definition = parseAddedElasticStiffness(definition,bendingType,options);

% Number of materials
numMat = 3;

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
reductionFactor = 1;

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
        case 'reduce'
            reductionFactor = varargin{i+1};
            i = i+1;
        otherwise
            error('Unknown Option');
    end
    i = i+1;
end

if isfield(obj,field)
    value = obj.(field);
    if checkForReplace && isequal(value,replaceValue)
        str = replaceWith;
    else
        str = sprintf(format,reductionFactor*value);
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
    case '2dx'
        nf2 = 'x';
    case '2dy'
        nf2 = 'y';
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
        case {'2dx','2dy'}
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
