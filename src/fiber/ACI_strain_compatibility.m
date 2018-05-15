classdef ACI_strain_compatibility < handle
    
    properties 
        strainAtExtremeConcreteFiber = -0.003;
        maxCompressiveStrength
        defaultTensileStrain = 0.005;
        AxesOrigin = 'AsDefined';
    end
    
    properties (SetAccess = private)
        fiberSection 
    end
    
    properties (Access = private)
        concreteBoundary_zPoint = zeros(0,1);
        concreteBoundary_yPoint = zeros(0,1);
        concreteBoundary_radius = zeros(0,1);
        steelBoundary_zPoint = zeros(0,1);
        steelBoundary_yPoint = zeros(0,1);
        steelBoundary_radius = zeros(0,1);
        materials   = cell(0,1);
        materialIDs = zeros(0,1);
    end

    methods
        %% Constructor
        function obj = ACI_strain_compatibility(fiberSection,AxesOrigin)
            if nargin > 0
                obj.fiberSection = fiberSection;
            end
            if nargin > 1
                obj.AxesOrigin = AxesOrigin;
            end
        end
               
        %% Define section
        function addConcreteBoundary(obj,zPoint,yPoint,radius)
            obj.concreteBoundary_zPoint = vertcat(obj.concreteBoundary_zPoint,zPoint);
            obj.concreteBoundary_yPoint = vertcat(obj.concreteBoundary_yPoint,yPoint);
            obj.concreteBoundary_radius = vertcat(obj.concreteBoundary_radius,radius);
        end
        function addSteelBoundary(obj,zPoint,yPoint,radius)
            obj.steelBoundary_zPoint = vertcat(obj.steelBoundary_zPoint,zPoint);
            obj.steelBoundary_yPoint = vertcat(obj.steelBoundary_yPoint,yPoint);
            obj.steelBoundary_radius = vertcat(obj.steelBoundary_radius,radius);
        end
        function addMaterial(obj,type,varargin)
            if isa(type,'ACI_strain_compatibility_material')
                mat = type;
            elseif ischar(type)
                switch lower(type)
                    case 'concrete'
                        mat = ACI_strain_compatibility_material_concrete(varargin{:});
                    case 'steel'
                        mat = ACI_strain_compatibility_material_steel(varargin{:});
                    otherwise
                        error('Unknown type: %s',type);
                end
            else
                error('Bad input');
            end
            assert(isempty(find(obj.materialIDs==mat.id,1)),'material id is not unique');
            obj.materials = vertcat(obj.materials,{mat});
            obj.materialIDs = vertcat(obj.materialIDs,mat.id);
        end
        function set.fiberSection(obj,fiberSection)
            assert(isa(fiberSection,'fiberSection'),...
                'fiberSection should be a fiberSection object');
            obj.fiberSection = fiberSection;
        end
        
        %% Computation
        function y = extremeConcreteCompressionFiber(obj,zPoint,yPoint,angle)
            assert(size(obj.concreteBoundary_radius,1) > 0 ,...
                'No concrete boundaries defined');
            
            a = sin(angle);
            b = -cos(angle);
            c = -sin(angle)*zPoint + cos(angle)*yPoint;
            
            y = a*obj.concreteBoundary_zPoint + b*obj.concreteBoundary_yPoint + c;
            y = min(y-obj.concreteBoundary_radius);
        end
        function y = extremeSteelTensionFiber(obj,zPoint,yPoint,angle)
            assert(size(obj.steelBoundary_radius,1) > 0 ,...
                'No steel boundaries defined');
            
            a = sin(angle);
            b = -cos(angle);
            c = -sin(angle)*zPoint + cos(angle)*yPoint;
            
            y = a*obj.steelBoundary_zPoint + b*obj.steelBoundary_yPoint + c;
            y = max(y+obj.steelBoundary_radius);
        end
        function et = extremeSteelTensileStrain(obj,zPoint,yPoint,angle)
            yc = obj.extremeConcreteCompressionFiber(zPoint,yPoint,angle);
            yt = obj.extremeSteelTensionFiber(zPoint,yPoint,angle);
            ec = obj.strainAtExtremeConcreteFiber;
            if yc < 0
                et = yt*(ec/yc);
            else
                et = obj.defaultTensileStrain;
            end
        end 
        function [d,maxDist,minDist] = sectionDepth(obj,angle)
            % Find points
            [~,~,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            dist = [-sin(angle) cos(angle)]*[z y]';
            maxDist = max(dist);
            minDist = min(dist);  
            d = maxDist-minDist;
            if nargout < 3
                clear('maxDist','minDist');
            end
        end
        function [P,Mz,My] = computePoint(obj,zPoint,yPoint,angle)
            mats = obj.fiberSection.matIDs;
            [mat,A,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            
            a = sin(angle);
            b = -cos(angle);
            c = -sin(angle)*zPoint + cos(angle)*yPoint;
            
            yECF = obj.extremeConcreteCompressionFiber(zPoint,yPoint,angle);
            if yECF < 0
                strain = obj.strainAtExtremeConcreteFiber/yECF*(a*z+b*y+c);
            else
                strain = obj.defaultTensileStrain*ones(size(mat));
            end
            
            stress = nan(size(mat));
            for i = 1:length(mats)
                % Find fibers of the material
                ind = find(mat==mats(i));
                
                % Find the constitutive relation
                iMat = find(obj.materialIDs==mats(i));
                assert(isscalar(iMat),'cannot find constitutive relation for material %i',mats(i))
                
                % Compute stress
                stress(ind) = obj.materials{iMat}.getStress(strain(ind));
            end
            P  = sum(stress.*A);
            Mz = sum(stress.*A.*-y);
            My = sum(stress.*A.*z);
            
            if ~isempty(obj.maxCompressiveStrength)
                P = max(P,obj.maxCompressiveStrength);
            end
        end      
        function [P,Mz,My] = computePoint_uniform(obj,axial_strain)
            mats = obj.fiberSection.matIDs;
            [mat,A,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            
            strain = axial_strain*ones(size(mat));            
            stress = nan(size(mat));
            
            for i = 1:length(mats)
                % Find fibers of the material
                ind = find(mat==mats(i));
                
                % Find the constitutive relation
                iMat = find(obj.materialIDs==mats(i));
                assert(isscalar(iMat),'cannot find constitutive relation for material %i',mats(i))
                
                % Compute stress
                stress(ind) = obj.materials{iMat}.getStress(strain(ind));
            end
            P  = sum(stress.*A);
            Mz = sum(stress.*A.*-y);
            My = sum(stress.*A.*z);
            
            if ~isempty(obj.maxCompressiveStrength)
                P = max(P,obj.maxCompressiveStrength);
            end
        end        
        function [P,Mz,My,et] = interactionSweep(obj,angle,numPoints)           
            [ymax,ymin] = obj.fiberSection.boundsAtAngle(angle,obj.AxesOrigin);
            d = ymax-ymin;
            points = [-Inf linspace(ymin-0.51*d,ymax+0.51*d,numPoints-2) Inf];
            
            % Compute interaction
            P  = zeros(2*numPoints,1);
            Mz = zeros(2*numPoints,1);
            My = zeros(2*numPoints,1);
            et = zeros(2*numPoints,1);
            for i = 1:numPoints
                if points(i) == Inf
                    [iP,iMz,iMy] = computePoint_uniform(obj,obj.defaultTensileStrain);
                    iet = obj.defaultTensileStrain;
                elseif points(i) == -Inf
                    [iP,iMz,iMy] = computePoint_uniform(obj,obj.strainAtExtremeConcreteFiber);
                    iet = obj.strainAtExtremeConcreteFiber;
                else
                    zPoint = -sin(angle)*points(i);
                    yPoint =  cos(angle)*points(i);
                    [iP,iMz,iMy] = computePoint(obj,zPoint,yPoint,angle);
                    iet = obj.extremeSteelTensileStrain(zPoint,yPoint,angle);
                end
                
                P(i)  = iP;
                Mz(i) = iMz;
                My(i) = iMy;
                et(i) = iet; 
                
                if points(i) == Inf
                    [iP,iMz,iMy] = computePoint_uniform(obj,obj.strainAtExtremeConcreteFiber);
                    iet = obj.strainAtExtremeConcreteFiber;
                elseif points(i) == -Inf
                    [iP,iMz,iMy] = computePoint_uniform(obj,obj.defaultTensileStrain);
                    iet = obj.defaultTensileStrain;
                else
                    zPoint = -sin(angle)*points(i);
                    yPoint =  cos(angle)*points(i);                
                    [iP,iMz,iMy] = computePoint(obj,zPoint,yPoint,angle+pi);
                    iet = obj.extremeSteelTensileStrain(zPoint,yPoint,angle+pi);
                end                

                P(numPoints+i)  = iP;
                Mz(numPoints+i) = iMz;
                My(numPoints+i) = iMy;
                et(numPoints+i) = iet;
            end
            
            if ~isempty(obj.maxCompressiveStrength)
                % Apply maximum strength cap
                P(P<obj.maxCompressiveStrength) = obj.maxCompressiveStrength;
                % Close the loop
                P  = vertcat(P,P(1));
                Mz = vertcat(Mz,Mz(1));
                My = vertcat(My,My(1));
                et = vertcat(et,et(1));
            end
            if nargout < 4
                clear('et');
            end
        end
        function [P,Mz,My] = interaction3d(obj,numPoints,numAngles)            
            if isempty(obj.maxCompressiveStrength)
                P  = zeros(2*numPoints,numAngles);
                Mz = zeros(2*numPoints,numAngles);
                My = zeros(2*numPoints,numAngles);
            else
                P  = zeros(2*numPoints+1,numAngles);
                Mz = zeros(2*numPoints+1,numAngles);
                My = zeros(2*numPoints+1,numAngles);
            end
            angles = linspace(0,pi,numAngles);
            
            for i = 1:numAngles
                [iP,iMz,iMy] = obj.interactionSweep(angles(i),numPoints);
                P(:,i) = iP;
                Mz(:,i) = iMz;
                My(:,i) = iMy;
            end
        end
        function results = interaction3d_2(obj,numP,numAngles,numPoints_calc,numAngles_calc)
            
            if nargin < 4
                numPoints_calc = 50;
            end
            if nargin < 5
                numAngles_calc = 2*numAngles;
            end
            
            angles_calc = linspace(0,2*pi,numAngles_calc+1);
            
            Mz_calc = nan(numP,numAngles_calc);
            My_calc = nan(numP,numAngles_calc);             
                        
            [~,~,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            for i = 1:numAngles_calc
                dist = [-sin(angles_calc(i)) cos(angles_calc(i))]*[z y]';
                maxDist = max(dist);
                minDist = min(dist);
                d = maxDist-minDist;
                points = [-Inf linspace(minDist-0.51*d,maxDist+0.51*d,numPoints_calc-2) Inf];
                
                iP  = nan(numPoints_calc,1);
                iMz = nan(numPoints_calc,1);
                iMy = nan(numPoints_calc,1);
                for j = 1:numPoints_calc
                    if points(j) == Inf
                        [jP,jMz,jMy] = computePoint_uniform(obj,obj.defaultTensileStrain);
                    elseif points(j) == -Inf
                        [jP,jMz,jMy] = computePoint_uniform(obj,obj.strainAtExtremeConcreteFiber);
                    else
                        zPoint = -sin(angles_calc(i))*points(j);
                        yPoint =  cos(angles_calc(i))*points(j);
                        [jP,jMz,jMy] = computePoint(obj,zPoint,yPoint,angles_calc(i));
                    end
                    iP(j)  = jP;
                    iMz(j) = jMz;
                    iMy(j) = jMy;
                end
                
                if i == 1
                    P_max  = iP(1);
                    Mz_max = iMz(1);
                    My_max = iMy(1);
                    P_min  = iP(end);
                    Mz_min = iMz(end);
                    My_min = iMy(end);
                    if isempty(obj.maxCompressiveStrength) || obj.maxCompressiveStrength >= P_max
                        P = linspace(P_min,P_max,numP);
                    else
                        P = [linspace(P_min,0.99*P_max,numP-1) P_max];
                    end
                end
                
                for j = 1:numP
                    [ind,x] = find_limit_point_in_vector(iP,P(j));
                    Mz_calc(j,i) = interpolate_vector(iMz,ind,x);
                    My_calc(j,i) = interpolate_vector(iMy,ind,x);
                end
            end

            results = struct;
            results.P  = nan(numP,numAngles);
            results.Mz = nan(numP,numAngles);
            results.My = nan(numP,numAngles);
            angles = linspace(0,2*pi,numAngles+1);
            angles = angles(1:(end-1));
            for i = 1:numP
                if i == 1
                    results.P(i,:)  = P_min;
                    results.Mz(i,:) = Mz_min;
                    results.My(i,:) = My_min;
                elseif i == numP
                    results.P(i,:)  = P_max;
                    results.Mz(i,:) = Mz_max;
                    results.My(i,:) = My_max;
                else
                    id = interactionDiagram2d(Mz_calc(i,:),My_calc(i,:));
                    d = id.radial_distance(angles);
                    results.P(i,:)  = P(i);
                    results.Mz(i,:) = d.*cos(angles);
                    results.My(i,:) = d.*sin(angles);
                end
            end
        end        
    end
end

