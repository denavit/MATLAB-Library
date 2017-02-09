classdef elastic_section_stiffness < handle
    
    properties
        AxesOrigin = 'AsDefined';
    end
    
    properties (SetAccess = private)
        fiberSection 
        materials_id = zeros(0,1);
        materials_E  = zeros(0,1);
        materials_t  = zeros(0,1);
        materials_c  = zeros(0,1);
    end

    methods
        %% Constructor
        function obj = elastic_section_stiffness(fiberSection,AxesOrigin)
            if nargin > 0
                obj.fiberSection = fiberSection;
            end
            if nargin > 1
                obj.AxesOrigin = AxesOrigin;
            end
        end
        
        %% Define section
        function addMaterial(obj,id,E,type)
            assert(sum(obj.materials_id==id)==0,'Material Already Defined');
            obj.materials_id = vertcat(obj.materials_id,id);
            obj.materials_E  = vertcat(obj.materials_E,E);
            if nargin < 4
                obj.materials_t = vertcat(obj.materials_t,true);
                obj.materials_c = vertcat(obj.materials_c,true);
            else
                switch lower(type)
                    case 'tensiononly'
                        obj.materials_t = vertcat(obj.materials_t,true);
                        obj.materials_c = vertcat(obj.materials_c,false);
                    case 'compressiononly'
                        obj.materials_t = vertcat(obj.materials_t,false);
                        obj.materials_c = vertcat(obj.materials_c,true);
                    otherwise
                        error('Unknown type: %s',type);
                end
            end
        end
        function set.fiberSection(obj,fiberSection)
            assert(isa(fiberSection,'fiberSection'),...
                'fiberSection should be a fiberSection object');
            obj.fiberSection = fiberSection;
        end
        
        %% Computation
        function [EIz,EIy,ez,ey] = momentOfInertia(obj,zPoint,yPoint,angle)
            [mat,A,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            
            a = sin(angle);
            b = -cos(angle);
            c = -sin(angle)*zPoint + cos(angle)*yPoint;

            strain  = a*z+b*y+c;
            [stress,tangent] = obj.computeStressAndTangent(mat,strain);
            
            P   = sum(stress.*A);
            Mz  = sum(stress.*A.*-y);
            My  = sum(stress.*A.*z);
            EIz = sum(tangent.*A.*y.^2);
            EIy = sum(tangent.*A.*z.^2);
            
            if nargout > 2
                ez = Mz/P;
                ey = My/P;
            end
        end      
        function [EIz,EIy,ez,ey] = momentOfInertiaSweep(obj,angle,numPoints)           
            % Find points
            [~,~,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            dist = [-sin(angle) cos(angle)]*[z y]';
            maxDist = max(dist);
            minDist = min(dist);
            points = linspace(minDist,maxDist,numPoints);
            points(1) = points(1)-1.0e-6*(maxDist-minDist);
            points(end) = points(end)+1.0e-6*(maxDist-minDist);

            % Compute interaction
            EIz = nan(2*numPoints,1);
            EIy = nan(2*numPoints,1);
            ez  = nan(2*numPoints,1);
            ey  = nan(2*numPoints,1);
            for i = 1:numPoints
                zPoint = -sin(angle)*points(i);
                yPoint =  cos(angle)*points(i);
                                
                [iEIz,iEIy,iez,iey] = obj.momentOfInertia(zPoint,yPoint,angle);
                EIz(i) = iEIz;
                EIy(i) = iEIy;
                ez(i)  = iez;
                ey(i)  = iey;
                
                [iEIz,iEIy,iez,iey] = obj.momentOfInertia(zPoint,yPoint,angle+pi);
                EIz(numPoints+i) = iEIz;
                EIy(numPoints+i) = iEIy;
                ez(numPoints+i)  = iez;
                ey(numPoints+i)  = iey;
            end
        end
        
        function [EIz,EIy] = momentOfInertiaAtZeroP(obj,angle)
            [mat,A,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            
            % Find strain distribution with zero axial load
            options = struct;
            options.Display = 'off';
            [c,~,exitflag] =  fsolve(...
                @(c)errorAxialLoad(obj,angle,c),0,options);
            if exitflag <= 0
                error('fsolve could not find solution');
            end
            
            % Compute stiffness
            a = sin(angle);
            b = -cos(angle);
            strain  = a*z+b*y+c;
            [~,tangent] = obj.computeStressAndTangent(mat,strain);
            EIz = sum(tangent.*A.*y.^2);
            EIy = sum(tangent.*A.*z.^2);
        end
        
        function [stress,tangent] = computeStressAndTangent(obj,mat,strain)
            mats = obj.fiberSection.matIDs;

            stress  = nan(size(strain));
            tangent = nan(size(strain));
            
            for i = 1:length(mats)
                % Find fibers of the material
                ind = find(mat==mats(i));
                
                % Find the constitutive relation
                iMat = find(obj.materials_id==mats(i));
                assert(isscalar(iMat),'cannot find constitutive relation for material %i',mats(i))
                
                % Compute stress
                iStrain  = strain(ind);
                iStress  = obj.materials_E(iMat)*iStrain;
                iTangent = repmat(obj.materials_E(iMat),size(iStrain));  
                if obj.materials_t(iMat) == false
                   % Set the stress and tangent to zero for fibers in tension
                   iStress(iStrain>0)  = 0;
                   iTangent(iStrain>0) = 0;                   
                end
                if obj.materials_c(iMat) == false
                   % Set the stress and tangent to zero for fibers in compression
                   iStress(iStrain<0)  = 0;
                   iTangent(iStrain<0) = 0;   
                end
                stress(ind)  = iStress;
                tangent(ind) = iTangent;                
            end            
        end
        
    end
end


function err = errorAxialLoad(obj,angle,c)
[mat,A,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
a = sin(angle);
b = -cos(angle);
strain  = a*z+b*y+c;
[stress,~] = obj.computeStressAndTangent(mat,strain);
err = sum(stress.*A);
end
