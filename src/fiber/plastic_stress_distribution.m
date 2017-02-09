classdef plastic_stress_distribution < handle
    
    properties
        AxesOrigin = 'AsDefined';
        fiberSection 
        materials_id = zeros(0,1);
        materials_Ft = zeros(0,1);
        materials_Fc = zeros(0,1);
    end

    methods
        %% Constructor
        function obj = plastic_stress_distribution(fiberSection,AxesOrigin)
            if nargin > 0
                obj.fiberSection = fiberSection;
            end
            if nargin > 1
                obj.AxesOrigin = AxesOrigin;
            end
        end
        
        %% Define section
        function addMaterial(obj,id,Ft,Fc)
            assert(sum(obj.materials_id==id)==0,'Material Already Defined');
            obj.materials_id = vertcat(obj.materials_id,id);
            obj.materials_Ft = vertcat(obj.materials_Ft,Ft);
            obj.materials_Fc = vertcat(obj.materials_Fc,Fc);
        end
        function set.fiberSection(obj,fiberSection)
            assert(isa(fiberSection,'fiberSection'),...
                'fiberSection should be a fiberSection object');
            obj.fiberSection = fiberSection;
        end
        
        %% Computation
        function [P,Mz,My] = computePoint(obj,zPoint,yPoint,angle)
            mats = obj.fiberSection.matIDs;
            [mat,A,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            
            a = sin(angle);
            b = -cos(angle);
            c = -sin(angle)*zPoint + cos(angle)*yPoint;
            P = 0;
            Mz = 0;
            My = 0;
            for i = 1:length(mats)
                ind = find(obj.materials_id==mats(i));
                assert(~isempty(ind),'Material %i is undefined',mats(i));
                Ft = obj.materials_Ft(ind);
                Fc = obj.materials_Fc(ind);
                tensFibers = find(mat(:,1) == mats(i) & a*z+b*y+c < 0);
                compFibers = find(mat(:,1) == mats(i) & a*z+b*y+c > 0);
                P  = P  + sum(Ft*A(tensFibers)) + sum(Fc*A(compFibers));
                Mz = Mz - sum(Ft*A(tensFibers).*y(tensFibers)) - ...
                    sum(Fc*A(compFibers).*y(compFibers));
                My = My + sum(Ft*A(tensFibers).*z(tensFibers)) + ...
                    sum(Fc*A(compFibers).*z(compFibers));
            end
        end      
        function [P,Mz,My] = interactionSweep(obj,angle,numPoints)           
            % Find points
            [~,~,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            dist = [-sin(angle) cos(angle)]*[z y]';
            maxDist = max(dist);
            minDist = min(dist);
            points = linspace(minDist,maxDist,numPoints);
            points(1) = points(1)-1.0e-6*(maxDist-minDist);
            points(end) = points(end)+1.0e-6*(maxDist-minDist);

            % Compute interaction
            P  = zeros(2*numPoints,1);
            Mz = zeros(2*numPoints,1);
            My = zeros(2*numPoints,1);
            for i = 1:numPoints
                zPoint = -sin(angle)*points(i);
                yPoint =  cos(angle)*points(i);
                                
                [iP,iMz,iMy] = computePoint(obj,zPoint,yPoint,angle);
                P(i)  = iP;
                Mz(i) = iMz;
                My(i) = iMy;
                
                [iP,iMz,iMy] = computePoint(obj,zPoint,yPoint,angle+pi);
                P(numPoints+i)  = iP;
                Mz(numPoints+i) = iMz;
                My(numPoints+i) = iMy;
            end
        end
        function [P,Mz,My] = interaction3d(obj,numPoints,numAngles)            
            P  = zeros(2*numPoints,numAngles);
            Mz = zeros(2*numPoints,numAngles);
            My = zeros(2*numPoints,numAngles);
            angles = linspace(0,pi,numAngles);
            
            for i = 1:numAngles
                [iP,iMz,iMy] = obj.interactionSweep(angles(i),numPoints);
                P(:,i) = iP;
                Mz(:,i) = iMz;
                My(:,i) = iMy;
            end
        end        
    end
end

