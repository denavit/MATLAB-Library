classdef RC < structural_shape
    
    properties
        fc
        Ec
        
        fy
        Es
        
        conc_cross_section
        reinforcement
        
        reinforcement_is_spiral = false;
        treat_reinforcement_as_point = true;
    end
    
    methods
        function obj = RC(fc,fy,conc_cross_section,reinforcement,units)
            obj.fc = fc;
            obj.fy = fy;
            obj.conc_cross_section = conc_cross_section;
            if iscell(reinforcement)
                obj.reinforcement = reinforcement;
            else
                obj.reinforcement = {reinforcement};
            end
            obj.units = units;
            
            switch obj.units
                case 'US'
                    obj.Es = 29000;
                    obj.Ec = 57*sqrt(obj.fc*1000);
                case 'SI'
                    obj.Es = unitConvert('pressure',29000,'ksi','MPa');
                    fc_psi = unitConvert('pressure',obj.fc,'MPa','psi');
                    Ec_ksi = 57*sqrt(fc_psi);
                    obj.Ec = unitConvert('pressure',Ec_ksi,'ksi','MPa');
                otherwise
                    warning('Design:Input','Unknown unit system, setting Es = 0 and Ec = 0')
                    obj.Es = 0;
                    obj.Ec = 0;
            end
        end
        function d = depth(obj,axis)
            d = obj.conc_cross_section.depth(axis);
        end
        function a = Ag(obj)
            a = obj.conc_cross_section.A;
        end
        function a = Ac(obj)
            a = obj.Ag-obj.Asr;
        end
        function a = Asr(obj)
            a = 0;
            for i = 1:length(obj.reinforcement)
                a = a + obj.reinforcement{i}.nb*obj.reinforcement{i}.Ab;
            end
        end
        function i = Ig(obj,axis)
            i = obj.conc_cross_section.I(axis);
        end
        function i = Ic(obj,axis)
            i = obj.Ig(axis)-obj.Isr(axis);
        end
        function i = Isr(obj,axis)
            i = 0;
            for j = 1:length(obj.reinforcement)
                i = i + obj.reinforcement{j}.I(axis);
            end
        end
        function pnco = Pnco(obj)
            error('Not yet implemented')
        end
        function pnc = Pnc(obj,axis)
            error('Not yet implemented')
        end
        function pnt = Pnt(obj)
            error('Not yet implemented')
        end
        function mno = Mno(obj,axis)
            error('Not yet implemented')
        end
        function [P,M] = sectionInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case 'aci'
                    scACI = obj.strainCompatibilityAciObject;
                    switch lower(axis)
                        case {'x','z'}
                            [P,M,~] = scACI.interactionSweep(0,50);
                        case 'y'
                            [P,~,M] = scACI.interactionSweep(pi/2,50);
                        otherwise
                            error('Unknown axis: %s',axis);
                    end
                case 'factoredaci'
                    scACI = obj.strainCompatibilityAciObject;
                    switch lower(axis)
                        case {'x','z'}
                            [P,M,~,et] = scACI.interactionSweep(0,50);
                        case 'y'
                            [P,~,M,et] = scACI.interactionSweep(pi/2,50);
                        otherwise
                            error('Unknown axis: %s',axis);
                    end
                    if obj.reinforcement_is_spiral
                        phi = ACI_phi('spiral',et,obj.Fy/obj.Es);
                    else
                        phi = ACI_phi('other',et,obj.Fy/obj.Es);
                    end
                    P = phi.*P;
                    M = phi.*M;
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [P,M] = beamColumnInteraction2d(obj,axis,type,quadrant)
            switch lower(strtok(type,'-'))
                case {'aci','factoredaci'}
                    [P,M] = sectionInteraction2d(obj,axis,type,quadrant);
                otherwise
                    error('Unknown type: %s',type);
            end
        end
        function [E,A,I] = sectionPropertiesForElasticAnalysis2d(obj,axis,type)
            error('Not yet implemented')
        end
        function [E,A,Iz,Iy,GJ] = sectionPropertiesForElasticAnalysis3d(obj,type)
            error('Not yet implemented')
        end
        function strain = longitudinalStrain2d(obj,axis,axialStrain,curvature,type)
            error('Not yet implemented')
        end
        function x = getSectionData(obj,type,axis)
            error('Not yet implemented')
        end
        function ratio = interactionCheck(obj,xi,P,Ms,Mw,Vs,Vw,T)
            error('Not yet implemented')
        end
        function pass_tf = proportioningCheck(obj,checkType,varagrin)
            error('Not yet implemented')
        end
        function description = sectionDescription(obj,flag)
            error('Not yet implemented')
        end
        function fs = fiberSectionObject(obj,idConc,idReinf)
            fs = fiberSection;
            obj.conc_cross_section.add_to_fiber_section(fs,idConc);
            for i = 1:length(obj.reinforcement)
                Ab    = obj.reinforcement{i}.Ab;
                [z,y] = obj.reinforcement{i}.coordinates;
                for j = 1:length(z)
                    if obj.treat_reinforcement_as_point
                        fs.addFiber( idConc,-Ab,z(j),y(j))
                        fs.addFiber(idReinf, Ab,z(j),y(j))
                    else
                        error('Not yet implemented');
                    end
                end
            end
        end        
        function psd = plasticStressDistributionObject(obj)
            error('Not yet implemented')
        end
        function scACI = strainCompatibilityAciObject(obj)
            idConc  = 1;
            idReinf = 2;
            fs = obj.fiberSectionObject(idConc,idReinf);
            scACI = ACI_strain_compatibility(fs);
            
            % Add Concrete Boundaries
            [x,y,r] = obj.conc_cross_section.boundary_points();
            for i = 1:length(x)
                scACI.addConcreteBoundary(x(i),y(i),r(i));
            end
            
            % Add Steel Boundaries
            for i = 1:length(obj.reinforcement)
                [z,y] = obj.reinforcement{i}.coordinates;
                for j = 1:length(z)
                    if obj.treat_reinforcement_as_point
                        scACI.addSteelBoundary(z(j),y(j),0); 
                    else 
                        error('Not yet implemented');
                    end
                end
            end
            
            scACI.maxCompressiveStrength = -0.85*(obj.Asr*obj.fy + 0.85*obj.Ac*obj.fc);
            
            scACI.addMaterial('concrete', idConc,obj.fc,obj.units);
            scACI.addMaterial(   'steel',idReinf,obj.fy,obj.Es);
        end           
        
        function plotSection(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            hold all
            obj.conc_cross_section.plotSection(lineWidth);
            for i = 1:length(obj.reinforcement)
                obj.reinforcement{i}.plotSection(lineWidth);
            end
            axis equal
        end
    end    
    methods (Static)
        function type = memberType()
            type = 'RC';
        end
        function tf = hasConcrete()
            tf = true;
        end
    end    
    
    
end
