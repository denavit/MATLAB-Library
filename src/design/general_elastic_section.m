classdef general_elastic_section < structural_shape
    % general_elastic_section
    %
    % Notes:
    %
    % REFERENCES:
    %
    
    properties
        A           % Cross sectional area
        Iz          % Moment of Inertia about the Z axis (strong axis)
        Iy          % Moment of Inertia about the Y axis (weak axis)
        E           % Modulus of Elasticity 
        GJ = 1e10   % Torsional Stiffness
    end
    
    methods
        %% Constructor
        function obj=general_elastic_section(E,A,Iz,Iy,GJ)
            % general_elastic_section(E,A,I)
            % general_elastic_section(E,A,Iz,Iy)
            % general_elastic_section(E,A,Iz,Iy,GJ)
            
            obj.E = E;
            obj.A = A;
            obj.Iz = Iz;
            if (nargin == 3)
                obj.Iy = Iz;
            end
            if (nargin >= 4)
                obj.Iy = Iy;
            end
            if (nargin >= 5)
                obj.GJ = GJ;
            end
        end
        
        function d = depth(obj,axis)
            d = [];            
        end
        
        %% Design Strengths
        function pnco = Pnco(obj)
            % Stub Column (L=0) Strength
            pnco = Inf;            
        end
        function pnc = Pnc(obj,axis)
            % Column Strength
            pnc = Inf; % @todo - euler load
        end
        function pnt = Pnt(obj)
            % Tensile Strength
            pnt = Inf;
        end
        function mno = Mno(obj,axis)
            % L=0 Moment Strength
            mno = Inf;
        end
        
        %% Design Checks
        function ratio = interactionCheck(obj,xi,P,Ms,Mw,Vs,Vw,T)           
            % Elastic Section, No Strength Limits
            ratio = 0;
        end
        function pass_tf = proportioningCheck(obj,checkType,varargin)
            switch lower(checkType)
                case 'slendernessratio_klr'
                    limit = varargin{1};
                    KLr_strong = obj.Kstrong*obj.L/sqrt(obj.Iz/obj.A);
                    KLr_weak = obj.Kweak*obj.L/sqrt(obj.Iy/obj.A);
                    pass_tf = max([KLr_strong KLr_weak]) < limit;
                otherwise
                    error('Unknown checkType');
            end
        end
        
        %% Interaction Strength
        function [P,M] = sectionInteraction2d(obj,axis,type,quadrant)
            error('Elastic Section, sectionInteraction2d not implemented');
        end
        function [P,M] = beamColumnInteraction2d(obj,axis,type,quadrant)
            error('Elastic Section, beamColumnInteraction2d not implemented');
        end
        
        %% Export and Information Functions
        function [E,A,I] = sectionPropertiesForElasticAnalysis2d(...
                obj,axis,type)
            switch lower(type)
                case {'gross','columnstrength'}
                    E = obj.E;
                    A = obj.A;
                    switch lower(axis)
                        case 'strong'
                            I = obj.Iz;
                        case 'weak'
                            I = obj.Iy;
                        otherwise
                            error('Bad axis');
                    end
                otherwise
                    error('Unknown type');
            end
        end
        function [E,A,Iz,Iy,GJ] = sectionPropertiesForElasticAnalysis3d(...
                obj,type)
            switch lower(type)
                case {'gross','columnstrength'}
                    E = obj.E;
                    A = obj.A;
                    Iz = obj.Iz;
                    Iy = obj.Iy;
                    GJ = obj.GJ;
                otherwise
                    error('Unknown type');
            end
        end
        function strain = longitudinalStrain2d(obj,axis,axialStrain,curvature,type)
            strain = nan(size(axialStrain));
            %error('Elastic Section, longitudinalStrain2d not implemented');
        end 
        function x = getSectionData(obj,type,axis)
            switch lower(type)
                otherwise
                    x = NaN;
            end
        end
        function description = sectionDescription(obj,flag)
            if nargin < 2
                flag = 1;
            end
            switch flag
                case 1
                    description = 'Elastic';
                case 2
                    description = sprintf('Elastic E = %g, A = %g, Iz = %g, Iy = %g',...
                        obj.E,obj.A,obj.Iz,obj.Iy);
                otherwise
                    error('Unknown flag');
            end
        end
        function fs = fiberSectionObject(obj)
            error('Elastic Section, No Shape for Fiber Section Object');
        end
        function psd = plasticStressDistributionObject(obj)
            error('Elastic Section, No Shape for Plastic Stress Distribution Object');
        end
        function plotSection(obj)
            warning('Elastic Section, No Shape for Plot');
        end        
    end
    
    methods (Static)
        function type = memberType()
            type = 'elastic';
        end   
        function tf = hasConcrete()
            tf = false;
        end 
    end
    
end
