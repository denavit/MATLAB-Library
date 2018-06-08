classdef ground_motion
    
    properties
        acc
        dt
    end
    
    methods
        function obj = ground_motion(acc,dt)
            % Constructor
            if isrow(acc)
                acc = acc';
            elseif ~iscolumn(acc)
                error('acc should be a numeric vector');
            end
            obj.acc = acc;
            obj.dt = dt; 
        end
        function npts = numPoints(obj)
            npts = length(obj.acc);
        end
        function t = time(obj)
            t = linspace(0,obj.numPoints*obj.dt,obj.numPoints)';
        end
        function pga = PGA(obj)
            pga = max(abs(obj.acc));
        end
        function pgv = PGV(obj)
            pgv = max(abs(cumtrapz(obj.time,obj.acc)));
        end
        function pgd = PGD(obj)
            pgd = max(abs(cumtrapz(obj.time,cumtrapz(obj.time,obj.acc))));
        end
        function [D,V,A] = elasticAnalysisSDOF(obj,Tn,zeta,D0,V0,gamma,beta)
            if nargin < 7
                beta = 1/6;
            end
            if nargin < 6
                gamma = 1/2;
            end
            if nargin < 5
                V0 = 0.0;
            end
            if nargin < 4
                D0 = 0.0;
            end
            
            D = zeros(obj.numPoints,1);
            V = zeros(obj.numPoints,1);
            A = zeros(obj.numPoints,1);
                               
            % Define SDOF system
            wn = 2*pi/Tn;
            m = 1;
            k = wn^2;
            c = 2*zeta*wn;
            f = -obj.acc;
            
            % Initial Step
            D(1) = V0;
            V(1) = D0;
            A(1) = (f(1) - c*V(1) - k*D(1))/m;
            
            % Remaining Steps
            khat = k + gamma/(beta*obj.dt)*c + 1/(beta*obj.dt^2)*m;
            a = m/(beta*obj.dt) + gamma*c/beta;
            b = 0.5*m/beta + obj.dt*(0.5*gamma/beta - 1)*c;
            for i=1:(obj.numPoints-1)
                dfhat = f(i+1)-f(i) + a*V(i) + b*A(i);
                dD = dfhat/khat;
                dV = gamma/(beta*obj.dt)*dD - gamma/beta*V(i) + ...
                    obj.dt*(1-gamma/(2*beta))*A(i);
                dA = 1/(beta*obj.dt^2)*dD - 1/(beta*obj.dt)*V(i) - ...
                    1/(2*beta)*A(i);
                D(i+1) = D(i) + dD;
                V(i+1) = V(i) + dV;
                A(i+1) = A(i) + dA;
            end
            A = A-f;
        end
        function [Sd,Sv,Sa] = elasticResponseSpectra(obj,Tn,zeta)
            numAnalyses = numel(Tn);
            Sd = zeros(size(Tn));
            for i = 1:numAnalyses
                [d,~,~] = obj.elasticAnalysisSDOF(Tn(i),zeta);
                Sd(i) = max(abs(d));
            end
            wn = 2*pi./Tn;
            Sv = wn.*Sd;
            Sa = wn.^2.*Sd;
        end
    end
end
