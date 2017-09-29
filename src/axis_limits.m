function axis_limits(type,varargin)

switch lower(type)
    case 'quadrant'
        quad = varargin{1};
        xl = xlim;
        yl = ylim;
        switch quad
            case 1
                xlim([0 xl(2)])
                ylim([0 yl(2)])
            case 2
                xlim([xl(1) 0])
                ylim([0 yl(2)])
            case 3
                xlim([xl(1) 0])
                ylim([yl(1) 0])
            case 4
                xlim([0 xl(2)])
                ylim([yl(1) 0])
            otherwise
                error('Unknown quadrant');
        end
    case 'marginx'
        switch length(varargin{1})
            case 1
                rm = varargin{1};
                lm = varargin{1};
            case 2
                rm = varargin{1}(1);
                lm = varargin{1}(2);
            otherwise
                error('Bad number of arguments')
        end
        [xmax,xmin,~,~] = extents();
        xlim([xmin-lm*(xmax-xmin) xmax+rm*(xmax-xmin)]);
    case 'marginy'
        switch length(varargin{1})
            case 1
                tm = varargin{1};
                bm = varargin{1};
            case 2
                tm = varargin{1}(1);
                bm = varargin{1}(2);
            otherwise
                error('Bad number of arguments')
        end
        [~,~,ymax,ymin] = extents();
        ylim([ymin-bm*(ymax-ymin) ymax+tm*(ymax-ymin)]);
    case 'margin'
        [tm,rm,bm,lm] = trbl(varargin{1});
        [xmax,xmin,ymax,ymin] = extents();
        xlim([xmin-lm*(xmax-xmin) xmax+rm*(xmax-xmin)]);
        ylim([ymin-bm*(ymax-ymin) ymax+tm*(ymax-ymin)]);
    case 'balancex'
        x = max(abs(xlim));
        xlim([-x x]);
    case 'balancey'
        y = max(abs(ylim));
        ylim([-y y]);
    case 'balance'
        figureTools.lim('balancex');
        figureTools.lim('balancey');
    otherwise
        error('Unknown type: %s',type);
end
end


function [t,r,b,l] = trbl(x)
switch length(x)
    case 1
        t = x;
        r = x;
        b = x;
        l = x;
    case 2
        t = x(1);
        r = x(2);
        b = x(1);
        l = x(2);
    case 3
        t = x(1);
        r = x(2);
        b = x(3);
        l = x(2);
    case 4
        t = x(1);
        r = x(2);
        b = x(3);
        l = x(4);
    otherwise
        error('Bad size')
end
end


function [xmax,xmin,ymax,ymin] = extents(ha)
if nargin < 1
    ha = gca;
end
xmax = -Inf;
xmin =  Inf;
ymax = -Inf;
ymin =  Inf;
h = get(ha,'Children');
if isempty(h)
    error('Nothing in the plot')
end
for i = 1:length(h)
    switch get(h(i),'Type')
        case {'line','hggroup','scatter'}
            xdata = get(h(i),'XData');
            ydata = get(h(i),'YData');
            xmax = max(xmax,max(xdata));
            xmin = min(xmin,min(xdata));
            ymax = max(ymax,max(ydata));
            ymin = min(ymin,min(ydata));
        case 'patch'
            xmax = max(xmax,max(h(i).Vertices(:,1)));
            xmin = min(xmin,min(h(i).Vertices(:,1)));
            ymax = max(ymax,max(h(i).Vertices(:,2)));
            ymin = min(ymin,min(h(i).Vertices(:,2)));
        otherwise
            warning('axis_limits:extents','Ignoring object of type: %s\n',get(h(i),'Type'));
    end
end
end