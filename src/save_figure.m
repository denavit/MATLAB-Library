function save_figure(h,filename,type,resolution)
    switch lower(type)
        case 'emf'
            %originalPosition = get(h,'Position');
            %paperSize = get(h,'PaperSize');
            %set(h,'Position',[1 1 paperSize(1)*72 paperSize(2)*72]);
            print(h,'-dmeta',sprintf('-r%i',resolution),filename);
            %set(h,'Position',originalPosition);
        case 'fig'
            saveas(h,filename);
        case 'tiff'
            print(h,'-dtiff',sprintf('-r%i',resolution),filename);
        case 'eps'
            print(h,'-depsc','-tiff',filename);
        case 'svg'
            print(h,'-dsvg',filename);                    
        case 'jpg'
            print(h,'-djpeg',sprintf('-r%i',resolution),filename);
        case 'png'
            print(h,'-dpng',sprintf('-r%i',resolution),filename);
        otherwise
            error('Unknown print type');
    end
end