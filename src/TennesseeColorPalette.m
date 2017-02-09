function rgb = TennesseeColorPalette(name)
switch lower(name)
    case 'orange'
        rgb = [255 130 0]/255;
    case 'white'
        rgb = [255 255 255]/255;
    case 'smokey'
        rgb = [88 89 91]/255;
    case 'valley'
        rgb = [0 116 111]/255;
    case 'torch'
        rgb = [230 89 51]/255;
    case 'globe'
        rgb = [0 108 147]/255;
    case 'limestone'
        rgb = [240 237 227]/255;
    case 'river'
        rgb = [81 124 150]/255;
    case 'leconte'
        rgb = [141 32 72]/255;
    case 'regalia'
        rgb = [117 74 126]/255;
    case 'sunsphere'
        rgb = [254 213 53]/255;
    case 'rock'
        rgb = [167 169 172]/255;
    case 'legacy'
        rgb = [87 149 132]/255;
    case 'summitt'
        rgb = [185 225 226]/255;
    case 'buckskin'
        rgb = [112 85 80]/255;
    case 'energy'
        rgb = [238 62 128]/255;
    case 'switchgrass'
        rgb = [171 193 120]/255;
    case 'fountain'
        rgb = [33 151 169]/255;
    case {'eureka','eureka!'}
        rgb = [235 234 100]/255;
    otherwise
            error('Unknown color: %s',name);
end
end