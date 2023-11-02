function [x_loc, y_loc] = geopix2loc(xy,R, type)
%calc location (lat, lon, x,y, etc) from pix points in geotif
%xy,:pix index in geotif {[x1,y1,x2,y2,...],[x5,y5,x6,y6,...],...}
%R: geotif spacial reference
    %get spacial reference
    switch type
        case 'plane'
            im_x0      = R.XWorldLimits;
            im_y0      = R.YWorldLimits;
            im_scale   = [R.CellExtentInWorldX, R.CellExtentInWorldY];
        case 'lonlat'
            im_x0      = R.LongitudeLimits;
            im_y0      = R.LatitudeLimits;
            im_scale   = [R.CellExtentInLongitude, R.CellExtentInLatitude];
    end

    
    xloc = cell(length(xy),1);
    yloc = cell(length(xy),1);
    if iscell(xy)
        for i = 1:size(xy,1)
            seg   = xy{i,1};
            x_loc{i,1} = [im_x0(1) + ((seg(1:2:end) - 1) * im_scale(1))];
            y_loc{i,1} = [im_y0(2) - ((seg(2:2:end) - 1) * im_scale(2))];
        end
    else
        seg   = xy;
        x_loc{1,1} = [im_x0(1) + ((seg(1:2:end) - 1) * im_scale(1))];
        y_loc{1,1} = [im_y0(2) - ((seg(2:2:end) - 1) * im_scale(2))];
    end
end