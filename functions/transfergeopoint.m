function [to_x, to_y] = transfergeopoint(xy, Rfrom, Rto)
    %xy:transfer points [x, y]
    %Rfrom: Spatial reference of geotif
    %Rto: Spatial reference of geotif

    %get spatial info
    from_x0      = Rfrom.XWorldLimits;
    from_y0      = Rfrom.YWorldLimits;
    from_scale   = [Rfrom.CellExtentInWorldX, Rfrom.CellExtentInWorldY];

    to_x0        = Rto.XWorldLimits;
    to_y0        = Rto.YWorldLimits;
    to_scale     = [Rto.CellExtentInWorldX, Rto.CellExtentInWorldY];

    %convert image location to plane projection
    pl_x = from_x0(1) + ((xy(:,1) - 0.5) * from_scale(1));%pl_x = from_x0(1) + ((xy(:,1) - 1) * from_scale(1));
    pl_y = from_y0(2) - ((xy(:,2) - 0.5) * from_scale(2));%pl_y = from_y0(2) - ((xy(:,2) - 1) * from_scale(2));

    %convert plane projection to dem location
    to_x = 1 + round(abs(pl_x - to_x0(1)) / to_scale(1));%to_x = 1 + round(abs(pl_x - to_x0(1)) / to_scale(1));
    to_y = 1 + round(abs(pl_y - to_y0(2)) / to_scale(2));%to_y = 1 + round(abs(pl_y - to_y0(2)) / to_scale(2));
end