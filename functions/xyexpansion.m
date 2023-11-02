function [out_x, out_y] = xyexpansion(x, y, ex)
    %seg:[x1,y1,x2,y2,...]
    %ex: expantion value from centroid
    %x = seg(1:2:end);
    %y = seg(2:2:end);

    c_x = median(x);
    c_y = median(y);
    
    [theta,rho] = cart2pol(x-c_x, y-c_y);


    [xe,ye] = pol2cart(theta,rho+ex);
    
    out_x = xe +c_x;
    out_y = ye +c_y;

    %output = seg;
    %output(1:2:end) = out_x;
    %output(2:2:end) = out_y;
end