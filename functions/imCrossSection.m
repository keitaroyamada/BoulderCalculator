function [output] = imCrossSection(im, cx, cy, theta)
    %cx: centroid
    %cy: centroid
    %theta: angle from 
    %output: [x, ,y val, xy diff]
    
    a = tan(deg2rad(rem(theta-90,360)));

    %x_pix = [1: size(im,2)]';
    %y_pix = round((x_pix * a) + (cy - a*cx))

    x_pix = [1:size(im,2)]';
    y_pix = ((x_pix * a) + (cy - a*cx));

    if y_pix(1)<1&y_pix(end)>size(im,1)||y_pix(end)<1&y_pix(1)>size(im,1)
        %nearly vertical
        y_pix = [1:size(im,1)]';
        x_pix = (y_pix - (cy - a*cx)) /a;
    end

    %check 
    x_idx = find(round(x_pix)>0&round(x_pix)<=size(im,2));
    y_idx = find(round(y_pix)>0&round(y_pix)<=size(im,1));
    xy_idx = intersect(x_idx,y_idx);
    x_pix = x_pix(xy_idx);
    y_pix = y_pix(xy_idx);

    output = zeros(length(x_pix),4);

    for n = 1:length(x_pix)
        output(n,1) = round(x_pix(n));
        output(n,2) = round(y_pix(n));
        output(n,3) = im(round(y_pix(n)), round(x_pix(n)),1);
        if isnan(output(n,3))
            output(n,3) = 0;
        end

        if n>1
            if output(n,3)~=0
                output(n,4) = sqrt((x_pix(n) - x_pix(n-1))^2 + (y_pix(n) - y_pix(n-1))^2);
            end
        else
            output(n,4) = 0;
        end
    end



end