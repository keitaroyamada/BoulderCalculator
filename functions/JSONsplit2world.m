function [jOut] = JSONsplit2world(jData, orthoName)
    %convert local xy pos to world xy pos
    jOut  = jData;
    im_info = jData.images(1);
    im_info.file_name = orthoName;
    %im_info.width = orthoSize(2);
    %im_info.height = orthoSize(1);
    jOut.images(2:end) = [];
    jOut.images(1) = im_info;

    for i = 1:height(jData.annotations)
        im_id = jData.annotations(i).image_id;
        im_name = [];
        im_size = [];
        for n = 1:height(jData.images)
            if jData.images(n).id == im_id
                im_name = jData.images(n).file_name;
                im_size = [jData.images(n).height,jData.images(n).width];
                break
            end
        end
        if isempty(im_name)==1
            disp('There is no image')
            return
        end
        name_parts = strsplit(im_name,{'_','.'});
        r = str2double(name_parts{2});
        c = str2double(name_parts{3});
        %disp(strcat('c:',num2str(c),'_r:',num2str(r)))
    
        r_min = 1 + (r-1)*im_size(1);
        c_min = 1 + (c-1)*im_size(2);
    
        %box
        jOut.annotations(i).bbox(1) = c_min + jData.annotations(i).bbox(1);
        jOut.annotations(i).bbox(2) = r_min + jData.annotations(i).bbox(2);
    
        %seg
        jOut.annotations(i).segmentation(1:2:end) = c_min + jData.annotations(i).segmentation(1:2:end);
        jOut.annotations(i).segmentation(2:2:end) = r_min + jData.annotations(i).segmentation(2:2:end);
    
        jOut.annotations(i).image_id = 1;
        %textprogress(i, height(jData.annotations))
    end
end