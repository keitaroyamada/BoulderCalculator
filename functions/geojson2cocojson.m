function varargout = geojson2cocojson(geoJson, ortho_info, varargin)
    
    %load data
    if ischar(ortho_info)
        [~, im_info] = readgeoraster(ortho_info);
        [~,im_name,~] = fileparts(ortho_info);
    else
        im_info = ortho_info;
        if size(varargin,1)>0
            im_name = varargin{1};
        else
            im_name = "geotiff"
        end
    end

    if ischar(geoJson)
        txt = fileread(geoJson);
        jData = jsondecode(txt);
    else
        jData = geoJson;
    end
    
    %initiarize
    outData = struct('licenses',[{struct('name','','id',0,'url','')}], ...
                     'info',{struct('contributor','','date_created','','description','','url','','version','','year','')},...
                     'categories',[], ...
                     'images',[], ...
                     'annotations',[]);
    
    categoryTable = [];
    for i=1:size(jData.features,1)
        if isempty(jData.features(i).properties.id)
            id = 0;
        else
            id = jData.features(i).properties.id;
        end
        categoryTable = [categoryTable; table(id, id, "none")];
    end
    categoryTable = unique(categoryTable);     
    categoryTable.Properties.VariableNames = {'id','name','supercategory'};
    
    imageTable      = [];
    imageTableNames      = {'id','width','height','file_name','license','flickr_url','coco_url','date_captured','window_size','window_step_rate','image_grid'};

    annotationTable = [];    
    annotationTableNames = {'id','image_id','category_id','segmentation','area','bbox','iscrowd','attributes'};
    
    categoryTableNames   = {'id','name','supercategory'};

    image_height = im_info.RasterSize(1);
    image_width  = im_info.RasterSize(2);
    image_name   = string(im_name);
    
    %main
    an_id = 1;
    for i=1:size(jData.features,1)
        %convert position fron world to image pix
        x = reshape(jData.features(i).geometry.coordinates(:,:,:,1),[],1);
        y = reshape(jData.features(i).geometry.coordinates(:,:,:,2),[],1);
        
        xr = (x - im_info.XWorldLimits(1)) / (im_info.XWorldLimits(2) - im_info.XWorldLimits(1));
        yr = (y - im_info.YWorldLimits(1)) / (im_info.YWorldLimits(2) - im_info.YWorldLimits(1));
    
        x_inim = 1 + image_width * xr;
        y_inim = 1 + image_height * (1-yr);
    
        %get images
        im_id = 1;
            
        imageTable = [imageTable;...
                      table(im_id, ...
                            image_width, ...
                            image_height, ...
                            image_name,...
                            0,...
                            "",...
                            "",...
                            "",...
                            [image_height, image_width],...
                            [1,1],...
                            [0 image_height 0 image_width])];
           
        %get category
        curr_cat_id = id;
        
        %get annotation
        segmentation = zeros(1, size(x,1)*2);
        segmentation(1:2:end) = x_inim;
        segmentation(2:2:end) = y_inim;
        seg_area =  round(polyarea(segmentation(1:2:end),segmentation(2:2:end)));
        if seg_area==0
            segmentation=[];
        end
    
        x0 = min(x);
        x1 = max(x);
        y0 = min(y);
        y1 = max(y);
        
        annotationTable = [annotationTable;... 
                           table(an_id,...                      %id
                           im_id,...                            %image_id
                           curr_cat_id ,...                     %category_id
                           {segmentation},...                   %segmentation
                           seg_area,...                         %area
                           {round([x0; y0; x1-x0; y1-y0])},...  %bbox
                           0,...                                %iscrowd
                           struct('occluded',logical(0)))];     %attributes
        an_id = an_id + 1;
    end
    
    %make json
    annotationTable.Properties.VariableNames = annotationTableNames;
    imageTable.Properties.VariableNames      = imageTableNames;
    categoryTable.Properties.VariableNames   = categoryTableNames;
    
    outData.images      = table2struct(imageTable);
    outData.annotations = table2struct(annotationTable);
    if height(categoryTable)==1
        outData.categories = {table2struct(categoryTable)};
    else
        outData.categories = table2struct(categoryTable);
    end
    
    outjData = jsonencode(outData);
    
    %add[] around segmentation points
    outjData = strrep(outjData,'],"area":',']],"area":');%Caution
    outjData = strrep(outjData,'"segmentation":[','"segmentation":[[');%Caution
    
    %add[] around licenses
    outjData = strrep(outjData,'"licenses":{','"licenses":[{');%Caution
    outjData = strrep(outjData,'},"info"','}],"info"');%Caution
    
    %add[] around categories
    outjData = strrep(outjData,'categories":{','categories":[{');%Caution
    outjData = strrep(outjData,'},"images"','}],"images"');%Caution
    
    if nargout==1
        varargout{1} = outData;
    else
        %write json
        fid = fopen(fullfile(p, '_annotations.coco.json'),'w');
        fprintf(fid, outjData);
        fclose(fid);
    end
    
end