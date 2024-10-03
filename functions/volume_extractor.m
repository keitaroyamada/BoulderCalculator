classdef volume_extractor < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access=public)
        save_dir;

        opts;

        im_name;
        im;
        im_info;
        dem;
        dem_info;

        json_indata;

        merged_raw;
        im_size_info;
        merged_world;

        height_threshold;
        test;
    end

    methods (Access=public)
        function obj = volume_extractor()
            %initiarise
            obj.save_dir = [];

            obj.opts = struct('output_order','large2small',...%[large2small, small2large,north2south,west2east]
                              'save_mat', true,...%save mat
                              'save_csv', true,...%save csv
                              'save_kml', true,... %save object kml of shape, major axis, minor axis
                              'save_each_image', false,...%save trimed object images
                              'save_each_3dmodel',false,...%save trimed object model for matlab
                              'save_image_grid',false ...%save analysied image data grid
                             );
            obj.height_threshold = [0.3, 20];%DEM hight threshold for vol
        end

        function [] = loadJsonForOrtho(obj, json_path)
            txt = fileread(json_path);
            obj.json_indata = jsondecode(txt);
            if isempty(obj.json_indata)
                disp('Json file is empty.')
            else
                disp('Json file loaded.')
            end

            %set default save dir
            [obj.save_dir,~,~] = fileparts(json_path);
        end

        function [] = loadJsonForSplit(obj, json_path)
            txt = fileread(json_path);
            obj.json_indata = JSONsplit2world(jsondecode(txt), obj.im_name);

            if isempty(obj.json_indata)
                disp('Json file is empty.')
            else
                disp('Json file loaded.')
            end

            %set default save dir
            [obj.save_dir,~,~] = fileparts(json_path);
        end

        function [] = loadOrtho(obj,ortho_path)
            [obj.im, obj.im_info]   = readgeoraster(ortho_path);
             [~,obj.im_name,~] = fileparts(ortho_path);
            disp('Ortho image loaded.')
        end

        function [] = loadDsm(obj,dsm_path)
            [obj.dem, obj.dem_info] = readgeoraster(dsm_path);
            disp('DSM loaded.')
        end

        function [] = mergeObjects(obj, th_orverlap, th_duplication, th_area)
            warning('off','all')
            %convert struct to table
            %get class
            json_temp = struct2table(obj.json_indata.annotations);
            [~,idx_categories] = ismember(json_temp.category_id, [obj.json_indata.categories.id]');
            ob_name = string({(obj.json_indata.categories(idx_categories,:).name)}');
            
            %get otherproperties
            ob_bdbox = ([obj.json_indata.annotations.bbox]');
            ob_seg = ({obj.json_indata.annotations.segmentation}');
            
            %get segdata
            ob_pgn = cell(size(ob_seg));
            ob_ct_x = zeros(size(ob_seg));
            ob_ct_y = zeros(size(ob_seg));
            ob_area = zeros(size(ob_seg));
            for i= 1:size(ob_seg,1)
                seg = ob_seg{i,1};
                ob = polyshape(seg(1:2:end), seg(2:2:end));
                ob = rmholes(ob);
                ob_pgn{i,1} = ob;
                [ob_ct_x(i,1), ob_ct_y(i,1)] = centroid(ob);
                ob_area(i,1) = area(ob);
            end

            ob_stats = [table(ob_name, ob_bdbox, ob_seg, ob_area, ob_ct_x, ob_ct_y), cell2table(ob_pgn)];
            it = struct2table(obj.json_indata.images,'AsArray',1);
            obj.im_size_info = [ones(height(it)), it.height, ones(height(it)), it.width];
        
            %initialize
            ob_stats = sortrows(ob_stats, 'ob_area', 'descend');%'as descend'
            ob_stats(find(isnan(ob_stats{:,5})),:) = [];%remove no shape lines
            ob_stats(find(ob_stats{:,4} <= th_area),:)=[];%area threshold

            merged_raw = [];
            
            N = height(ob_stats);
            wc = 0;
            while height(ob_stats)>0
                wc = wc+1;
                i = 1;
                
                nr_idx = knnsearch(ob_stats{:,5:6},ob_stats{i,5:6},'K',30);%search centroid
                C   = intersect(ob_stats{i,7}, ob_stats{nr_idx,7});
        
                baseAreaRate1 = area(C)./area(ob_stats{i,7});
                baseAreaRate2 = area(C)./area(ob_stats{nr_idx,7});
                
                idx_base = find(any([baseAreaRate1>th_orverlap, baseAreaRate2>th_orverlap],2));%check iou
        
                idx = idx_base;
        
                if numel(idx)>=th_duplication
                    %if objects overlap 
                    U = union(ob_stats{nr_idx(idx),7});
                    temp = ob_stats(nr_idx(idx(1)),:);
                    temp.ob_bdbox        = [min(U.Vertices(:,1)),...
                                             min(U.Vertices(:,2)),...
                                             max(U.Vertices(:,1)) - min(U.Vertices(:,1)),...
                                             max(U.Vertices(:,2)) - min(U.Vertices(:,2))];
                    temp.ob_seg          = zeros(1,numel(U.Vertices));
                    temp.ob_seg(1:2:end) = U.Vertices(:,1);
                    temp.ob_seg(2:2:end) = U.Vertices(:,2);
                    temp.ob_seg          = {temp.ob_seg};
                    temp.ob_area         = area(U);
                    [temp.ob_ct_x, temp.ob_ct_y] = centroid(U);
                    temp_pgn              = U;
                    temp.duplicate        = numel(idx);
                    
                    %add into new list
                    merged_raw = [merged_raw; temp];
                    %delete old list
                    ob_stats(nr_idx(idx),:) = [];
                else
                    %delete old list
                    ob_stats(nr_idx(idx),:) = [];
                end
            end
            warning('on')

            %sort
            switch obj.opts.output_order
                case 'large2small'
                    merged_raw = merged_raw;
                case'small2large'
                    merged_raw = sortrows(merged_raw,'ob_area','ascend');
                case'north2south'
                    merged_raw = sortrows(merged_raw,'ob_ct_y','ascend');
                case'south2north'
                    merged_raw = sortrows(merged_raw,'ob_ct_y','descend');
                case'west2east'
                    merged_raw = sortrows(merged_raw,'ob_ct_x','ascend');
                case'east2west'
                    merged_raw = sortrows(merged_raw,'ob_ct_x','descend');
                otherwise
                    merged_raw = merged_raw;
            end

            obj.merged_raw = merged_raw;
            disp('Objects merged.')
        end

        function [] = extractVolume(obj)
            
            %make save directories
            mkdir(fullfile(obj.save_dir));

            if obj.opts.save_each_image == true
                mkdir(fullfile(obj.save_dir, 'object images'))
            end
            if obj.opts.save_each_3dmodel == true
                mkdir(fullfile(obj.save_dir, 'object 3dmodels'))
            end
        
            %get plane position
            disp(strcat('Calc object volume from DEM'))
            dem_x0     = obj.dem_info.XWorldLimits;
            dem_y0     = obj.dem_info.YWorldLimits;
            dem_scale  = [obj.dem_info.CellExtentInWorldX, obj.dem_info.CellExtentInWorldY];
            
            im_x0      = obj.im_info.XWorldLimits;
            im_y0      = obj.im_info.YWorldLimits;
            im_scale   = [obj.im_info.CellExtentInWorldX, obj.im_info.CellExtentInWorldY];
            
            merged_world    = [];
            for i=1:height(obj.merged_raw)
                textprogress(i, height(obj.merged_raw))
                description = "";

                %make trimed ORTHO
                %get obj location&expantion
                ex = 5; % expanding size
                seg  = obj.merged_raw.ob_seg{i};
                [im_xe, im_ye] = xyexpansion(seg(1:2:end), seg(2:2:end), ex);%expanding object shape
                
                im_xmin = round(min(im_xe));
                im_xmax = round(max(im_xe));
                im_ymin = round(min(im_ye));
                im_ymax = round(max(im_ye));
        
                if isnan(im_xmin)==1||isnan(im_xmax)==1||isnan(im_ymin)==1||isnan(im_ymax)==1
                    continue
                end
        
                %find polygon in the target area
                im_frame      = polyshape(im_xe, im_ye);
                poly_in_frame = intersect(obj.merged_raw.ob_pgn, im_frame);
                idx_in        = find(area(poly_in_frame)~=0);
        
                %get object centre in plane
                [ob_ct_plx, ob_ct_ply] = geopix2loc([obj.merged_raw.ob_ct_x(i), obj.merged_raw.ob_ct_y(i)], obj.im_info, 'plane');
            
                %get object seg in plane 
                [ob_seg_plx, ob_seg_ply] = geopix2loc(seg, obj.im_info, 'plane');
        
                %transfer from xy location in ortho to in dem 
                [dem_xe, dem_ye] = transfergeopoint([im_xe', im_ye'], obj.im_info, obj.dem_info);
            
                %get boxboundaries
                dem_xmin = min(dem_xe,[],1,"omitnan");
                dem_xmax = max(dem_xe);
                dem_ymin = min(dem_ye);
                dem_ymax = max(dem_ye);
                
                %get trim area
                dem_pad = round(2 / obj.dem_info.CellExtentInWorldX);
                dem_xfrom   = dem_xmin-dem_pad;
                dem_xto     = dem_xmax+dem_pad;
                dem_yfrom   = dem_ymin-dem_pad;
                dem_yto     = dem_ymax+dem_pad;
                
                if isnan(dem_xfrom)==1||isnan(dem_xto)==1||isnan(dem_yfrom)==1||isnan(dem_yto)==1
                    continue
                end
            
                if dem_yto>obj.dem_info.RasterSize(1)||dem_xto>obj.dem_info.RasterSize(2)
                    dem_xto   = obj.dem_info.RasterSize(2);
                    dem_yto   = obj.dem_info.RasterSize(1);
                elseif dem_xfrom<1||dem_yfrom<1
                    dem_xfrom = 1;
                    dem_yfrom = 1;
                end
            
                if dem_yfrom>=dem_yto||dem_xfrom>=dem_xto
                    continue
                end
            
                dem_trim = obj.dem(dem_yfrom:dem_yto, dem_xfrom:dem_xto);
                dem_trim(dem_trim<-100)=NaN;%remove too low meth

                %noise reduction
                dem_retrim     = 0;
                dem_smooth     = 20;%unit:cm
                dem_smooth_pix = round(dem_smooth/(obj.dem_info.CellExtentInWorldX*100));

                dem_trim_rd = medfilt2(im2double(dem_trim),[dem_smooth_pix, dem_smooth_pix],'symmetric');%omit noise
                dem_trim_rd = dem_trim_rd(1+dem_retrim:end-dem_retrim, 1+dem_retrim:end-dem_retrim,:);%trim edges
        
                if numel(find(dem_trim_rd==0))>0
                    continue
                end
        
                %make BWimage
                im_shape = zeros(size(dem_trim_rd));%mask base;
                [dem_xtemp, dem_ytemp] = transfergeopoint([seg(1:2:end)', seg(2:2:end)'], obj.im_info, obj.dem_info);
                seg_temp = seg;
                seg_temp(1:2:end) = dem_xtemp - (dem_xmin-dem_pad+dem_retrim);
                seg_temp(2:2:end) = dem_ytemp - (dem_ymin-dem_pad+dem_retrim);
                im_shape = insertShape(im_shape, 'FilledPolygon', seg_temp, 'Color','white','Opacity',1);
                rprops = regionprops(im_shape(:,:,1),{'MajorAxisLength','MinorAxisLength','Orientation'});
        
                if size(rprops,1)==0
                    continue
                end
        
                MajorAxis   = rprops.MajorAxisLength * mean(dem_scale);
                MinorAxis   = rprops.MinorAxisLength * mean(dem_scale);
                Orientation = rprops.Orientation;
        
                %make object mask(expaned shape)
                mask = zeros(size(dem_trim_rd));%mask base 
                mask_seg  = obj.merged_raw.ob_seg{i};
                mask_seg(1:2:end) = dem_xe - (dem_xmin-dem_pad+dem_retrim);
                mask_seg(2:2:end) = dem_ye - (dem_ymin-dem_pad+dem_retrim);
                mask = insertShape(mask, 'FilledPolygon', mask_seg, 'Color','white','Opacity',1);
                mask(mask==0)=nan; 
            
                %make basement mask(expaned shape)
                mask_pad = round(obj.dem_info.CellExtentInWorldX * 100 * (1.5));
                mask_basement = ones(size(dem_trim_rd));%mask base
                for m=1:numel(idx_in)
                    mask_seg  = obj.merged_raw.ob_seg{idx_in(m)};
                    [im_xtemp, im_ytemp]   = xyexpansion(mask_seg(1:2:end), mask_seg(2:2:end), mask_pad);%礫形状の膨張
                    [dem_xtemp, dem_ytemp] = transfergeopoint([im_xtemp', im_ytemp'], obj.im_info, obj.dem_info);
                    if isnan(dem_xtemp)==1|isnan(dem_ytemp)
                        continue
                    end
                    mask_seg(1:2:end) = dem_xtemp - (dem_xmin-dem_pad);
                    mask_seg(2:2:end) = dem_ytemp - (dem_ymin-dem_pad);
            
                    mask_basement = insertShape(mask_basement, 'FilledPolygon', mask_seg, 'Color','black','Opacity',1);
                end
                mask_basement(mask_basement~=1)=nan;

                %calc height of basement
                fit_type='poly33';
                basement = dem_trim_rd.*mask_basement(:,:,1);
                [X,Y]    = meshgrid([1:size(basement,2)],[1:size(basement,1)]);
                idx      = find(not(isnan(basement)));
                if length(idx)<10
                    continue
                end
                [B, gof] = fit([reshape(X(idx),[],1),reshape(Y(idx),[],1)],reshape(double(basement(idx)),[],1),fit_type);
                fit_rmse = gof.rmse;
        
                %height_basement = B.p00;
                height_basement = B(X,Y);
            
                %calc relative height
                height_obj = dem_trim_rd - height_basement;
           
                %calc height
                mask_vol = imcomplement(mask_basement);
                mask_vol(isnan(mask_vol))=1;
                dem_obj = height_obj.*mask_vol(:,:,1);
                dem_obj(dem_obj<0)=0;
        
                base_obj = height_basement.*mask_vol(:,:,1);
                base_obj(base_obj==0)=nan;
        
                %check fitness & calc height
                height_obj2 = dem_obj;
                height_obj2(height_obj2==0)=nan;
                meanh_obj  = median(median(height_obj2,'omitnan'),'omitnan');
                maxh_obj   = double(max(max(height_obj2)));
                if fit_rmse>1
                    base_height_obj = nan;
                    description = [description, "Height is set as Nan because it exceeds the fitting rmse."];
                else
                    base_height_obj = median(median(base_obj,'omitnan'),'omitnan');
                end
            
                isboulder = 0;
                if maxh_obj>obj.height_threshold(1)&maxh_obj<=obj.height_threshold(2)
                    vol = sum(dem_obj,"all","omitnan")*dem_scale(1)*dem_scale(2);
                    if vol <= 0
                        dsm_vol_m3 = 0;
                    else
                        dsm_vol_m3 = vol;
                        isboulder = 1;
                    end
                else
                    dsm_vol_m3 = nan;
                    description = [description, "Vdsm is set as Nan because it exceeds the height thresholds."];
                end
        
                %------------------------------------------------------------------
                if obj.opts.save_each_image == true
                    %save trim images
                    [im_xfrom, im_yfrom] = transfergeopoint([dem_xfrom+dem_retrim, dem_yfrom+dem_retrim], obj.dem_info, obj.im_info);
                    [im_xto, im_yto]     = transfergeopoint([dem_xto-dem_retrim, dem_yto-dem_retrim], obj.dem_info, obj.im_info);
                    [im_xseg, im_yseg]   = transfergeopoint([[round(seg_temp(1:2:end)),round(seg_temp(1))]', [round(seg_temp(2:2:end)),round(seg_temp(2))]'], obj.dem_info, obj.im_info);
        
                    %im_trim = im(im_ymin:im_ymax, im_xmin:im_xmax, 1:3);
                    im_trim = obj.im(im_yfrom:im_yto, im_xfrom:im_xto, 1:3);
        
                    f=figure('visible','off');
                    imshow(im_trim)
                    hold on
                    plot([seg(1:2:end),seg(1)]'-im_xfrom, [seg(2:2:end),seg(2)]'-im_yfrom,'w',"LineWidth",2)%seg based on original resolution
                    hold on
                    plot(im_xseg, im_yseg,'r',"LineWidth",2)%seg based on DEM resolution
                    
                    exportgraphics(f,fullfile(obj.save_dir,'object images', strcat(obj.im_name, num2str(i,'_obj%03d'),'.jpg')))
                    close(f)
                end
                
                if obj.opts.save_each_3dmodel
                    %save object
                    meshColR = mask_basement(:,:,1);
                    meshColG = mask_basement(:,:,1);
                    meshColB = mask_basement(:,:,1);
        
                    %backgroundcolor
                    meshColR(meshColR==1) = 0.00;
                    meshColG(meshColG==1) = 0.00;
                    meshColB(meshColB==1) = 1.00;
        
                    %gravel color
                    meshColR(isnan(meshColR))=1.00;
                    meshColG(isnan(meshColG))=0.65;
                    meshColB(isnan(meshColB))=0.00;
        
                    meshCol = cat(3, meshColR,meshColG,meshColB);
            
                    f=figure('visible','on');%if set off, figure does not open automatically.
                    a=mesh(dem_trim_rd);
                    a.CData = meshCol;
                    pbaspect([1, 1, 0.5])
                    hold on
                    plot(B)
                    xt=[-10:1:15];
                    xticklabels(string(num2str(xt')))
                    xticks(xt/im_scale(1))
                    yticklabels(string(num2str(xt')))
                    yticks(xt/im_scale(1))
                    xlabel('W-E (m)')
                    ylabel('N-S (m)')
                    zlabel('Altitude (m)')
                    view(-45,50)
                    set(gca, 'YDir', 'reverse')
        
                    saveas(f,fullfile(obj.save_dir,'object 3dmodels',strcat(obj.im_name, num2str(i,'_obj%03d'),'.fig')))
                    %exportgraphics(f,fullfile(save_dir,'object mesh',strcat(im_name, num2str(i,'_obj%03d'),'.png')))
                    close(f)
                end

                %calc cross sectional area 
                major_theta = Orientation;%(Orientation)/180*pi;
                minor_theta = Orientation + 90; %(Orientation+90)/180*pi;

                [dsm_cx_pix, dsm_cy_pix] = transfergeopoint([obj.merged_raw.ob_ct_x(i), obj.merged_raw.ob_ct_y(i)], obj.im_info, obj.dem_info);

                trim_dsm_cx_pix = dsm_cx_pix - dem_xfrom;
                trim_dsm_cy_pix = dsm_cy_pix - dem_yfrom;
                %obj.test = {dem_obj, trim_dsm_cx_pix,trim_dsm_cy_pix,major_theta,minor_theta};

                major_cross = imCrossSection(dem_obj, trim_dsm_cx_pix, trim_dsm_cy_pix, major_theta);
                minor_cross = imCrossSection(dem_obj, trim_dsm_cx_pix, trim_dsm_cy_pix, minor_theta);

                major_cross_area = sum(major_cross(:,3) .* major_cross(:,4)) * mean([dem_scale(1),dem_scale(2)]);
                minor_cross_area = sum(minor_cross(:,3) .* minor_cross(:,4)) * mean([dem_scale(1),dem_scale(2)]);
                
                %make data for summary properties
                ob_area_m2  = obj.merged_raw.ob_area(i) * im_scale(1)*im_scale(2);
                ob_id       = string(num2str(i,'ob_%03d'));
                ob_name     = string(obj.merged_raw.ob_name{i});
                crs         = obj.im_info.ProjectedCRS.Name;
                abc_vol_m3  = MajorAxis.*MinorAxis.*maxh_obj;
                eli_vol_m3  = 4/3*pi*abc_vol_m3/8;
                calc_id = i;
                if fit_rmse>1
                    isboulder = 0;
                end
        
                [ob_ct_lat, ob_ct_lon] = projinv(obj.im_info.ProjectedCRS, ob_ct_plx{:}, ob_ct_ply{:});
                merged_world = [merged_world; table(calc_id, ob_id, ob_name, isboulder,... %1:4
                                                    ob_ct_lat, ob_ct_lon, ob_ct_plx, ob_ct_ply, base_height_obj,...%5:9
                                                    MajorAxis, MinorAxis, maxh_obj, meanh_obj, Orientation,...%10:14
                                                    ob_area_m2, major_cross_area, minor_cross_area,...%15:17
                                                    abc_vol_m3, eli_vol_m3, dsm_vol_m3,...%18:20
                                                    fit_rmse, ob_seg_plx, ob_seg_ply, crs, {description})];%21:25
            end
            %==============================================================
            %save result
            obj.merged_world = merged_world;
            merge_raw = obj.merged_raw;
            im_info   = obj.im_info;
            dem_info  = obj.dem_info;
            
            if  obj.opts.save_mat == true
                %save mat file
                save(fullfile(obj.save_dir, 'detectron2_results.mat'),'merged_world', 'merge_raw', 'im_info', 'dem_info');
            end

            if obj.opts.save_csv == true
                %save csv
                n = [1:4, 5:6,9, 21,7:8, 10:12,14, 15:17, 18:20, 21, 25];
                out_data = merged_world(:,n);

                writetable(out_data, fullfile(obj.save_dir,strcat(obj.im_name, '_results.csv')))
            end

            if obj.opts.save_kml==true
                %save kml
                save_kml_with_nandata = true;
                if save_kml_with_nandata==false
                    idx = find(not(isnan(merged_world.dsm_vol_m3)));
                else
                    idx = [1:height(merged_world)]';
                end
        
                ob_name  = strings(size(idx,1),1);   
                seg_lat   = cell(size(idx,1),1);
                seg_lon   = cell(size(idx,1),1);
                seg_elv   = cell(size(idx,1),1);
                major_lat = cell(size(idx,1),1);
                major_lon = cell(size(idx,1),1);
                major_elv = cell(size(idx,1),1);
                minor_lat = cell(size(idx,1),1);
                minor_lon = cell(size(idx,1),1);
                minor_elv = cell(size(idx,1),1);
                for n=1:size(idx,1)
                    i   = idx(n);
                    x   = merged_world.ob_seg_plx{i};
                    y   = merged_world.ob_seg_ply{i};
                    [lat, lon]   = projinv(obj.im_info.ProjectedCRS, x, y);
                    elv = zeros(1,size(lat,2));
            
                    %get object name
                    ob_name{n,1} = merged_world.ob_id{i};
            
                    %get shape
                    seg_lat{n,1} = [lat,lat(end)];
                    seg_lon{n,1} = [lon,lon(end)];
                    seg_elv{n,1} = [elv,elv(end)];
            
                    %get axies
                    major_theta        = (merged_world.Orientation(i)/180)*pi;
                    minor_theta        = (merged_world.Orientation(i)+90)/180*pi;
                    major_len          = merged_world.MajorAxis(i)/2;
                    minor_len          = merged_world.MinorAxis(i)/2;
                    [win_x, win_y]     = pol2cart([major_theta;major_theta;minor_theta;minor_theta], [major_len;-major_len;minor_len;-minor_len]);
                    ob     = polyshape(x, y); 
                    [cx,cy] = centroid(ob);
        
                    [win_lat, win_lon] = projinv(obj.im_info.ProjectedCRS, cx+win_x, cy+win_y);
                
                    major_lat{n,1} = [win_lat(1:2)'];
                    major_lon{n,1} = [win_lon(1:2)'];
                    major_elv{n,1} = [0 0];
                    minor_lat{n,1} = [win_lat(3:4)'];
                    minor_lon{n,1} = [win_lon(3:4)'];
                    minor_elv{n,1} = [0 0];
                end
                segDataTable = table(ob_name, seg_lat, seg_lon, seg_elv);
                segDataTable.Properties.VariableNames = {'name','latitude','longitude','elevation'};
                majorDataTable = table(ob_name, major_lat, major_lon, major_elv);
                majorDataTable.Properties.VariableNames = {'name','latitude','longitude','elevation'};
                minorDataTable = table(ob_name, minor_lat, minor_lon, minor_elv);
                minorDataTable.Properties.VariableNames = {'name','latitude','longitude','elevation'};
                
                saveLineKML(fullfile(obj.save_dir,strcat(obj.im_name, 'Shapes.kml')), segDataTable);
                saveLineKML(fullfile(obj.save_dir,strcat(obj.im_name, 'MajorAxis.kml')), majorDataTable);
                saveLineKML(fullfile(obj.save_dir,strcat(obj.im_name, 'MinorAxis.kml')), minorDataTable);
            end

            if obj.opts.save_image_grid
                    %get grid data
                    grid_data = obj.json_indata.images.image_grid;

                    %save analised image grid
                    grid_name   = strings(height(grid_data),1);  
                    grid_lat    = cell(height(grid_data),1);
                    grid_lon    = cell(height(grid_data),1);
                    grid_elv    = cell(height(grid_data),1);
                    for n = 1:height(grid_data)
                        y0 = grid_data(n,1) + 1;
                        y1 = grid_data(n,2) + 1;
                        x0 = grid_data(n,3) + 1;
                        x1 = grid_data(n,4) + 1;

                        %x = [x0, x1, x1, x0, x0]';
                        %y = [y0, y0, y1, y1, y0]';
                        xy = [x0, y0, x1, y0, x1, y1, x0, y1, x0, y0];

                        %get object seg in plane 
                        [plx, ply] = geopix2loc(xy, obj.im_info, 'plane');

                        [lat, lon]   = projinv(obj.im_info.ProjectedCRS, plx{:}, ply{:});
                        grid_lat{n,1} = [lat];
                        grid_lon{n,1} = [lon];
                        grid_elv{n,1} = [repmat(50,5,1)];
                        grid_name{n,1}= num2str(n,'%03d');

                    end

                    gridDataTable = table(grid_name, grid_lat, grid_lon, grid_elv);
                    gridDataTable.Properties.VariableNames = {'name','latitude','longitude','elevation'};
                    saveLineKML(fullfile(obj.save_dir,strcat(obj.im_name, 'Grid.kml')), gridDataTable);
            end

        end
    end
end