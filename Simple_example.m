clear

%make instance
exr = volume_extractor;

%set data path
[ortho_n, ortho_p] = uigetfile('*.tif', 'Chose Planar Ortho image') 
[dsm_n, dsm_p]     = uigetfile(fullfile(ortho_p, '*.tif'), 'Chose Planar DSM image')
[json_n, json_p]   = uigetfile(fullfile(ortho_p, '*.json'), 'Chose Annotation Json data') 

%load ortho image
exr.loadOrtho(fullfile(ortho_p, ortho_n));

%load dsm
exr.loadDsm(fullfile(dsm_p, dsm_n));

%load overall annotation (COCO format)
exr.loadJsonForOrtho(fullfile(json_p, json_n));

%load split annotation (COCO format)
%exr.loadJsonForSplit('full path of json annotation of split images');

%load overall annotation (QGIS Geojson)
%[geojson_n, geojson_p] = uigetfile('*.geojson', 'Chose Annotation QGIS GeoJson data') 
%exr.json_indata = geojson2cocojson(fullfile(geojson_p, geojson_n), exr.im_info, exr.im_name);


%merge overlapped objects in the annotation
overlap_area_threshold = 0.15;
overlap_num_threshold = 1;
area_threshold = 0;
exr.opts.output_order = 'north2south';
exr.mergeObjects(overlap_area_threshold, overlap_num_threshold, area_threshold);

%set export settings
exr.save_dir = fullfile(geojson_p, '_outputs');
exr.opts = struct('save_mat', true,...%save mat
                  'save_csv', true,...%save csv
                  'save_kml', true,... %save object kml of shape, major axis, minor axis
                  'save_each_image', false,...%save trimed object images
                  'save_each_3dmodel',false,...%save trimed object model for matlab
                  'save_image_grid',true...
                  );

%extract & export
exr.extractVolume();



