clear

%make instance
exr = volume_extractor;

%load ortho image
exr.loadOrtho('full path of ortho image')

%load dsm
exr.loadDsm('full path of dsm image')

%load annotation
exr.loadJsonForOrtho('full path of json annotation of ortho image');
%exr.loadJsonForSplit('full path of json annotation of split images');

%merge overlapped objects in the annotation
overlap_area_threshold = 0.15;
overlap_num_threshold = 1;
area_threshold = 0;
exr.mergeObjects(overlap_area_threshold, overlap_num_threshold, area_threshold)

%set export settings
exr.save_dir = fullfile('full path of save directory')
exr.opts = struct('save_mat', true,...%save mat
                  'save_csv', true,...%save csv
                  'save_kml', true,... %save object kml of shape, major axis, minor axis
                  'save_each_image', false,...%save trimed object images
                  'save_each_3dmodel',false,...%save trimed object model for matlab
                  'save_image_grid',true...
                  );

%extract & export
exr.extractVolume()





