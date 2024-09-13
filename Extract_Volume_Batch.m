clear

%batch get process
json_dir = 'C:\Users\';
json_dir = uigetdir(json_dir,'Chose json data folder');
json_list= dir(fullfile(json_dir,'*.json'))
im_dir   = uigetdir(json_dir, 'Chose ortho & dsm image folder');

for json_n =1:size(json_list,1)
    %make instance
    exr = volume_extractor;

    %load ortho and dsm
    [~,im_name_temp,~] = fileparts(json_list(json_n).name);
    im_name  = strcat(im_name_temp,'.tif');
    dem_name = strrep(im_name,'ortho','dsm');
    disp(strcat("[",num2str(json_n),"/",num2str(size(json_list,1)),"]:",im_name));
    
    %load ortho image
    exr.loadOrtho(fullfile(im_dir, im_name))
    
    %load dsm
    exr.loadDsm(fullfile(im_dir, dem_name))

    %load annotation
    exr.loadJsonForOrtho(fullfile(json_dir, json_list(json_n).name));
    %exr.loadJsonForSplit(fullfile(json_dir, json_list(json_n).name));
    
    %merge overlapped objects in the annotation
    overlap_area_threshold = 0.15;
    overlap_num_threshold = 1;
    area_threshold = 0;
    exr.mergeObjects(overlap_area_threshold, overlap_num_threshold, area_threshold)
    
    %set export settings
    exr.save_dir = fullfile(json_dir,'outputs',im_name);
    exr.opts = struct('save_mat', true,...%save mat
                      'save_csv', true,...%save csv
                      'save_kml', true,... %save object kml of shape, major axis, minor axis
                      'save_each_image', true,...%save trimed object images
                      'save_each_3dmodel',true ...%save trimed object model for matlab
                      );
    
    %extract & export
    exr.extractVolume()
end



