%% Match field UTM points to volume_extractor boulder polygons.
% The original ob_id is preserved. The source MAT file name is written
% separately as matched_mat_file.

clear; clc;

%% Config
work_dir = 'C:\working\folder\path';

excel_file = fullfile(work_dir, 'source_position.xlsx');
mat_files = dir(fullfile(work_dir, '*detectron2_results.mat'));

% Excel column assignment. Row 1 is ignored because it is empty.
excel_first_data_row = 2;
excel_input_id_col = 1;   % A column
excel_gps_id_col   = 2;   % B column
excel_utm_x_col    = 4;   % D column
excel_utm_y_col    = 5;   % E column

include_boundary = true;
multiple_match_rule = 'nearest_centroid';

output_csv = fullfile(work_dir, 'matched_boulder_positions.csv');
output_kml = fullfile(work_dir, 'matched_boulder_positions.kml');
write_kml = true;

%% Read Excel points
point_raw = readcell(excel_file);
point_raw = point_raw(excel_first_data_row:end, :);

input_id = string(point_raw(:, excel_input_id_col));
gps_id = string(point_raw(:, excel_gps_id_col));
las_file_name = string(point_raw(:, 3));
utm_x = cellfun(@double, point_raw(:, excel_utm_x_col));
utm_y = cellfun(@double, point_raw(:, excel_utm_y_col));

valid_point = isfinite(utm_x) & isfinite(utm_y);

%% Read MAT boulder polygons
boulder_table = table();
kml_projected_crs = [];

for f = 1:numel(mat_files)
    mat_file_name = mat_files(f).name;
    mat_path = fullfile(mat_files(f).folder, mat_file_name);
    S = load(mat_path, 'merged_world', 'im_info');

    if ~isfield(S, 'merged_world')
        warning('%s does not contain merged_world. Skipped.', mat_file_name);
        continue
    end

    if isempty(kml_projected_crs) && isfield(S, 'im_info')
        if isobject(S.im_info) && isprop(S.im_info, 'ProjectedCRS')
            kml_projected_crs = S.im_info.ProjectedCRS;
        elseif isstruct(S.im_info) && isfield(S.im_info, 'ProjectedCRS')
            kml_projected_crs = S.im_info.ProjectedCRS;
        end
    end

    merged_world = S.merged_world;
    if isempty(merged_world)
        continue
    end

    for i = 1:height(merged_world)
        seg_x = merged_world.ob_seg_plx{i};
        seg_y = merged_world.ob_seg_ply{i};

        seg_x = seg_x(:);
        seg_y = seg_y(:);
        valid_xy = isfinite(seg_x) & isfinite(seg_y);
        seg_x = seg_x(valid_xy);
        seg_y = seg_y(valid_xy);

        if numel(seg_x) < 3
            continue
        end

        ct_x = merged_world.ob_ct_plx{i};
        ct_y = merged_world.ob_ct_ply{i};
        if iscell(ct_x), ct_x = ct_x{:}; end
        if iscell(ct_y), ct_y = ct_y{:}; end
        ct_x = double(ct_x(1));
        ct_y = double(ct_y(1));

        ob_id = merged_world.ob_id{i};
        if iscell(ob_id), ob_id = ob_id{:}; end
        ob_id = string(ob_id);

        new_row = table( ...
            string(mat_file_name), ob_id, {seg_x}, {seg_y}, ct_x, ct_y, ...
            'VariableNames', {'mat_file', 'ob_id', 'seg_x', 'seg_y', 'centroid_x', 'centroid_y'});
        boulder_table = [boulder_table; new_row]; %#ok<AGROW>
    end
end

if isempty(boulder_table)
    error('No valid boulder polygons were loaded from MAT files.');
end

%% Match points to polygons
n_points = numel(utm_x);

match_status = strings(n_points, 1);
matched_mat_file = strings(n_points, 1);
matched_ob_id = strings(n_points, 1);
matched_centroid_x = nan(n_points, 1);
matched_centroid_y = nan(n_points, 1);
distance_to_centroid = nan(n_points, 1);
num_candidates = zeros(n_points, 1);

for p = 1:n_points
    if ~valid_point(p)
        match_status(p) = "invalid_coordinate";
        continue
    end

    candidate_idx = [];
    candidate_dist = [];

    for b = 1:height(boulder_table)
        bx = boulder_table.seg_x{b};
        by = boulder_table.seg_y{b};

        [in, on] = inpolygon(utm_x(p), utm_y(p), bx, by);
        if include_boundary
            is_match = in || on;
        else
            is_match = in;
        end

        if is_match
            candidate_idx(end+1, 1) = b; %#ok<SAGROW>
            dx = utm_x(p) - boulder_table.centroid_x(b);
            dy = utm_y(p) - boulder_table.centroid_y(b);
            candidate_dist(end+1, 1) = hypot(dx, dy); %#ok<SAGROW>
        end
    end

    num_candidates(p) = numel(candidate_idx);

    if isempty(candidate_idx)
        match_status(p) = "unmatched";
        continue
    end

    if numel(candidate_idx) == 1
        selected_idx = candidate_idx(1);
        match_status(p) = "matched";
    else
        switch multiple_match_rule
            case 'nearest_centroid'
                [~, min_idx] = min(candidate_dist);
                selected_idx = candidate_idx(min_idx);
                match_status(p) = "multiple_candidates";
            otherwise
                error('Unsupported multiple_match_rule: %s', multiple_match_rule);
        end
    end

    matched_mat_file(p) = boulder_table.mat_file(selected_idx);
    matched_ob_id(p) = boulder_table.ob_id(selected_idx);
    matched_centroid_x(p) = boulder_table.centroid_x(selected_idx);
    matched_centroid_y(p) = boulder_table.centroid_y(selected_idx);
    distance_to_centroid(p) = hypot(utm_x(p) - matched_centroid_x(p), utm_y(p) - matched_centroid_y(p));
end

%% Save CSV
result_table = table( ...
    input_id, gps_id, las_file_name, utm_x, utm_y, ...
    match_status, matched_mat_file, matched_ob_id, ...
    matched_centroid_x, matched_centroid_y, distance_to_centroid, num_candidates, ...
    'VariableNames', { ...
        'ID', 'GPS', 'Las_file_name', 'X_231116ver', 'Y_231116ver', ...
        'match_status', 'matched_mat_file', 'matched_ob_id', ...
        'matched_centroid_x', 'matched_centroid_y', 'distance_to_centroid', 'num_candidates'});

writetable(result_table, output_csv);

%% Save KML
if write_kml
    if isempty(kml_projected_crs)
        warning('KML was not written because im_info.ProjectedCRS is not available.');
    else
        point_lat = nan(size(utm_x));
        point_lon = nan(size(utm_y));
        [point_lat(valid_point), point_lon(valid_point)] = ...
            projinv(kml_projected_crs, utm_x(valid_point), utm_y(valid_point));
        matched_valid = isfinite(matched_centroid_x) & isfinite(matched_centroid_y);
        centroid_lat = nan(size(matched_centroid_x));
        centroid_lon = nan(size(matched_centroid_y));
        [centroid_lat(matched_valid), centroid_lon(matched_valid)] = ...
            projinv(kml_projected_crs, matched_centroid_x(matched_valid), matched_centroid_y(matched_valid));

        fid = fopen(output_kml, 'w');
        if fid < 0
            error('Failed to open KML output: %s', output_kml);
        end

        fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
        fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n');
        fprintf(fid, '<name>matched_boulder_positions</name>\n');

        fprintf(fid, '<Folder><name>input_points</name>\n');
        for p = 1:n_points
            if ~valid_point(p)
                continue
            end
            fprintf(fid, '<Placemark>\n');
            fprintf(fid, '<name>input_%s</name>\n', char(string(input_id(p))));
            fprintf(fid, '<ExtendedData>\n');
            fprintf(fid, '<Data name="gps_id"><value>%s</value></Data>\n', char(string(gps_id(p))));
            fprintf(fid, '<Data name="match_status"><value>%s</value></Data>\n', char(match_status(p)));
            fprintf(fid, '<Data name="matched_mat_file"><value>%s</value></Data>\n', char(matched_mat_file(p)));
            fprintf(fid, '<Data name="matched_ob_id"><value>%s</value></Data>\n', char(matched_ob_id(p)));
            fprintf(fid, '<Data name="num_candidates"><value>%d</value></Data>\n', num_candidates(p));
            fprintf(fid, '</ExtendedData>\n');
            fprintf(fid, '<Point><coordinates>%.12f,%.12f,0</coordinates></Point>\n', point_lon(p), point_lat(p));
            fprintf(fid, '</Placemark>\n');
        end
        fprintf(fid, '</Folder>\n');

        fprintf(fid, '<Folder><name>matched_centroids</name>\n');
        for p = 1:n_points
            if ~matched_valid(p)
                continue
            end
            fprintf(fid, '<Placemark>\n');
            fprintf(fid, '<name>%s</name>\n', char(matched_ob_id(p)));
            fprintf(fid, '<ExtendedData>\n');
            fprintf(fid, '<Data name="input_id"><value>%s</value></Data>\n', char(string(input_id(p))));
            fprintf(fid, '<Data name="gps_id"><value>%s</value></Data>\n', char(string(gps_id(p))));
            fprintf(fid, '<Data name="matched_mat_file"><value>%s</value></Data>\n', char(matched_mat_file(p)));
            fprintf(fid, '<Data name="distance_to_centroid"><value>%.6f</value></Data>\n', distance_to_centroid(p));
            fprintf(fid, '</ExtendedData>\n');
            fprintf(fid, '<Point><coordinates>%.12f,%.12f,0</coordinates></Point>\n', centroid_lon(p), centroid_lat(p));
            fprintf(fid, '</Placemark>\n');
        end
        fprintf(fid, '</Folder>\n');
        fprintf(fid, '</Document>\n</kml>\n');
        fclose(fid);
    end
end

%% Summary
disp("Position matching finished.");
disp(strcat("Loaded boulders: ", string(height(boulder_table))));
disp(strcat("Input points: ", string(n_points)));
disp(strcat("Matched points: ", string(sum(match_status ~= "unmatched"))));
disp(strcat("Unmatched points: ", string(sum(match_status == "unmatched"))));
disp(strcat("CSV: ", output_csv));
if write_kml
    disp(strcat("KML: ", output_kml));
end
