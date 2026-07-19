function [] = saveLineKML(save_path, dataTable)
% save_path: Output KML file path.
% dataTable: Table with name, latitude, longitude, elevation,
%            and optional attribute columns.

requiredNames = {'name','latitude','longitude','elevation'};
tableNames = dataTable.Properties.VariableNames;

missingNames = setdiff(requiredNames, tableNames, 'stable');
if ~isempty(missingNames)
    error('Missing required columns: %s', strjoin(missingNames, ', '));
end

attributeNames = setdiff(tableNames, requiredNames, 'stable');

docNode = com.mathworks.xml.XMLUtils.createDocument('kml');
kml = docNode.getDocumentElement;
kml.setAttribute('xmlns','http://www.opengis.net/kml/2.2');

Document = docNode.createElement('Document');
kml.appendChild(Document);

DocName = docNode.createElement('name');
DocName.appendChild(docNode.createTextNode('object'));
Document.appendChild(DocName);

Folder = docNode.createElement('Folder');
Document.appendChild(Folder);

DirName = docNode.createElement('name');
DirName.appendChild(docNode.createTextNode('shape lines'));
Folder.appendChild(DirName);

for i = 1:size(dataTable,1)

    data = dataTable(i,:);

    lat = data.latitude;
    lon = data.longitude;
    elv = data.elevation;

    if iscell(lat)
        lat = lat{1};
    end
    if iscell(lon)
        lon = lon{1};
    end
    if iscell(elv)
        elv = elv{1};
    end

    if numel(lat) ~= numel(lon) || numel(lat) ~= numel(elv)
        error('Coordinate sizes do not match at row %d.', i);
    end

    Placemark = docNode.createElement('Placemark');
    Folder.appendChild(Placemark);

    Snippet = docNode.createElement('Snippet');
    Snippet.setAttribute('maxLines','0');
    Placemark.appendChild(Snippet);

    description = docNode.createElement('description');
    Placemark.appendChild(description);

    nameValue = data.name;
    if iscell(nameValue)
        nameValue = nameValue{1};
    end

    ObjName = docNode.createElement('name');
    ObjName.appendChild(docNode.createTextNode(char(string(nameValue))));
    Placemark.appendChild(ObjName);

    if ~isempty(attributeNames)

        ExtendedData = docNode.createElement('ExtendedData');
        Placemark.appendChild(ExtendedData);

        for j = 1:numel(attributeNames)

            attributeName = attributeNames{j};
            attributeValue = data.(attributeName);

            if iscell(attributeValue)
                attributeValue = attributeValue{1};
            end

            if isempty(attributeValue)
                attributeText = '';
            elseif isnumeric(attributeValue) || islogical(attributeValue)
                attributeText = char(strjoin(string(attributeValue(:)), ','));
            else
                attributeText = char(strjoin(string(attributeValue(:)), ','));
            end

            DataNode = docNode.createElement('Data');
            DataNode.setAttribute('name', attributeName);
            ExtendedData.appendChild(DataNode);

            ValueNode = docNode.createElement('value');
            ValueNode.appendChild(docNode.createTextNode(attributeText));
            DataNode.appendChild(ValueNode);
        end
    end

    style = docNode.createElement('Style');
    Placemark.appendChild(style);

    LineStyle = docNode.createElement('LineStyle');
    style.appendChild(LineStyle);

    color = docNode.createElement('color');
    color.appendChild(docNode.createTextNode('ffffffff'));
    LineStyle.appendChild(color);

    LineString = docNode.createElement('LineString');
    Placemark.appendChild(LineString);

    segData = "";

    for n = 1:numel(lat)
        segData = segData + sprintf( ...
            '%.12f,%.12f,%.2f ', ...
            lon(n), lat(n), elv(n));
    end

    segData = segData + sprintf( ...
        '%.12f,%.12f,%.2f', ...
        lon(1), lat(1), elv(1));

    coordinates = docNode.createElement('coordinates');
    coordinates.appendChild(docNode.createTextNode(char(segData)));
    LineString.appendChild(coordinates);
end

xmlwrite(save_path, docNode);

end