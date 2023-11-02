function [] = saveLineKML(save_path, dataTable)
%dataTable: name, latitude, longitude, elevation, Color(string)
docNode = com.mathworks.xml.XMLUtils.createDocument('kml');
kml =  docNode.getDocumentElement;
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
                lat  = data.latitude{:};
                lon  = data.longitude{:};
                elv  = data.elevation{:};

                Placemark = docNode.createElement('Placemark');
                Folder.appendChild(Placemark);

                    Snippet = docNode.createElement('Snippet');
                    Snippet.setAttribute('maxLines','0');
                    Placemark.appendChild(Snippet);

                    description = docNode.createElement('description');
                    Placemark.appendChild(description);

                    ObjName = docNode.createElement('name');
                    ObjName.appendChild(docNode.createTextNode(data.name));
                    Placemark.appendChild(ObjName);

                    style = docNode.createElement('style');
                    Placemark.appendChild(style);
                        
                        LineStyle = docNode.createElement('LineStyle');
                        style.appendChild(LineStyle);

                            color = docNode.createElement('color');
                            color.appendChild(docNode.createTextNode('ffffffff'));
                            LineStyle.appendChild(color);
            
                    
                    LineString = docNode.createElement('LineString');
                    Placemark.appendChild(LineString);
                        
                        segData = "";
                        for n = 1:size(lat,2)
                            segData = strcat(segData, num2str(lon(n),'%.12f'),',', num2str(lat(n),'%.12f'), ',', num2str(elv(n),'%.2f')," ");
                        end
                        segData = strcat(segData, num2str(lon(1),'%.12f'),',',num2str(lat(1),'%.12f'), ',', num2str(elv(1),'%.2f'));

                        coordinates = docNode.createElement('coordinates');
                        coordinates.appendChild(docNode.createTextNode(segData));
                        LineString.appendChild(coordinates);
            end

xmlwrite(save_path, docNode);

end