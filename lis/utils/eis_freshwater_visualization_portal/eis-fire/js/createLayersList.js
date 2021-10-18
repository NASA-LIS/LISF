var http = require('http'),
    concat = require('concat-stream'),
    parseString = require('xml2js').parseString,
    fs = require('fs');
    
var maskIDs = ['MODIS_Terra_Data_No_Data', 'MODIS_Aqua_Data_No_Data'];

http.get('http://map1.vis.earthdata.nasa.gov/wmts-webmerc/1.0.0/WMTSCapabilities.xml', function(response) {
    response.pipe(concat(function(buffer) {
        parseString(buffer.toString(), function(err, res) {
            var layerHash = {};
            var maskHash = {};
            var xmlLayers = res.Capabilities.Contents[0].Layer;
            xmlLayers.forEach(function(layer) {
                var tileMatrixSet = layer.TileMatrixSetLink[0].TileMatrixSet[0],
                    zoom = parseInt(tileMatrixSet.match(/Level(\d+)/)[1]),
                    title = layer['ows:Identifier'][0];
                var info = {
                    title: title,
                    template: layer.ResourceURL[0].$.template,
                    zoom: zoom,
                    date: !!layer.Dimension
                };
                
                if (maskIDs.indexOf(title) !== -1) {
                    maskHash[title] = info;
                } else {
                    layerHash[title] = info;
                }
            });
            
            fs.writeFileSync('GIBSMetadata.js', 
                'L.GIBS_LAYERS = ' + JSON.stringify(layerHash) + 
                ';L.GIBS_MASKS = ' + JSON.stringify(maskHash) + ';'
            );
            
            console.log('Layers: ', xmlLayers.length - 2);
            console.log('Masks: ', 2);
        });
    }));
});