//New layer class for Leaflet-GeoMixer plugin (https://github.com/ScanEx/Leaflet-GeoMixer)
//Meta properties:
//  * gibs-layername Imagery Layer ID from https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+Available+Imagery+Products
//  * gibs-transparent If present, data mask layer will be used to hide black no-data imagery background
(function() {
    
var GIBSProxyLayer = function() {};

GIBSProxyLayer.prototype.initFromDescription = function(layerDescription) {
    var props = layerDescription.properties;
    
    if (!props.MetaProperties['gibs-layername']) {
         return new L.gmx.DummyLayer(props);
    }
    
    var layerName = props.MetaProperties['gibs-layername'].Value,
        isTransparent = !!props.MetaProperties['gibs-transparent'];

    try {
        var layer = new L.GIBSLayer(layerName, {transparent: isTransparent});
    } catch(e) {
        return new L.gmx.DummyLayer(props);
    }

    layer.getGmxProperties = function() {
        return props;
    }

    layer.setDateInterval = function(dateBegin, dateEnd) {

        //if dateEnd is exactly midnight, we show previous day iff dateBegin is not equal to dateEnd (next day otherwise)
        var shift = 0;
        if (dateEnd.valueOf() % (3600*1000*24) === 0) {
            shift = dateBegin < dateEnd ? 1 : 0;
        }

        this.setDate(new Date(dateEnd - shift));
        return this;
    }

    return layer;
}

if (L.gmx && L.gmx.addLayerClass) {
    L.gmx.addLayerClass('GIBS', GIBSProxyLayer);
}

window.gmxCore && gmxCore.addModule('GIBSVirtualLayer', {layerClass: GIBSProxyLayer}, {
    init: function(module, path) {
        if (!L.GIBSLayer) {
            return $.when(
                gmxCore.loadScript(path + 'GIBSLayer.js'),
                gmxCore.loadScript(path + 'GIBSMetadata.js')
            )
        }
    }
});

})();