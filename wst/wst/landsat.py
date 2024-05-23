import ee
from wst.wst import model

class Image():
    '''Get the water surface temperature for an image'''


    def __init__(self,image):

        self.image = image
       

    def from_image_id(image_id):

        # prep image
        image_sr = ee.Image(model.prep_image(image_id))

        # properties
        _id = ee.Image(image_id).get('system:id')
        _index = ee.Image(image_id).get('system:index')
        _time_start = ee.Image(image_id).get('system:time_start')
        _properties = {
            'system:index': _index,
            'system:time_start':_time_start,
            'image_id':_id
        }

        # calculating ndvi
        ndvi = model.calc_ndvi(image_sr)

        # calculating ndwi
        ndwi = model.calc_ndwi(image_sr)

        # calculating wst
        wst_image = model.calc_wst(ee.Image(image_id),ndvi)

        return  wst_image.updateMask(ndwi.gte(0)).rename('wst').set(_properties)