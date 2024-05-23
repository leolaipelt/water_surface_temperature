import ee
 
def cloud_mask_sr_l457(landsat_image):
    """Cloud mask (Landsat 4/5/7)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 66  (01000010)
    water value = 68   (01000100)

    References
    ----------

    """
    quality = landsat_image.select('pixel_qa')
    c01 = quality.eq(66)  # Clear (01000010)
    c02 = quality.eq(68)  # Water (01000100)
    mask = c01.Or(c02)

    return mask


def cloud_mask_sr_l8(landsat_image):
    """Cloud mask (Landsat 8)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 322  (00101000010)
    water value = 324   (00101000100)

    References
    ----------

    """
    quality = landsat_image.select('pixel_qa')
    c01 = quality.eq(322)
    c02 = quality.eq(324)
    c03 = quality.eq(1346)  # (10101000010)
    mask = c01.Or(c02).Or(c03)

    return mask


def cloud_mask_C2_l457(landsat_image):
    """Cloud mask (Landsat 4/5/7)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 5440  (0001010101000000)
    water value = 5504   (0001010110000000)

    References
    ----------

    """
    quality = landsat_image.select('QA_PIXEL')
    c01 = quality.eq(5440)
    c02 = quality.eq(5504)
    mask = c01.Or(c02)

    return mask


def cloud_mask_C2_l89(landsat_image):
    """Cloud mask (Landsat 8/9)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 21824  (0001010101000000)
    water value = 21952   (0001010110000000)

    References
    ----------

    """
    quality = landsat_image.select('QA_PIXEL')
    c01 = quality.eq(21824)
    c02 = quality.eq(21952)
    mask = c01.Or(c02)

    return mask

def calc_ndwi(landsat_image):
    """Normalized difference water index

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    return ee.Image(landsat_image).normalizedDifference(['green', 'nir'])\
        .rename('ndwi').unmask(0)

def prep_image(sr_image):
        
    sr_image = ee.Image(sr_image)

    # Use the SPACECRAFT_ID property identify each Landsat type
    spacecraft_id = ee.String(sr_image.get('SPACECRAFT_ID'))

    # Rename bands to generic names
    input_bands = ee.Dictionary({
        'LANDSAT_4': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4',
                        'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
        'LANDSAT_5': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4',
                        'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
        'LANDSAT_7': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4',
                        'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
        'LANDSAT_8': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5',
                        'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
        'LANDSAT_9': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5',
                        'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
    })
    output_bands = ee.Dictionary({
        'LANDSAT_4': ['blue', 'green', 'red', 'nir',
                        'swir1', 'swir2', 'lst', 'QA_PIXEL'],
        'LANDSAT_5': ['blue', 'green', 'red', 'nir',
                        'swir1', 'swir2', 'lst', 'QA_PIXEL'],
        'LANDSAT_7': ['blue', 'green', 'red', 'nir',
                        'swir1', 'swir2', 'lst', 'QA_PIXEL'],
        'LANDSAT_8': ['ultra_blue', 'blue', 'green', 'red', 'nir',
                        'swir1', 'swir2', 'lst', 'QA_PIXEL'],
        'LANDSAT_9': ['ultra_blue', 'blue', 'green', 'red', 'nir',
                        'swir1', 'swir2', 'lst', 'QA_PIXEL'],
    })
    scalars = ee.Dictionary({
        'LANDSAT_4': [0.0000275, 0.0000275, 0.0000275, 0.0000275,
                        0.0000275, 0.0000275, 0.00341802, 1],
        'LANDSAT_5': [0.0000275, 0.0000275, 0.0000275, 0.0000275,
                        0.0000275, 0.0000275, 0.00341802, 1],
        'LANDSAT_7': [0.0000275, 0.0000275, 0.0000275, 0.0000275,
                        0.0000275, 0.0000275, 0.00341802, 1],
        'LANDSAT_8': [0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275,
                        0.0000275, 0.0000275, 0.00341802, 1],
        'LANDSAT_9': [0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275,
                        0.0000275, 0.0000275, 0.00341802, 1],
    })
    offsets = ee.Dictionary({
        'LANDSAT_4': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0],
        'LANDSAT_5': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0],
        'LANDSAT_7': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0],
        'LANDSAT_8': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0],
        'LANDSAT_9': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0],
    })

    result = sr_image\
        .select(input_bands.get(spacecraft_id), output_bands.get(spacecraft_id))\
        .multiply(ee.Image.constant(ee.List(scalars.get(spacecraft_id))))\
        .add(ee.Image.constant(ee.List(offsets.get(spacecraft_id))))
    
    cloud_mask = ee.Algorithms.If(
        ee.List(['LANDSAT_8', 'LANDSAT_9']).contains(spacecraft_id),
        cloud_mask_C2_l89(sr_image),
        cloud_mask_C2_l457(sr_image),
    )

    return result

def prep_collection(start_date,end_date,point,cloud_cover):
  
  bad_data_filter = ee.Filter.date(ee.Date.fromYMD(2023,10,29), ee.Date.fromYMD(2023,10,30))

  l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filter(bad_data_filter.Not())
  l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  
  #var col_merge = l5.merge(l8).merge(l9);
  col_merge = l5.merge(l7).merge(l8).merge(l9)
  
  input_col = col_merge.filterDate(start_date, end_date)\
                        .filterBounds(point)\
                        .filterMetadata('CLOUD_COVER_LAND', 'less_than',cloud_cover)
                                    
  return input_col


def calc_ndvi(landsat_image):
    """Normalized difference vegetation index

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    return ee.Image(landsat_image).normalizedDifference(['nir', 'red'])\
        .rename(['NDVI']) #.unmask(0)


def func_emissivity(ndvi):
   
   
    pv = ndvi.expression('((ndvi - 0.2) / 0.3) ** 2', {'ndvi': ndvi})
  
    #Assuming typical Soil Emissivity of 0.97 and Veg Emissivity of 0.99
    #   and shape Factor mean value of 0.553
    de = pv.expression('(1 - 0.97) * (1 - Pv) * (0.55 * 0.99)', {'Pv': pv})
    RangeEmiss = de.expression('(0.99 * Pv) + (0.97 * (1 - Pv)) + dE', {'Pv': pv, 'dE': de})

    return ndvi.where(ndvi.lt(0), 0.985)\
        .where(ndvi.gte(0).And(ndvi.lt(0.2)), 0.977)\
        .where(ndvi.gt(0.5), 0.99)\
        .where(ndvi.gte(0.2).And(ndvi.lte(0.5)), RangeEmiss)\
        .clamp(0.977, 0.99)\
        .rename(['emissivity'])
    

def calc_wst(img,ndvi):
  
    # Get NDVI
    #ndvi =  calc_ndvi(img).clamp(0,1)

    collection_id = ee.Dictionary({
    'LANDSAT_5':  ee.ImageCollection('LANDSAT/LT05/C02/T1'),
    'LANDSAT_7':  ee.ImageCollection('LANDSAT/LE07/C02/T1_RT'),
    'LANDSAT_8': ee.ImageCollection('LANDSAT/LC08/C02/T1_RT'),
    'LANDSAT_9':  ee.ImageCollection('LANDSAT/LC09/C02/T1'),
    })

    ls_product_id = img.get('LANDSAT_PRODUCT_ID')
    spacecraft_id = img.get('SPACECRAFT_ID')
    system_index = img.get("LANDSAT_SCENE_ID")

    col_rt  = ee.ImageCollection(collection_id.get(spacecraft_id))\
                    .filterMetadata('LANDSAT_SCENE_ID','equals',system_index)\
                    .first()

    coef_dict_multi = ee.Dictionary({

    'LANDSAT_5':'RADIANCE_MULT_BAND_6',
    'LANDSAT_7':'RADIANCE_MULT_BAND_6_VCID_1',
    'LANDSAT_8':'RADIANCE_MULT_BAND_10',
    'LANDSAT_9':'RADIANCE_MULT_BAND_10',
    })

    coef_dict_add = ee.Dictionary({

    'LANDSAT_5':'RADIANCE_ADD_BAND_6',
    'LANDSAT_7':'RADIANCE_ADD_BAND_6_VCID_1',
    'LANDSAT_8':'RADIANCE_ADD_BAND_10',
    'LANDSAT_9':'RADIANCE_ADD_BAND_10',
    })

    coef_k1 = ee.Dictionary({
    'LANDSAT_5':607.76,
    'LANDSAT_7':666.09,
    'LANDSAT_8':774.8853,
    'LANDSAT_9':774.8853,
    })

    coef_k2 = ee.Dictionary({
    'LANDSAT_5':1260.56,
    'LANDSAT_7':1282.71,
    'LANDSAT_8':1321.0789,
    'LANDSAT_9':1321.0789,
    })

    band_lst = ee.Dictionary({
    'LANDSAT_5':'B6',
    'LANDSAT_7':'B6_VCID_1',
    'LANDSAT_8':'B10',
    'LANDSAT_9':'B10',
    })


    fc_Landsat = ndvi.multiply(1.0)\
            .subtract(0.15)\
            .divide(0.65)\
                .clamp(0,1.0)
          
    toa = img.expression(
    'ML*Qcal+AL',{
    'ML':ee.Number(col_rt.get(coef_dict_multi.get(spacecraft_id))),
    'Qcal': col_rt.select([band_lst.get(ee.String(spacecraft_id))]),
    'AL': ee.Number(col_rt.get(coef_dict_add.get(spacecraft_id)))
    }
    )

    bt = ee.Algorithms.Landsat.TOA(col_rt)
    bt = bt.select([band_lst.get(ee.String(spacecraft_id))])
    bt = bt.rename('BT')

    # get maximum NDVI
    max = ee.Number(ndvi.reduceRegion(
    reducer=ee.Reducer.max(),
    geometry=  ndvi.geometry(),
    scale= 30,
    bestEffort= True,
    maxPixels= 1e7
    ).values().get(0))

    # get minimum NDVI
    min = ee.Number(ndvi.reduceRegion(
    reducer= ee.Reducer.min(),
    geometry= ndvi.geometry(),
    scale= 30,
    bestEffort= True,
    maxPixels= 1e7
    ).values().get(0))

    fvc = ndvi.expression(
    '((ndvi - ndvi_min)/(ndvi_max - ndvi))**2',{
    'ndvi': ndvi,
    'ndvi_min': min,
    'ndvi_max':  max,
    })

    emiss = fvc.expression(
    '0.004*fvc+0.986',{
    'fvc': fvc
    }
    )

    '''lst = fvc.expression(
    '(bt / 1) + 10.895 * (bt / 14380) * log(emiss)',{
    'bt':bt,
    'emiss':emiss
    }
    ).set('system:time_start',img.get('system:time_start'))

    lst = fvc.expression(
    '(bt / (1 + (0.00115 * bt / 1.4388) * log(emiss))) ',{
    'bt':bt,
    'emiss':emiss
    }
    ).set('system:time_start',img.get('system:time_start'));
    '''  
    # narrow band transmissivity of air
    tnb = 0.866 

    # path radiance 
    rp = 0.91  

    # narrow band clear sky downward thermal radiation   
    rsky = 1.32   

    k1 = ee.Number(coef_k1.get(spacecraft_id))
    k2 = ee.Number(coef_k2.get(spacecraft_id))
    
    emissivity_img = func_emissivity(ndvi)

    thermal_rad_toa = bt.expression(
    'k1 / (exp(k2 / bt) - 1)',
    {'bt': bt, 'k1': k1, 'k2': k2}
    )

    rc = thermal_rad_toa.expression(
    '((thermal_rad_toa - rp) / tnb) - ((1 - emiss) * rsky)',
    {
        'thermal_rad_toa': thermal_rad_toa,
        'emiss': emissivity_img,
        'rp': 0.91, 'tnb': 0.866, 'rsky': 1.32,
    }
    )
    wst = rc.expression(
    'k2 / log(emiss * k1 / rc + 1)',
    {'emiss': emissivity_img, 'rc': rc, 'k1': k1, 'k2': k2}
    )

    return wst.rename('LST_C2')
  
  
  
