# Water surface temperature - Landsat
=================

<img src="https://github.com/leolaipelt/surfacewatertemperature/blob/main/images/fig_capa.png" width="360">

## Estimating Water surface temperature using Landsat in Google Earth Engine platform.

* The algorithm is developed to estimate water surface temperature for Landsat 5, 7, 8 and 9 based on the Collection 2 product.

## Input

* The following collections are compatible:

| Image Collections IDs  |
| :--------------------: |
| LANDSAT/LC09/C02/T1_L2 |
| LANDSAT/LC08/C02/T1_L2 |
| LANDSAT/LE07/C02/T1_L2 |
| LANDSAT/LT05/C02/T1_L2 |

## Functions

### Image.from_image_id()

* Compute water surface temperature for a single image.
* Allow to obtain water surface temperature with image collections by mapping over Landsat collections.

### Example

	from wst import wst

	ls_id = 'LANDSAT/LC08/C02/T1_L2/LC08_001062_20230906'
	
    lst = wst.Image.from_image_id(ls_id)

## Examples Notebooks

Examples of how to use the algorithm is detailed in *examples* folder:

[wst example.](https://github.com/leolaipelt/surfacewatertemperature/tree/main/examples "Example")

## Installation

	pip install wst

### Depedencies

 * `earthengine-api` <https://github.com/google/earthengine-api>`

## References

