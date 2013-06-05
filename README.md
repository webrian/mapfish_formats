mapfish_formats
===============

mapfish_formats extends the MapFish protocol with custom output formats.

How it works
------------
Class `FormatsProtocol` in module `mapfish_formats.protocol` subclasses MapFish's `Protocol` class and adds a new
parameter `format` to set the output format.

Currently available:
* `format=geojson` returns GeoJSON objects (MapFish default output format)
* `format=ext` returns JSON objects that can be consumed in [Ext JS JsonStores](http://docs.sencha.com/extjs/3.4.0/#!/api/Ext.data.JsonStore)
* `format=xls` returns an Excel spreadsheet with a data sheet and metadata sheet
* `format=shp` returns a compressed Shapefile
* `format=hist` returns a histogram and barplot respectively as PNG image
  * Additional parameters `width` and `height` to set the image width and height

Prerequisites
-------------
* [MapFish](http://www.mapfish.org) 2.2
* [xlwt](http://www.python-excel.org/) for spreadsheet support
* [Python Shapefile Library](https://code.google.com/p/pyshp/) for Shapefile support
* [R](http://www.r-project.org/) and [RPy2](http://rpy.sourceforge.net/rpy2.html) for histograms and barplots
* [Matplotlib](http://matplotlib.org/) for histograms and barplots as an alternative renderer to the R renderer
