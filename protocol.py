#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__ = "Adrian Weber, Centre for Development and Environment, University of Bern"
__date__ = "$Apr 29, 2013 6:55:21 AM$"

from geoalchemy import functions
import geojson
from mapfish.protocol import *
import rpy2.rinterface as rinterface
import rpy2.robjects as robjects
import shapefile
from shapely.wkb import loads
import simplejson as json
from sqlalchemy import func
from tempfile import NamedTemporaryFile
try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO
from zipfile import ZIP_DEFLATED
from zipfile import ZipFile


# Map of EPSG codes to write the .prj files
# taken from spatialreference.org
epsg_code = {
4326: 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]',
32648: 'PROJCS["WGS_1984_UTM_Zone_48N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["Meter",1]]'
}

class FormatsProtocol(Protocol):

    def read(self, request, filter=None, id=None, format='geojson', ** kwargs):
        """
        Build a query based on the filter or the idenfier, send the query
        to the database, and return a Feature or a FeatureCollection.
        """

        if format == 'geojson':

            ret = None
            if id is not None:
                o = self.Session.query(self.mapped_class).get(id)
                if o is None:
                    abort(404)
                ret = self._filter_attrs(o.toFeature(), request)
            else:
                objs = self._query(request, filter)
                ret = FeatureCollection(
                                        [self._filter_attrs(o.toFeature(), request) \
                                        for o in objs])

            return geojson.dumps(ret)

        if format == 'ext':
            if id is not None:
                query = self.Session.query(self.mapped_class).get(id)
            else:
                query = self._query(request, filter, False)
            return self._read_ext(request, query, filter=filter, name_mapping=kwargs['name_mapping'])

        if format == 'hist':
            if filter is None:
                filter = create_default_filter(request, self.mapped_class)
            return self._plot_histogram(request, self.Session.query(self.mapped_class).filter(filter), kwargs['categories'])

        if format == 'shp':

            epsg = kwargs.get("epsg", 4326)

            if filter is None:
                filter = create_default_filter(request, self.mapped_class)

            mapped_attributes = []
            mapped_attributes.append(functions.wkb(functions.transform(self.mapped_class.geometry_column(), epsg)).label("geometry_column"))
            for attr in request.params.get("attrs").split(","):
                mapped_attributes.append(getattr(self.mapped_class, attr))

            return self._read_shp(request, self.Session.query(* mapped_attributes).filter(filter), epsg=epsg)

    def _read_ext(self, request, query, filter=None, name_mapping=None):
        """
        A format suitable for Ext json stores.
        """

        output = {}
        output['totalResults'] = self.count(request, filter)

        attrs = request.params['attrs'].split(',')

        output['metaData'] = {'totalProperty': 'totalResults', 'root': 'rows'}
        output['metaData']['fields'] = []

        output['rows'] = []

        for i in query.all():
            row = {}
            for k in attrs:
                value = getattr(i, k)
                if name_mapping is not None and k in name_mapping:
                    k = name_mapping[k]

                row[k] = value


                # Check if metadata of current attribute is already declared
                is_declared = False
                for f in output['metaData']['fields']:
                    try:
                        if f['name'] == k:
                            is_declared = True
                    except TypeError:
                        pass

                if not is_declared:
                    metadataField = {}
                    metadataField['name'] = k
                    if isinstance(value, (float)):
                        metadataField['type'] = 'float'
                    elif isinstance(value, (int)):
                        metadataField['type'] = 'int'
                    else:
                        metadataField['type'] = 'string'
                    output['metaData']['fields'].append(metadataField)

            # Append the current row to the rows
            output['rows'].append(row)

        return json.dumps(output)


    def _plot_histogram(self, request, query, categories=None):

        # Get the requested attribute
        attr = request.params.get('attrs').split(",")[0]

        defaultSide = 480

        try:
            height = int(request.params.get("height", defaultSide))
        except ValueError:
            height = defaultSide
        try:
            width = int(request.params.get("width", defaultSide))
        except ValueError:
            width = defaultSide

        rbreaks = 10
        if 'breaks' in request.params:
            rbreaks = int(request.params['breaks'])

        rinterface.initr()

        r = robjects.r
        r.library('grDevices')

        mappedAttribute = getattr(self.mapped_class, attr)

        # Create a temporary file
        file = NamedTemporaryFile()
        r.png(file.name, width=width, height=height)
        r.par(bg="#F0F0F0", mar=robjects.FloatVector([2.6, 4.1, 3.1, 1.1]))

        bar_color = "#BEBEBE"

        # If categories is not none, then the current attribute has categories and
        # we want to draw a barplot instead of a histogram
        if categories is not None:

            names = []
            v = []
            #for i in self.Session.query(func.count(mappedAttribute)).filter(mappedAttribute.in_(keys)).group_by(mappedAttribute):
            for a, count in query.from_self(mappedAttribute, func.count(mappedAttribute)).filter(mappedAttribute.in_(categories.keys())).group_by(mappedAttribute):
                v.append(int(count))
                names.append(categories[unicode(a)].encode('UTF-8'))
            x = robjects.IntVector(v)

            x.names = robjects.StrVector(names)
            r.barplot(x, col=bar_color, xlab=str(), ylab=str(), main=str()) #, ** {"names.arg": robjects.StrVector(names)})
            #r.par(bg="#F0F0F0", mar=robjects.FloatVector([1.5, 1.5, 1.5, 1.5]))
            #r.pie(x, labels=robjects.StrVector(names), clockwise=True)

        # Handle quantitative data
        else:

            v = []
            for i in query.all():
                v.append(getattr(i, attr))

            x = robjects.FloatVector(v)

            r.hist(x, col=bar_color, breaks=rbreaks, ylab=str(), xlab=str(), main=str())

        # Finish drawing
        r('dev.off()')

        f = open(file.name, 'r')

        return f

    def _read_shp(self, request, query, **kwargs):

        requested_attrs = request.params.get("attrs").split(",")

        w = shapefile.Writer(shapefile.POLYGON)
        w.autoBalance = 1

        # Get the first feature to guess the datatype
        first_record = query.first()

        # Loop all requested attributes
        for attr in requested_attrs:

            first_value = getattr(first_record, attr)

            # Guess the datatype
            if isinstance(first_value, int):
                w.field(str(attr), 'N', 40)
            elif isinstance(first_value, float):
                w.field(str(attr), 'N', 40, 10)
            else:
                w.field(str(attr), 'C', 40)

        # Now query all features
        for i in query.all():

            # Create geometry from self.mapped_class
            #g = loads(str(getattr(i, 'geometry_column').geom_wkb))
            # Create geometry from AsBinary query
            g = loads(str(getattr(i, 'geometry_column')))

            if g.geom_type == "Polygon":

                ring_list = []

                point_list = []

                for j in g.exterior.coords:
                    point_list.append([j[0], j[1]])

                ring_list.append(point_list)

                for interior in g.interiors:

                    point_list = []

                    for k in interior.coords:
                        point_list.append([k[0], k[1]])

                    ring_list.append(point_list)

                w.poly(shapeType=shapefile.POLYGON, parts=ring_list)

                values = []
                for v in requested_attrs:
                    try:
                        values.append(str(getattr(i, v)))
                    except UnicodeEncodeError:
                        values.append(str(getattr(i, v).encode("UTF-8")))
                
                w.record(* values)

        # Create the required files and fill them
        shp = StringIO()
        shx = StringIO()
        dbf = StringIO()
        cpg = StringIO()
        prj = StringIO()
        w.saveShp(shp)
        w.saveShx(shx)
        w.saveDbf(dbf)
        cpg.write("UTF-8")
        prj.write(epsg_code[kwargs.get("epsg", 4326)])

        # Create a memory file-like deflated zip file
        s = StringIO()
        f = ZipFile(s, 'w', ZIP_DEFLATED)
        f.writestr("data.shp", shp.getvalue())
        f.writestr("data.dbf", dbf.getvalue())
        f.writestr("data.shx", shx.getvalue())
        f.writestr("data.cpg", cpg.getvalue())
        f.writestr("data.prj", prj.getvalue())

        # Close the zip file
        f.close()

        # And return the content
        return s