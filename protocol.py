#
# mapnik_formats
# Copyright (C) 2013 Centre for Development and Environment, University of Bern
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = "Adrian Weber, Centre for Development and Environment, University of Bern"
__date__ = "$Apr 29, 2013 6:55:21 AM$"

from geoalchemy import functions
import geojson
from mapfish.protocol import *
import shapefile
from shapely.wkb import loads
import simplejson as json
from sqlalchemy import func
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
from zipfile import ZIP_DEFLATED
from zipfile import ZipFile
import xlwt
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


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

            if id is not None:
                query = self.Session.query(self.mapped_class).get(id)
            else:
                query = self._query(request, filter, False)
            return self._read_geojson(request, query, filter=filter, name_mapping=kwargs.get("name_mapping", {}))

        if format == 'ext':
            if id is not None:
                query = self.Session.query(self.mapped_class).get(id)
            else:
                query = self._query(request, filter, False)
            return self._read_ext(request, query, filter=filter, name_mapping=kwargs['name_mapping'])

        if format == 'hist':
            if filter is None:
                filter = create_default_filter(request, self.mapped_class)

            renderer = self._plot_histogram_matplotlib
            if "renderer" in request.params and request.params.get("renderer") == "r":
                renderer = self._plot_histogram_r
            return renderer(request, self.Session.query(self.mapped_class).filter(filter), ** kwargs)

        if format == 'xls':

            query = self._query(request, filter, False)
            return self._read_xls(request, query, filter=filter, ** kwargs)

        if format == 'shp':

            epsg = kwargs.get("epsg", 4326)

            metadata = kwargs.get("metadata", None)

            if filter is None:
                filter = create_default_filter(request, self.mapped_class)

            mapped_attributes = []
            mapped_attributes.append(functions.wkb(functions.transform(self.mapped_class.geometry_column(), epsg)).label("geometry_column"))
            for attr in request.params.get("attrs").split(","):
                mapped_attributes.append(getattr(self.mapped_class, attr))

            return self._read_shp(request, self.Session.query(* mapped_attributes).filter(filter), epsg=epsg, ** kwargs)

    def _read_geojson(self, request, query, filter=None, name_mapping={}):
        """
        Output GeoJSON format
        """

        ret = None
        if query.count() == 1:
            o = query.first()
            if o is None:
                abort(404)
            ret = self._filter_attrs(o.toFeature(), request)
        else:
            ret = FeatureCollection(
                                    [self._filter_attrs(o.toFeature(), request) \
                                    for o in query.all()])

        def _replace_mapping(feature):
            """
            Replacing code keys with better understandable keys
            """
            for p in feature.properties.keys():
                if p in name_mapping and name_mapping[p] != p:
                    feature.properties[name_mapping[p]] = feature.properties[p]
                    del feature.properties[p]
            return feature

        if hasattr(ret, "properties"):
            ret = _replace_mapping(ret)
        else:
            ret = FeatureCollection([_replace_mapping(feature) for feature in ret.features])

        return geojson.dumps(ret)

    def _read_ext(self, request, query, filter=None, name_mapping=None):
        """
        A format suitable for Ext json stores.
        """

        output = {}
        output['totalResults'] = self.count(request, filter)

        attrs = request.params['attrs'].split(',')
        mappedAttributes = [getattr(self.mapped_class, i) for i in attrs]

        output['metaData'] = {'totalProperty': 'totalResults', 'root': 'rows'}
        output['metaData']['fields'] = []

        output['rows'] = []

        for i in query.from_self(*mappedAttributes).all():
            row = {}
            for k in attrs:
                value = getattr(i, k)
                #if name_mapping is not None and k in name_mapping:
                #    k = name_mapping[k]

                row[k] = value


                # Check if metadata of current attribute is already declared
                is_declared = False
                for f in output['metaData']['fields']:
                    try:
                        if f['mapping'] == k:
                            is_declared = True
                    except TypeError:
                        pass

                if not is_declared:
                    metadataField = {}
                    metadataField['mapping'] = k
                    if name_mapping is not None and k in name_mapping:
                        metadataField['name'] = name_mapping[k]
                    else:
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

    def _plot_histogram_matplotlib(self, request, query, ** kwargs):
        """
        Alternative implementation with matplotlib instead of R
        """

        # Get the first requested attribute
        attr = request.params.get('attrs').split(",")[0]

        # Set a default value for the image size in pixel
        defaultSide = 480.0
        # Set the dpi
        dpi = 96.0

        try:
            height = float(request.params.get("height", defaultSide))
        except ValueError:
            height = defaultSide
        try:
            width = float(request.params.get("width", defaultSide))
        except ValueError:
            width = defaultSide

        mappedAttribute = getattr(self.mapped_class, attr)
        
        fig = plt.figure(figsize=(width / dpi, height / dpi))
        ax = fig.add_subplot(111)

        # Set fontProperties
        fontProperties = FontProperties(family="sans-serif", size='x-small')

        # Set smaller fonts
        [i.set_fontproperties(fontProperties) for i in ax.get_yticklabels()]
        [j.set_fontproperties(fontProperties) for j in ax.get_xticklabels()]
        
        # If categories is not none, then the current attribute has categories and
        # we want to draw a barplot instead of a histogram
        if kwargs["categories"] is not None:
            
            categories = kwargs["categories"]

            v = []
            names = []
            for a, count in query.from_self(mappedAttribute, func.count(mappedAttribute)).filter(mappedAttribute.in_(categories.keys())).group_by(mappedAttribute):
                v.append(int(count))
                names.append(categories[unicode(a)].encode('UTF-8'))

            N = len(v)

            ind = range(N)

            # the histogram of the data
            ax.bar(np.array(ind) + 0.1, v, width=0.8, color=kwargs.get("color"))

            # hist uses np.histogram under the hood to create 'n' and 'bins'.
            # np.histogram returns the bin edges, so there will be 50 probability
            # density values in n, 51 bin edges in bins and 50 patches.  To get
            # everything lined up, we'll compute the bin centers
            ax.set_xticks(np.arange(len(v)) + 0.5)
            ax.set_xticklabels(names)

        else:

            # Get the number of distinct values
            distinct_value = query.from_self(mappedAttribute).distinct(mappedAttribute).count()

            # Limit the breaks to 100
            if distinct_value > 100:
                distinct_value = 100
            # In case of less distinct values, limit the number of breaks
            elif distinct_value > 20 and distinct_value < 100:
                distinct_value = int(distinct_value / 2)

            v = [i for i, in query.from_self(mappedAttribute).all()]

            # the histogram of the data
            n, bins, patches = ax.hist(v, bins=distinct_value, facecolor=kwargs.get("color"), alpha=0.75)

            """
            n, bins = np.histogram(v, distinct_value)

            # get the corners of the rectangles for the histogram
            left = np.array(bins[:-1])
            right = np.array(bins[1:])
            bottom = np.zeros(len(left))
            top = bottom + n
            nrects = len(left)

            nverts = nrects*(1+3+1)
            verts = np.zeros((nverts, 2))
            codes = np.ones(nverts, int) * matplotlib.path.Path.LINETO
            codes[0::5] = matplotlib.path.Path.MOVETO
            codes[4::5] = matplotlib.path.Path.CLOSEPOLY
            verts[0::5,0] = left
            verts[0::5,1] = bottom
            verts[1::5,0] = left
            verts[1::5,1] = top
            verts[2::5,0] = right
            verts[2::5,1] = top
            verts[3::5,0] = right
            verts[3::5,1] = bottom

            barpath = matplotlib.path.Path(verts, codes)
            patch = matplotlib.patches.PathPatch(barpath, facecolor=kwargs.get("color"), alpha=0.75)
            ax.add_patch(patch)

            ax.set_xlim(0, right[-1])
            ax.set_ylim(bottom.min(), top.max())
            """

            # hist uses np.histogram under the hood to create 'n' and 'bins'.
            # np.histogram returns the bin edges, so there will be 50 probability
            # density values in n, 51 bin edges in bins and 50 patches.  To get
            # everything lined up, we'll compute the bin centers

        if 'xlabel' in kwargs:
            ax.set_xlabel(kwargs['xlabel'], fontproperties=fontProperties)
        if 'ylabel' in kwargs:
            ax.set_ylabel(kwargs['ylabel'], fontproperties=fontProperties)


        ax.grid(True)

        if "filename" in kwargs:
            file = open(kwargs["filename"], 'wb')
            fig.savefig(kwargs["filename"], dpi=dpi, format="png")
        else:
            file = StringIO()
            fig.savefig(file, dpi=dpi, format="png")

        file.seek(0)  # rewind the data

        return file

    """
    @deprecated
    def _plot_histogram_r(self, request, query, ** kwargs):

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

        mappedAttribute = getattr(self.mapped_class, attr)

        rinterface.initr()

        r = robjects.r
        r.library('grDevices')

        if 'filename' in kwargs:
            file = open(kwargs['filename'], 'wb')
        else:
            # Create a temporary file
            file = NamedTemporaryFile()
        r.png(file.name, width=width, height=height)
        r.par(bg="#F0F0F0", mar=robjects.FloatVector([2.6, 4.1, 3.1, 1.1]))

        bar_color = "#BEBEBE"

        # If categories is not none, then the current attribute has categories and
        # we want to draw a barplot instead of a histogram
        if kwargs["categories"] is not None:

            categories = kwargs["categories"]

            names = []
            v = []
            for a, count in query.from_self(mappedAttribute, func.count(mappedAttribute)).filter(mappedAttribute.in_(categories.keys())).group_by(mappedAttribute):
                v.append(int(count))
                names.append(categories[unicode(a)].encode('UTF-8'))
            
            x = robjects.IntVector(v)
            x.names = robjects.StrVector(names)

            r.barplot(x, col=bar_color, xlab=str(), ylab=str(), main=str(), ** {"names.arg": robjects.StrVector(names)})
            #r.par(bg="#F0F0F0", mar=robjects.FloatVector([1.5, 1.5, 1.5, 1.5]))
            #r.pie(x, labels=robjects.StrVector(names), clockwise=True)

        # Handle quantitative data
        else:

            # Get the number of distinct values
            distinct_value = query.from_self(mappedAttribute).distinct(mappedAttribute).count()

            # Limit the breaks to 100
            if distinct_value > 100:
                distinct_value = 100
            # In case of less distinct values, limit the number of breaks
            elif distinct_value > 20 and distinct_value < 100:
                distinct_value = int(distinct_value / 2)

            rbreaks = int(request.params.get("breaks", distinct_value))

            v = [i for i, in query.from_self(mappedAttribute).all()]

            x = robjects.FloatVector(v)

            r.hist(x, col=bar_color, breaks=rbreaks, ylab=str(), xlab=str(), main=str())

        # Finish drawing
        r('dev.off()')

        f = open(file.name, 'r')

        return f
    """


    def _read_xls(self, request, query, ** kwargs):

        requested_attrs = request.params["attrs"].split(",")

        mappedAttributes = [getattr(self.mapped_class, i) for i in requested_attrs]

        workbook = xlwt.Workbook(encoding='utf-8')
        sheet = workbook.add_sheet("data")

        default_header = xlwt.easyxf('font: bold true; borders: bottom THIN;')
        
        row = 0
        column = 0
        for a in requested_attrs:
            sheet.write(row, column, a, default_header)
            column += 1

        row += 1

        for i in query.from_self(*mappedAttributes).all():
            column = 0
            #for a in requested_attrs:
            for attribute in i:
                sheet.write(row, column, attribute)
                column += 1

            row += 1

        if "metadata" in kwargs:

            self._write_metadata(workbook, kwargs["metadata"])
            # Write the workbook to a file-like object
            xls = StringIO()
            # Save the workbook to the memory object
            workbook.save(xls)

        if "filename" in kwargs:
            s = open(kwargs["filename"], "wb")
        else:
            # Create a file-like object
            s = StringIO()
        # Save the workbook to the memory object
        workbook.save(s)
        return s

    def _read_shp(self, request, query, ** kwargs):

        requested_attrs = request.params["attrs"].split(",")

        # Get the first feature to guess the datatype
        first_record = query.first()

        # Create geometry from AsBinary query
        first_geom = loads(str(getattr(first_record, 'geometry_column')))

        w = shapefile.Writer(shapefile.POLYGON)
        if first_geom.geom_type == "Point":
            w = shapefile.Writer(shapefile.POINT)
        elif first_geom.geom_type == "LineString":
            w = shapefile.Writer(shapefile.POLYLINE)
        w.autoBalance = 1

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

            # Handle point geometries
            if g.geom_type == "Point":

                w.point(g.coords[0][0], g.coords[0][1])

            # Handle linestring geometries
            elif g.geom_type == "LineString":

                point_list = []

                for p in g.coords:
                    point_list.append([p[0], p[1]])

                w.line(parts=[point_list])

            # Handle polygon geometries
            elif g.geom_type == "Polygon":

                ring_list = []

                exterior_ring = [[j[0], j[1]] for j in g.exterior.coords]

                ring_list.append(exterior_ring)

                for interior in g.interiors:

                    interior_ring = [[k[0], k[1]] for k in interior.coords]

                    ring_list.append(interior_ring)

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
        if "filename" in kwargs:
            s = open(kwargs["filename"], "wb")
        else:
            # Create a file-like object
            s = StringIO()
        f = ZipFile(s, 'w', ZIP_DEFLATED)
        f.writestr("data.shp", shp.getvalue())
        f.writestr("data.dbf", dbf.getvalue())
        f.writestr("data.shx", shx.getvalue())
        f.writestr("data.cpg", cpg.getvalue())
        f.writestr("data.prj", prj.getvalue())

        if "metadata" in kwargs:
            wb = xlwt.Workbook(encoding='utf-8')

            self._write_metadata(wb, kwargs["metadata"])
            # Write the workbook to a file-like object
            xls = StringIO()
            # Save the workbook to the memory object
            wb.save(xls)

            f.writestr("metadata.xls", xls.getvalue())

        # Close the zip file
        f.close()

        # And return the content
        return s

    def _write_metadata(self, workbook, metadata):

        sheet = workbook.add_sheet("metadata")

        row = 0

        # Create a style that draws a bottom line
        bottomMediumStlye = xlwt.easyxf('font: bold true; borders: bottom THIN;')

        # Write the column headers
        column = 0
        for h in metadata.get_headers():
            sheet.write(row, column, h, bottomMediumStlye)
            column += 1
            
        row += 1

        # Write the variable metadata
        for r in metadata.get_rows():
            column = 0
            for c in r:
                sheet.write(row, column, c)
                column += 1

            row += 1

        # Write contact address
        if metadata.get_address() is not None:
            # One row as space
            row += 1
            sheet.write(row, 0, "*************************************************************")
            row += 1
            sheet.write(row, 0, "*")
            sheet.write(row, 1, "Points of Contact", xlwt.easyxf('font: bold true;'))
            row += 1
            sheet.write(row, 0, "*************************************************************")

        for a in metadata.get_address():
            row += 1
            column = 0
            for c in a:
                sheet.write(row, 0, "*")
                sheet.write(row, 1, c)
                row += 1
            sheet.write(row, 0, "*************************************************************")


    def _get_order_by(self, request):
        """ Return an SA order_by """
        column_name = None
        raw_value = None
        if 'sort' in request.params:
            raw_value = request.params['sort']
            column_name = raw_value.split("#")[0]
        elif 'order_by' in request.params:
            raw_value = request.params['order_by']
            column_name = raw_value.split("#")[0]

        if column_name and column_name in self.mapped_class.__table__.c:
            column = self.mapped_class.__table__.c[column_name]
            if 'dir' in request.params and request.params['dir'].upper() == 'DESC':
                return desc(column)
            elif len(raw_value.split("#")) > 1:
                dir = raw_value.split("#")[1]
                if dir.upper() == "DESC":
                    return desc(column)
                else:
                    return asc(column)
            else:
                return asc(column)

        else:
            return None
