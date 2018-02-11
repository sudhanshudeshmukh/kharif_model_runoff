from __future__ import division
from qgis.gui import QgsMapToolEmitPoint
import qgis.gui
from qgis.core import QgsSpatialIndex, QgsPoint, QgsRectangle, QgsRaster, QgsRasterLayer, QgsVectorLayer, QgsFeature, QgsField, QgsMapLayerRegistry, QgsSymbolV2, QgsRendererRangeV2, QgsGraduatedSymbolRendererV2, QgsVectorFileWriter
from qgis.analysis import  QgsGeometryAnalyzer
from PyQt4.QtGui import *
from PyQt4.QtCore import QVariant
from PyQt4.QtCore import QFileInfo
import csv
import os
import time
import processing
import sys
import shutil
from collections import OrderedDict
from math import exp
from math import log
from constants_dicts_lookups import *
from qgis.gui import QgsMapTool
from PyQt4.QtGui import QCursor, QPixmap
from PyQt4.QtCore import Qt

class VectorLayer:

	def __init__(self, qgsLayer, name=''):
		self.qgsLayer = qgsLayer
		self.name = name
		self.feature_dict = {f.id(): f for f in qgsLayer.getFeatures()}
		self.index = QgsSpatialIndex(qgsLayer.getFeatures())


class RunOffCalculator(QgsMapTool):
    
    def __init__(self, canvas, path, stream_layer, drainage_direction_layer_path, output_grid_points, zones_layer):
        
        super(QgsMapTool, self).__init__(canvas)
        self.canvas = canvas
        self.cursor = QCursor(Qt.CrossCursor)
        self.stream_layer = VectorLayer(stream_layer)
        self.drainage_direction_layer_path = drainage_direction_layer_path
        self.drainage_direction_layer = QgsRasterLayer(drainage_direction_layer_path, 'Drainage_Direction')
        self.path = path
        self.watershed_vector = None
        self.watershed_raster = None
        self.xmin =  self.drainage_direction_layer.extent().xMinimum()
    	self.xmax = self.drainage_direction_layer.extent().xMaximum()
    	self.ymin = self.drainage_direction_layer.extent().yMinimum()
    	self.ymax = self.drainage_direction_layer.extent().yMaximum()
    	self.output_grid_points = output_grid_points
    	self.zones_layer = zones_layer


    def activate(self):
        self.canvas.setCursor(self.cursor)
    
    def screenToLayerCoords(self, screenPos, layer):
        
        transform = self.canvas.getCoordinateTransform()
        canvasPoint = QgsMapToPixel.toMapCoordinates(   transform, 
                                                        screenPos.x(), 
                                                        screenPos.y() )
        
        # Transform if required
        layerEPSG = layer.crs().authid()
        projectEPSG = self.canvas.mapRenderer().destinationCrs().authid()
        if layerEPSG != projectEPSG:
            renderer = self.canvas.mapRenderer()
            layerPoint = renderer.mapToLayerCoordinates(    layer, 
                                                            canvasPoint )
        else:
            layerPoint = canvasPoint
        
        # Convert this point (QgsPoint) to a QgsGeometry
        return QgsGeometry.fromPoint(layerPoint)


    def canvasReleaseEvent(self, mouseEvent):
        """ 
        Each time the mouse is clicked on the map canvas, perform 
        the following tasks:
            Loop through all visible vector layers and for each:
                Ensure no features are selected
                Determine the distance of the closes feature in the layer to the mouse click
                Keep track of the layer id and id of the closest feature
            Select the id of the closes feature 
        """
        
        # Determine the location of the click in real-world coords
        if (self.watershed_raster != None):
        	QgsMapLayerRegistry.instance().removeMapLayer( self.watershed_raster.id() )
        if (self.watershed_vector != None):
        	QgsMapLayerRegistry.instance().removeMapLayer( self.watershed_vector.id() )
        watershed_path= self.path+'/TempWatersheds'
        if os.path.exists(watershed_path):
        	shutil.rmtree(watershed_path)
        if not os.path.exists(watershed_path):
		    os.makedirs(watershed_path)

        point = mouseEvent.pos() 
        x = point.x()
        y = point.y()
        point = self.canvas.getCoordinateTransform().toMapCoordinates(x,y)
        print "x",point.x()
        print "y",point.y()
        nearestIds = self.stream_layer.index.nearestNeighbor(point,1)
        stream_seg = nearestIds[0]
        points_stream_seg = self.stream_layer.feature_dict[stream_seg].geometry().asPolyline()
        min_dist = float('inf')
        for pt in points_stream_seg:
        	if (pt.distance(point[0],point[1]) < min_dist):
        		feat =pt
        		min_dist = pt.distance(point[0],point[1])
        print feat[0],feat[1]
        basin_raster_path = watershed_path+'/watershed_rast.tif'
        basin_vector_path = watershed_path+'/watershed_vector.shp'
        print "Creating Raster layer for Basin"
        processing.runalg('grass7:r.water.outlet',  self.drainage_direction_layer_path , "%f,%f"%(feat[0] , feat[1]) ,"%f,%f,%f,%f"% (self.xmin, self.xmax, self.ymin, self.ymax),0.0,basin_raster_path)
        print "Loading basin Raster"
        fileInfo = QFileInfo(basin_raster_path)
        baseName = fileInfo.baseName()
        self.watershed_raster = QgsRasterLayer(basin_raster_path, baseName)
        if not self.watershed_raster.isValid():
        	print "Basin Raster Layer failed to load!"
        	sys.exit()
        else:
        	QgsMapLayerRegistry.instance().addMapLayer(self.watershed_raster)

        print "Creating Vector layer for Basin"
        processing.runalg("grass7:r.to.vect",basin_raster_path,2,True,"%f,%f,%f,%f"% (self.xmin, self.xmax, self.ymin, self.ymax),0,basin_vector_path)
        self.watershed_vector = QgsVectorLayer(basin_vector_path,'watershed_vector' , 'ogr')
        if not self.watershed_vector.isValid():
        	print "Basin vector Layer failed to load!"
        	sys.exit()
        else:
        	QgsMapLayerRegistry.instance().addMapLayer(self.watershed_vector)

        filtered_points = []
        for pts in self.output_grid_points:
        	for polygon in self.watershed_vector.getFeatures():
        		if (polygon.geometry().contains(pts.qgsPoint)):
        			filtered_points.append(pts)
        print (len(filtered_points))
        area = 0
        runoff_in_mm = sum([p.budget.runoff_monsoon_end for p in filtered_points]) / len(filtered_points)

        for zone_id in self.zones_layer.feature_dict :
        	for polygon in self.watershed_vector.getFeatures():
        		area += polygon.geometry().intersection(self.zones_layer.feature_dict[zone_id].geometry()).area() 
        print "runoff_in_mm", runoff_in_mm
        print "Area", area/10000
        print "run off ", (runoff_in_mm/1000.0 * area) / 1000 
        strn = "Run off in mm : "+str(runoff_in_mm)+"mm\nArea of watershed intersecting with Zones: "+str(area/10000) + "hectare\n Run of in TCM: "+str((runoff_in_mm/1000.0 * area) / 1000 )
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setWindowTitle("Runoff Details:")
        msg.setText(strn)
        retval = msg.exec_()