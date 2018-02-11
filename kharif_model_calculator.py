from __future__ import division
from qgis.gui import QgsMapToolEmitPoint
import qgis.gui
from qgis.core import QgsSpatialIndex, QgsPoint, QgsRectangle, QgsRaster, QgsVectorLayer, QgsFeature, QgsField, QgsMapLayerRegistry, QgsSymbolV2, QgsRendererRangeV2, QgsGraduatedSymbolRendererV2, QgsVectorFileWriter
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

printed = False

BOUNDARY_LABEL = 'Zones'
SOIL_LABEL = 'Soil'
LULC_LABEL = 'Land-Use-Land-Cover'
SLOPE_LABEL = 'Slope'
CADASTRAL_LABEL = 'Cadastral'

CALCULATE_FOR_LULC_TYPES = ['agriculture', 'fallow land']

class Budget:

	def __init__(self):
		self.sm, self.runoff, self.infil, self.AET, self.GW_rech = [],[],[],[],[]

	def summarize(self, PET_sum, PET_sum_cropend, start_date_index, end_date_index, crop_end_index):
		self.sm_crop_end = self.sm[crop_end_index]
		self.sm_monsoon_end = self.sm[123]
		self.runoff_crop_end = sum(self.runoff[start_date_index:crop_end_index+1])
		self.runoff_monsoon_end = sum(self.runoff[start_date_index:123])
		self.infil_crop_end = sum(self.infil[start_date_index:crop_end_index+1])
		self.infil_monsoon_end = sum(self.infil[start_date_index:123])
		self.AET_crop_end = sum(self.AET[start_date_index:crop_end_index+1])
		self.AET_monsoon_end = sum(self.AET[start_date_index:123])
		self.GW_rech_crop_end = sum(self.GW_rech[start_date_index:crop_end_index+1])
		self.GW_rech_monsoon_end = sum(self.GW_rech[start_date_index:123])
		self.PET_minus_AET_monsoon_end = PET_sum - self.AET_monsoon_end
		self.PET_minus_AET_post_monsoon= (PET_sum_cropend - self.AET_crop_end)-self.PET_minus_AET_monsoon_end
		self.PET_minus_AET_crop_end= (PET_sum_cropend - self.AET_crop_end)


class Point:

	def __init__(self, qgsPoint):
		self.qgsPoint = qgsPoint
		self.container_polygons = {}
		self.slope = None
		self.budget = Budget()
		self.crop_name = None

	def run_model(self, rain, PET, PET_sum, PET_sum_cropend, start_date_index, end_date_index,crop_end_index):
		self.setup_for_daily_computations()
		PET = PET[self.crop_name]
		PET_sum = PET_sum[self.crop_name]
		PET_sum_cropend = PET_sum_cropend[self.crop_name]
		self.SM1_fraction = self.layer2_moisture = self.WP

		self.not_small_KS = []
		# self.not_small_SM = []
		for day in range (365):
			self.primary_runoff(day, rain)
			self.aet(day, PET)
			self.percolation_below_root_zone(day)
			self.secondary_runoff(day)
			self.percolation_to_GW(day)
		# print len(self.not_small_KS)#, len(self.not_small_SM)
		# if not printed:
		# 	print PET[183:]
		# 	printed = True
		self.budget.summarize(PET_sum, PET_sum_cropend, start_date_index, end_date_index,crop_end_index)
		self.budget.sm_crop_end -= self.WP_depth	# requirement expressed by users
		self.budget.sm_monsoon_end -= self.WP_depth	# requirement expressed by users

	def setup_for_daily_computations(self):
		"""
		"""
		poly_soil = self.container_polygons[SOIL_LABEL]
		texture = poly_soil[TEX].lower()
		Ksat = round(float(dict_SoilContent[texture][7]),4)
		self.Sat = round(float(dict_SoilContent[texture][6]),4)
		self.WP = round(float(dict_SoilContent[texture][4]),4)
		self.FC = round(float(dict_SoilContent[texture][5]),4)
		depth_value = float(dict_SoilDep[poly_soil[Depth].lower()])

		lu_Type = dict_lulc[self.container_polygons[LULC_LABEL][Desc].lower()]

		HSG =  dict_SoilContent[texture][0]
		cn_val = int(dict_RO[lu_Type][HSG])

		Sat_depth = self.Sat * depth_value*1000
		self.WP_depth = self.WP*depth_value*1000
		FC_depth = self.FC*depth_value*1000
		ROOT_LEVEL = dict_crop_root[self.crop_name]
		if(depth_value <= ROOT_LEVEL): #thin soil layer
			self.SM1 = depth_value - 0.01;	self.SM2 = 0.01
		else:
			self.SM1 = ROOT_LEVEL;	self.SM2 = depth_value - ROOT_LEVEL

		cn_s = cn_val
		cn3 = cn_s *exp(0.00673*(100-cn_s))
		if (self.slope > 5.0):
			cn_s = (((cn3-cn_val)/float(3))*(1-2*exp(-13.86*self.slope * 0.01))) + cn_val
		cn1_s = cn_s - 20*(100-cn_s)/float(100-cn_s+exp(2.533-0.0636*(100-cn_s)))
		cn3_s = cn_s *exp(0.00673*(100-cn_s))

		self.Smax = 25.4 * (1000/float(cn1_s) - 10)
		S3 = 25.4 * (1000/float(cn3_s) - 10)
		self.W2 = (log((FC_depth- self.WP_depth)/(1-float(S3/self.Smax)) - (FC_depth - self.WP_depth )) - log ((Sat_depth - self.WP_depth)/(1-2.54/self.Smax) - (Sat_depth - self.WP_depth)))/((Sat_depth- self.WP_depth) - (FC_depth - self.WP_depth))
		self.W1 = log((FC_depth- self.WP_depth)/(1- S3/self.Smax) - (FC_depth - self.WP_depth)) + self.W2 * (FC_depth -self.WP_depth)

		TT_perc = (Sat_depth- FC_depth)/Ksat	#SWAT equation 2:3.2.4
		self.daily_perc_factor = 1 - exp(-24 / TT_perc)	#SWAT equation 2:3.2.3

		self.depletion_factor = dict_crop[self.crop_name][1]

	def primary_runoff(self, day, rain):
		"""
		Retention parameter 'S_swat' using SWAT equation 2:1.1.6
		Curve Number for the day 'Cn_swat' using SWAT equation 2:1.1.11
		Initial abstractions (surface storage,interception and infiltration prior to runoff)
			'Ia_swat' derived approximately as recommended by SWAT
		Primary Runoff 'Swat_RO' using SWAT equation 2:1.1.1
		"""
		self.budget.sm.append((self.SM1_fraction * self.SM1 + self.layer2_moisture * self.SM2) * 1000)
		# if not printed and self.sm[-1] > 100:    self.
		self.SW = self.budget.sm[-1] - self.WP_depth
		S_swat = self.Smax*(1 - self.SW/(self.SW + exp(self.W1 - self.W2 * self.SW)))

		Cn_swat = 25400/float(S_swat+254)
		Ia_swat = 0.2 * S_swat
		#~ print 'len(rain), day : ', len(rain), day
		if(rain[day] > Ia_swat):
			self.budget.runoff.append(((rain[day]-Ia_swat)**2)/(rain[day] + 0.8*S_swat))
		else:
			self.budget.runoff.append(0)
		self.budget.infil.append(rain[day] - self.budget.runoff[day])
		assert len(self.budget.runoff) == day+1, (self.budget.runoff, day)
		assert len(self.budget.infil) == day+1

	def aet(self, day, PET):
		"""
		Water Stress Coefficient 'KS' using FAO Irrigation and Drainage Paper 56, page 167 and
			page 169 equation 84
		Actual Evapotranspiration 'AET' using FAO Irrigation and Drainage Paper 56, page 6 and 
			page 161 equation 81
		"""
		global printed
		if (self.SM1_fraction < self.WP):
			KS= 0
		elif (self.SM1_fraction > (self.FC *(1- self.depletion_factor) + self.depletion_factor * self.WP)):
			KS = 1
		else :
			KS = (self.SM1_fraction - self.WP)/(self.FC - self.WP) /(1- self.depletion_factor)
		if day > 182 and KS > 0.1:	self.not_small_KS.append(KS)
		self.budget.AET.append( KS * PET[day] )

	def percolation_below_root_zone(self, day):
		"""
		Calculate soil moisture (fraction) 'SM1_before' as the one after infiltration and (then) AET occur,
		but before percolation starts below root-zone. Percolation below root-zone starts only if
		'SM1_before' is more than field capacity and the soil below root-zone is not saturated,i.e.
		'layer2_moisture' is less than saturation. When precolation occurs it is derived as
		the minimum of the maximum possible percolation (using SWAT equation 2:3.2.3) and
		the amount available in the root-zone for percolation.
		"""
		self.SM1_before = (self.SM1_fraction*self.SM1 +((self.budget.infil[day]-self.budget.AET[day])/float(1000)))/self.SM1
		if (self.SM1_before < self.FC):
			self.R_to_second_layer =0
		elif (self.layer2_moisture < self.Sat) :
			self.R_to_second_layer = min((self.Sat - self.layer2_moisture) * self.SM2 * 1000,
										 (self.SM1_before - self.FC) * self.SM1 * 1000 * self.daily_perc_factor)
		else :
			self.R_to_second_layer = 0
		self.SM2_before = (self.layer2_moisture*self.SM2*1000 + self.R_to_second_layer)/self.SM2/1000

	def secondary_runoff(self, day):
		"""
		
		"""
		if (((self.SM1_before*self.SM1 - self.R_to_second_layer/1000)/self.SM1) > self.Sat):
			sec_run_off= (((self.SM1_before*self.SM1 - self.R_to_second_layer/1000)/self.SM1) - self.Sat) *self.SM1*1000
		else:
			sec_run_off = 0
		self.SM1_fraction = min((self.SM1_before*self.SM1*1000 - self.R_to_second_layer)/self.SM1/1000,self.Sat)

	def percolation_to_GW(self, day):
		"""
		
		"""
		self.budget.GW_rech.append(max((self.SM2_before - self.FC)*self.SM2*self.daily_perc_factor*1000,0))
		self.layer2_moisture = min(((self.SM2_before*self.SM2*1000- self.budget.GW_rech[day])/self.SM2/1000),self.Sat)


class VectorLayer:

	def __init__(self, qgsLayer, name=''):
		self.qgsLayer = qgsLayer
		self.name = name
		self.feature_dict = {f.id(): f for f in qgsLayer.getFeatures()}
		self.index = QgsSpatialIndex(qgsLayer.getFeatures())

	def get_polygon_containing_point(self, point):
		intersector_ids = self.index.intersects( QgsRectangle( point.qgsPoint, point.qgsPoint ) )
		for intersector_id in intersector_ids:
			if self.name == CADASTRAL_LABEL and intersector_id not in self.feature_dict:
				continue
			polygon = self.feature_dict[intersector_id]
			if (polygon.geometry().contains(point.qgsPoint)):
				return polygon
		return None


class KharifModelCalculator:
	"""
	The actual algorithm for calculating results of the Kharif Model
	"""

	def __init__(self, path, zones_layer, soil_layer, lulc_layer, cadastral_layer, slope_layer, rainfall_csv_path):
		self.zones_layer = VectorLayer(zones_layer, BOUNDARY_LABEL)
		self.soil_layer = VectorLayer(soil_layer, SOIL_LABEL)
		self.lulc_layer = VectorLayer(lulc_layer, LULC_LABEL)
		self.cadastral_layer = VectorLayer(cadastral_layer, CADASTRAL_LABEL)
		zone_polygon_ids = self.zones_layer.feature_dict.keys()
		self.zone_points_dict = dict(zip(zone_polygon_ids, [[]	for i in range(len(zone_polygon_ids))]))
		cadastral_polygon_ids = self.cadastral_layer.feature_dict.keys()
		#~ self.cadastral_points_dict = dict(zip(cadastral_polygon_ids, [[]	for i in range(len(cadastral_polygon_ids))]))
		#~ print 'zone_points_dict : ', self.zone_points_dict

		self.slope_layer = slope_layer

		# Working Directory path
		self.path = path

		self.path_et = path + '/ET0_file.csv'

		#~ rainfall_csv_path = path + '/rainfall.csv'
		rainfall_csv = open(rainfall_csv_path)
		self.rain = [float(row["Rainfall"]) for row in csv.DictReader(rainfall_csv)]
		print 'len(rain) = ', len(self.rain)

	def pet_calculation(self, crop_name, sowing_threshold):
		test_csv = open(self.path_et)
		a = [float(row["ET0"]) for row in csv.DictReader(test_csv)]
		Kc=[]
		stage=[]
		print (crop_name)
		self.PET={}

		#computes PET array for crop
		def pet_crop(crop):
			Kc_crop = []
			for i in range (len(dict_crop[crop][0])):
				stage = dict_crop[crop][0][i][0]
				Kc = dict_crop[crop][0][i][1]
				Kc_crop.extend([Kc]*stage)
			Kc_crop = [0]*self.initial_dry+ Kc_crop
			Kc_len=len(Kc_crop)
			if(Kc_len<= self.duration):
				Kc_crop = Kc_crop + [0]*(self.duration-Kc_len)
			elif(Kc_len>self.duration):
				Kc_crop = Kc_crop[0:self.duration]
			self.PET[crop]= [self.et0[i]*Kc_crop[i] for i in range (0,self.duration)]

		#computes initial period where total rainfall is below the sowing-threshold for sowing date
		def compute_initial_period():
			rain_sum = 0
			for i in range (0,len(self.rain)):
				if (rain_sum < sowing_threshold):
					rain_sum += self.rain[i]
				else :
					break
			return i
		self.initial_dry = compute_initial_period() if sowing_threshold != 0 else 0
		self.duration=sum (i[0] for i in dict_crop[crop_name][0])
		self.duration+=self.initial_dry
		self.duration=max(self.duration,365)
		#Computation of ET0 from June to May
		self.et0=[]
		for i in range (0,len(a)):
			if (i in [0,3,5,10]):
				c = [a[i]]*30
				self.et0.extend(c)
			elif (i in [9]):
				c = [a[i]]*28
				self.et0.extend(c)
			else:
				c = [a[i]]*31
				self.et0.extend(c)
		self.et0=self.et0[0:self.duration]

		if (self.duration>183):
			self.rain = self.rain + [0]*(self.duration-len(self.rain))
		pet_crop(crop_name)
		pet_crop('deciduous - dense crop')
		pet_crop('deciduous open crop')
		pet_crop('scrub dense crop')
		pet_crop('scrub forest crop')
		pet_crop('scrub open crop')

	def filter_out_cadastral_plots_outside_boundary(self):
		#~ QgsGeometryAnalyzer().dissolve(self.zones_layer.qgsLayer, 'temp.shp')
		#~ dissolved_zones_layer = QgsVectorLayer('temp.shp', 'dissolved boundary', 'ogr')
		filtered_feature_dict = {}
		for polygon_id in self.cadastral_layer.feature_dict:
			for feature in self.zones_layer.qgsLayer.getFeatures():
				if self.cadastral_layer.feature_dict[polygon_id].geometry().intersects(feature.geometry()):
					filtered_feature_dict[polygon_id] = self.cadastral_layer.feature_dict[polygon_id]
					break
		self.cadastral_layer.feature_dict = filtered_feature_dict

	def generate_output_points_grid(self, input_points_filename=None):
		if input_points_filename is None:
			xminB =  self.zones_layer.qgsLayer.extent().xMinimum()
			xmaxB = self.zones_layer.qgsLayer.extent().xMaximum()
			yminB = self.zones_layer.qgsLayer.extent().yMinimum()
			ymaxB = self.zones_layer.qgsLayer.extent().yMaximum()
			print 'boundary min, max : ' , xminB, xmaxB, yminB, ymaxB
			def frange(start,end,step):
				i = start
				while i<=end :
					yield i
					i = i+step
			x_List = [x for x in frange(xminB,xmaxB,STEP)]
			y_List = [x for x in frange(yminB,ymaxB,STEP)]
			print len(x_List)
			print len (y_List)
		else:
			# Read from file
			pass
		output_points = [Point(QgsPoint(x,y))	for x in x_List	for y in y_List]
		'''print "No of points: ", len (output_points)
		csvwrite = open(self.path + "/tt.csv",'w+b')
		writer = csv.writer(csvwrite)
		writer.writerow(['X', 'Y'])
		writer.writerow([xminB,yminB])
		writer.writerow([xminB,ymaxB])
		writer.writerow([xmaxB,yminB])
		writer.writerow([xmaxB,ymaxB])		
		csvwrite.close()'''

		return output_points

	def filter_out_points_outside_boundary(self):
		filtered_points = []
		'''csvwrite = open(self.path + "/tt1.csv",'w+b')
		writer = csv.writer(csvwrite)
		writer.writerow(['X', 'Y'])'''

		for point in self.output_grid_points:
			polygon = self.zones_layer.get_polygon_containing_point(point)
			if polygon is not None:
				point.container_polygons[BOUNDARY_LABEL] = polygon
				self.zone_points_dict[polygon.id()].append(point)
				filtered_points.append(point)
			'''else:
				writer.writerow([point.qgsPoint.x(), point.qgsPoint.y()])
		csvwrite.close()'''

		self.output_grid_points = filtered_points
		'''print(len(filtered_points))
		for zone_id in self.zone_points_dict:
			print str(self.zones_layer.feature_dict[zone_id]['Zone_name'])+ ' '+ str(len(self.zone_points_dict[zone_id]))'''


	def generate_output_points_for_cadastral_plots(self):
		#~ cadastral_points_dict = {};
		output_cadastral_points = []
		for polygon_id in self.cadastral_layer.feature_dict:
			qgsPoint = self.cadastral_layer.feature_dict[polygon_id].geometry().centroid().asPoint()
			polygon_geom = self.cadastral_layer.feature_dict[polygon_id].geometry()
			if polygon_geom.contains(qgsPoint):
				point = Point(qgsPoint)
			else:
				for feature in self.zones_layer.qgsLayer.getFeatures():
					if polygon_geom.intersects(feature.geometry()):
						polygon_intersection_some_zone = polygon_geom.intersection(feature.geometry())
						point = Point(polygon_intersection_some_zone.pointOnSurface().asPoint())
						break
			point.container_polygons[CADASTRAL_LABEL] = self.cadastral_layer.feature_dict[polygon_id]
			#~ cadastral_points_dict[polygon_id] = point
			output_cadastral_points.append(point)
		#~ return cadastral_points_dict, output_cadastral_points
		return output_cadastral_points

	def set_container_polygon_of_points_for_layers(self, points, polygon_vector_layers):
		for layer in polygon_vector_layers:
			for p in points:
				p.container_polygons[layer.name] = layer.get_polygon_containing_point(p)
				#~ if p.container_polygons[layer.name] is None:	raise Exception(layer.name, p.qgsPoint.x(), p.qgsPoint.y())

	def set_slope_at_points(self, points):
		for point in points:
			identify_results = self.slope_layer.dataProvider().identify(point.qgsPoint, QgsRaster.IdentifyFormatValue)
			if identify_results.isValid():    point.slope = identify_results.results()[1]
			else:                             point.slope = None

	def set_crop_at_points(self, points,crop_name):
		for point in points:
			if (dict_lulc[point.container_polygons[LULC_LABEL][Desc].lower()] in ['agriculture', 'fallow land']):
				point.crop_name =crop_name
			elif(dict_lulc[point.container_polygons[LULC_LABEL][Desc].lower()] in ['deciduous - dense','deciduous open','scrub dense','scrub forest','scrub open']):
				point.crop_name =other_LU_crops[dict_lulc[point.container_polygons[LULC_LABEL][Desc].lower()]]


	def filter_out_points_with_incomplete_data(self, points):
		log_file = open(os.path.join(self.path, 'log'), 'a')
		log_file.write(time.ctime(time.time()) + '\n')
		filtered_points = []
		for point in points:
			if (None not in [
				point.container_polygons[SOIL_LABEL],
				point.container_polygons[LULC_LABEL],
				point.slope
			]):
				filtered_points.append(point)
			else:
				if point.container_polygons[SOIL_LABEL] is None:
					log_file.write('Soil polygon could not be obtained for point at: x = '
								   + str(point.qgsPoint.x()) + ', y = ' + str(point.qgsPoint.y()))
				if point.container_polygons[LULC_LABEL] is None:
					log_file.write('LULC polygon could not be obtained for point at: x = '
								   + str(point.qgsPoint.x()) + ', y = ' + str(point.qgsPoint.y()))
				if point.slope is None:
					log_file.write('Slope could not be obtained for point at: x = '
								   + str(point.qgsPoint.x()) + ', y = ' + str(point.qgsPoint.y()))
		log_file.close()
		return filtered_points

	def output_point_results_to_csv(self, pointwise_output_csv_filename):
		csvwrite = open(self.path + pointwise_output_csv_filename,'w+b')
		writer = csv.writer(csvwrite)
		writer.writerow(['X', 'Y','Monsoon PET-AET','Crop duration PET-AET' ,'Soil Moisture','Infiltration', 'Runoff', 'GW Recharge'])
		for point in self.output_grid_points:
			if not point.container_polygons[BOUNDARY_LABEL]:	continue
			if point.budget.sm_crop_end > point.budget.infil_monsoon_end:
				print point.container_polygons[SOIL_LABEL][Depth]
			writer.writerow([point.qgsPoint.x(), point.qgsPoint.y(), point.budget.PET_minus_AET_monsoon_end, point.budget.PET_minus_AET_crop_end, point.budget.sm_crop_end, point.budget.infil_monsoon_end, point.budget.runoff_monsoon_end, point.budget.GW_rech_monsoon_end])
		csvwrite.close()

	def compute_zonewise_budget(self):
		self.zonewise_budgets = OrderedDict()
		for zone_id in self.zone_points_dict:
			if (self.zones_layer.qgsLayer.fieldNameIndex('Zone_name') != -1):
				zone_name = self.zones_layer.feature_dict[zone_id]['Zone_name']
			else:
				zone_name = zone_id

			zone_points = self.zone_points_dict[zone_id]
			no_of_zone_points = len(zone_points)
			if no_of_zone_points == 0:	continue
			self.zonewise_budgets[zone_name] = {'agricultural': OrderedDict(), 'non-agricultural': OrderedDict()}
			all_agricultural_points = filter(lambda p:	dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()] in ['agriculture', 'fallow land'], zone_points)
			non_agricultural_points_dict = {lulc_type: filter(lambda p:	dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()] == lulc_type, zone_points)	for lulc_type in self.lulc_types if lulc_type not in ['agriculture', 'fallow land', 'water']}
			no_of_agricultural_type_points = len(all_agricultural_points)

			if(no_of_agricultural_type_points!=0):
				no_of_soil_type_points = {}
				for soil_type in self.soil_types:
					soil_type_points = filter(lambda p:	p.container_polygons[SOIL_LABEL][TEX].lower() == soil_type, all_agricultural_points)
					no_of_soil_type_points[soil_type] = len(soil_type_points)
					if no_of_soil_type_points[soil_type] == 0:	continue

					zb = self.zonewise_budgets[zone_name]['agricultural'][soil_type] = Budget()
					zb.sm_crop_end = sum([p.budget.sm_crop_end	for p in soil_type_points]) / no_of_soil_type_points[soil_type]
					zb.sm_monsoon_end = sum([p.budget.sm_monsoon_end for p in soil_type_points])/no_of_soil_type_points[soil_type]
					zb.runoff_monsoon_end = sum([p.budget.runoff_monsoon_end	for p in soil_type_points]) / no_of_soil_type_points[soil_type]
					zb.infil_monsoon_end = sum([p.budget.infil_monsoon_end	for p in soil_type_points]) / no_of_soil_type_points[soil_type]
					zb.AET_crop_end = sum([p.budget.AET_crop_end	for p in soil_type_points]) / no_of_soil_type_points[soil_type]
					zb.AET_monsoon_end = sum([p.budget.AET_monsoon_end	for p in soil_type_points]) / no_of_soil_type_points[soil_type]
					zb.GW_rech_monsoon_end = sum([p.budget.GW_rech_monsoon_end	for p in soil_type_points]) / no_of_soil_type_points[soil_type]
					zb.PET_minus_AET_monsoon_end = sum([p.budget.PET_minus_AET_monsoon_end	for p in soil_type_points]) / no_of_soil_type_points[soil_type]
					zb.PET_minus_AET_crop_end = sum([p.budget.PET_minus_AET_crop_end	for p in soil_type_points]) / no_of_soil_type_points[soil_type]
				zb = self.zonewise_budgets[zone_name]['agricultural']['Agricultural Total'] = Budget()
				zb.sm_crop_end = sum([p.budget.sm_crop_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.sm_monsoon_end = sum([p.budget.sm_monsoon_end for p in all_agricultural_points] ) /no_of_agricultural_type_points
				zb.runoff_monsoon_end = sum([p.budget.runoff_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.infil_monsoon_end = sum([p.budget.infil_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.AET_crop_end = sum([p.budget.AET_crop_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.AET_monsoon_end = sum([p.budget.AET_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.GW_rech_monsoon_end = sum([p.budget.GW_rech_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.PET_minus_AET_monsoon_end = sum([p.budget.PET_minus_AET_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.PET_minus_AET_crop_end = sum([p.budget.PET_minus_AET_crop_end	for p in all_agricultural_points]) / no_of_agricultural_type_points


	def output_zonewise_budget_to_csv_agri(self, zonewise_budget_csv_filename, PET_sum_cropend, PET_sum,rain_sum):
		csvwrite = open(self.path + zonewise_budget_csv_filename,'wb')
		writer = csv.writer(csvwrite)
		writer.writerow([''] + ['zone-'+str(ID)+'-'+t	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Rainfall'] + [rain_sum	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Runoff in Monsoon'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].runoff_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Infiltration in Monsoon'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].infil_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Soil Moisture Crop end'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].sm_crop_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Soil Moisture Monsoon end'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].sm_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['GW Recharge in Monsoon'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].GW_rech_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['AET Crop End'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].AET_crop_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['PET Crop End'] + [PET_sum_cropend[self.crop_name] if (ag_or_non_ag == 'agricultural' or ag_or_non_ag == 'zone' ) else PET_sum_cropend[other_LU_crops[t]] for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag] ])
		writer.writerow(['AET Monsoon End'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].AET_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['PET Monsoon End'] + [PET_sum[self.crop_name] if (ag_or_non_ag == 'agricultural' or ag_or_non_ag == 'zone' ) else PET_sum_cropend[other_LU_crops[t]] for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag] ])
		writer.writerow(['Monsoon Deficit(PET-AET)'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].PET_minus_AET_monsoon_end	for ID in self.zonewise_budgets		for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Crop duration Deficit(PET-AET)'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].PET_minus_AET_crop_end	for ID in self.zonewise_budgets		for ag_or_non_ag in ['agricultural']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		csvwrite.close()

	def compute_zonewise_budget_LU(self):
		self.zonewise_budgets = OrderedDict()
		for zone_id in self.zone_points_dict:
			if(self.zones_layer.qgsLayer.fieldNameIndex('Zone_name') != -1):
				zone_name = self.zones_layer.feature_dict[zone_id]['Zone_name']
			else:
				zone_name = zone_id
			zone_points = self.zone_points_dict[zone_id]
			no_of_zone_points = len(zone_points)
			if no_of_zone_points == 0:	continue
			self.zonewise_budgets[zone_name] = {'agricultural': OrderedDict(), 'non-agricultural': OrderedDict()}
			all_agricultural_points = filter(lambda p:	dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()] in ['agriculture', 'fallow land'], zone_points)
			non_agricultural_points_dict = {lulc_type: filter(lambda p:	dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()] == lulc_type, zone_points)	for lulc_type in self.lulc_types if lulc_type not in ['agriculture', 'fallow land', 'water']}

			no_of_agricultural_type_points = len(all_agricultural_points)
			if(no_of_agricultural_type_points!=0):
				zb = self.zonewise_budgets[zone_name]['agricultural']['Agricultural Total'] = Budget()
				zb.sm_crop_end = sum([p.budget.sm_crop_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.sm_monsoon_end = sum([p.budget.sm_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.runoff_monsoon_end = sum([p.budget.runoff_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.infil_monsoon_end = sum([p.budget.infil_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.AET_crop_end = sum([p.budget.AET_crop_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.AET_monsoon_end = sum([p.budget.AET_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.GW_rech_monsoon_end = sum([p.budget.GW_rech_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.PET_minus_AET_monsoon_end = sum([p.budget.PET_minus_AET_monsoon_end	for p in all_agricultural_points]) / no_of_agricultural_type_points
				zb.PET_minus_AET_crop_end = sum([p.budget.PET_minus_AET_crop_end	for p in all_agricultural_points]) / no_of_agricultural_type_points

			no_of_non_ag_lulc_type_points = {}
			for lulc_type in non_agricultural_points_dict:
				lulc_type_points = non_agricultural_points_dict[lulc_type]
				no_of_non_ag_lulc_type_points[lulc_type] = len(lulc_type_points)
				if no_of_non_ag_lulc_type_points[lulc_type] == 0:	continue

				zb = self.zonewise_budgets[zone_name]['non-agricultural'][lulc_type] = Budget()
				zb.sm_crop_end = sum([p.budget.sm_crop_end	for p in lulc_type_points]) / no_of_non_ag_lulc_type_points[lulc_type]
				zb.sm_monsoon_end = sum([p.budget.sm_monsoon_end	for p in lulc_type_points]) / no_of_non_ag_lulc_type_points[lulc_type]
				zb.runoff_monsoon_end = sum([p.budget.runoff_monsoon_end	for p in lulc_type_points]) / no_of_non_ag_lulc_type_points[lulc_type]
				zb.infil_monsoon_end = sum([p.budget.infil_monsoon_end	for p in lulc_type_points]) / no_of_non_ag_lulc_type_points[lulc_type]
				zb.AET_crop_end = sum([p.budget.AET_crop_end	for p in lulc_type_points]) / no_of_non_ag_lulc_type_points[lulc_type]
				zb.AET_monsoon_end = sum([p.budget.AET_monsoon_end	for p in lulc_type_points]) / no_of_non_ag_lulc_type_points[lulc_type]
				zb.GW_rech_monsoon_end = sum([p.budget.GW_rech_monsoon_end	for p in lulc_type_points]) / no_of_non_ag_lulc_type_points[lulc_type]
				zb.PET_minus_AET_monsoon_end = sum([p.budget.PET_minus_AET_monsoon_end	for p in lulc_type_points]) / no_of_non_ag_lulc_type_points[lulc_type]
				zb.PET_minus_AET_crop_end = sum([p.budget.PET_minus_AET_crop_end	for p in lulc_type_points]) / no_of_non_ag_lulc_type_points[lulc_type]

			zb = Budget()
			zb.sm_crop_end = sum([p.budget.sm_crop_end	for p in zone_points]) / no_of_zone_points
			zb.sm_monsoon_end = sum([p.budget.sm_monsoon_end	for p in zone_points]) / no_of_zone_points
			zb.runoff_monsoon_end = sum([p.budget.runoff_monsoon_end	for p in zone_points]) / no_of_zone_points
			zb.infil_monsoon_end = sum([p.budget.infil_monsoon_end	for p in zone_points]) / no_of_zone_points
			zb.AET_crop_end = sum([p.budget.AET_crop_end	for p in zone_points]) / no_of_zone_points
			zb.AET_monsoon_end = sum([p.budget.AET_monsoon_end	for p in zone_points]) / no_of_zone_points
			zb.GW_rech_monsoon_end = sum([p.budget.GW_rech_monsoon_end	for p in zone_points]) / no_of_zone_points
			zb.PET_minus_AET_monsoon_end = sum([p.budget.PET_minus_AET_monsoon_end	for p in zone_points]) / no_of_zone_points
			zb.PET_minus_AET_crop_end = sum([p.budget.PET_minus_AET_crop_end	for p in zone_points]) / no_of_zone_points
			self.zonewise_budgets[zone_name]['zone'] = {'overall':zb}  # dict {'overall':zb} assigned instead of simple zb for convenience in iterating with ag and non-ag

	def compute_zonewise_budget_areawise(self):
		self.zonewise_budgets = OrderedDict()
		for zone_id in self.zone_points_dict:
			if (self.zones_layer.qgsLayer.fieldNameIndex('Zone_name') != -1):
				zone_name = self.zones_layer.feature_dict[zone_id]['Zone_name']
			else:
				zone_name = zone_id

			zone_points = self.zone_points_dict[zone_id]
			if len(zone_points) == 0:	continue
			rain_in_mm = sum(self.rain[:183])
			gw_rech_in_mm = sum([p.budget.GW_rech_monsoon_end for p in zone_points]) / len(zone_points)
			runoff_in_mm = sum([p.budget.runoff_monsoon_end for p in zone_points]) / len(zone_points)
			ag_area_total = non_ag_area_total = 0
			for ID, polygon in self.lulc_layer.feature_dict.items():
				if dict_lulc[polygon[Desc].lower()] not in ['water', 'habitation']:
					intersection_area = polygon.geometry().intersection(self.zones_layer.feature_dict[zone_id].geometry()).area()
				if dict_lulc[polygon[Desc].lower()] in ['agriculture', 'fallow land']:
					ag_area_total += intersection_area
				elif dict_lulc[polygon[Desc].lower()] not in ['water', 'habitation']:
					non_ag_area_total += intersection_area
			zone_points_ag = [p	for p in zone_points
							  	if dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()]
							  		in ['agriculture', 'fallow land']]
			zone_points_non_ag = [p for p in zone_points
									if dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()]
							  			not in ['agriculture', 'fallow land', 'water', 'habitation']]
			self.zonewise_budgets[zone_name] = {}

			if(len(zone_points_ag)!=0):
				sm_in_mm = sum([p.budget.sm_crop_end	for p in zone_points_ag])/len(zone_points_ag)
				sm_monsoon_in_mm = sum([p.budget.sm_monsoon_end	for p in zone_points_ag])/len(zone_points_ag)
				deficit_in_mm = sum([p.budget.PET_minus_AET_crop_end for p in zone_points_ag]) / len(zone_points_ag)
				monsoon_deficit_in_mm = sum([p.budget.PET_minus_AET_monsoon_end for p in zone_points_ag]) / len(zone_points_ag)
			else:
				sm_in_mm = 0
				sm_monsoon_in_mm = 0
				deficit_in_mm = 0
				monsoon_deficit_in_mm = 0

			self.zonewise_budgets[zone_name]['ag_area'] = ag_area_total/10000.0
			self.zonewise_budgets[zone_name]['non_ag_area'] = non_ag_area_total/10000.0
			self.zonewise_budgets[zone_name]['rain_mm'] = rain_in_mm
			self.zonewise_budgets[zone_name]['rain_TCM'] = (rain_in_mm/1000.0 * (ag_area_total + non_ag_area_total)) / 1000
			self.zonewise_budgets[zone_name]['gw_rech'] = (gw_rech_in_mm/1000.0 * (ag_area_total + non_ag_area_total)) / 1000
			self.zonewise_budgets[zone_name]['runoff'] = (runoff_in_mm/1000.0 * (ag_area_total + non_ag_area_total)) / 1000
			self.zonewise_budgets[zone_name]['sm'] = sm_in_mm
			self.zonewise_budgets[zone_name]['sm_monsoon'] = sm_monsoon_in_mm
			self.zonewise_budgets[zone_name]['monsoon_deficit'] = monsoon_deficit_in_mm
			self.zonewise_budgets[zone_name]['deficit'] = deficit_in_mm

	def output_zonewise_budget_areawise_to_csv(self, zonewise_budget_areawise_csv_filename):
		csvwrite = open(self.path + zonewise_budget_areawise_csv_filename, 'wb')
		writer = csv.writer(csvwrite)
		writer.writerow([''] + ['zone-' + str(ID)	for ID in self.zonewise_budgets])
		writer.writerow(['Ag. Area'] + [self.zonewise_budgets[ID]['ag_area'] for ID in self.zonewise_budgets])
		writer.writerow(['Non-Ag. Area'] + [self.zonewise_budgets[ID]['non_ag_area'] for ID in self.zonewise_budgets])
		writer.writerow(['Rainfall in mm'] + [self.zonewise_budgets[ID]['rain_mm'] for ID in self.zonewise_budgets])
		writer.writerow(['Rainfall in TCM'] + [self.zonewise_budgets[ID]['rain_TCM'] for ID in self.zonewise_budgets])
		writer.writerow(['GW Recharge'] + [self.zonewise_budgets[ID]['gw_rech'] for ID in self.zonewise_budgets])
		writer.writerow(['Run-off'] + [self.zonewise_budgets[ID]['runoff'] for ID in self.zonewise_budgets])
		writer.writerow(['Usable SM'] + [self.zonewise_budgets[ID]['sm'] for ID in self.zonewise_budgets])
		writer.writerow(['Usable SM after Monsoon'] + [self.zonewise_budgets[ID]['sm_monsoon'] for ID in self.zonewise_budgets])
		writer.writerow(['Monsoon Deficit'] + [self.zonewise_budgets[ID]['monsoon_deficit'] for ID in self.zonewise_budgets])
		writer.writerow(['Deficit'] + [self.zonewise_budgets[ID]['deficit'] for ID in self.zonewise_budgets])
		csvwrite.close()

	def output_zonewise_budget_to_csv(self, zonewise_budget_csv_filename, PET_sum_cropend, PET_sum,rain_sum):
		csvwrite = open(self.path + zonewise_budget_csv_filename,'wb')
		writer = csv.writer(csvwrite)
		writer.writerow([''] + ['zone-'+str(ID)+'-'+t	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Rainfall'] + [rain_sum	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Runoff in Monsoon'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].runoff_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Infilitration in Monsoon'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].infil_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Soil Moisture Crop end'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].sm_crop_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Soil Moisture Monsoon end'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].sm_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['GW Recharge in Monsoon'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].GW_rech_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['AET Crop End'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].AET_crop_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['PET Crop End'] + [PET_sum_cropend[self.crop_name] if (ag_or_non_ag == 'agricultural' or ag_or_non_ag == 'zone' ) else PET_sum_cropend[other_LU_crops[t]] for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag] ])
		writer.writerow(['AET Monsoon End'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].AET_monsoon_end	for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['PET Monsoon End'] + [PET_sum[self.crop_name] if (ag_or_non_ag == 'agricultural' or ag_or_non_ag == 'zone' ) else PET_sum_cropend[other_LU_crops[t]] for ID in self.zonewise_budgets	for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag] ])
		writer.writerow(['Monsoon Deficit(PET-AET)'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].PET_minus_AET_monsoon_end	for ID in self.zonewise_budgets		for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		writer.writerow(['Crop duration Deficit(PET-AET)'] + [self.zonewise_budgets[ID][ag_or_non_ag][t].PET_minus_AET_crop_end	for ID in self.zonewise_budgets		for ag_or_non_ag in ['agricultural', 'non-agricultural', 'zone']	for t in self.zonewise_budgets[ID][ag_or_non_ag]])
		csvwrite.close()

	def compute_and_output_cadastral_vulnerability_to_csv(self, cadastral_vulnerability_csv_filename):
		plot_vulnerability_dict = {
			p.container_polygons[CADASTRAL_LABEL].id(): (p.budget.PET_minus_AET_crop_end,p.budget.PET_minus_AET_monsoon_end)
				for p in self.output_cadastral_points if dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()] in ['agriculture', 'fallow land']
		}
		sorted_keys = sorted(plot_vulnerability_dict.keys(), key=lambda ID:	plot_vulnerability_dict[ID], reverse=True)
		csvwrite = open(self.path + cadastral_vulnerability_csv_filename,'w+b')
		writer = csv.writer(csvwrite)
		writer.writerow(['Plot ID', 'Crop end Vulnerability', 'Crop end Deficit Waterings', 'Monsoon end Vulnerability', 'Crop end Deficit Waterings'])
		for key in sorted_keys:	writer.writerow([key, '{0:.2f}'.format(plot_vulnerability_dict[key][0]), round((plot_vulnerability_dict[key][0])/50), '{0:.2f}'.format(plot_vulnerability_dict[key][1]), round((plot_vulnerability_dict[key][1])/50)])
		csvwrite.close()

	def compute_and_display_cadastral_vulnerability(self):
		self.cadastral_points_per_plot = {}
		for p in (self.output_grid_points+self.output_cadastral_points):
			if p.container_polygons[CADASTRAL_LABEL] is None:	continue
			if (dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()]
                    not in ['agriculture', 'fallow land']): continue
			if p.container_polygons[CADASTRAL_LABEL].id() in self.cadastral_points_per_plot:
				self.cadastral_points_per_plot[p.container_polygons[CADASTRAL_LABEL].id()].append(p.budget.PET_minus_AET_crop_end)
			else:
				self.cadastral_points_per_plot[p.container_polygons[CADASTRAL_LABEL].id()] = [p.budget.PET_minus_AET_crop_end]
		for k,v in self.cadastral_points_per_plot.items():
			if len(v) > 0:  self.cadastral_points_per_plot[k] = sum(v)/len(v)
			else:           del self.cadastral_points_per_plot[k]

		#	Create duplicate cadastral layer in memory
		memory_cadastral_layer = QgsVectorLayer('Polygon?crs=epsg:32643', 'Cadastral Level Vulnerability', 'memory')
		memory_cadastral_layer.startEditing()
		memory_cadastral_layer.dataProvider().addAttributes(self.cadastral_layer.qgsLayer.dataProvider().fields().toList())
		memory_cadastral_layer.updateFields()
		dict_new_feature_id_to_old_feature_id = {}
		for old_plot_id in self.cadastral_points_per_plot:
			result, output_features = memory_cadastral_layer.dataProvider().addFeatures([self.cadastral_layer.feature_dict[old_plot_id]])
			dict_new_feature_id_to_old_feature_id[output_features[0].id()] = old_plot_id
		memory_cadastral_layer.dataProvider().addAttributes([QgsField('Deficit', QVariant.Double)])
		memory_cadastral_layer.updateFields()
		for new_feature in memory_cadastral_layer.getFeatures():
			new_feature['Deficit'] = self.cadastral_points_per_plot[dict_new_feature_id_to_old_feature_id[new_feature.id()]]
			memory_cadastral_layer.updateFeature(new_feature)
		memory_cadastral_layer.commitChanges()

		#	Graduated Rendering
		graduated_symbol_renderer_range_list = []
		ET_D_max = max(self.cadastral_points_per_plot.values())
		opacity = 1
		geometry_type = memory_cadastral_layer.geometryType()
		intervals_count = CADASTRAL_VULNERABILITY_DISPLAY_COLOUR_INTERVALS_COUNT
		dict_interval_colour = {0: [200, 200, 255], 1: [255, 0, 255], 2: [255, 165, 0], 3: [255, 0, 0]}
		for i in range(intervals_count):
			interval_min = 0 if i == 0 else (ET_D_max * float(i)/CADASTRAL_VULNERABILITY_DISPLAY_COLOUR_INTERVALS_COUNT) + 0.01
			interval_max = ET_D_max * float(i+1)/CADASTRAL_VULNERABILITY_DISPLAY_COLOUR_INTERVALS_COUNT
			label = "{0:.2f} - {1:.2f}".format(interval_min, interval_max)
			colour = QColor(*dict_interval_colour[i])
			symbol = QgsSymbolV2.defaultSymbol(geometry_type)
			symbol.setColor(colour)
			symbol.setAlpha(opacity)
			interval_range = QgsRendererRangeV2(interval_min, interval_max, symbol, label)
			graduated_symbol_renderer_range_list.append(interval_range)
		renderer = QgsGraduatedSymbolRendererV2('', graduated_symbol_renderer_range_list)
		renderer.setMode(QgsGraduatedSymbolRendererV2.EqualInterval)
		renderer.setClassAttribute('Deficit')
		memory_cadastral_layer.setRendererV2(renderer)

		QgsMapLayerRegistry.instance().addMapLayer(memory_cadastral_layer)
		memory_cadastral_layer.setCustomProperty('labeling', 'pal')
		memory_cadastral_layer.setCustomProperty('labeling/enabled', 'true')
		memory_cadastral_layer.setCustomProperty('labeling/fieldName', 'Number')
		memory_cadastral_layer.setCustomProperty('labeling/fontSize', '10')
		memory_cadastral_layer.setCustomProperty('labeling/placement', '0')
		QgsVectorFileWriter.writeAsVectorFormat(memory_cadastral_layer, self.path+'/kharif_cadastral_level_vulnerability.shp', "utf-8", None, "ESRI Shapefile")

	def calculate(self,
					crop_name,
					pointwise_output_csv_filename,
					zonewise_budget_csv_filename,
					zonewise_budget_csv_filename_LU,
			  		zonewise_budget_areawise_csv_filename,
					cadastral_vulnerability_csv_filename,
					sowing_threshold,
					start_date_index=0,
					end_date_index=182,
					input_points_filename=None
				):

		start_time = time.time()
		self.crop_name = crop_name
		self.pet_calculation(crop_name.lower(), sowing_threshold)
		PET_sum = dict((crop,sum(pet_values[start_date_index:123]) )  for crop,pet_values in self.PET.items())
		PET_sum_cropend = dict((crop,sum(pet_values[start_date_index:self.duration]) )  for crop,pet_values in self.PET.items())
		rain_sum = sum(self.rain[start_date_index:end_date_index+1])

		self.soil_types = dict_SoilContent.keys();	self.soil_types.remove('soil type')
		self.lulc_types = dict_RO.keys()

		self.output_grid_points = self.generate_output_points_grid(input_points_filename)
		self.filter_out_points_outside_boundary()
		self.set_container_polygon_of_points_for_layers(self.output_grid_points, [self.soil_layer, self.lulc_layer, self.cadastral_layer])
		self.set_slope_at_points(self.output_grid_points)
		self.output_grid_points = self.filter_out_points_with_incomplete_data(self.output_grid_points)
		for zone_id in self.zone_points_dict:
			self.zone_points_dict[zone_id] = self.filter_out_points_with_incomplete_data(self.zone_points_dict[zone_id])
		self.output_grid_points = filter(lambda p:	dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()] not in ['habitation', 'water'], self.output_grid_points)
		print 'Number of grid points to process : ', len(self.output_grid_points)
		for zone_id in self.zone_points_dict:
			self.zone_points_dict[zone_id] = filter(lambda p:	dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()] not in ['habitation', 'water'], self.zone_points_dict[zone_id])
		self.set_crop_at_points(self.output_grid_points,crop_name)

		i = 0
		
		for point in self.output_grid_points:
			point.run_model(self.rain, self.PET, PET_sum, PET_sum_cropend, start_date_index, end_date_index,self.duration-1)
			if point.budget.sm_crop_end > 100:	i += 1
		print 'no. of points with sm > 100 : ', i
		self.output_point_results_to_csv(pointwise_output_csv_filename)
		self.compute_zonewise_budget()
		self.output_zonewise_budget_to_csv_agri(zonewise_budget_csv_filename, PET_sum_cropend,PET_sum, rain_sum)
		self.compute_zonewise_budget_LU()
		self.output_zonewise_budget_to_csv(zonewise_budget_csv_filename_LU, PET_sum_cropend,PET_sum, rain_sum)
		self.compute_zonewise_budget_areawise()
		self.output_zonewise_budget_areawise_to_csv(zonewise_budget_areawise_csv_filename)

		self.filter_out_cadastral_plots_outside_boundary()
		#~ self.cadastral_points_dict,
		self.output_cadastral_points = self.generate_output_points_for_cadastral_plots()
		#~ print 'self.cadastral_layer.feature_dict : ', self.cadastral_layer.feature_dict
		self.set_container_polygon_of_points_for_layers(self.output_cadastral_points, [self.soil_layer, self.lulc_layer, self.cadastral_layer])
		self.set_slope_at_points(self.output_cadastral_points)
		self.output_cadastral_points = self.filter_out_points_with_incomplete_data(self.output_cadastral_points)
		self.output_cadastral_points = filter(lambda p:	dict_lulc[p.container_polygons[LULC_LABEL][Desc].lower()] not in ['habitation', 'water'], self.output_cadastral_points)
		print 'Number of cadastral points to process : ', len(self.output_cadastral_points)
		self.set_crop_at_points(self.output_cadastral_points,crop_name)
		for point in self.output_cadastral_points:
			point.run_model(self.rain, self.PET, PET_sum, PET_sum_cropend, start_date_index, end_date_index,self.duration-1)
		self.compute_and_output_cadastral_vulnerability_to_csv(cadastral_vulnerability_csv_filename)
		self.compute_and_display_cadastral_vulnerability()

		print("--- %s seconds ---" % (time.time() - start_time))
		print("done")
