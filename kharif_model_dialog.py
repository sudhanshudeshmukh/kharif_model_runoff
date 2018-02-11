
"""
/***************************************************************************
 KharifModelDialog
                                 A QGIS plugin
 Generates kharif season vulnerability map
                             -------------------
        begin                : 2017-11-18
        git sha              : $Format:%H$
        copyright            : (C) 2017 by IITB
        email                : sohoni@cse.iitb.ac.in
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os

from PyQt4 import QtGui, uic
#~ from PyQt4.QtCore import QString

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'kharif_model_dialog_base.ui'))


class KharifModelDialog(QtGui.QDialog, FORM_CLASS):
	def __init__(self, parent=None, crops=[]):
		"""Constructor."""
		super(KharifModelDialog, self).__init__(parent)
		# Set up the user interface from Designer.
		# After setupUI you can access any designer object by doing
		# self.<objectname>, and you can use autoconnect slots - see
		# http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
		# #widgets-and-dialogs-with-auto-connect
		self.setupUi(self)
		
		self.crop_combo_box.addItems(crops)
		
		self.last_path = ''
		
		self.folder_path_browse.clicked.connect(lambda : self.on_browse(self.folder_path, 'Folder containing the data-set', folder=True))
		self.zones_layer_browse.clicked.connect(lambda : self.on_browse(self.zones_layer_filename, 'Zones Vector Layer', 'Shapefiles (*.shp)'))
		self.soil_layer_browse.clicked.connect(lambda : self.on_browse(self.soil_layer_filename, 'Soil-cover Vector Layer', 'Shapefiles (*.shp)'))
		self.lulc_layer_browse.clicked.connect(lambda : self.on_browse(self.lulc_layer_filename, 'Land-use-land-cover Vector Layer', 'Shapefiles (*.shp)'))
		self.cadastral_layer_browse.clicked.connect(lambda : self.on_browse(self.cadastral_layer_filename, 'Cadastral Map Vector Layer', 'Shapefiles (*.shp)'))
		self.slope_layer_browse.clicked.connect(lambda : self.on_browse(self.slope_layer_filename, 'Slope Raster Layer', 'TIFF files (*.tif *.tiff)'))
		self.drainage_layer_browse.clicked.connect(lambda : self.on_browse(self.drainage_layer_filename, 'Drainage Vector Layer', 'Shapefiles (*.shp)'))
		self.rainfall_csv_browse.clicked.connect(lambda : self.on_browse(self.rainfall_csv_filename, 'Daily Rainfall CSV File', 'CSV files (*.csv)'))
		self.save_image_browse.clicked.connect(lambda : self.on_browse(self.save_image_filename, 'Save As Image In Folder', 'PNG files (*.png)', folder=True, save=True))
		
		self.colour_code_interval_points = [0, 100]
		self.colour_code_intervals_split_button.clicked.connect(self.on_split)
		self.colour_code_intervals_merge_button.clicked.connect(self.on_merge)
		self.colour_code_intervals_list_widget.addItem('0-100')
	
	def on_browse(self, lineEdit, caption, fltr='', folder=False, save=False):
		if folder:
			if save:
				path = QtGui.QFileDialog.getSaveFileName(self, caption, self.last_path, '.png')
			else:
				path = QtGui.QFileDialog.getExistingDirectory(self, caption, self.last_path)
				self.last_path = path
				self.autofill(path)
		else:
			path = QtGui.QFileDialog.getOpenFileName(self, caption, self.last_path, fltr)
		lineEdit.setText(path)
		if not self.folder_path.text():	self.last_path = os.path.dirname(path)
	
	def autofill(self, path):
		inputfiles_lineEdit_dict = {
			'Zones.shp': self.zones_layer_filename,
			'Soil.shp': self.soil_layer_filename,
			'LULC.shp': self.lulc_layer_filename,
			'Cadastral.shp': self.cadastral_layer_filename,
			'Slope.tif': self.slope_layer_filename,
			'Drainage.shp': self.drainage_layer_filename,
			'Rainfall.csv': self.rainfall_csv_filename
		}
		for inputfile in inputfiles_lineEdit_dict:
			if os.path.exists(path + '/' + inputfile):
				inputfiles_lineEdit_dict[inputfile].setText(path + '/' + inputfile)
	
	def on_split(self):
		split_at = self.colour_code_intervals_split_value_spin_box.value()
		if split_at not in self.colour_code_interval_points:
			i = 0
			while i < len(self.colour_code_interval_points) and self.colour_code_interval_points[i] < split_at:
				i += 1
			self.colour_code_interval_points.insert(i, split_at)
			self.colour_code_intervals_list_widget.takeItem(i-1)
			self.colour_code_intervals_list_widget.insertItem(i-1, str(self.colour_code_interval_points[i-1])+'-'+str(self.colour_code_interval_points[i]))
			self.colour_code_intervals_list_widget.insertItem(i, str(self.colour_code_interval_points[i])+'-'+str(self.colour_code_interval_points[i+1]))

	def on_merge(self):
		selection = self.colour_code_intervals_list_widget.currentRow()
		if selection > 0:
			self.colour_code_intervals_list_widget.takeItem(selection-1)
			self.colour_code_intervals_list_widget.takeItem(selection-1)
			#~ print self.colour_code_interval_points, selection
			del self.colour_code_interval_points[selection]
			self.colour_code_intervals_list_widget.insertItem(selection-1, str(self.colour_code_interval_points[selection-1])+'-'+str(self.colour_code_interval_points[selection]))
