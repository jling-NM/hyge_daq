#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2017 josef ling <jling@mrn.org>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
import sys
import os
import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
from lib.plotarea import PlotArea
from lib.plotarea import PlotClearDepth
from lib.experiment import Experiment
from lib.experiment import SensorRotation


# Run program even without NI driver for daq
try:
    import PyDAQmx as pmx
except ImportError as e:
    print(e)
    print("Continuing without DAQ driver")



# sensor map from id to radianConvFactor
SENSOR_MAP = {"0390": 0.0481, "0392": 0.0464}


class GUI(QtWidgets.QMainWindow):
        
    def __init__(self):
        super().__init__()
        
        # init a default experiment
        self.experiment = Experiment()
                
        # get app settings
        self.read_app_settings()
        
        self.main_frame = QtWidgets.QWidget(self)
        self.plot_area = None
        
        self.init_gui()


    def init_gui(self):
        """
        QT GUI
        """
        # create menus
        exitAction = QtWidgets.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtWidgets.qApp.quit)

        waitTraceAction = QtWidgets.QAction('Wait For &Trigger', self)
        waitTraceAction.setShortcut('Ctrl+T')
        waitTraceAction.setStatusTip('Waits for trigger and collects data')
        waitTraceAction.triggered.connect(self.get_experiment_params)
        
        clearTraceAction = QtWidgets.QAction('&Clear Traces', self)
        clearTraceAction.setStatusTip('Clear all traces from plot')
        clearTraceAction.triggered.connect(self.confirm_clear)

        openFileAction = QtWidgets.QAction('&Open Trace File', self)
        openFileAction.setShortcut('Ctrl+O')
        openFileAction.setStatusTip('Load Trace Data File')
        openFileAction.triggered.connect(self.load_trace)

        sensorSettingsAction = QtWidgets.QAction('Sensor Settings', self)
        sensorSettingsAction.setStatusTip('Update Sensor Settings')
        sensorSettingsAction.triggered.connect(self.get_sensor_params)


        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(waitTraceAction)
        fileMenu.addAction(openFileAction)
        fileMenu.addAction(clearTraceAction)
        fileMenu.addAction(sensorSettingsAction)
        fileMenu.addAction(exitAction)

        # options menu
        optMenu = menubar.addMenu('&Options')
        
        # overlay option submenu
        self.overlayMenu = optMenu.addMenu('New Trace:')
        
        # group so options are exclusive
        ag = QtWidgets.QActionGroup(self.overlayMenu, exclusive=True)
        # add menu items
        a = ag.addAction(QtWidgets.QAction('Overlays Existing', self.overlayMenu, checkable=True))
        a.setData(PlotClearDepth.NON)
        if self.plot_clear_level == PlotClearDepth.NON:
            a.setChecked(True)
        self.overlayMenu.addAction(a)
        
        a = ag.addAction(QtWidgets.QAction('Replaces Top', self.overlayMenu, checkable=True))
        a.setData(PlotClearDepth.TOP)
        if self.plot_clear_level == PlotClearDepth.TOP:
            a.setChecked(True)        
        self.overlayMenu.addAction(a)
        
        a = ag.addAction(QtWidgets.QAction('Replaces All', self.overlayMenu, checkable=True))
        a.setData(PlotClearDepth.ALL)
        if self.plot_clear_level == PlotClearDepth.ALL:
            a.setChecked(True)        
        self.overlayMenu.addAction(a)        
        self.overlayMenu.triggered.connect(self.overlayMenu_changed)

        # plot filter submenu
        self.plotFilterMenu = optMenu.addMenu('Filter Trace Plot:')
        
        # group so options are exclusive
        ag = QtWidgets.QActionGroup(self.plotFilterMenu, exclusive=True)
        # add menu items
        a = ag.addAction(QtWidgets.QAction('Yes', self.plotFilterMenu, checkable=True))
        a.setData(True)
        if self.plot_filter_data:
            a.setChecked(True)
        self.plotFilterMenu.addAction(a)
        
        a = ag.addAction(QtWidgets.QAction('No', self.plotFilterMenu, checkable=True))
        a.setData(False)
        if not self.plot_filter_data:
            a.setChecked(True)
        self.plotFilterMenu.addAction(a)
        self.plotFilterMenu.triggered.connect(self.plotFilterMenu_changed)

        
        
        # plot filter submenu
        self.plotAnnotationMenu = optMenu.addMenu('Annotate Plot:')
        
        # group so options are exclusive
        ag = QtWidgets.QActionGroup(self.plotAnnotationMenu, exclusive=True)
        # add menu items
        a = ag.addAction(QtWidgets.QAction('Yes', self.plotAnnotationMenu, checkable=True))
        a.setData(True)
        if self.plot_annotate:
            a.setChecked(True)
        self.plotAnnotationMenu.addAction(a)
        
        a = ag.addAction(QtWidgets.QAction('No', self.plotAnnotationMenu, checkable=True))
        a.setData(False)
        if not self.plot_annotate:
            a.setChecked(True)
        self.plotAnnotationMenu.addAction(a)
        self.plotAnnotationMenu.triggered.connect(self.plotAnnotationMenu_changed)


        
        # create status bar
        self.statusBar()

        self.setWindowTitle('HYGE DAQ')
        self.setWindowIcon(QtGui.QIcon('rc/appicon.png'))
        self.statusBar().showMessage('Ready')

        self.plot_area = PlotArea(self.main_frame)
        self.main_frame.setLayout(self.plot_area)
        self.setCentralWidget(self.main_frame)


    def get_experiment_params(self):
        """
        Query experiment variables from user; 
        update experiment instance;
        chain to wait for trigger
        """
        def frmAccept():
            
            # make sure we have data
            if (len(inFldId.text()) != 0) and (len(inFldPsiLoad.text()) != 0) and (len(inFldPsiSet.text()) != 0):
                
                # start new experiment
                self.experiment = Experiment()
                self.experiment.subjectId = str(inFldId.text())
                self.experiment.PsiLoad = str(inFldPsiLoad.text())
                self.experiment.PsiSet = str(inFldPsiSet.text())
                
                # close user input dialog
                experimentDlg.close()
                
                # go to trigger wait
                self.wait_trigger_collection()
                                
            
        def frmCancel():
            experimentDlg.close()


        # see if we can collect an experiment first; perhaps better place to put this?
        if 'PyDAQmx' not in sys.modules:
            self.display_msg("No NI Driver", "You cannot collect data.", "Data collection requires PyDAQmx")
            
            
        experimentDlg = QtWidgets.QDialog(self)
        experimentDlg.setWindowTitle("Experiment Settings")
        experimentDlg.setWindowModality(QtCore.Qt.ApplicationModal)
        
        # input form
        inLblId = QtWidgets.QLabel("Subject ID:")
        inFldId = QtWidgets.QLineEdit()
        #inFldId.setInputMask('M99999999')   # input validation
        inFldId.setText(self.experiment.subjectId)
        #inFldId.setCursorPosition(len(self.experiment.subjectId)) # input validation, cursor position
 
        inLblPsiLoad = QtWidgets.QLabel("Load PSI:")
        inFldPsiLoad = QtWidgets.QLineEdit()
        inFldPsiLoad.setValidator(QtGui.QIntValidator())
        inFldPsiLoad.setMaxLength(4)
        inFldPsiLoad.setText(self.experiment.PsiLoad)
        
        inLblPsiSet = QtWidgets.QLabel("Set PSI:")
        inFldPsiSet = QtWidgets.QLineEdit()
        inFldPsiSet.setValidator(QtGui.QIntValidator())
        inFldPsiSet.setMaxLength(3)
        inFldPsiSet.setText(self.experiment.PsiSet)
        
        btnBox = QtWidgets.QDialogButtonBox()
        btnBox.addButton("Continue", QtWidgets.QDialogButtonBox.AcceptRole)
        btnBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel)
        btnBox.accepted.connect(frmAccept)
        btnBox.rejected.connect(frmCancel)
        btnBox.centerButtons()
        
        dlgLayout = QtWidgets.QFormLayout()
        dlgLayout.addRow(inLblId, inFldId)
        dlgLayout.addRow(inLblPsiSet, inFldPsiSet)
        dlgLayout.addRow(inLblPsiLoad, inFldPsiLoad)
        dlgLayout.addRow(btnBox)
        
        experimentDlg.setLayout(dlgLayout)
        experimentDlg.exec_()



    def get_sensor_params(self):
        """
        Query sensor variables from user; 
        update gui instance;
        """
        def frmAccept():
            
            # update parameters
            self.sensor_channel1_id = cbSensor1Id.currentText()
            self.sensor_channel1_rotation = cbSensor1Rot.currentData()
            self.sensor_channel2_id = cbSensor2Id.currentText()
            self.sensor_channel2_rotation = cbSensor2Rot.currentData()
            # close user input dialog
            sensorDlg.close()
                                
            
        def frmCancel():
            sensorDlg.close()

        sensorDlg = QtWidgets.QDialog(self)
        sensorDlg.setWindowTitle("Sensor Settings")
        sensorDlg.setWindowModality(QtCore.Qt.ApplicationModal)
        
        # input form
        lblSensor1 = QtWidgets.QLabel("Sensor 1")
        lblSensor1Id = QtWidgets.QLabel("ID:")
        cbSensor1Id = QtWidgets.QComboBox()
        
        for i, item in enumerate(SENSOR_MAP):
            cbSensor1Id.addItem(item)
            if self.sensor_channel1_id == item:
                cbSensor1Id.setCurrentIndex(i)
            
        lblSensor1Rot = QtWidgets.QLabel("Rotation:")
        cbSensor1Rot = QtWidgets.QComboBox()
        cbSensor1Rot.addItem("CounterClockwise", SensorRotation.COUNTERCLOCKWISE)            
        cbSensor1Rot.addItem("Clockwise", SensorRotation.CLOCKWISE)
        if self.sensor_channel1_rotation == SensorRotation.COUNTERCLOCKWISE:
            cbSensor1Rot.setCurrentIndex(0)
        else:
            cbSensor1Rot.setCurrentIndex(1)        
            
        lblSensor2 = QtWidgets.QLabel("Sensor 2")
        lblSensor2Id = QtWidgets.QLabel("ID:")
        cbSensor2Id = QtWidgets.QComboBox()
        
        for i, item in enumerate(SENSOR_MAP):
            cbSensor2Id.addItem(item)
            if self.sensor_channel2_id == item:
                cbSensor2Id.setCurrentIndex(i)
            
        lblSensor2Rot = QtWidgets.QLabel("Rotation:")
        cbSensor2Rot = QtWidgets.QComboBox()
        cbSensor2Rot.addItem("CounterClockwise", SensorRotation.COUNTERCLOCKWISE)
        cbSensor2Rot.addItem("Clockwise", SensorRotation.CLOCKWISE)
        if self.sensor_channel2_rotation == SensorRotation.COUNTERCLOCKWISE:
            cbSensor2Rot.setCurrentIndex(0)
        else:
            cbSensor2Rot.setCurrentIndex(1)
            
        
        btnBox = QtWidgets.QDialogButtonBox()
        btnBox.addButton("Continue", QtWidgets.QDialogButtonBox.AcceptRole)
        btnBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel)
        btnBox.accepted.connect(frmAccept)
        btnBox.rejected.connect(frmCancel)
        btnBox.centerButtons()
        
        dlgLayout = QtWidgets.QFormLayout()
        dlgLayout.addRow(lblSensor1)
        dlgLayout.addRow(lblSensor1Id, cbSensor1Id)
        dlgLayout.addRow(lblSensor1Rot, cbSensor1Rot)
        dlgLayout.addRow(QtWidgets.QLabel(""))
        dlgLayout.addRow(lblSensor2)
        dlgLayout.addRow(lblSensor2Id, cbSensor2Id)
        dlgLayout.addRow(lblSensor2Rot, cbSensor2Rot)
        dlgLayout.addRow(QtWidgets.QLabel(""))
        dlgLayout.addRow(btnBox)
        
        sensorDlg.setLayout(dlgLayout)
        sensorDlg.exec_()
        
        
    def display_msg(self, txt, iTxt, dTxt):
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText(txt)
        msg.setInformativeText(iTxt)
        msg.setWindowTitle('HYGE DAQ')
        msg.setDetailedText(dTxt)
        msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msg.exec_()
        

    def confirm_clear(self, event):
        """
        Handle request to clear all traces and text from plot
        """
        reply = QtWidgets.QMessageBox.question(self, 'Confirm', "Clear plot traces?",
                                               QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                               QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            self.plot_area.clear_plot(PlotClearDepth.ALL)



    def overlayMenu_changed(self): 
        """ 
        Connect to Option menu for changing how plot overlay works
        """ 
        for action in self.overlayMenu.actions():
            if action.isChecked():
                self.plot_clear_level = action.data()
            
        

    def plotFilterMenu_changed(self): 
        """ 
        Connect to Option menu for sensor
        """ 
        for action in self.plotFilterMenu.actions():
            if action.isChecked():
                self.plot_filter_data = action.data()

    def plotAnnotationMenu_changed(self): 
        """ 
        Connect to Option menu for sensor
        """ 
        for action in self.plotAnnotationMenu.actions():
            if action.isChecked():
                self.plot_annotate = action.data()
                
                
    def wait_trigger_collection(self):
        """ 
        init wait for data            
        
        1) Create channel reader
        2) scale voltage to account for voltage divider. Now working with scaled voltage like what is loaded from hyg file
        
        3) Load scaled voltage channel0 data starting after time offset
        4) calculate baseline shift in small noise window  
        5) subtract baseline shift from voltage data and convert to rad/s
        6) get summary parameters and plot all

        """
        
        # clear plot before waiting
        self.plot_area.clear_plot(self.plot_clear_level)
        
        # update experiment for correct header settings in case user changed them prior to capturing data
        self.experiment.channel1_radianConvFactor = SENSOR_MAP[self.sensor_channel1_id]
        self.experiment.channel2_radianConvFactor = SENSOR_MAP[self.sensor_channel2_id]

    
        try:
            
            # init task
            analog_input = pmx.Task()
            # define task datatype
            read = pmx.int32()
            
            #single channel
                #data = np.zeros((self.params.readNumSamples,), dtype=np.float64)
            # two channels - init storate
            sensor_data_raw = np.zeros((self.experiment.readNumSamples*2,), dtype=np.float64)
            
            try:
                # single channel
                    #analog_input.CreateAIVoltageChan("Dev1/ai0:1","Sensor1",pmx.DAQmx_Val_RSE,-10.0,10.0,pmx.DAQmx_Val_Volts,None)
                    #analog_input.CfgSampClkTiming("",self.params.readRate,pmx.DAQmx_Val_Rising,pmx.DAQmx_Val_FiniteSamps,self.params.readNumSamples)
                
                # two channel - init channel
                analog_input.CreateAIVoltageChan("Dev1/ai0:1","",pmx.DAQmx_Val_RSE,-10.0,10.0,pmx.DAQmx_Val_Volts,None)
                            
                # init sampling clock
                analog_input.CfgSampClkTiming("", self.experiment.readRate, pmx.DAQmx_Val_Rising,pmx.DAQmx_Val_FiniteSamps, self.experiment.readNumSamples)
                
                # digital edge detect to wait for trigger signal
                # TRIGGER                
                analog_input.CfgDigEdgeStartTrig("/Dev1/PFI1",pmx.DAQmx_Val_Rising) 
                
                # DAQmx Start task and therefore wait...
                analog_input.StartTask()
                
                # message to user in status bar...
                self.statusBar().showMessage('Waiting for Trigger...')
                
                # DAQmx Read
                    # single channel
                        #analog_input.ReadAnalogF64(self.params.readNumSamples,pmx.DAQmx_Val_WaitInfinitely,pmx.DAQmx_Val_GroupByChannel,sensor_data_raw,self.params.readNumSamples,pmx.byref(read), None)        
                    # two channel
                        # group data by channel
                            #analog_input.ReadAnalogF64(25000,10.0,pmx.DAQmx_Val_GroupByChannel,sensor_data_raw,50000,pmx.byref(read), None)        
                # two channel interleaved data
                analog_input.ReadAnalogF64(self.experiment.readNumSamples,pmx.DAQmx_Val_WaitInfinitely,pmx.DAQmx_Val_GroupByScanNumber,sensor_data_raw,self.experiment.readNumSamples*2,pmx.byref(read), None)
                
            except pmx.DAQError as err:
                #raise( "\n\nDAQmx Error:\n{}\n\n".format(err) )
                raise(err)
                
            finally:
                if analog_input:
                    # task clean up to return to GUI
                    analog_input.StopTask()
                    analog_input.ClearTask()
            
            



            # if the data is just noise from not having anything hooked to the device, we just plot it without all the bells and whistles of the summary calculations
            # this check could be moved into getSummaryData()
            if np.max(abs(sensor_data_raw)) > 0.01:
                
                # flip the sign according to sensor mount location so data is positive
                if self.sensor_channel1_rotation == SensorRotation.COUNTERCLOCKWISE:
                    voltage_scalar = self.experiment.voltageScalingFactor
                elif self.sensor_channel1_rotation == SensorRotation.CLOCKWISE:
                    voltage_scalar = -1 * self.experiment.voltageScalingFactor
                
                # message to user
                self.statusBar().showMessage('Saving...', 2000)
                # scale voltage to account for voltage divider
                sensor_data_voltage_scaled = sensor_data_raw * voltage_scalar           
                ### interleave data ###                
                # 1) get only the first channel 'channel 0'
                self.experiment.voltage_data = sensor_data_voltage_scaled[0::2]
                # first save data to file in case of a program error later. This slows plotting but protects data.
                self.save_trace()
                                
                # plot
                self.plot_area.plot(self.experiment, self.plot_filter_data, self.plot_annotate)
                
                del sensor_data_raw
                
            else:
                
                # NO GOOD DATA; just plot what you got
                channel0_data = sensor_data_raw[(self.experiment.triggerOffset * 2)::2]
                self.plot_area.plot(self.experiment, channel0_data)
                del sensor_data_raw

                
        except Exception as e:
            self.display_msg("Error:", "Waiting for Trigger", str(e))
            return


    def load_trace(self):
        """
        Read hyge data file and display in plot
        
        1) Read Header
        2) Load scaled voltage channel0 data starting after time offset
        3) plot
        """
        try:
            
            fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file',
                                                             self.experiment.lastDataPath, "Trace Files (*.hyg)")
            if fname:

                # update experiment parameters with header from file being loaded
                self.experiment = Experiment.load(fname)
                
                # if the data is just noise from not having anything hooked to the device,
                # we just plot it without all the bells and whistles of the summary calculations
                # this just avoids throwing errors on summary calculations for invalid data
                # and could be removed from the code
                if np.max(abs(self.experiment.voltage_data)) > 0.01:

                    # clear the plot
                    self.plot_area.clear_plot(self.plot_clear_level)
                
                    # plot data
                    self.plot_area.plot(self.experiment, self.plot_filter_data, self.plot_annotate)
                    self.statusBar().showMessage('Ready')

                else:
                    # plot null data
                    self.plot_area.plot(self.experiment, self.experiment.voltage_data)
                    self.statusBar().showMessage('Ready')

                
        except Exception as e:
            self.display_msg("Error:", "Loading Trace file", str(e))
            return


    def save_trace(self):
        """ 
        write out trace data in voltage 
        """
        self.experiment.save(os.path.join(os.getcwd(), "data", self.experiment.file_name))
        
        
    def read_app_settings(self):
        """
        read app session settings
        """
        self.settings = QtCore.QSettings()
        
        # sensor settings
        self.sensor_channel1_id = self.settings.value('sensor_channel1_id', list(SENSOR_MAP.keys())[0])
        self.sensor_channel2_id = self.settings.value('sensor_channel2_id', list(SENSOR_MAP.keys())[1])
        self.sensor_channel1_rotation = self.settings.value('sensor_channel1_rotation', SensorRotation.COUNTERCLOCKWISE)
        self.sensor_channel2_rotation = self.settings.value('sensor_channel2_rotation', SensorRotation.CLOCKWISE)
        
        # plot settings
        # default plot overlay setting
        self.plot_clear_level = self.settings.value('plot_clear_level', PlotClearDepth.TOP)
        # default plot filter
        self.plot_filter_data = self.settings.value('plot_filter_data', False, type=bool)
        # default plot annotation
        self.plot_annotate = self.settings.value('plot_annotate', False, type=bool)
        
        # change to script directory
        script_home = os.path.dirname(os.path.abspath(__file__))
        os.chdir(script_home)
        # point experiment output to the data subdirectory
        self.experiment.lastDataPath = self.settings.value('lastDataPath', os.path.join(script_home, 'data'))
        
                
    def save_app_settings(self):
        """
        save app session settings
        """
        # update settings
        self.settings.setValue('plot_filter_data', self.plot_filter_data)
        self.settings.setValue('plot_annotate', self.plot_annotate)
        self.settings.setValue('plot_clear_level', self.plot_clear_level)
        self.settings.setValue('sensor_channel1_id', self.sensor_channel1_id)
        self.settings.setValue('sensor_channel2_id', self.sensor_channel2_id)
        self.settings.setValue('sensor_channel1_rotation', self.sensor_channel1_rotation)
        self.settings.setValue('sensor_channel2_rotation', self.sensor_channel2_rotation)
        self.settings.setValue('lastDataPath', self.experiment.lastDataPath)
        # this writes to native storage
        del self.settings
        
    @QtCore.pyqtSlot()
    def close(self):
        """
        Anything you want to do before we close gui?
        :return: 
        """
        self.save_app_settings()
        super().close()

def main():

    try:
        app = QtWidgets.QApplication(sys.argv)
        # set application values once here for QSettings
        app.setOrganizationName("MayerLab")
        app.setApplicationName("HygeDAQ")
        
        gui = GUI()
        # exit menu action calls app quit. Notify mainwindow gui that we are quiting.
        app.aboutToQuit.connect(gui.close)
        gui.showFullScreen()
        gui.showMaximized()
        
        

    except Exception as e:
        print(e)

    finally:
        sys.exit(app.exec_())



if __name__ == '__main__':
    main()
