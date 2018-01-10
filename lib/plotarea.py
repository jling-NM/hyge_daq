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
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
from enum import Enum
from matplotlib.widgets import SpanSelector
#import matplotlib.markers


class PlotClearDepth(Enum):
    """
    plot clearing option
    """
    NON = 0 # don't clear anything; overlay
    ALL = 1 # clear everything
    TOP = 2 # clear the last one if more than one
    
    

class PlotArea(QtWidgets.QVBoxLayout):

    """
    Plots an Experiment
    """
    def __init__(self, parent=None):
        super(PlotArea, self).__init__()

        self.fig = Figure((5.0, 4.0), dpi=100)
        self.fig.set_tight_layout(True)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(parent)

        # Create the navigation toolbar, tied to the canvas
        mpl_toolbar = NavigationToolbar(self.canvas, parent)

        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)

        self.axes.set_xlabel('Time(ms)')
        # set xlim to 1.0 to 2.0 seconds post trigger
#        self.axes.set_xlim(0, 12500) # display: 500ms
        self.axes.set_xlim(0, 25000) # display: 1000ms
        start, end = self.axes.get_xlim()
#        self.axes.xaxis.set_ticks(np.arange(start, end, 500)) # display: 500ms
        self.axes.xaxis.set_ticks(np.arange(start, end, 1000)) # display: 1000ms
#        self.axes.xaxis.set_ticklabels(np.arange(0, 500, 20)) # display: 500ms
        self.axes.xaxis.set_ticklabels(np.arange(0, 1000, 40)) # display: 1000ms

        self.axes.set_ylabel('Velocity (rad/s)')
        self.axes.set_ylim([-150.0, 350.0])
        start, end = self.axes.get_ylim()
        self.axes.yaxis.set_ticks(np.arange(start, end+50, 50))
        self.axes.yaxis.set_ticklabels(np.arange(start, end+50, 50))
        # add gridlines to plot
        self.axes.grid()

        self.addWidget(self.canvas)
        self.addWidget(mpl_toolbar)

        # initialize summarybox and provide reference
        self.summaryBox = None
        self.selected_xspan = None
        self.positive_phase_area = None
        self.axes.leg = None
        
        # init plot history
        self.reset_history()

        # store any annotations for easier clearing
        self.annotation_list = []
        
        # keep track of this for entire plot area
        self.display_annotations = False
            
            
    def plot(self, experiment, plot_filter_data, plot_annotate):
        """
        Plot experiment radian data
        """

        # update this setting
        self.display_annotations = plot_annotate
        
        
        # see if we have a baseline shift and remove
        baselineShift = np.median(experiment.voltage_data[0:int(experiment.readRate)]) 
        if np.abs(baselineShift) < 1.0:
           experiment.voltage_data = experiment.voltage_data - baselineShift


        # scale data from volts to rad/s
        radian_data = experiment.get_radian_data()
        
        # filter data
        filtered_radian_data = experiment.filter_radian_data(radian_data, filter_type='binomial')


        # get summary metrics from radian data
        try:
            
            # generate summary data for display
            summary_data = self.getSummaryVals(experiment, filtered_radian_data)
    
            # TEXT BOX for summary parameters using FULL data set instead of the shortened plotting data set
            self.draw_summary_box(experiment, summary_data)

        except Exception as e:
            print(e)


        # which to plot and truncate
        if plot_filter_data:
            plot_data = filtered_radian_data[experiment.triggerOffset:]
        else:
            plot_data = radian_data[experiment.triggerOffset:]

        # resample plotting data if necessary to match current underlay
        if self.current_sample_rate == 0.0:
            self.current_sample_rate = int(experiment.readRate)            
        elif int(experiment.readRate) != self.current_sample_rate:
            plot_data = experiment.resampled_data(plot_data, self.current_sample_rate)

            
        # following selection of plotting data, resampling center plots on peaks in data to be displayed
        shift_plot_pos = 0
        if self.underlay_peak_index == 0:
            self.underlay_peak_index = np.argmax(plot_data)
        else:
            shift_plot_pos = np.argmax(plot_data) - self.underlay_peak_index
            if shift_plot_pos < 0:
                # zeropad front, clip from end
                plot_data = np.concatenate((np.zeros(np.abs(shift_plot_pos)), plot_data[0:shift_plot_pos]))
            else:
                # cut from front, zeropad end
                plot_data = np.concatenate((plot_data[shift_plot_pos:], np.zeros(np.abs(shift_plot_pos))))            
        
        
        
        self.axes.set_xlim(0, self.current_sample_rate) # display: 1000ms
        start, end = self.axes.get_xlim()
        # create 25(includes 0) ticks for 1 sec
        self.axes.xaxis.set_ticks(np.arange(start, end, int(self.current_sample_rate/25))) # display: 1000ms
        # label those 25 ticks for 1 sec
        self.axes.xaxis.set_ticklabels(np.arange(0, 1000, 40)) # display: 1000ms        
                          
        # line plot of data
        self.axes.plot(plot_data, label=experiment.get_label())
        
        
        
        # shall we add annotations?
        if self.display_annotations:
            
            peak_indx = (summary_data['peak_ind']+(-shift_plot_pos))
            rise_start_indx = (summary_data['rise_start_index']+(-shift_plot_pos))
            rise_end_indx = (summary_data['rise_end_index']+(-shift_plot_pos))
            
            #marker=matplotlib.markers.CARETUP, 
            self.axes.plot(
                    [rise_start_indx, peak_indx, rise_end_indx], 
                    [plot_data[rise_start_indx], plot_data[peak_indx], plot_data[rise_end_indx]], 
                    's', 
                    marker='.',
                    markersize='5',
                    color="red")
            #self.axes.scatter([summary_data['rise_start_index'], summary_data['rise_end_index']], [plot_data[summary_data['rise_start_index']], plot_data[summary_data['rise_end_index']]])
    
            self.annotation_list.append(self.axes.annotate(
                    'peak', 
                    xy=(peak_indx, plot_data[peak_indx]),
                    size=10,
                    xycoords='data', 
                    xytext=(-50,10), 
                    textcoords='offset points', 
                    arrowprops=dict(arrowstyle="->", connectionstyle="angle, angleA=0, angleB=90, rad=10")
                )
            )
            
            self.annotation_list.append(self.axes.annotate(
                'rise start', 
                xy=(rise_start_indx, plot_data[rise_start_indx]), 
                size=10,
                xycoords='data', 
                xytext=(-50,30), 
                textcoords='offset points', 
                arrowprops=dict(arrowstyle="->")
                )
            )
            
            self.annotation_list.append(self.axes.annotate(
                'rise end', 
                xy=(rise_end_indx, plot_data[rise_end_indx]), 
                size=10,
                xycoords='data', 
                xytext=(10,30), 
                textcoords='offset points', 
                arrowprops=dict(arrowstyle="->")
                )
            )
              
        
        
        
        
        # update previous titles
        self.axes.set_title(experiment.get_label())
        if len(self.axes.lines) > 1:
            self.axes.leg = self.axes.legend(loc=3, fancybox=True, framealpha=0.5)

        # feature for selecting an x range and then firing "onselect" method
        self.span = SpanSelector(self.axes, self.on_select_xspan, 'horizontal', useblit=True,
                            rectprops=dict(alpha=0.5, facecolor='red'))


        # positive_phase_area indicator
        # only show on single plot without annotations
        if len(self.axes.lines) == 1:
            self.positive_phase_area = self.axes.fill_between(np.arange(summary_data['rise_start_index'],
                                                              summary_data['rise_end_index']), 0,
                                                              plot_data[summary_data['rise_start_index']:summary_data['rise_end_index']],
                                                              facecolor='beige', alpha=0.50)

        # refresh canvas so plot is updated
        self.canvas.draw()

        del filtered_radian_data, radian_data


    def on_select_xspan(self, xmin, xmax):
        """
        Called after user selects an x range in plot
        just visual right now, doesn't recalculate anything
        :param xmin: 
        :param xmax: 
        :return: 
        """
        if self.selected_xspan is not None:
            self.selected_xspan.remove()
            self.selected_xspan = None

        self.selected_xspan = self.axes.axvspan(xmin, xmax, facecolor='b', alpha=0.05)


    def clear_plot(self, clear_depth=0):
        """ clear the plot """
        
        # remove plot lines
        if clear_depth == PlotClearDepth.ALL:
            # remove all plots
            del self.axes.lines[:]
            # reset history when removing all lines
            self.reset_history()
            
        elif clear_depth == PlotClearDepth.TOP:
            # remove last plot if multiple      
            if self.display_annotations:
                if (len(self.axes.lines) > 3):
                    del self.axes.lines[2:] # [2:] to remove line but leave markers, [1:] to remove line and markers
            else:
                if (len(self.axes.lines) > 1):
                    del self.axes.lines[-1]


                    
        # remove title
        self.axes.set_title("")          
            
        # remove text box and reset reference
        if self.summaryBox is not None:
            self.summaryBox.remove()
            self.summaryBox = None

        # remove xspan and reset reference
        if self.selected_xspan is not None:
            self.selected_xspan.remove()
            self.selected_xspan = None

        # remove positive phase and reset reference
        if self.positive_phase_area is not None:
            self.positive_phase_area.remove()
            self.positive_phase_area = None

        # clean up legend if it is initialized
        if self.axes.leg is not None:
            self.axes.leg.remove()
            self.axes.leg = None  
            
        # remove annotations
        for ann in self.annotation_list:
            ann.remove()
        self.annotation_list = []
                
                
            # refresh canvas
        self.canvas.draw()



    def reset_history(self):
        """ clear plot history """
        # init plot history
        self.underlay_peak_index = 0
        self.current_sample_rate = 0    
        
        
        
    def getSummaryVals(self, experiment, filtered_radian_data):
        """
        Calculate some summary values about our trace to overlay on the plot
        
        noisefloor : max positive value of first second of data (1 sec before fire) after removing baseline
        peak_ind : position in data with max value after subtracting baseline and after smoothing
        trough_ind : position in data with minimum value after subtracting baseline and after smoothing calculated from peak_ind to end of data
        rise_start_index : first position in data moving left from 100 ms left of peak to the right that exceeds noise floor
        rise_end_index : first position from peak_ind moving left to right that exceeds noise floor
        peak_accel_ind : position or maximum acceleration calculated from filtered data converted to acceleration(derivative or velocity) left to right from rise_start_index through peak_ind
        peak_accel : acceleration value at peak_accel_ind
        peak_decel_ind :position or minimum acceleration calculated from filtered data converted to acceleration(derivative or velocity) left to right from peak_ind through trough_ind]
        peak_decel : acceleration value at peak_accel_ind
        delta_t : cumulative time in milleseconds from rise_end_index through rise_start_index which is function of experiment.readRate
        time_to_peak : cumulative time in milleseconds from peak_ind through rise_start_index which is function of experiment.readRate
        decel_time : cumulative time in milleseconds from peak_ind through rise_end_index which is function of experiment.readRate
        peak_vel : value at peak_ind but from UNSMOOTHED radian data
        excursionAngleRadians : sum of filtered radian data from rise_start_index to rise_end_index multiplied by time of each step
        excursionAngleDegrees : np.degrees(excursionAngleRadians)
        excursionTravelDistanceInches : excursionAngleRadians converted to length based on experiment.swingArmToPivotInches
        
        """
            
        
        # get noise floor from radian data 1 sec waiting for fire
        noise_floor = np.max(filtered_radian_data[0:int(experiment.readRate)]) * 3
        
        #print("noisefloor:{}".format(noise_floor))
        #print("triggeroffset:{}".format(experiment.triggerOffset))
                    
        # don't need noise floor and want to return proper indices for plotting so truncate now
        # easier to report back proper indices for plotting
        filtered_radian_data = filtered_radian_data[experiment.triggerOffset:]
        
        # Find peak index
        peak_ind = np.argmax(filtered_radian_data)
        trough_ind = np.argmin(filtered_radian_data[peak_ind:]) + peak_ind
        rise_start_index = 0
        rise_end_index = len(filtered_radian_data)

        # Go back 100 ms from peak
        search_start = peak_ind - int(experiment.readRate / 10)



        # Search for rise
        for x in np.arange(search_start, peak_ind, 1):
            if filtered_radian_data[x] > noise_floor:
                rise_start_index = x
                break
        
        # rise endpeak_decel
        for x in np.arange(peak_ind, peak_ind + int(experiment.readRate / 10), 1):
            if filtered_radian_data[x] < noise_floor:
                rise_end_index = x
                break


        # enable this to see the filtered data overlayed on raw data
        #self.axes.plot(filtered_data, label="filtered")
        
        #print("noise_floor:{} start:{} end:{}".format(noise_floor, rise_start_index, rise_end_index))
        
        # calculate derivative, derivative for acceleration
        acceleration = np.gradient(filtered_radian_data, 1. / experiment.readRate)
        # peak acceleration
        peak_accel_ind = np.argmax(acceleration[rise_start_index:peak_ind]) + rise_start_index
        peak_accel = acceleration[peak_accel_ind]

        # peak deceleration
        peak_decel_ind = np.argmin(acceleration[peak_ind:trough_ind]) + peak_ind
        peak_decel = acceleration[peak_decel_ind]

        delta_t = (rise_end_index - rise_start_index) / (experiment.readRate / 1000)  # (self.experiment.readRate/1000) is samples/millesecond
        time_to_peak = (peak_ind - rise_start_index) / (experiment.readRate / 1000)
        decel_time = (rise_end_index - peak_ind) / (experiment.readRate / 1000)
        peak_vel = filtered_radian_data[peak_ind]

        ###
        # calculate angle of excursion(degrees) and linear distance travelled by swing arm
        # from rise_start to rise_end segment
        ###
        ##        excursionAngleRadians = np.sum(filtered_data[rise_start_index:rise_end_index]) * (1.0/self.experiment.readRate)
        excursionAngleRadians = np.sum(filtered_radian_data[rise_start_index:rise_end_index]) * (1.0 / experiment.readRate)
        excursionAngleDegrees = np.degrees(excursionAngleRadians)
        # excursionTravelDistanceCm = self.experiment.swingArmToPivotCm * 2 * np.sin(excursionAngleRadians/2)
        excursionTravelDistanceInches = experiment.swingArmToPivotInches * 2 * np.sin(excursionAngleRadians / 2)

        ###
        # final resting positionenvironment
        # from first non-zero velocity to final final non-zero velocity
        ###
        # calculate linear excursion distance
        #restAngleRadians = np.sum(filtered_data[rise_start_index:rest_index]) * (1.0 / experiment.readRate)
        #restPistonDistanceInches = experiment.swingArmToPivotInches * 2 * np.sin(restAngleRadians / 2)

        ###
        # find point of impact
        ###
        #impact_ind = np.argmax(np.diff(filtered_data[peak_ind:rise_end_index]) < -0.3) + peak_ind

        ###
        # distance from fire to impact
        ###
        #impactAngleRadians = np.sum(filtered_data[rise_start_index:impact_ind]) * (1.0 / experiment.readRate)
        #impactTravelDistanceInches = experiment.swingArmToPivotInches * 2 * np.sin(impactAngleRadians / 2)
        
        summary_data = {'peak_ind': peak_ind, 'trough_ind': trough_ind, 'peak_accel_ind': peak_accel_ind,
                'peak_decel_ind': peak_decel_ind, 'delta_t': delta_t,
                'time_to_peak': time_to_peak, 'decel_time': decel_time, 'peak_vel': peak_vel, 'peak_accel': peak_accel,
                'peak_decel': peak_decel,
                'rise_start_index': rise_start_index, 'rise_end_index': rise_end_index,
                'excursionAngleDegrees': excursionAngleDegrees,
                'maxPistonTravelDistanceInches': excursionTravelDistanceInches }

        return summary_data



    def draw_summary_box(self, experiment, summary):
        """
        Add a text box to plot
        """
        summary_txt = "Set: {} $psi$  Load: {} $psi$\nPeak velocity: {:.4f} $rad/s$\nTime to peak: {:.2f} $ms$\nDeceleration Time: {:.2f} $ms$\nDelta t: {:.2f} $ms$\nMax acceleration: {:.2f} $rad/s^2$\nMax deceleration: {:.2f} $rad/s^2$\nPP Excursion Angle: {:.2f}$^\circ$\nPP Excursion: {:.2f} $inches$".format(
            experiment.PsiSet, experiment.PsiLoad, summary['peak_vel'], summary['time_to_peak'],
            summary['decel_time'], summary['delta_t'], summary['peak_accel'], summary['peak_decel'],
            summary['excursionAngleDegrees'], summary['maxPistonTravelDistanceInches'])

        if self.summaryBox is not None:
            self.summaryBox.remove()
            self.summaryBox = None

        from matplotlib.offsetbox import AnchoredText
        anchored_text = AnchoredText(summary_txt, loc=1, prop=dict(family='sans-serif', size=10, linespacing=1.5))
        anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        anchored_text.patch.set_facecolor('beige')
        anchored_text.patch.set_alpha(0.95)
        self.summaryBox = self.axes.add_artist(anchored_text)
