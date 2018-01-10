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

import os
import datetime
from enum import Enum
import numpy as np
from pandas import read_csv
import scipy.signal

class SensorRotation(Enum):
    """
    An enumeration of possible sensor mount rotations
    """
    COUNTERCLOCKWISE = 1
    CLOCKWISE = 2


class Experiment:
    """
    An instance of data collection either from the sensor or a data file.
    Settings specific to that instance and its data.
    """

    """ default settings """
    readRate                = 25000.0 # The number of samples per channel per second
    readNumSamples          = 50000   # total number of samples. self.readRate * number of channels
    triggerOffset           = 25000   # sample space from trigger to fire
    swingArmToPivotInches   = 2.3622  # distance from linkage bolt to the swing arm bolt in inches. Used to calculate linear travel.
    voltageScalingFactor    = 1.5     # to scale voltage divider back up

    channel1_radianConvFactor        = 1.0 # sensor scaling factor (this is for one of the sensors and each is different)
    channel2_radianConvFactor        = 1.0
    channel1_sensorRotation          = SensorRotation.COUNTERCLOCKWISE
    channel2_sensorRotation          = SensorRotation.CLOCKWISE
    
    voltage_data            = np.array([], dtype=float)

    


    def __init__(self):
        self._date = datetime.datetime.today().strftime("%Y%m%d")
        self._time = datetime.datetime.today().strftime("%H%M%S")
        self._subject_id = ""  # default subject id
        self.PsiLoad = ""  # default Load psi setting
        self.PsiSet = ""  # default Set psi setting
        self.lastDataPath = ''
        self.file_name = self.subjectId + "_" + self.date + "_" + self.time + ".hyg"

    @property
    def subjectId(self):
        return self._subject_id

    @subjectId.setter
    def subjectId(self, value):
        self._subject_id = value
        self.file_name = self._subject_id + "_" + self._date + "_" + self._time + ".hyg"
        
    @property
    def date(self):
        return self._date

    @date.setter
    def date(self, value):
        self._date = value
        self.file_name = self.subjectId + "_" + self._date + "_" + self._time + ".hyg"

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        self._time = value
        self.file_name = self.subjectId + "_" + self._date + "_" + self._time + ".hyg"  

    def set_timestamp(self):
        import datetime
        self._date = datetime.datetime.today().strftime("%Y%m%d")
        self._time = datetime.datetime.today().strftime("%H%M%S")
        self.file_name = self.subjectId + "_" + self._date + "_" + self._time + ".hyg"

    def get_header(self):
        """
        Provide consistent header for experiments
        :return: 
        """
        return "Date:{}\nTime:{}\nSubjectID:{}\nPSI Load:{}\nPSI Set:{}\nSampling Rate(Hz):{}\nChannel1 SensorScalar:{}\nChannel2 SensorScalar:{}\nSamples Per Channel:{}\nY_Unit_Label:{}\nX_Dimension:{}\nChannel Order:AI0,AI1".format(
           self._date, self.time, self.subjectId, self.PsiLoad,
           self.PsiSet, self.readRate, self.channel1_radianConvFactor, self.channel2_radianConvFactor, self.readNumSamples, "Volts", "Time")

    def get_label(self):
        """
        name for experiment
        """
        # return self.subjectId + "_" + self.date + "_" + self.time
        # didn't do that incase someone change file name
        return self.file_name
    
    def get_radian_data(self):
        """
        """
        return (self.voltage_data / self.channel1_radianConvFactor)
    
        
    def filter_radian_data(self, radian_data, filter_type='binomial'):
        """
        """
        if filter_type == 'binomial':
        
            # binomial filter used in UPENN code
            #
            # create Pascal coefficients
            a = b = [0.5, 0.5]
            # create a length 19 binamial filter
            for x in range(17):
                b = np.convolve(a, b, mode='full')
            # filter data
            return np.convolve(radian_data, b, mode='same')
        
        elif filter_type == 'Savitsky-Golay':
            #
            # Savitsky-Golay filter. Perhaps too smooth and not what UPENN did
            #
            # filtered_data = scipy.signal.savgol_filter(data, window_length=19, polyorder=3)
            #            
            return scipy.signal.savgol_filter(radian_data, window_length=19, polyorder=3)

        else:
            return radian_data
        

    def resampled_data(self, data, new_rate):
        
        # first downsample the trigger offset setting to compensate for data downsampling to match underlay
        #self.triggerOffset = int(self.triggerOffset/(int(self.readRate)/new_rate))
        # update read rate of experment to match underlay
        #self.readRate = float(new_rate)
        # downsample the data
        return scipy.signal.resample(data, new_rate)
        
        
        
    @classmethod
    def load(cls, data_file_path):

        experiment = Experiment()

        # parse header
        with open(data_file_path, 'r') as hdr:
            for line in hdr:
                if line.startswith("#"):
                    # use some fields from header
                    fld_name, fld_value = line[2:].strip().split(":")

                    if fld_name == "Date":
                        experiment.date = fld_value

                    if fld_name == "Time":
                        experiment.time = fld_value

                    if fld_name == "SubjectID":
                        experiment.subjectId = fld_value

                    if fld_name == "PSI Load":
                        experiment.PsiLoad = fld_value

                    if fld_name == "PSI Set":
                        experiment.PsiSet = fld_value

                    if fld_name == "Sampling Rate(Hz)":
                        experiment.readRate = float(fld_value)
                        experiment.triggerOffset = int(float(experiment.readRate))
                        
                    if fld_name == "Channel1 SensorScalar":
                        experiment.channel1_radianConvFactor = float(fld_value)
                    # backwards compatibility channel1
                    if fld_name == "SensorScalar":
                        experiment.channel1_radianConvFactor = float(fld_value)  
                        
                    if fld_name == "Channel2 SensorScalar":
                        experiment.channel2_radianConvFactor = float(fld_value)
                        
                else:
                    break

        experiment.lastDataPath = os.path.sep.join(str(data_file_path).split('/')[0:-1])
        experiment.file_name = str(data_file_path).split('/')[-1]

        # load data but skip the time we waited for the fire
        # load only channel one ", usecols=[0]"
        # do NOT use skiplines in loadtxt() to truncate
        # switched from np.loadtxt() to pandas.read_csv() for performance
        experiment.voltage_data = np.squeeze(read_csv(str(data_file_path),
                                                      delimiter=",", header=None,
                                                      usecols=[0], comment='#',
                                                      dtype='float').values)

        return experiment


    def save(self, path):
        """ 
        write out trace data in voltage 
        """

        if len(self.voltage_data.shape) > 1:
            # two channel interleaved
            data_unleaved = np.array([self.voltage_data[0::2], self.voltage_data[1::2]]).transpose()
            # datetime stamp experiment here
            self.set_timestamp()
            np.savetxt(path,
                       data_unleaved, fmt='%.11f', delimiter=',',
                       header=self.get_header())            
        else:
            # datetime stamp experiment here
            self.set_timestamp()
            np.savetxt(path,
                       self.voltage_data, fmt='%.11f', delimiter=',',
                       header=self.get_header())  
            

        
