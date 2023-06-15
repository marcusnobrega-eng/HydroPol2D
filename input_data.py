###################################################################
#                                                                 #
#                        Produced by:                             #
#                 Marcus Nobrega Junior gomes                     #
#              (marcusnobrega.engcivil@gmail.com)                 #
#                               &                                 #
#                 Luis Miguel Castillo Rápalo                     #
#                   (luis.castillo@unah.hn)                       #
#                    September 2021                               #
#                                                                 #
#       Last update : 7 July, 2021                                #
###################################################################

import configparser
import numpy as np
import pandas as pd


class GeneralData:
    def __init__(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)
        self.resolution = float(config['GENERAL']['resolution'])
        self.slope_min = float(config['GENERAL']['slope_min'])
        self.alfa = float(config['GENERAL']['alfa'])
        self.depth_wse = float(config['GENERAL']['depth_wse'])
        self.ADD = float(config['GENERAL']['ADD'])
        self.slope_alfa = float(config['GENERAL']['slope_alfa'])
        self.alfa_min = float(config['GENERAL']['alfa_min'])
        self.v_thresold = float(config['GENERAL']['v_threshold'])
        
        self.routing_time = float(config['ROUTING']['routing_time'])
        self.time_step_model = float(config['ROUTING']['time_step_model'])/60 #it is already in minutes????????w
        self.min_time_step = float(config['ROUTING']['min_time_step'])
        self.max_time_step = float(config['ROUTING']['max_time_step'])
        self.time_step_increments = float(config['ROUTING']['time_step_increments'])
        self.time_step_change = float(config['ROUTING']['time_step_change'])
        self.time_step_matrices = float(config['ROUTING']['time_step_matrices'])

        self.factor_cells = float(config['TOLERANCE']['factor_cells'])
        self.flow_tolerance = float(config['TOLERANCE']['flow_tolerance'])
        self.depth_tolerance = float(config['TOLERANCE']['depth_tolerance'])
        
        self.outlet_type = float(config['OUTLET']['outlet_type'])
        self.slope_outlet = float(config['OUTLET']['slope_outlet'])
        self.x_outlet_begin = float(config['OUTLET']['x_outlet_begin'])
        self.x_outlet_end = float(config['OUTLET']['x_outlet_end'])
        self.y_outlet_begin = float(config['OUTLET']['y_outlet_begin'])
        self.y_outlet_end = float(config['OUTLET']['y_outlet_end'])
        
        self.record_time_maps = float(config['RECORD']['record_time_maps'])
        self.record_time_hydrographs = float(config['RECORD']['record_time_hydrographs'])

        self.x_min = float(config['BOUNDARIES']['x_min'])
        self.x_max = float(config['BOUNDARIES']['x_max'])
        self.y_min = float(config['BOUNDARIES']['y_min'])
        self.y_max = float(config['BOUNDARIES']['y_max'])
        self.xllupcorner = float(config['BOUNDARIES']['xllupcorner'])
        self.yllupcorner = float(config['BOUNDARIES']['yllupcorner'])
        
        self.flag_rainfall = float(config['FLAGS']['flag_rainfall'])
        self.flag_inflow = float(config['FLAGS']['flag_inflow'])
        self.flag_abstraction = float(config['FLAGS']['flag_abstraction'])
        self.flagwse = float(config['FLAGS']['flagwse'])
        self.flag_waterbalance = float(config['FLAGS']['flag_waterbalance'])
        self.flag_waterquality = float(config['FLAGS']['flag_waterquality'])
        self.flag_timestep = float(config['FLAGS']['flag_timestep'])
        self.flag_warmup = float(config['FLAGS']['flag_warmup'])
        self.flag_initial_buildup = float(config['FLAGS']['flag_initial_buildup'])
        self.flag_wq_model = float(config['FLAGS']['flag_wq_model'])
        self.flag_infiltration = float(config['FLAGS']['flag_infiltration'])
        self.flag_critical = float(config['FLAGS']['flag_critical'])
        self.pol_min = float(config['FLAGS']['pol_min'])
        self.flag_floodmaps = float(config['FLAGS']['flag_floodmaps'])
        self.flag_riskmaps = float(config['FLAGS']['flag_riskmaps'])
        self.flag_gifs = float(config['FLAGS']['flag_gifs'])
        self.flag_txt = float(config['FLAGS']['flag_txt'])
        self.flag_vel = float(config['FLAGS']['flag_vel'])
        self.flag_risk_model = float(config['FLAGS']['flag_risk_model'])

        self.g = float(config['HUMAN_RISK']['g'])
        self.f = float(config['HUMAN_RISK']['f'])
        self.Fc = float(config['HUMAN_RISK']['Fc'])
        self.Cc = float(config['HUMAN_RISK']['Cc'])
        self.m_c_m = float(config['HUMAN_RISK']['m_c_m'])
        self.y_c_m = float(config['HUMAN_RISK']['y_c_m'])
        self.D_c_m = float(config['HUMAN_RISK']['Di_c_m'])
        self.d_c_m = float(config['HUMAN_RISK']['d_c_m'])
        self.m_t_m = float(config['HUMAN_RISK']['m_t_m'])
        self.y_t_m = float(config['HUMAN_RISK']['y_t_m'])
        self.D_t_m = float(config['HUMAN_RISK']['Di_t_m'])
        self.d_t_m = float(config['HUMAN_RISK']['d_t_m'])
        self.m_a_m = float(config['HUMAN_RISK']['m_a_m'])
        self.y_a_m = float(config['HUMAN_RISK']['y_a_m'])
        self.D_a_m = float(config['HUMAN_RISK']['Di_a_m'])
        self.d_a_m = float(config['HUMAN_RISK']['d_a_m'])
        self.m_o_m = float(config['HUMAN_RISK']['m_o_m'])
        self.y_o_m = float(config['HUMAN_RISK']['y_o_m'])
        self.D_o_m = float(config['HUMAN_RISK']['Di_o_m'])
        self.d_o_m = float(config['HUMAN_RISK']['d_o_m'])
        self.m_c_f = float(config['HUMAN_RISK']['m_c_f'])
        self.y_c_f = float(config['HUMAN_RISK']['y_c_f'])
        self.D_c_f = float(config['HUMAN_RISK']['Di_c_f'])
        self.d_c_f = float(config['HUMAN_RISK']['d_c_f'])
        self.m_t_f = float(config['HUMAN_RISK']['m_t_f'])
        self.y_t_f = float(config['HUMAN_RISK']['y_t_f'])
        self.D_t_f = float(config['HUMAN_RISK']['Di_t_f'])
        self.d_t_f = float(config['HUMAN_RISK']['d_t_f'])
        self.m_a_f = float(config['HUMAN_RISK']['m_a_f'])
        self.y_a_f = float(config['HUMAN_RISK']['y_a_f'])
        self.D_a_f = float(config['HUMAN_RISK']['Di_a_f'])
        self.d_a_f = float(config['HUMAN_RISK']['d_a_f'])
        self.m_o_f = float(config['HUMAN_RISK']['m_o_f'])
        self.y_o_f = float(config['HUMAN_RISK']['y_o_f'])
        self.D_o_f = float(config['HUMAN_RISK']['Di_o_f'])
        self.d_o_f = float(config['HUMAN_RISK']['d_o_f'])

class LULCData: #this should be simplified to n number of lulc classes, an iterative procedure!!** here!!**
    def __init__(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)
        self.imp1_index = float(config['Imp1']['index'])
        self.imp1_n = float(config['Imp1']['n'])
        self.imp1_ksat = float(config['Imp1']['ksat'])
        self.imp1_d0 = float(config['Imp1']['d0'])
        self.imp1_h0 = float(config['Imp1']['h0'])
        self.imp1_psi = float(config['Imp1']['psi'])
        self.imp1_I0 = float(config['Imp1']['I0'])
        self.imp1_tetasat = float(config['Imp1']['tetasat'])
        self.imp1_tetai = float(config['Imp1']['tetai'])
        self.imp1_C1 = float(config['Imp1']['C1'])
        self.imp1_C2 = float(config['Imp1']['C2'])
        self.imp1_C3 = float(config['Imp1']['C3'])
        self.imp1_C4 = float(config['Imp1']['C4'])
        
        self.per2_index = float(config['Per2']['index'])
        self.per2_n = float(config['Per2']['n'])
        self.per2_ksat = float(config['Per2']['ksat'])
        self.per2_d0 = float(config['Per2']['d0'])
        self.per2_h0 = float(config['Per2']['h0'])
        self.per2_psi = float(config['Per2']['psi'])
        self.per2_I0 = float(config['Per2']['I0'])
        self.per2_tetasat = float(config['Per2']['tetasat'])
        self.per2_tetai = float(config['Per2']['tetai'])
        self.per2_C1 = float(config['Per2']['C1'])
        self.per2_C2 = float(config['Per2']['C2'])
        self.per2_C3 = float(config['Per2']['C3'])
        self.per2_C4 = float(config['Per2']['C4'])
        
        self.per3_index = float(config['Per3']['index'])
        self.per3_n = float(config['Per3']['n'])
        self.per3_ksat = float(config['Per3']['ksat'])
        self.per3_d0 = float(config['Per3']['d0'])
        self.per3_h0 = float(config['Per3']['h0'])
        self.per3_psi = float(config['Per3']['psi'])
        self.per3_I0 = float(config['Per3']['I0'])
        self.per3_tetasat = float(config['Per3']['tetasat'])
        self.per3_tetai = float(config['Per3']['tetai'])
        self.per3_C1 = float(config['Per3']['C1'])
        self.per3_C2 = float(config['Per3']['C2'])
        self.per3_C3 = float(config['Per3']['C3'])
        self.per3_C4 = float(config['Per3']['C4'])
        
        self.per4_index = float(config['Per4']['index'])
        self.per4_n = float(config['Per4']['n'])
        self.per4_ksat = float(config['Per4']['ksat'])
        self.per4_d0 = float(config['Per4']['d0'])
        self.per4_h0 = float(config['Per4']['h0'])
        self.per4_psi = float(config['Per4']['psi'])
        self.per4_I0 = float(config['Per4']['I0'])
        self.per4_tetasat = float(config['Per4']['tetasat'])
        self.per4_tetai = float(config['Per4']['tetai'])
        self.per4_C1 = float(config['Per4']['C1'])
        self.per4_C2 = float(config['Per4']['C2'])
        self.per4_C3 = float(config['Per4']['C3'])
        self.per4_C4 = float(config['Per4']['C4'])
        
        self.per5_index = float(config['Per5']['index'])
        self.per5_n = float(config['Per5']['n'])
        self.per5_ksat = float(config['Per5']['ksat'])
        self.per5_d0 = float(config['Per5']['d0'])
        self.per5_h0 = float(config['Per5']['h0'])
        self.per5_psi = float(config['Per5']['psi'])
        self.per5_I0 = float(config['Per5']['I0'])
        self.per5_tetasat = float(config['Per5']['tetasat'])
        self.per5_tetai = float(config['Per5']['tetai'])
        self.per5_C1 = float(config['Per5']['C1'])
        self.per5_C2 = float(config['Per5']['C2'])
        self.per5_C3 = float(config['Per5']['C3'])
        self.per5_C4 = float(config['Per5']['C4'])
        
        self.per6_index = float(config['Per6']['index'])
        self.per6_n = float(config['Per6']['n'])
        self.per6_ksat = float(config['Per6']['ksat'])
        self.per6_d0 = float(config['Per6']['d0'])
        self.per6_h0 = float(config['Per6']['h0'])
        self.per6_psi = float(config['Per6']['psi'])
        self.per6_I0 = float(config['Per6']['I0'])
        self.per6_tetasat = float(config['Per6']['tetasat'])
        self.per6_tetai = float(config['Per6']['tetai'])
        self.per6_C1 = float(config['Per6']['C1'])
        self.per6_C2 = float(config['Per6']['C2'])
        self.per6_C3 = float(config['Per6']['C3'])
        self.per6_C4 = float(config['Per6']['C4'])

    def variable_to_matrix(self): ###same here, iterative procedure depending the number of LULC classes **!!here!!**
        data = np.array([[self.imp1_n, self.imp1_ksat, self.imp1_d0, self.imp1_h0, self.imp1_psi,
                          self.imp1_I0, self.imp1_tetasat, self.imp1_tetai, self.imp1_C1, self.imp1_C2, 
                          self.imp1_C3, self.imp1_C4], 
                         [self.per2_n, self.per2_ksat, self.per2_d0, self.per2_h0, self.per2_psi,
                          self.per2_I0, self.per2_tetasat, self.per2_tetai, self.per2_C1, self.per2_C2, 
                          self.per2_C3, self.per2_C4],
                         [self.per3_n, self.per3_ksat, self.per3_d0, self.per3_h0, self.per3_psi,
                          self.per3_I0, self.per3_tetasat, self.per3_tetai, self.per3_C1, self.per3_C2, 
                          self.per3_C3, self.per3_C4],
                         [self.per4_n, self.per4_ksat, self.per4_d0, self.per4_h0, self.per4_psi,
                          self.per4_I0, self.per4_tetasat, self.per4_tetai, self.per4_C1, self.per4_C2, 
                          self.per4_C3, self.per4_C4], 
                         [self.per5_n, self.per5_ksat, self.per5_d0, self.per5_h0, self.per5_psi,
                          self.per5_I0, self.per5_tetasat, self.per5_tetai, self.per5_C1, self.per5_C2, 
                          self.per5_C3, self.per5_C4], 
                         [self.per6_n, self.per6_ksat, self.per6_d0, self.per6_h0, self.per6_psi,
                          self.per6_I0, self.per6_tetasat, self.per6_tetai, self.per6_C1, self.per6_C2, 
                          self.per6_C3, self.per6_C4]])
        return data

    def lulc_index_array(self):
        lulc_index = np.array([self.imp1_index, self.per2_index, self.per3_index, self.per4_index,
                               self.per5_index, self.per6_index])
        return lulc_index

class rainfalls_number:
    def __init__(self, precipitation_file):
        df = pd.read_csv(precipitation_file, sep=',')
        self.rains = int(len(df.columns)/2)

class PrecipitationData:  #here to version no2, as matrix for distribuited rainfall, loop for iterative list of rasters
    def __init__(self, precipitation_file, rain):
        df = pd.read_csv(precipitation_file, sep=',')
        self.intensity_rainfall = df.iloc[:, (rain*2)+1]
        self.time_rainfall = df.iloc[:, (rain*2)+0]
        self.time_step_rainfall = df.iloc[2, (rain*2)+0] - df.iloc[1, (rain*2)+0]
        self.rainfall_duration = df.iloc[2, (rain*2)+0] - df.iloc[1, (rain*2)+0]

class InflowData:
    def __init__(self, inflow_file, resolution):
        df = pd.read_csv(inflow_file, sep=',') #check the separator of 176
        self.inflow_discharge = df.iloc[:, 1].dropna()
        self.time_inflow = df.iloc[:, 0]. dropna()
        self.time_step_inflow = df.iloc[1, 0] - df.iloc[0, 0]
        self.inlet_cells = df.iloc[:, 2:4].dropna()
        self.n_inlets = len(self.inlet_cells)
        if self.n_inlets == 0:
            raise Exception('Please, insert the inlet point(s) in the Inflow_Hydrograph.csv file')

        self.area = self.area_calculation(resolution)
        self.inflow_hydrograph_rate = self.inflow_rate()

    def area_calculation(self, resolution):
        area = resolution ** 2 * self.n_inlets
        return area

    def inflow_rate(self):
        inflow_hydrograph_rate = self.inflow_discharge / self.area * 1000 * 3600 # to (mm/h)
        return inflow_hydrograph_rate



