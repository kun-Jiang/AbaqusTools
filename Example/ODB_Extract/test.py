# -*- coding: utf-8 -*-
from odbAccess import *
from textRepr import prettyPrint as prettyprint
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import threading
from multiprocessing import Pool, pool
from joblib import Parallel, delayed

def logwrite(message):
    with open('log.txt','a') as log_file:
        log_file.write(message+'\n')
        print(message)
        log_file.close()

class Get_field_output():
#class Get_field_output(np.ndarray):    
    # def __new__(cls, current_frame, output_option):
    #     output = current_frame.fieldOutputs[output_option]
    #     field_output_values = output.values
    #     obj = np.asarray(field_output_values).view(cls)
    #     obj.current_frame = current_frame
    #     obj.output_option = output_option
    #     return obj
    
    def __init__(self,current_frame, output_option, field_vector_length):
        self.current_frame = current_frame
        self.output_option = output_option
        self.field_vector_length = field_vector_length
        # Identify the length of the field vector, that is, the number of components

    def get_values_array(self):
        self.field_output = self.current_frame.fieldOutputs[self.output_option]
        # prettyprint(self.field_output)
        self.field_output_values = self.field_output.values
        if 'UVARM' in self.output_option:
            # Values of UVARM and SDV are stored in the intergration points, so the length of the array is 4 times of the number of elements
            # The first column is the element label
            # The second column is the value in the first intergration point
            # The third column is the value in the second intergration point
            # The fourth column is the value in the third intergration point
            # The fifth column is the value in the fourth intergration point
            field = np.zeros((len(self.field_output_values)/4, self.field_vector_length))
        else:
            # Values of EVOL are stored in the element, so the length of the array is the number of elements
            # The first column is the element label
            # The second column is the value of scale field
            field = np.zeros((len(self.field_output_values),self.field_vector_length))
        count = 0
        for scale_field in self.field_output_values:
            if count == 0:
                # As usual, the number of the first element is 1, 
                # if the number of the first element is not 1, there is a offset
                element_offset = scale_field.elementLabel -1
                count += 1
            if self.output_option == 'EVOL':
                # The EVOL variable is defined in the dummy element, so there is a offset in the element number
                element_index = scale_field.elementLabel-element_offset-1
                field[element_index,0] = scale_field.elementLabel
                field[element_index,1] = scale_field.data
            elif 'UVARM' in self.output_option:
                # The UVARM variable is also defined in the dummy element, so there is a offset in the element number
                # In the case of ith element and jth intergration point
                # The index of the element in the array is 
                #               elementLable - element_offset-1 (counting from 0 in python)
                element_index = scale_field.elementLabel-element_offset-1
                intergration_point_sequence = scale_field.integrationPoint
                # print(intergration_point_sequence)
                # print(np.shape(field))
                field[element_index,0] = scale_field.elementLabel
                field[element_index,intergration_point_sequence] = scale_field.data
        # print(scale_field.elementLabel-1)
        return field
    
class plot_x_y_curve():
    def __init__(self,data_vector,x_label,y_label,title,file_name):
        self.data_x = [vector[0] for vector in data_vector]
        self.data_y = [vector[1] for vector in data_vector]
        self.label_x = x_label
        self.label_y = y_label
        self.title_name = title
        self.export_file_name = file_name
        self.plot()
    def plot(self):
        # Plot the curve
        logwrite('Plot the curve')
        fig=plt.figure(num=1,figsize=(6,6))
        plt.plot(self.data_x,self.data_y,'r-',linewidth=1)
        plt.title(self.title_name,fontsize=16, weight='bold', fontname='sans-serif')
        plt.xlabel(self.label_x,fontsize=16, weight='bold', fontname='sans-serif')
        plt.ylabel(self.label_y,fontsize=16, weight='bold', fontname='sans-serif')
        # plt.grid(True)
        plt.show()

def process_frame(frame):
    current_time = frame.frameValue
    # ********************************************************************************************
    # Calculate the weight by ratio of element volume to total volume
    # ********************************************************************************************
    # Get the values of element volume for every element
    # output_option = 'EVOL'
    field_vector_length = 2
    field_array = Get_field_output(frame, 'EVOL', field_vector_length).get_values_array()
    # Get the column that contains element volume and calculate the total volume
    EVOL_Total = np.sum(field_array[:,1])
    # Get the weight of every element
    Element_Weight = field_array[:,1]/EVOL_Total
    # Write the weight of every element into a txt file
    with open('EVOL.txt','w') as file:
        for row in field_array:
            file.write(str(row)+'\n')

    # ********************************************************************************************
    # Calculate the average value of user-defined field output
    # ********************************************************************************************
    # Get the values of user-defined field output for every element
    output_option = 'UVARM2'
    field_vector_length = 5
    field_array = Get_field_output(frame, output_option, field_vector_length).get_values_array()
    Concentration = field_array[:,1:]
    Concentration_Element_Weighted = np.sum(Concentration[:,0] + 
                                            Concentration[:,1] + 
                                            Concentration[:,2] + 
                                            Concentration[:,3])/4
    # Write the Concentration data into a txt file
    with open('Concentration.txt','w') as file:
        for row in field_array:
            file.write(str(row)+'\n')
    
    # ********************************************************************************************
    # Calculate the average value of user-defined field output
    # ********************************************************************************************    
    Concentration_Total_Weighted = np.sum(Concentration_Element_Weighted*Element_Weight)
    
    Data_plot[frame.incrementNumber,0] = current_time
    Data_plot[frame.incrementNumber,1] = Concentration_Total_Weighted
    frame_count += 1
    # finish_percent = float(frame_count/frame_number*100)
    print("Progress: Frame:%d" % frame_count)

    return Data_plot

if __name__ == '__main__':
    # Check if the log file exists
    if os.path.exists('log.txt'):
        os.remove('log.txt')
    logwrite('*'*70 + '\n' +
             '*' + '{0:^68}'.format('Extracting Data from ODB') + '*' + '\n' +
             '*'*70)
    # *********************************************************************************************************************
    # Set the parameters
    # *********************************************************************************************************************
    odb_file_name = 'Model_Chem.odb'
    output_option = 'UVARM2'
    logwrite('{0:<25s}'.format('ODB file name') + ':' + odb_file_name)
    logwrite('{0:<25s}'.format('Output field parameter') + ':' + output_option + '\n' +
             '-'*60)
    try:
        odb = openOdb(odb_file_name,readOnly=True)
        logwrite("Open the ODB file successfully")
    except:
        logwrite("There is no ODB file named " + odb_file_name)
    # Extract the data in the specified step
    step = odb.steps['DIffusion']
    if output_option == 'EVOL':
        field_vector_length = 2
    if 'SDV' in output_option:
        field_vector_length = 2
    if 'UVARM' in output_option:
        field_vector_length = 5
    if 'U' in output_option:
        # 识别别是二维单元还是三维单元，二维单元只有两个分量，三维单元有三个分量，将这些结果都输出一下
        field_vector_length = 2
    if 'S' in output_option:
        field_vector_length = 3
        
    # Extract the data in the every frame
    Data_plot = np.zeros([len(step.frames),field_vector_length])
    # print(len(step.frames))
    # 把关键字输出一下
    # print(step.frames[0].fieldOutputs.keys())
    frame_count = 0
    frame_number = len(step.frames)
    logwrite('Extract the data from every frame...')
    # 创建一个Pool对象，使用所有可用的CPU核心
    pool = Pool()

    # 对于每个时间步，调用process_frame函数处理数据
    Data_plot = Parallel(n_jobs=-1)(delayed(process_frame)(frame) for frame in step.frames)
    try:
        plot_x_y_curve(data_vector=Data_plot, 
                    x_label='Time', y_label='Displacement', title='Displacement-Time Curve', 
                    file_name='Displacement-Time Curve')
    except:
        logwrite('Plot the curve failed')
    


        # displacement_last=current_frame.fieldOutputs['U']
        # displacementValues_last=displacement_last.values
        
        # state_dependent_variable=current_frame.fieldOutputs['SDV2']
        # state_dependent_variable_values=state_dependent_variable.values
        
        # SDV1 = []
        # for scale in state_dependent_variable_values:
        #     # print(scale)
        #     # os._exit(0)
        #     SDV1.append([scale.elementLabel,scale.data])
        
        # DISP = []
        # for v in displacementValues_last:
        #     DISP.append([v.nodeLabel,v.data[0],v.data[1]])
        #     # DISP.append([v.data[0],v.data[1]])
        # try:
        #     plot_x_y_curve(data_vector=DISP, 
        #                 x_label='Time', y_label='Displacement', title='Displacement-Time Curve', 
        #                 file_name='Displacement-Time Curve')
        # except:
        #     logwrite('Plot the curve failed')


        # with open('field.txt', 'w') as f:
        #     for i in DISP:
        #         i=(','.join(str(j) for j in i))
        #         f.write(str(i) + '\n')
        # odb.close
        # logwrite("Write the txt file successfully")    