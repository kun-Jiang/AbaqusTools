# -*- coding: utf-8 -*-
from odbAccess import *
from textRepr import prettyPrint as prettyprint
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import threading

def logwrite(message):
    with open('log.txt','a') as log_file:
        log_file.write(message+'\n')
        print(message)
        log_file.close()

def save_as_file(data_vector, file_name):
    file_path = os.path.join(extracted_data_folder, '%s.csv'% file_name)
    [rows_num, columns_num] = np.shape(data_vector)
    with open(file_path,'w') as data_file:
        for row in range(rows_num):
            row_data = []
            for column in range(columns_num):
                row_data.append(float(data_vector[row,column]))
            row_line = ','.join(str(x) for x in row_data) + '\n'
            data_file.write(row_line)

class Run_in_new_thread():
    def __init__(self, func):
        self.func = func
        # Create new thread to run the function
        threading.Thread(target=func).start()

class Get_field_output():
    
    def __init__(self, current_frame, output_option, output_suboption=None, field_vector_length=None):
        self.current_frame = current_frame
        self.output_option = output_option
        self.output_suboption = output_suboption
        self.field_vector_length = field_vector_length
        # Identify the length of the field vector, that is, the number of components

    def get_values_array(self):
        self.field_output = self.current_frame.fieldOutputs[self.output_option]
        # Get the type of the field output
        field_type = self.field_output.type
        self.field_output_values = self.field_output.values
        # **********************************************************************
        # Prepare the array to store the field values, set its size in advance
        # **********************************************************************
        if self.output_option == 'EVOL':
            # Values of EVOL are stored in the element, so the length of the array is the number of elements
            # The first column is the element label
            # The second column is the value of scale field
            field = np.zeros((len(self.field_output_values),self.field_vector_length))
        if 'UVARM' in self.output_option or 'SDV' in self.output_option:
            # Values of UVARM and SDV are stored in the intergration points, so the length of the array is 4 times of the number of elements
            # The first column is the element label
            # The second column is the value in the first intergration point
            # The third column is the value in the second intergration point
            # The fourth column is the value in the third intergration point
            # The fifth column is the value in the fourth intergration point
            # ......
            field = np.zeros((len(self.field_output_values)/8, self.field_vector_length))
        if self.output_option == 'S':
            field = np.zeros((len(self.field_output_values),self.field_vector_length))
        if self.output_option == 'U':
            field = np.zeros((len(self.field_output_values),self.field_vector_length))
            pass

        # **********************************************************************
        # Obtain the output field values and store them in the array
        # There are kinds of output field, such as UVARM, EVOL, SDV, etc. 
        # Their data structure are different which need to be treated differently
        # **********************************************************************
        count = 0
        for field_value in self.field_output_values:
            # prettyprint(field_value)
            if count == 0:
                # As usual, the number of the first element is 1, 
                # if the number of the first element is not 1, there is maybe a offset
                if field_value.elementLabel != None:
                    # Some field outputs are defined on the node, so there isn't element label
                    element_offset = field_value.elementLabel -1
                count += 1
            # **************************************************************************************
            # Output the volume of the element
            # **************************************************************************************
            if self.output_option == 'EVOL':
                # The EVOL variable is defined in the dummy element, so there is a offset in the element number
                element_index = field_value.elementLabel-element_offset-1
                field[element_index,0] = field_value.elementLabel
                field[element_index,1] = field_value.data
            # **************************************************************************************
            # Output the stress of the element
            # **************************************************************************************
            elif self.output_option == 'S':
                # prettyprint(field_value.mises)
                # print(len(self.field_output_values))
                
                element_index = field_value.elementLabel
                intergration_point_sequence = field_value.integrationPoint
                # Remove the temp element
                if element_index > len(self.field_output_values):
                    continue
                field[element_index,0] = field_value.elementLabel
                # Output the stress
                # prettyprint(odb.rootAssembly.elementSets)
                # os._exit(0)
                if self.output_suboption == None:
                    # If the output_suboption is not specified, default output the von mises stress
                    # field[element_index,intergration_point_sequence] = self.field_output.getSubset(region=odb.rootAssembly.elementSets['Ele_Mech']).getScalarField(componentLabel='Mises')
                    field[element_index,intergration_point_sequence] = field_value.mises
                elif self.output_suboption == 'S11':
                    field[element_index,intergration_point_sequence] = field_value.data[0]
                elif self.output_suboption == 'S22':
                    field[element_index,intergration_point_sequence] = field_value.data[1]
                elif self.output_suboption == 'S33':
                    field[element_index,intergration_point_sequence] = field_value.data[2]   
            # **************************************************************************************
            # Output the displacement of the element
            # **************************************************************************************
            elif self.output_option == 'U':
                # prettyprint(field_value.mises)
                
                element_index = field_value.nodeLabel
                intergration_point_sequence = field_value.integrationPoint
                # Remove the temp element
                if element_index > len(self.field_output_values):
                    continue
                field[element_index,0] = field_value.nodeLabel
                # Output the displacement
                if self.output_suboption == None:
                    # If the output_suboption is not specified, default output the magnitude of the displacement
                    field[element_index,intergration_point_sequence] = field_value.magnitude
                elif self.output_suboption == 'U1':
                    field[element_index,intergration_point_sequence] = field_value.data[0]
                elif self.output_suboption == 'U2':
                    field[element_index,intergration_point_sequence] = field_value.data[1]
            # **************************************************************************************
            # Output the user defined output variable of the element
            # **************************************************************************************
            elif 'UVARM' in self.output_option:
                # 添加集合的选择
                
                # The UVARM variable is defined in the user defined elements, so there is maybe an offset in the element number
                # In the case of ith element and jth intergration point
                # The index of the element in the array is 
                #               elementLable - element_offset-1 (counting from 0 in python)
                element_index = field_value.elementLabel-element_offset-1
                intergration_point_sequence = field_value.integrationPoint
                # print(intergration_point_sequence)
                # print(np.shape(field))
                field[element_index,0] = field_value.elementLabel
                field[element_index,intergration_point_sequence] = field_value.data
            # **************************************************************************************
            # Output the state dependent variable of the element
            # **************************************************************************************
            elif 'SDV' in self.output_option:
                # The SDV variable is defined in normal elements, so there isn't offset in the element number
                element_index = field_value.elementLabel-1
                intergration_point_sequence = field_value.integrationPoint
                field[element_index,0] = field_value.elementLabel
                field[element_index,intergration_point_sequence] = field_value.data
        # print(field_value.elementLabel-1)
        return field
    
class output_field_process():
    def __init__(self, frame, output_option=None, field_array=None):
        self.frame = frame
        self.output_option = output_option
        self.field_array = field_array
    def get_element_weight(self):
        field_array = Get_field_output(self.frame, 'EVOL', field_vector_length=2).get_values_array()
        # Get the column that contains element volume and calculate the total volume
        EVOL_Total = np.sum(field_array[:,1])
        # Get the weight of every element
        self.Element_Weight = field_array[:,1]/EVOL_Total
        return self.Element_Weight
    
    def get_element_volume(self):
        field_array = Get_field_output(self.frame, 'EVOL', field_vector_length=2).get_values_array()
        # Get the column that contains element volume and calculate the total volume
        EVOL_Individual = field_array[:,1]
        EVOL_Total = np.sum(EVOL_Individual)
        return EVOL_Total
    
    def get_uvarm_value(self):
        field_array = Get_field_output(self.frame, self.output_option, field_vector_length=9).get_values_array()
        UVARM = field_array[:,1:]
        [ rows_num, columns_num] = np.shape(UVARM)
        UVARM_Element_Weighted = np.zeros([rows_num,1])
        field_vector_length = 9
        intergration_point_Weight = [1.0/(field_vector_length-1) for i in range(field_vector_length-1)]
        # print(columns_num)
        # print(intergration_point_Weight)
        for row in range(rows_num):
            for column in range(columns_num):
                # Summarize the UVARM of every integration point in element by the intergration point weight
                UVARM_Element_Weighted[row] = UVARM[row,column]*intergration_point_Weight[column] + UVARM_Element_Weighted[row]

        # Write the user defined output field data into a txt file
        UVARM_file_folder = os.path.join(extracted_data_folder,'UVARM')
        if os.path.exists(UVARM_file_folder) == False:
            os.mkdir(UVARM_file_folder)
        with open(UVARM_file_folder+'\%s_%s.txt'%(self.output_option, frame_count),'w') as file:
            for row in field_array:
                file.write(str(row)+'\n')
        Element_Weight = self.get_element_weight()

        UVARM_Total_Weighted = np.sum(UVARM_Element_Weighted[:,0]*Element_Weight)
        return UVARM_Total_Weighted
    
    def get_stress_value(self):
        field_array = Get_field_output(self.frame, self.output_option, field_vector_length=9).get_values_array()
        Stress = field_array[:,1:]
        [ rows_num, columns_num] = np.shape(Stress)
        Stress_Element_Weighted = np.zeros([rows_num,1])
        intergration_point_Weight = [1./4,1./4,1./4,1./4]
        # print(Stress)
        # os._exit(0)
        for row in range(rows_num):
            for column in range(columns_num):
                # print(column)
                # Summarize the UVARM of every integration point in element by the intergration point weight
                Stress_Element_Weighted[row] = Stress[row,column]*intergration_point_Weight[column] + Stress_Element_Weighted[row]
                # print(Stress_Element_Weighted)
                # os._exit(0)

        # Write the user defined output field data into a txt file
        Stress_file_folder = os.path.join(extracted_data_folder,'Stress')
        if os.path.exists(Stress_file_folder) == False:
            os.mkdir(Stress_file_folder)
        with open(Stress_file_folder+'\%s_%s.txt'%(self.output_option, frame_count),'w') as file:
            for row in field_array:
                file.write(str(row)+'\n')
        Element_Weight = self.get_element_weight()
        Stress_Total_Weighted = np.sum(Stress_Element_Weighted*Element_Weight)
        return Stress_Total_Weighted
    
    def get_displacement_value(self):
        field_array = Get_field_output(self.frame, self.output_option, output_suboption='U1', field_vector_length=3).get_values_array()
        Displacement = field_array[:,1:]
        [ rows_num, columns_num] = np.shape(Displacement)
        Displacement_Total_Weighted = np.sum(Displacement[:,1])/float(rows_num)
        # Write the displacement data into a txt file
        Displacement_file_folder = os.path.join(extracted_data_folder,'Displacement')
        if os.path.exists(Displacement_file_folder) == False:
            os.mkdir(Displacement_file_folder)
        with open(Displacement_file_folder+'\%s_%s.txt'%(self.output_option, frame_count),'w') as file:
            for row in field_array:
                file.write(str(row)+'\n')
        return Displacement_Total_Weighted

    
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
        # Plot the curve of each calculation result and then close the figure
        fig_sub = plt.figure(num=1,figsize=(6,6))
        plt.plot(self.data_x,self.data_y,linewidth=1)
        plt.title(self.title_name,fontsize=16, weight='bold', fontname='sans-serif')
        plt.xlabel(self.label_x,fontsize=16, weight='bold', fontname='sans-serif')
        plt.ylabel(self.label_y,fontsize=16, weight='bold', fontname='sans-serif')
        # Save the figure in each sub folder
        fig_path = os.path.join(extracted_data_folder,'%s.png'%self.export_file_name)
        plt.savefig(fig_path,dpi=300,bbox_inches='tight')
        plt.close(fig_sub)
        # Plot all curves in one figure and save the figure in the primary folder
        fig_Total = plt.figure(num=2,figsize=(6,6))
        plt.plot(self.data_x,self.data_y,linewidth=1,label=sub_folder)
        plt.title(self.title_name,fontsize=16, weight='bold', fontname='sans-serif')
        plt.xlabel(self.label_x,fontsize=16, weight='bold', fontname='sans-serif')
        plt.ylabel(self.label_y,fontsize=16, weight='bold', fontname='sans-serif')
        plt.legend()
        # plt.grid(True)
        # plt.show()
        # plt.show(block=True)
        plt.savefig('%s_Total.png'%self.export_file_name,dpi=300,bbox_inches='tight')
        logwrite('Plot the curve successfully!')


def extract_odb(working_directory):
    global extracted_data_folder
    # os.chdir(working_directory)
    extracted_data_folder = os.path.join(working_directory,'Extracted Results')
    if os.path.exists(extracted_data_folder) == False:
        os.mkdir(extracted_data_folder)
    # Check if the log file exists
    if os.path.exists('log.txt'):
        os.remove('log.txt')
    logwrite('*'*70 + '\n' +
             '*' + '{0:^68}'.format('Extracting Data from ODB') + '*' + '\n' +
             '*'*70)
    # *********************************************************************************************************************
    # Set the parameters
    # *********************************************************************************************************************
    # Parameters interpretation:
    # ----------------------------------------------------------------------------------------------------------------
    # odb_file_name: The name of the ODB file
    # output_option: The name of output field parameter
    #           Example: UVARM2, S, U, etc.
    # output_suboption: The name of output field subparameter
    #           Example: S: S11, S22, S33, S12, S13, S23 (If not specified, the Von Mises stress will be extracted)
    #                    U: U1, U2, U3
    # specified_step: The relative sequence number of the step
    #           If you want to extract the data in a specified step, please set the relative sequence number of the step
    #           Example: specified_step = 1 means the first step
    #                    specified_step = None(default) means all steps
    odb_file_name = 'Mod_MultiC.odb'
    output_option = 'UVARM3'
    output_suboption = None
    specified_step = 1
    # *********************************************************************************************************************
    # Display the information of ODB file and output field parameter
    logwrite('{0:<25s}'.format('ODB file name') + ':' + odb_file_name)
    logwrite('{0:<25s}'.format('Output field parameter') + ':' + output_option + '\n' +
             '-'*60)
    # Open the ODB file
    try:
        odb_file_path = os.path.join(working_directory,odb_file_name)
        odb = openOdb(odb_file_path,readOnly=True)
        logwrite("Open the ODB file successfully")
    except:
        logwrite("There is no ODB file named " + odb_file_name)
    # Extract the data in the specified step
    # prettyprint(len(odb.steps.keys()))
    # os._exit(0)
    # --------------------------------------------------------------------------
    # If the specified_step is None, the data in all steps will be extracted
    # We can obtain the step name by odb.steps.keys()
    # This code is to get string of step name
    step_num = len(odb.steps.keys())
    step_name = []
    if specified_step == None:
        step_num = len(odb.steps.keys())
        for i in range(step_num):
            # Add the step name into the list
            step_name.append(odb.steps.keys()[i])
    else:
        step_num = 1
        # Add the specified step name into the list
        step_name.append(odb.steps.keys()[specified_step-1])
    # --------------------------------------------------------------------------
    # *********************************************************************************************************************
    # Extract the data in the every step
    # *********************************************************************************************************************
    # Define a golobal counter for the frames
    global frame_count
    frame_count = 0
    for step_sequence in range(step_num):
        step = odb.steps[step_name[step_sequence]]
        if output_option == 'EVOL':
            field_vector_length = 2
        if 'SDV' in output_option:
            field_vector_length = 2
        if 'UVARM' in output_option:
            field_vector_length = 9
        if output_option == 'U':
            # 识别别是二维单元还是三维单元，二维单元只有两个分量，三维单元有三个分量，将这些结果都输出一下
            field_vector_length = 2
        if output_option == 'S':
            field_vector_length = 3
        # ----------------------------------------------------------------------------------------------------------------
        # Extract the data in the every frame
        # ----------------------------------------------------------------------------------------------------------------
        # Prepare the data matrix for plotting
        Data_plot = np.zeros([len(step.frames),2])
        # print(len(step.frames))
        # 把关键字输出一下
        # print(step.frames[0].fieldOutputs.keys())
        # os._exit(0)

        frame_number = len(step.frames)
        logwrite('Extract the data from every frame...')
        for frame in step.frames:
            
            current_time = frame.frameValue
            if output_option == 'EVOL':
                EVOL_Total = output_field_process(frame,output_option).get_element_volume()
                Variable_Total_Weighted = EVOL_Total
            elif 'UVARM' in output_option:
                UVARM_Total_Weighted = output_field_process(frame,output_option).get_uvarm_value()
                Variable_Total_Weighted = UVARM_Total_Weighted
            elif output_option == 'S':
                Stress_Total_Weighted = output_field_process(frame,output_option).get_stress_value()
                Variable_Total_Weighted = Stress_Total_Weighted
            elif output_option == 'U':
                Displacement_Total_Weighted = output_field_process(frame,output_option).get_displacement_value()
                Variable_Total_Weighted = Displacement_Total_Weighted

            Data_plot[frame_count,0] = current_time
            Data_plot[frame_count,1] = Variable_Total_Weighted
            frame_count += 1
            logwrite("Progress: Frame:%d" % frame_count)
            # print("Progress: [{0}] {1:.2f}%".format("="*(frame_count + 1), finish_percent))
        # print(Data_plot)
        # Plot the curve of field output variable vs. time
        try:
            plot_x_y_curve(data_vector=Data_plot, 
                        x_label='Time', y_label=output_option, title='%s-Time Curve'%output_option, 
                        file_name='%s-Time Curve'%output_option)
            save_as_file(data_vector=Data_plot, file_name='%s-Time'%output_option)
        except:
            [rows_num, columns_num] = np.shape(Data_plot)
            if columns_num != 2:
                logwrite('The data_vector array should have two columns')
            logwrite('Plot the curve failed')

if __name__ == '__main__':
    # *********************************************************************************************************************
    #                                      Using tutorial
    # *********************************************************************************************************************
    # Putting this program in the primary directory of the calculation, and using the following command to run it
    # Command line: 
    #               abaqus script=ODB_Extract_Batch.py
    # ---------------------------------------------------------------------------------------------------------------------
    # The program will extract the data in the every sub-folder
    # The directory structure is as follows:
    # Primary directory (working_directory)
    #   |---sub-folder1------------------------------------|
    #       |---Calculation1.inp
    #       |---Calculation1.odb
    #       |---Extracted Results(Folder)------------------|
    #           |---Calculation1-Time Curve.png
    #           |---Calculation1-Time.csv
    #   |---sub-folder2------------------------------------|
    #       |---Calculation2.inp
    #       |---Calculation2.odb
    #       |---Extracted Results(Folder)------------------|
    #           |---Calculation2-Time Curve.png
    #           |---Calculation2-Time.csv
    #   ...
    global sub_folder
    Calculation_folders = os.listdir(os.getcwd())
    os.chdir(os.getcwd())
    for sub_folder in Calculation_folders:
        if os.path.isdir(sub_folder) == True:
            working_directory = os.path.join(os.getcwd(),sub_folder)
            print(working_directory)
            extract_odb(working_directory)
            