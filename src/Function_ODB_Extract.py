# -*- coding: utf-8 -*-
from odbAccess import *
from textRepr import prettyPrint as prettyprint
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import threading
import argparse

def logwrite(message):
    with open('extract_odb.log','a') as log_file:
        log_file.write(message+'\n')
        print(message)
        log_file.close()
        
def Remove_file(folder, extension):
    for file in os.listdir(folder):
        if file.endswith(extension):
            logwrite('Remove file: %s'%file)
            os.remove(os.path.join(folder, file))
            
def save_as_file(data_vector, file_name):
    file_path = os.path.join(extracted_data_folder, '%s.dat'% file_name)
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
        # prettyprint(self.current_frame.fieldOutputs.keys())
        # os._exit(0)
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
            # print(len(self.field_output_values))
            field = np.zeros((len(self.field_output_values),self.field_vector_length))
        if 'UVARM' in self.output_option or 'SDV' in self.output_option:
            # Values of UVARM and SDV are stored in the intergration points, 
            # so the length of the array is 4 (node quantity per element. e.g. 4 for CPS4) times of the number of elements
            # The first column is the element label
            # The second column is the value in the first intergration point
            # The third column is the value in the second intergration point
            # The fourth column is the value in the third intergration point
            # The fifth column is the value in the fourth intergration point
            # ......
            field = np.zeros((len(self.field_output_values)/4, self.field_vector_length))
        # --------------------------------Stress field--------------------------------
        if self.output_option in ['S','S11','S22','S33','S12','S13','S23','SEQV','Mises']:
            field = np.zeros((len(self.field_output_values)/4,self.field_vector_length))
        # -----------------------------Displacement field-----------------------------
        if self.output_option in ['U','U1','U2','U3']:
            field = np.zeros((len(self.field_output_values)/4,self.field_vector_length))
            pass
        if self.output_option in ['LE','E']:
            field = np.zeros((len(self.field_output_values)/4,self.field_vector_length))
        # **********************************************************************
        # Obtain the output field values and store them in the array
        # There are kinds of output field, such as UVARM, EVOL, SDV, etc. 
        # Their data structure are different which need to be treated differently
        # **********************************************************************
        count = 0
        for field_value in self.field_output_values:
            # prettyprint(field_value)
            # Check the element number of the first element (count == 0)
            if count == 0:
                # As usual, the number of the first element is 1, 
                # if the number of the first element is not 1, there exist a offset
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
                
                element_index = field_value.elementLabel - 1
                intergration_point_sequence = field_value.integrationPoint
                # Pass by the temporary element
                if element_index > len(self.field_output_values):
                    continue
                # print(element_index)
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
            elif self.output_option in ['LE','E']:
                # prettyprint(field_value)
                # os._exit(0)
                element_index = field_value.elementLabel-element_offset-1
                intergration_point_sequence = field_value.integrationPoint
                field[element_index,0] = field_value.elementLabel
                # prettyprint(field_value.data)
                # os._exit(0)
                # field_value.data[0] is LE11
                # field_value.data[1] is LE22
                # field_value.data[2] is LE33
                # field_value.data[3] is LE12
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
                element_index = field_value.elementLabel-element_offset-1
                intergration_point_sequence = field_value.integrationPoint
                # prettyprint(intergration_point_sequence)
                field[element_index,0] = field_value.elementLabel
                field[element_index,intergration_point_sequence] = field_value.data
        # print(field_value.elementLabel-1)
        return field
    
class output_field_access():
    def __init__(self, frame, output_option, output_suboption=None, method='gauss_integration'):
        self.frame = frame
        self.output_option = output_option
        self.output_suboption = output_suboption
        # There are two methods to calculate the field value:
        # 1. gauss_integration: calculate the area/volume integral of the field value (accurate but slow)
        # 2. volume_average   : calculate the volume average of the field value (fast but not accurate)
        self.method = method
        
        self.dimension = 2
        self.element_type = 'isoparametric'
        self.coords = self.getNodalCoordinates()
        self.EleConnectivity = self.getElementsConnectivity()
        # print(self.EleConnectivity)
    def getNodalCoordinates(self):
        nodes = self.frame.fieldOutputs[self.output_option].values[0].instance.nodes
        coords = np.zeros((len(nodes),3))
        for node in nodes:
            coords[node.label-1] = node.coordinates
        np.savetxt('coords.txt',coords,fmt='%f')
        return coords
    
    def getElementsConnectivity(self):
        elements = self.frame.fieldOutputs[self.output_option].values[0].instance.elements
        if self.dimension == 2:
            if self.element_type == 'isoparametric':
                connectivity = np.zeros((len(elements),4))
        for element in elements:
            connectivity[element.label-1] = element.connectivity
        connectivity = connectivity.astype(int)
        np.savetxt('connectivity.txt',connectivity,fmt='%d')
        return connectivity
    
    def gauss_integ(self, field_values):
        if self.dimension == 2:
            if self.element_type == 'isoparametric':
                gauss_point = np.array([
                    [-1/np.sqrt(3),-1/np.sqrt(3)],
                    [1/np.sqrt(3),-1/np.sqrt(3)],
                    [1/np.sqrt(3),1/np.sqrt(3)],
                    [-1/np.sqrt(3),1/np.sqrt(3)]
                ])
                gauss_weight = np.array([1,1,1,1])
        
        
        coords = self.coords
        # !!!! Need to be modified
        # The connectivity of the element is hard coded here, because the first and second layer
        # of the element connectivity are repeated, so only the third layer is used
        connectivity = self.EleConnectivity[np.shape(self.EleConnectivity)[0]/3*2:,:]
        field_integ_total = 0
        field_integ_element = np.zeros(len(connectivity))
        for i,nodes in enumerate(connectivity):
            node_coords = np.zeros((4,2))
            for k,node in enumerate(nodes):
                node_coords[k] = coords[node-1,:2]
            for j in range(len(gauss_point)):
                integ = 0
                xi, eta = gauss_point[j]
                weight = gauss_weight[j]
                N, dN_dxi = self.shape_function(xi, eta)
                # J = np.dot(dN_dxi[:,j],node_coords)
                J = np.zeros((2,2))
                for l in range(len(gauss_point)):
                    J[0, 0] += dN_dxi[0, l] * node_coords[l, 0]
                    J[0, 1] += dN_dxi[0, l] * node_coords[l, 1]
                    J[1, 0] += dN_dxi[1, l] * node_coords[l, 0]
                    J[1, 1] += dN_dxi[1, l] * node_coords[l, 1]
                detJ = np.linalg.det(J)
                integ = field_values[i,j] * detJ * weight
                field_integ_total += integ

                field_integ_element[i] += integ
        return field_integ_element, field_integ_total
    
    def shape_function(self, xi, eta):
        if self.dimension == 2:
            if self.element_type == 'isoparametric':
                N1 = 0.25 * (1 - xi) * (1 - eta)
                N2 = 0.25 * (1 + xi) * (1 - eta)
                N3 = 0.25 * (1 + xi) * (1 + eta)
                N4 = 0.25 * (1 - xi) * (1 + eta)
                
                dN_dxi = np.array([
                    [-0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta)],
                    [-0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)]
                ])
        return np.array([N1,N2,N3,N4]), dN_dxi
    
    def get_element_weight(self):
        field_array = Get_field_output(self.frame, 'EVOL', field_vector_length=2).get_values_array()
        # Get the column that contains element volume and calculate the total volume
        EVOL_Total = np.sum(field_array[:,1])
        # Get the weight of every element
        self.Element_Weight = field_array[:,1:2]/EVOL_Total
        return self.Element_Weight
    
    def get_element_volume(self):
        field_array = Get_field_output(self.frame, 'EVOL', field_vector_length=2).get_values_array()
        # Get the column that contains element volume and calculate the total volume
        EVOL_Individual = field_array[:,1]
        EVOL_file_path = os.path.join(extracted_data_folder,'Element_Volume.txt')
        with open(EVOL_file_path,'w') as EVOL_file:
            EVOL_file.write('Element number'+'\t'+'Element volume'+'\n')
            for i in range(len(EVOL_Individual)):
                EVOL_file.write(str(int(field_array[i,0]))+'\t'+str(field_array[i,1])+'\n')
            
        #     for i in range(len(EVOL_Individual)):
        #         f.write(str(EVOL_Individual[i])+'\n')
        EVOL_Total = np.sum(EVOL_Individual)
        return EVOL_Total
    
    def get_uvarm_value(self):
        field_array = Get_field_output(self.frame, self.output_option, field_vector_length=5).get_values_array()
        UVARM = field_array[:,1:]
        [ rows_num, columns_num] = np.shape(UVARM)
        UVARM_Element_Weighted = np.zeros([rows_num,1])
        intergration_point_Weight = [[1./4],[1./4],[1./4],[1./4]]
        # Weight the UVARM of every integration point in element by the intergration point weight
        UVARM_Element_Weighted = np.dot(UVARM,intergration_point_Weight)
        # Write the user defined output field data into a txt file
        # UVARM_file_folder = os.path.join(extracted_data_folder,'UVARM')
        # if os.path.exists(UVARM_file_folder) == False:
        #     os.mkdir(UVARM_file_folder)
        # with open(UVARM_file_folder+'\%s_%s.txt'%(self.output_option, frame_count),'w') as file:
        #     for row in field_array:
        #         file.write(str(row)+'\n')
        Element_Weight = self.get_element_weight()
        UVARM_Total_Weighted = np.sum(UVARM_Element_Weighted*Element_Weight)
        return UVARM_Total_Weighted

    def get_sdv_value(self):
        field_array = Get_field_output(self.frame, self.output_option, field_vector_length=5).get_values_array()
        SDV = field_array[:,1:]
        if self.method == 'gauss_integration':
            SDV_Element_Weighted, SDV_total = self.gauss_integ(SDV)
        elif self.method == 'volume_average':
            [ rows_num, columns_num] = np.shape(SDV)
            SDV_Element_Weighted = np.zeros([rows_num,1])
            intergration_point_Weight = [[1./4],[1./4],[1./4],[1./4]]
            # Weight the SDV of every integration point in element by the intergration point weight
            SDV_Element_Weighted = np.dot(SDV,intergration_point_Weight)
            Element_Weight = self.get_element_weight()
            SDV_total = np.sum(SDV_Element_Weighted*Element_Weight)
        
        # Write the solution dependent variable output field data into file
        SDV_file_folder = os.path.join(extracted_data_folder,'SDV')
        if os.path.exists(SDV_file_folder) == False:
            os.mkdir(SDV_file_folder)
        # Write the data for each frame to the document
        # np.savetxt(SDV_file_folder+'\%s_%s.txt'%(self.output_option, frame_count), 
        #     SDV_Element_Weighted, fmt='%-4.5e', delimiter=',')
        return SDV_total

    def get_stress_value(self):
        field_array = Get_field_output(self.frame, self.output_option, field_vector_length=5).get_values_array()
        Stress = field_array[:,1:]
        [rows_num, columns_num] = np.shape(Stress)
        Stress_Element_Weighted = np.zeros([rows_num,1])
        intergration_point_Weight = [[1./4],[1./4],[1./4],[1./4]]
        # Weight the stress of every integration point in element by the intergration point weight
        Stress_Element_Weighted = np.dot(Stress,intergration_point_Weight)
        # Write the user defined output field data into a txt file
        Stress_file_folder = os.path.join(extracted_data_folder,'Stress')
        if os.path.exists(Stress_file_folder) == False:
            os.mkdir(Stress_file_folder)
        # Write the data for each frame to the document
        # with open(Stress_file_folder+'\%s_%s.txt'%(self.output_option, frame_count),'w') as file:
        #     for row in field_array:
        #         file.write(str(row)+'\n')
        Element_Weight = self.get_element_weight()
        Stress_Total_Weighted = np.sum(Stress_Element_Weighted*Element_Weight)
        return Stress_Total_Weighted
    
    def get_strain_value(self):
        field_array = Get_field_output(self.frame, self.output_option, field_vector_length=5).get_values_array()
        Strain = field_array[:,1:]
        [rows_num, columns_num] = np.shape(Strain)
        Strain_Element_Weighted = np.zeros([rows_num,1])
        intergration_point_Weight = [[1./4],[1./4],[1./4],[1./4]]
        # Weight the strain of every integration point in element by the intergration point weight
        Strain_Element_Weighted = np.dot(Strain,intergration_point_Weight)
        # Write the user defined output field data into a txt file
        Strain_file_folder = os.path.join(extracted_data_folder,'Strain')
        if os.path.exists(Strain_file_folder) == False:
            os.mkdir(Strain_file_folder)
        with open(Strain_file_folder+'\%s_%s.txt'%(self.output_option, frame_count),'w') as file:
            for row in field_array:
                file.write(str(row)+'\n')
        Element_Weight = self.get_element_weight()
        Strain_Total_Weighted = np.sum(Strain_Element_Weighted*Element_Weight)
        return Strain_Total_Weighted
    
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
    def __init__(self,xmin, xmax, ymin, ymax,data_vector,x_label,y_label,title,file_name):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
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
        if (self.xmax - self.xmax) != 0:
            plt.xlim(self.xmin, self.xmax)  # 设置 x 轴范围
        if (self.ymax - self.ymin) != 0:
            plt.ylim(self.ymin, self.ymax)  # 设置 y 轴范围 
        # Save the figure in each sub folder
        fig_path = os.path.join(extracted_data_folder,'%s.png'%self.export_file_name)
        plt.savefig(fig_path,dpi=300,bbox_inches='tight')
        plt.close(fig)
        plt.show()
        plt.show(block=True)

def Output_filed_result(frame,output_option, output_suboption=None):
    if output_option == 'EVOL':
        EVOL_Total = output_field_access(frame,output_option).get_element_volume()
        Variable_Total_Weighted = EVOL_Total
    elif 'UVARM' in output_option:
        UVARM_Total_Weighted = output_field_access(frame,output_option).get_uvarm_value()
        Variable_Total_Weighted = UVARM_Total_Weighted
    elif 'SDV' in output_option:
        SDV_Total_Weighted = output_field_access(frame,output_option).get_sdv_value()
        Variable_Total_Weighted = SDV_Total_Weighted
    elif output_option in ['S','S11','S22','S33','S12','S13','S23','SEQV','Mises']:
        Stress_Total_Weighted = output_field_access(frame,output_option).get_stress_value()
        Variable_Total_Weighted = Stress_Total_Weighted
    elif output_option in ['U','U1','U2','U3']:
        Displacement_Total_Weighted = output_field_access(frame,output_option).get_displacement_value()
        Variable_Total_Weighted = Displacement_Total_Weighted
    elif output_option in ['LE', 'E']:
        Linear_Strain_Total_Weighted = output_field_access(frame,output_option).get_strain_value()
        Variable_Total_Weighted = Linear_Strain_Total_Weighted
    return Variable_Total_Weighted

def extract_odb(working_directory, odb_file_name, output_option, output_suboption, dimension, 
                specified_step, specified_frame,
                xmin, xmax, ymin, ymax):
    global extracted_data_folder
    os.chdir(working_directory)

    extracted_data_folder = os.path.join(working_directory,'Extracted Results')
    if os.path.exists(extracted_data_folder) == False:
        os.mkdir(extracted_data_folder)
    # os.chdir(extracted_data_folder)

    logwrite('*'*70 + '\n' +
             '*' + '{0:^68}'.format('Extracting Data from ODB') + '*' + '\n' +
             '*'*70)
    files_list = os.listdir(working_directory)
    Remove_file(working_directory, 'lck')
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
    # odb_file_name = 'Model_RVE.odb'
    # output_option = 'UVARM1'
    # output_suboption = None
    # specified_step = 1
    # specified_frame = None
    # *********************************************************************************************************************
    # Display the information of ODB file and output field parameter
    logwrite('{0:<25s}'.format('ODB file name') + ':' + odb_file_name)
    logwrite('{0:<25s}'.format('Output field parameter') + ':' + output_option + '\t' + output_suboption + '\n' +
             '-'*60)
    # Open the ODB file
    try:
        odb_file_path = os.path.join(working_directory,odb_file_name)
        logwrite('Odb file path: ' + odb_file_path)
        odb = openOdb(odb_file_path,readOnly=True)
        logwrite("Open the ODB file successfully")
    except:
        try:
            upgrade_odb_file_name = 'upgraded_' + odb_file_name
            upgrade_odb_path = os.path.join(working_directory,upgrade_odb_file_name)
            logwrite("Fail to open the ODB file, try to upgrade the ODB file" + '\n' +
                     'New odb file path: ' + upgrade_odb_path)
            if os.path.exists(upgrade_odb_path):
                logwrite('The upgraded Odb file already exists.')
            else:
                upgradeOdb(existingOdbPath=odb_file_path, upgradedOdbPath=upgrade_odb_path)
            logwrite("Upgrade the ODB file successfully")
            odb_file_path = upgrade_odb_path
            odb = openOdb(odb_file_path, readOnly=True)
        except:
            logwrite('Fail to open the Odb file, please check it.')
    # Extract the data in the specified step
    # prettyprint(len(odb.steps.keys()))
    # os._exit(0)
    # --------------------------------------------------------------------------
    # If the specified_step is None, the data in all steps will be extracted
    # We can obtain the step name by odb.steps.keys()
    # This code is to get string of step name
    step_num = len(odb.steps.keys())
    step_name = []
    if specified_step == 0:
        step_num = len(odb.steps.keys())
        for i in range(step_num):
            # Add the step name into the list
            step_name.append(odb.steps.keys()[i])
        logwrite('Extract the data from every step...')
    else:
        step_num = 1
        # Add the specified step name into the list
        if specified_step == -1:
            logwrite('Extract the data from the last setp...')
            step_name.append(odb.steps.keys()[-1])
        else:
            logwrite('Extract the data from the ' + str(specified_frame) + 'th step...')
            step_name.append(odb.steps.keys()[specified_step-1])
    # *********************************************************************************************************************
    # Extract the data in the every step
    # *********************************************************************************************************************
    # Prepare the data matrix for plotting
    frame_total = 0
    for step_sequence in range(step_num):
        step = odb.steps[step_name[step_sequence]]
        frame_total = frame_total + len(step.frames)
    Data_plot = np.zeros([frame_total,2])
    current_time = np.zeros([frame_total,1])
    print("There are %s frames in total" % frame_total)
    # Prepare parameters related to output option
    if output_option == 'EVOL':
        field_vector_length = 2
    if 'SDV' in output_option:
        field_vector_length = 2
        output_option = output_option + output_suboption
    if 'UVARM' in output_option:
        field_vector_length = 5
        output_option = output_option + output_suboption
    if output_option in ['U','U1','U2','U3']:
        # 识别别是二维单元还是三维单元，二维单元只有两个分量，三维单元有三个分量，将这些结果都输出一下
        field_vector_length = 2
    if output_option in ['S','S11','S22','S33','S12','S13','S23','SEQV','Mises']:
        field_vector_length = 3
    if output_option in ['LE']:
        field_vector_length = 3
    # ------------------------------------------------------------------
    # Define a golobal counter for the frames
    global frame_count
    frame_count = 0
    for step_sequence in range(step_num):
        logwrite('{0:-^70}'.format('Step: %s'%step_name[step_sequence]))
        step = odb.steps[step_name[step_sequence]]
        # ----------------------------------------------------------------------------------------------------------------
        # Extract the data in the every frame
        # ----------------------------------------------------------------------------------------------------------------
        # print(len(step.frames))
        # 把关键字输出一下
        # print(step.frames[0].fieldOutputs.keys())
        # os._exit(0)
        frame_number = len(step.frames)
        if specified_frame == 0:
            logwrite('Extract the data from every frame...')
            for frame in step.frames:
                
                current_time[step_sequence] = frame.frameValue
                Variable_Total_Weighted=Output_filed_result(frame,output_option)
                # Data_plot[frame.incrementNumber,0] = current_time
                # Data_plot[frame.incrementNumber,1] = Variable_Total_Weighted
                time_total = np.sum(current_time)
                # # Remove the first frame due to the initial value of time is zero
                # if time_total == 0.0:
                #     continue
                try:
                    Data_plot[frame_count,0] = time_total
                    Data_plot[frame_count,1] = Variable_Total_Weighted
                except:
                    print('There is something wrong with the data matrix for plotting')
                frame_count += 1
                logwrite("Progress: Frame:%d\t%.2f s" % (frame_count,time_total) + '\t' + "Equi. value:%f"%Variable_Total_Weighted)
                # print("Progress: [{0}] {1:.2f}%".format("="*(frame_count + 1), finish_percent))
            # print(Data_plot)
            # Plot the curve of field output variable vs. time
            try:
                plot_x_y_curve(xmin, xmax, ymin, ymax, data_vector=Data_plot, 
                            x_label='Time', y_label=output_option, title='%s-Time Curve'%output_option, 
                            file_name=odb_file_name+'-%s-Time Curve'%output_option)
                save_as_file(data_vector=Data_plot, file_name=odb_file_name+'%s-Time'%output_option)
            except:
                [rows_num, columns_num] = np.shape(Data_plot)
                if columns_num != 2:
                    logwrite('The data_vector array should have two columns')
                logwrite('Plot the curve failed')
        else:
            # Get field output of specified frame
            if specified_frame == -1:
                logwrite('Extract the data from the last frame...')
                Variable_Total_Weighted = Output_filed_result(step.frames[-1],output_option)
                logwrite('The weighted %s of the specified frame is %f'%(output_option,Variable_Total_Weighted))
            else:
                logwrite('Extract the data from the ' + str(specified_frame) + 'th frame...')
                Variable_Total_Weighted = Output_filed_result(step.frames[specified_frame-1],output_option)
                logwrite('The weighted %s of the specified frame is %f'%(output_option,Variable_Total_Weighted))
    # *********************************************************************************************************************
    # The extraction is finished
    logwrite('{0:<25s}'.format('ODB file name') + ':' + odb_file_path + '\n' +
             '{0:<25s}'.format('Output field parameter') + ':' + output_option + '\t' + output_suboption + '\n' +
             '{0:*^70}'.format('Finished'))
          
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The script of spliting the ABAQUS inp file')
    parser.add_argument('-f', '--odb_file', type=str, metavar='', required=True, help='odb file path')
    parser.add_argument('-o', '--output', type=str, metavar='', required=True, help='output option')
    parser.add_argument('-so', '--suboption', type=str, metavar='', required=False, help='suboption')
    parser.add_argument('-d', '--dimension', type=int, metavar='', required=True, help='dimension of the model')
    parser.add_argument('-s', '--step', type=int, metavar='', required=False, help='specified step')
    parser.add_argument('-r', '--frame', type=int, metavar='', required=False, help='specified frame')
    parser.add_argument('-x1', '--xmin', type=float, metavar='', required=False, help='minimum value of x axis')
    parser.add_argument('-x2', '--xmax', type=float, metavar='', required=False, help='maximum value of x axis')
    parser.add_argument('-y1', '--ymin', type=float, metavar='', required=False, help='minimum value of y axis')
    parser.add_argument('-y2', '--ymax', type=float, metavar='', required=False, help='maximum value of y axis')
    args = parser.parse_args()
    # *********************************************************************************************************************
    odb_file_path = args.odb_file
    output_option = args.output
    output_suboption = args.suboption
    dimension = args.dimension
    specified_step = args.step
    specified_frame = args.frame
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    [working_directory, odb_file]= os.path.split(odb_file_path)
    # Check if the log file exists
    # logfile_path = os.path.join(extracted_data_folder,'extract_odb.log')
    # if os.path.exists(logfile_path):
    #     os.remove(logfile_path)
    extract_odb(working_directory, odb_file, output_option, output_suboption, dimension, 
                specified_step, specified_frame,
                xmin, xmax, ymin, ymax)
  
