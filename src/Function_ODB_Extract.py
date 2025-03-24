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
            
class Run_in_new_thread():
    def __init__(self, func):
        self.func = func
        # Create new thread to run the function
        threading.Thread(target=func).start()

class Get_field_output():
    
    def __init__(self, current_frame, output_option, output_index=None, field_vector_length=None):
        self.current_frame = current_frame
        self.output_option = output_option
        self.output_index = output_index
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
        if 'NT' in self.output_option:
            field = np.zeros(len(self.field_output_values))
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
                # prettyprint(field_value.elementLabel)
                if field_value.elementLabel != None:
                    # Some field outputs are defined on the node, so there isn't element label
                    element_offset = field_value.elementLabel -1
                else:
                    element_offset = 0
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
                if self.output_index == None:
                    # If the output_index is not specified, default output the von mises stress
                    # field[element_index,intergration_point_sequence] = self.field_output.getSubset(region=odb.rootAssembly.elementSets['Ele_Mech']).getScalarField(componentLabel='Mises')
                    field[element_index,intergration_point_sequence] = field_value.mises
                elif self.output_index == 'S11':
                    field[element_index,intergration_point_sequence] = field_value.data[0]
                elif self.output_index == 'S22':
                    field[element_index,intergration_point_sequence] = field_value.data[1]
                elif self.output_index == 'S33':
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
                if self.output_index == None:
                    # If the output_index is not specified, default output the magnitude of the displacement
                    field[element_index,intergration_point_sequence] = field_value.magnitude
                elif self.output_index == 'U1':
                    field[element_index,intergration_point_sequence] = field_value.data[0]
                elif self.output_index == 'U2':
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
            # **************************************************************************************
            # Output the nodal temperature of the element
            # **************************************************************************************
            elif 'NT' in self.output_option:
                node_label = field_value.nodeLabel - 1
                field[node_label] = field_value.data
                
        # print(field_value.elementLabel-1)
        return field
    
class output_field_access():
    def __init__(self, frame, output_option, dimension, output_index=None, method='gauss_integration'):
        self.frame = frame
        self.output_option = output_option
        self.output_index = output_index
        # There are two methods to calculate the field value:
        # 1. gauss_integration: calculate the area/volume integral of the field value (accurate but slow)
        # 2. volume_average   : calculate the volume average of the field value (fast but not accurate)
        self.method = method
        
        self.dimension = dimension
        self.element_type = 'isoparametric'
        self.coords = self.getNodalCoordinates()
        # print(self.EleConnectivity)
    def getNodalCoordinates(self):
        nodes = self.frame.fieldOutputs[self.output_option].values[0].instance.nodes
        coords = np.zeros((len(nodes),3))
        for node in nodes:
            coords[node.label-1] = node.coordinates
        file_path = os.path.join(extracted_data_folder,'Nodal_Coordinates.txt')
        np.savetxt(file_path,coords,fmt='%f')
        return coords
    
    def getElementsConnectivity(self):
        elements = self.frame.fieldOutputs[self.output_option].values[0].instance.elements
        if self.dimension == 2:
            if self.element_type == 'isoparametric':
                connectivity = np.zeros((len(elements),4))
        elif self.dimension == 3:
            if self.element_type == 'isoparametric':
                connectivity = np.zeros((len(elements),8))
        for element in elements:
            connectivity[element.label-1] = element.connectivity
        connectivity = connectivity.astype(int)
        file_path = os.path.join(extracted_data_folder,'Element_Connectivity.txt')
        np.savetxt(file_path, connectivity,fmt='%d')
        return connectivity
    
    def gauss_integ(self, field_values):
        self.EleConnectivity = self.getElementsConnectivity()
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
        field_array = Get_field_output(self.frame, self.output_option, output_index='U1', field_vector_length=3).get_values_array()
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
    
    def get_nt_value(self):
        Nodal_Temperature = Get_field_output(self.frame, self.output_option, field_vector_length=1).get_values_array()
        # Save the nodal temperature data into a .dat file
        NT_file_folder = os.path.join(extracted_data_folder,'Nodal_Temperature')
        if os.path.exists(NT_file_folder) == False:
            os.mkdir(NT_file_folder)
        np.savetxt(NT_file_folder + '\%s_%s.dat'%(self.output_option, frame_count), Nodal_Temperature, fmt='%-4.5e', delimiter=',')
        
        return float('nan')

    
class plot_x_y_curve():
    def __init__(self,xmin, xmax, ymin, ymax,data_array,x_label,y_label,title,file_name):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.data_x = data_array[:,0]
        self.data_y = data_array[:,1]
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

def Output_filed_result(frame,output_option, dimension):
    if output_option == 'EVOL':
        EVOL_Total = output_field_access(frame,output_option).get_element_volume()
        return EVOL_Total
    elif 'UVARM' in output_option:
        UVARM_Total_Weighted = output_field_access(frame,output_option).get_uvarm_value()
        return UVARM_Total_Weighted
    elif 'SDV' in output_option:
        SDV_Total_Weighted = output_field_access(frame,output_option).get_sdv_value()
        return SDV_Total_Weighted
    elif output_option in ['S','S11','S22','S33','S12','S13','S23','SEQV','Mises']:
        Stress_Total_Weighted = output_field_access(frame,output_option).get_stress_value()
        return Stress_Total_Weighted
    elif output_option in ['U','U1','U2','U3']:
        Displacement_Total_Weighted = output_field_access(frame,output_option).get_displacement_value()
        return Displacement_Total_Weighted
    elif output_option in ['LE', 'E']:
        Linear_Strain_Total_Weighted = output_field_access(frame,output_option).get_strain_value()
        return Linear_Strain_Total_Weighted
    elif 'NT' in output_option:
        Nodal_Temperature_Total_Weighted = output_field_access(frame, output_option, dimension).get_nt_value()
        return Nodal_Temperature_Total_Weighted


def extract_odb(working_directory, odb_file_name, output_option, output_index, dimension, 
                specified_step, specified_frame,
                xmin, xmax, ymin, ymax):
    global extracted_data_folder
    os.chdir(working_directory)
    # Create the folder to store the extracted data
    extracted_data_folder = os.path.join(working_directory,'Extracted Results')
    if os.path.exists(extracted_data_folder) == False:
        os.mkdir(extracted_data_folder)
    # *********************************************************************************************************************
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
    # output_index: The index of output field parameter
    #           Example: S: S11, S22, S33, S12, S13, S23 (If not specified, the Von Mises stress will be extracted)
    #                    U: U1, U2, U3
    # specified_step: The relative sequence number of the step
    #           If you want to extract the data in a specified step, please set the relative sequence number of the step
    #           Example: specified_step = 1 means the first step
    #                    specified_step = None(default) means all steps
    # odb_file_name = 'Model_RVE.odb'
    # output_option = 'UVARM1'
    # output_index = None
    # specified_step = 1
    # specified_frame = None
    # *********************************************************************************************************************
    # Display the information of ODB file and output field parameter
    logwrite('{0:<25s}'.format('ODB file name') + ':' + odb_file_name)
    if output_index != '0':
        output_option = output_option + output_index
    logwrite('{0:<25s}'.format('Output field parameter') + ':' + output_option + '\n' +
             '-'*60)
    # Open the ODB file
    try:
        odb_file_path = os.path.join(working_directory,odb_file_name)
        logwrite('Odb file path: ' + odb_file_path)
        odb = openOdb(odb_file_path, readOnly=True)
        logwrite("Open the ODB file successfully")
    except:      
        try:
            # The ODB file may be created by the old version of Abaqus, so we need to upgrade the ODB file
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
    # If the specified_step = 0, extract the data in every step
    #                       = -1, extract the data in the last step
    #                       = n, extract the data in the nth step
    # We can obtain the step name by odb.steps.keys()
    # This code is to get string of step name
    if specified_step == 0:
        steps_name = odb.steps.keys()
        step_num = len(steps_name)
        logwrite('Extract the data from every step...')
    else:
        step_num = 1
        # Add the specified step name into the list
        if specified_step == -1:
            steps_name = [odb.steps.keys()[-1]]
            logwrite('Extract the data from the last setp...')
        else:
            steps_name = [odb.steps.keys()[specified_step-1]]
            logwrite('Extract the data from the ' + str(specified_frame) + 'th step...')
    # *********************************************************************************************************************
    # Extract the data in the every step
    # *********************************************************************************************************************
    # Prepare the data matrix for plotting
    frame_total = 0
    for name in steps_name:
        step = odb.steps[name]
        frame_total = frame_total + len(step.frames)
    
    print("There are %s frames in total" % frame_total)
    # Prepare parameters related to output option
    # ------------------------------------------------------------------
    # Define a golobal counter for the frames
    global frame_count
    frame_count = 0
    step_time = 0
    for i, name in enumerate(steps_name):
        logwrite('{0:-^70}'.format('Step: %s'%name))
        step = odb.steps[name]
        keywords = step.frames[0].fieldOutputs.keys()
        logwrite('Available field outputs: %s' % keywords)
        # ----------------------------------------------------------------------------------------------------------------
        # Extract the data in the every frame
        # ----------------------------------------------------------------------------------------------------------------
        if specified_frame == 0:
            Data_plot = np.zeros((frame_total,2))
            logwrite('Extract the data from every frame...')
            for frame in step.frames:
                current_time = frame.frameValue + step_time
                # Extract the data corresponding to the output option
                Weighted_Variable = Output_filed_result(frame, output_option, dimension)
                logwrite("Progress: Frame:%d\t%.2f s" %(frame_count,current_time))
                
                Data_plot[frame_count,0] = current_time
                Data_plot[frame_count,1] = Weighted_Variable
                
                logwrite('The weighted %s is %.4f' %(output_option, Weighted_Variable))
                frame_count += 1
            step_time = current_time
        # ----------------------------------------------------------------------------------------------------------------
        # Extract the data in the last frame
        # ----------------------------------------------------------------------------------------------------------------
        elif specified_frame == -1:
            Data_plot = None
            logwrite('Extract the data from the last frame...')
            Weighted_Variable = Output_filed_result(step.frames[-1],output_option)
            if Weighted_Variable != None:
                logwrite('The weighted %s is %.4f' %(output_option, Weighted_Variable))
        # ----------------------------------------------------------------------------------------------------------------
        # Extract the data in the specified frame
        # ----------------------------------------------------------------------------------------------------------------
        else:
            Data_plot = None
            logwrite('Extract the data from the ' + str(specified_frame) + 'th frame...')
            Weighted_Variable = Output_filed_result(step.frames[specified_frame-1], output_option)
            if Weighted_Variable != None:
                logwrite('The weighted %s is %.4f' %(output_option, Weighted_Variable))
    # *********************************************************************************************************************
    # Plot the curve of the output field parameter w.r.t. time
    # *********************************************************************************************************************
    if Data_plot.any():
        plot_x_y_curve(xmin, xmax, ymin, ymax, data_array=Data_plot, 
                    x_label='Time', y_label=output_option, title='%s-Time Curve'%output_option, 
                    file_name=odb_file_name+'-%s-Time Curve'%output_option)
        np.savetxt(os.path.join(extracted_data_folder,odb_file_name+'%s-Time'%output_option),Data_plot,fmt='%-4.5e',delimiter=',')

    # *********************************************************************************************************************
    # The extraction is finished
    logwrite('{0:<25s}'.format('ODB file name') + ':' + odb_file_path + '\n' +
             '{0:<25s}'.format('Output field parameter') + ':' + output_option + '\n' +
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
    output_index = args.suboption
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
    extract_odb(working_directory, odb_file, output_option, output_index, dimension, 
                specified_step, specified_frame,
                xmin, xmax, ymin, ymax)
  
