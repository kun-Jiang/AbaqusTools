# -*- coding: utf-8 -*-
from odbAccess import *
from textRepr import prettyPrint as prettyprint
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import threading
import argparse
from collections import defaultdict

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
        # Get the type of the field output
        field_type = self.field_output.type
        self.field_output_values = self.field_output.values
        # Read the first element to get the offset of the element number
        first_element_num = self.field_output_values[0].elementLabel
        # Read the first node to get the offset of the node number
        first_node_num = self.field_output_values[0].nodeLabel
        elemental_type = False
        nodal_type = False
        if first_element_num != None:
            element_num_offset = first_element_num -1
            elemental_type = True
        else:
            element_num_offset = 0
        if first_node_num != None:
            node_num_offset = first_node_num - 1
            nodal_type = True
        else:
            node_num_offset = 0
        
        print(self.output_option)
        if elemental_type:
            # prettyprint(self.current_frame)
            logwrite('Elemental field output detected.\n' +
                     'Element number offset: %d' % element_num_offset)
            output_field = self.element_data_access()
        elif nodal_type:
            # prettyprint(self.current_frame)
            logwrite('Nodal field output detected.\n' +
                     'Node number offset: %d' % node_num_offset)
            output_field = self.nodal_data_access()
        return output_field
        os._exit(0)
        prettyprint(first_element_num)
        prettyprint(first_node_num)
        prettyprint(elemental_type)
        prettyprint(nodal_type)
        prettyprint(field_type)
        # field = self.field_output.getSubset(position=ELEMENT_NODAL).values
        field = self.field_output.getSubset(position=UNDEFINED_POSITION).values
        # 准备数据结构：每个 elementLabel 对应多个积分点数据
        element_data = defaultdict(list)
        element_type_map = {}
        os._exit(0)
        
        for i in field:  # field 是 frame.fieldOutputs['S'].values 这类数据
            eid = i.elementLabel
            element_data[eid].append(i.data)
            element_type_map[eid] = i.baseElementType  # 也可以记录单元类型

        # 将每个单元的数据堆叠成 numpy 数组，并按类型分组
        type_dict = defaultdict(list)

        for eid, gp_data in element_data.items():
            gauss_array = np.vstack(gp_data)  # shape: n_gauss × var_dim（如 8×6）
            element_type = element_type_map[eid]
            type_dict[element_type].append({
                "elementLabel": eid,
                "data": gauss_array
            })
            if eid == 2:
                # print(type_dict)
                # print(type_dict['C3D8T'])
                # for d in type_dict['C3D8T']:
                #     print(d)
                eid_list = []
                tensor_list = []
                print(type_dict.keys())
                for elem in type_dict[type_dict.keys()[0]]:
                    eid_list.append(elem['elementLabel'])
                    tensor_list.append(elem['data'])
                print(eid_list)
                print(tensor_list)
                eids = np.array(eid_list)  # shape: (n_elements,)
                tensor = np.stack(tensor_list, axis=0)  # shape: (n_elements, n_gauss, var_dim)
                print('-----')
                print(eids)
                print(tensor)
                print(np.shape(tensor))
                np.savez('tensor_data.npz', tensor=tensor, eids=eids)
                os._exit(0)
        print(type_dict)
        print(len(field))
        # prettyprint(field)
        prettyprint(field[0])
        prettyprint(self.field_output.locations[0].position)
        prettyprint(self.field_output.locations[0].sectionPoints)
        # 统计唯一单元
        element_labels =  set(loc.elementLabel for loc  in self.field_output.locations)
        print(len(element_labels))
        # print(f"Number of unique elements (from 'S'): {len(element_labels)}")

        # prettyprint(self.field_output.getSubset(position=ELEMENT_NODAL).values)
        # prettyprint(self.field_output_values[0])
        # prettyprint(self.field_output_values[0].data)
        # output_array = self.data_access(field)
        os._exit(0)
        return output_array
    
        
        
        
        
        
        
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
        if self.output_option in ['GRADT','GRADT1','GRADT2','GRADT3']:
            field = np.zeros((len(self.field_output_values)/4,self.field_vector_length))
        if self.output_option in ['HFL','HFL1','HFL2','HFL3']:
            field = np.zeros((len(self.field_output_values)/4, self.field_vector_length))
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
                    element_num_offset = field_value.elementLabel -1
                else:
                    element_num_offset = 0
                    
                if field_value.nodeLabel != None:
                    node_num_offset = field_value.nodeLabel - 1
                else:
                    node_num_offset = 0
                count += 1
            # **************************************************************************************
            # Output the volume of the element
            # **************************************************************************************
            if self.output_option == 'EVOL':
                # The EVOL variable is defined in the dummy element, so there is a offset in the element number
                element_index = field_value.elementLabel-element_num_offset-1
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
                element_index = field_value.elementLabel-element_num_offset-1
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
                #               elementLable - element_num_offset-1 (counting from 0 in python)
                element_index = field_value.elementLabel-element_num_offset-1
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
                element_index = field_value.elementLabel-element_num_offset-1
                intergration_point_sequence = field_value.integrationPoint
                # prettyprint(intergration_point_sequence)
                field[element_index,0] = field_value.elementLabel
                field[element_index,intergration_point_sequence] = field_value.data
            # **************************************************************************************
            # Output the nodal temperature of the element
            # **************************************************************************************
            elif 'NT' in self.output_option:
                node_label = field_value.nodeLabel - node_num_offset - 1
                field[node_label] = field_value.data
            
            # **************************************************************************************
            # Output the heat flux of the element
            # **************************************************************************************
            elif 'HFL' in self.output_option:
                # print(len(self.field_output_values)/4)
                # prettyprint(self.current_frame)
                # prettyprint(field_value)
                # os._exit(0)
                element_index = field_value.elementLabel-element_num_offset-1
                inpt = field_value.integrationPoint
                field[element_index,0] = field_value.elementLabel
                if self.output_index == None:
                    # If the output_index is not specified, default output the magnitude of the heat flux
                    field[element_index,inpt] = field_value.magnitude
            # **************************************************************************************
            # Output the temperature gradient of the element
            # **************************************************************************************
            elif 'GRADT' in self.output_option:
                element_index = field_value.elementLabel-element_num_offset-1
                inpt = field_value.integrationPoint
                field[element_index,0] = field_value.elementLabel
                if self.output_index == None:
                    # If the output_index is not specified, default output the magnitude of the temperature gradient
                    field[element_index,inpt] = field_value.magnitude
        # print(field_value.elementLabel-1)
        return field
    
    def element_data_access(self):
        # This function is used to access the field data associated with the output variable (e.g., EVOL, UVARM, SDV, etc.)
        data_type = self.field_output.type
        ele_types = self.field_output.baseElementTypes
        componentLabel = self.field_output.componentLabels
        variable_name = self.field_output.name
        logwrite('Output variable type: %s' % data_type)
        # print(ele_types)
        # prettyprint(componentLabel)
        # print(self.field_output.locations[0])
        logwrite('Output variable name: %s' % variable_name)
        # Parallel processing of the field data
        database = self.field_output.getSubset(position=UNDEFINED_POSITION).values
        field_data_list = []
        if str(data_type) == "SCALAR":
            for elemental_data in database:
                element_num = elemental_data.elementLabel
                field_value = elemental_data.data
                field_data_list.append([element_num, field_value])
        else:
            for elemental_data in database:
                # prettyprint(elemental_data)
                # print(getattr(elemental_data, 'data'))
                # os._exit(0)
                element_num = elemental_data.elementLabel
                # field_value = elemental_data.data
                field_value = elemental_data.magnitude
                if field_value != None:
                    field_data_list.append([element_num, field_value])
                else:
                    if self.output_index == None:
                        field_vector = elemental_data.data
                        field_data_list.append(np.concatenate((np.array([element_num]), field_vector)))
                    else:
                        field_value = getattr(elemental_data, self.output_index)
                        field_data_list.append([element_num, field_value])
                # print(field_data_list)
                # if element_num == 2:
                #     os._exit(0)
        # Convert the list to a numpy array
        field_data_array = np.array(field_data_list)
        # print(field_data_array[:10,:])
        # logwrite('Successed\n' +
        #          'The shape of the field data array is: %s' % str(np.shape(field_data_array)))
        # os._exit(0)
        return field_data_array
        
    def nodal_data_access(self):
        # This function is used to access the nodal field data associated with the output variable (e.g., U, NT, HFL, etc.)
        data_type = self.field_output.type
        ele_types = self.field_output.baseElementTypes
        componentLabel = self.field_output.componentLabels
        variable_name = self.field_output.name
        logwrite('Output variable type: %s' % data_type)
        # print(ele_types)
        # prettyprint(componentLabel)
        # print(self.field_output.locations[0])
        logwrite('Output variable name: %s' % variable_name)
        # Parallel processing of the field data
        database = self.field_output.getSubset(position=UNDEFINED_POSITION).values
        field_data_list = []
        if str(data_type) == "SCALAR":
            for nodal_data in database:
                # prettyprint(nodal_data)
                # os._exit(0)
                node_num = nodal_data.nodeLabel
                field_value = nodal_data.data
                field_data_list.append([node_num, field_value])
        else:
            for nodal_data in database:
                node_num = nodal_data.nodeLabel
                field_value = nodal_data.data
                field_data_list.append(np.concatenate((np.array([node_num]), field_value)))
        # Convert the list to a numpy array
        field_data_array = np.array(field_data_list)
        return field_data_array

class output_field_access():
    def __init__(self, frame, output_option, dimension=2, output_index=None, method='gauss_integration'):
        self.frame = frame
        self.output_option = output_option
        self.output_index = output_index
        # There are two methods to calculate the field value:
        # 1. gauss_integration: calculate the area/volume integral of the field value (accurate but slow)
        # 2. volume_average   : calculate the volume average of the field value (fast but not accurate)
        self.method = method
        
        self.dimension = dimension
        self.element_type = 'isoparametric'
        if self.dimension == 2:
            self.node_output_length = 1
            self.element_output_length = 5
            self.nodes_per_element = 4
            self.degree_of_freedom = 2
        elif self.dimension == 3:
            self.node_output_length = 1
            self.element_output_length = 9
            self.nodes_per_element = 8
            self.degree_of_freedom = 3
        # self.coords = self.getNodalCoordinates()
        # print(self.EleConnectivity)
    def getNodalCoordinates(self):
        nodal_coords_folder = os.path.join(extracted_data_folder, 'Nodal_coords')
        if not os.path.exists(nodal_coords_folder):
            os.mkdir(nodal_coords_folder)
        nodes = self.frame.fieldOutputs[self.output_option].values[0].instance.nodes
        current_coords = np.zeros((len(nodes), self.dimension))
        file_path = os.path.join(extracted_data_folder, 'Nodal_coords' ,'Nodal_Coordinates_{}.dat'.format(frame_count))
        initial_file_path = os.path.join(extracted_data_folder, 'Nodal_coords' ,'Nodal_Coordinates_0.dat')
        if frame_count == 0:
            origin_coords = np.zeros((len(nodes), self.dimension))
            for node in nodes:
                origin_coords[node.label-1] = node.coordinates[:self.dimension]
        else:
            origin_coords = np.loadtxt(initial_file_path)
        disp_fields = self.frame.fieldOutputs['U'].values
        # for node, disp in zip(nodes, disp_fields):
        for disp in disp_fields:
            current_coords[disp.nodeLabel-1] = origin_coords[disp.nodeLabel-1,:self.dimension] + disp.data
        np.savetxt(file_path, current_coords,fmt='%f')
        return current_coords
    
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
        file_path = os.path.join(extracted_data_folder,'Element_Connectivity.dat')
        np.savetxt(file_path, connectivity,fmt='%d')
        return connectivity
    
    def shape_function(self, xi, eta, zeta):
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
        elif self.dimension == 3:
            if self.element_type == 'isoparametric':
                N1 = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta)
                N2 = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta)
                N3 = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta)
                N4 = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta)
                N5 = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta)
                N6 = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta)
                N7 = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta)
                N8 = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta)
                
                dN_dxi = np.array([
                    # dN/dξ (ξ-derivative)
                    [
                        -0.125*(1-eta)*(1-zeta),   # N1
                        0.125*(1-eta)*(1-zeta),   # N2
                        0.125*(1+eta)*(1-zeta),   # N3
                        -0.125*(1+eta)*(1-zeta),   # N4
                        -0.125*(1-eta)*(1+zeta),   # N5
                        0.125*(1-eta)*(1+zeta),   # N6
                        0.125*(1+eta)*(1+zeta),   # N7
                        -0.125*(1+eta)*(1+zeta)    # N8
                    ],
                    # dN/dη (η-derivative)
                    [
                        -0.125*(1-xi)*(1-zeta),    # N1
                        -0.125*(1+xi)*(1-zeta),    # N2
                        0.125*(1+xi)*(1-zeta),    # N3
                        0.125*(1-xi)*(1-zeta),    # N4
                        -0.125*(1-xi)*(1+zeta),    # N5
                        -0.125*(1+xi)*(1+zeta),    # N6
                        0.125*(1+xi)*(1+zeta),    # N7
                        0.125*(1-xi)*(1+zeta)     # N8
                    ],
                    # dN/dζ (ζ-derivative) - 修正部分
                    [
                        -0.125*(1-xi)*(1-eta),    # N1
                        -0.125*(1+xi)*(1-eta),     # N2
                        -0.125*(1+xi)*(1+eta),     # N3
                        -0.125*(1-xi)*(1+eta),     # N4
                        0.125*(1-xi)*(1-eta),     # N5
                        0.125*(1+xi)*(1-eta),     # N6
                        0.125*(1+xi)*(1+eta),     # N7
                        0.125*(1-xi)*(1+eta)      # N8
                    ]
                ])
            return np.array([N1,N2,N3,N4,N5,N6,N7,N8]), dN_dxi

    def gauss_integ(self, field_values):
        self.EleConnectivity = self.getElementsConnectivity()
        if self.dimension == 2:
            if self.element_type == 'isoparametric':
                gauss_point = np.array([
                    [-1/np.sqrt(3),-1/np.sqrt(3), 0],
                    [1/np.sqrt(3),-1/np.sqrt(3), 0],
                    [1/np.sqrt(3),1/np.sqrt(3), 0],
                    [-1/np.sqrt(3),1/np.sqrt(3), 0]
                ])
                gauss_weight = np.array([1,1,1,1])
        elif self.dimension == 3:
            if self.element_type == 'isoparametric':
                gauss_point = np.array([
                    [-1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)],
                    [1/np.sqrt(3),-1/np.sqrt(3),-1/np.sqrt(3)],
                    [1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)],
                    [-1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)],
                    [-1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)],
                    [1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)],
                    [1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)],
                    [-1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]
                ])
                gauss_weight = np.array([
                    1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0
                ])
        
        coords = self.getNodalCoordinates()
        # # !!!! Need to be modified
        # # The connectivity of the element is hard coded here, because the first and second layer
        # # of the element connectivity are repeated, so only the third layer is used
        # connectivity = self.EleConnectivity[np.shape(self.EleConnectivity)[0]/3*2:,:]
        connectivity = self.EleConnectivity
        field_integ_total = 0
        field_integ_element = np.zeros(len(connectivity))
        for i,nodes in enumerate(connectivity):
            node_coords = np.zeros((self.nodes_per_element, self.degree_of_freedom))
            integ = 0
            for k,node in enumerate(nodes):
                node_coords[k] = coords[node-1,:self.degree_of_freedom]
            for j in range(len(gauss_point)):
                
                xi, eta, zeta = gauss_point[j]
                weight = gauss_weight[j]
                N, dN_dxi = self.shape_function(xi, eta, zeta)
                J = np.dot(dN_dxi,node_coords)
                # J = np.zeros((2,2))
                # for l in range(len(gauss_point)):
                #     J[0, 0] += dN_dxi[0, l] * node_coords[l, 0]
                #     J[0, 1] += dN_dxi[0, l] * node_coords[l, 1]
                #     J[1, 0] += dN_dxi[1, l] * node_coords[l, 0]
                #     J[1, 1] += dN_dxi[1, l] * node_coords[l, 1]
                detJ = np.linalg.det(J)
                integ += field_values[i,j] * detJ * weight
                # print(detJ)
                # print(weight)
                # print(dN_dxi)
                # os._exit(0)
            field_integ_total += integ

            field_integ_element[i] += integ
        return field_integ_element, field_integ_total
    
    
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
        print(EVOL_Individual.shape)
        print('Element volue: %f'%EVOL_Total)
        
        # os._exit(0)
        return EVOL_Total
    
    def get_integration_point_weight(self):
        field_array = Get_field_output(self.frame, 'IVOL', field_vector_length=2).get_values_array()
        self.integration_point_weight = field_array[:,1]
        return self.integration_point_weight

    def get_integration_point_volume(self):
        field_array = Get_field_output(self.frame, 'IVOL', field_vector_length=2).get_values_array()
        IVOL_Total = np.sum(field_array[:,1])
        print(field_array.shape)
        return IVOL_Total
    
    def get_uvarm_value(self):
        field_array = Get_field_output(self.frame, self.output_option, field_vector_length=self.element_output_length).get_values_array()
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
        # with open(UVARM_file_folder+'\%s_%s.dat'%(self.output_option, frame_count),'w') as file:
        #     for row in field_array:
        #         file.write(str(row)+'\n')
        Element_Weight = self.get_element_weight()
        UVARM_Total_Weighted = np.sum(UVARM_Element_Weighted*Element_Weight)
        return UVARM_Total_Weighted

    def get_sdv_value(self):
        SDV_value_at_nodes = Get_field_output(self.frame, self.output_option, field_vector_length=self.element_output_length).get_values_array()
        SDV_field = SDV_value_at_nodes[:,1:].reshape(len(SDV_value_at_nodes)/self.nodes_per_element, self.nodes_per_element)
        if self.method == 'gauss_integration':
            SDV_integ, SDV_total = self.gauss_integ(SDV_field)
            IVOL_Total = self.get_integration_point_volume()
            SDV_Average = SDV_total / IVOL_Total
            print('Total integrated %s: %f' %(self.output_option, SDV_total))
            print('Total integration point volume: %f' %IVOL_Total)
        elif self.method == 'volume_average':
            [ rows_num, columns_num] = np.shape(SDV_field)
            SDV_Element_Weighted = np.zeros([rows_num,1])
            intergration_point_Weight = [[1./4],[1./4],[1./4],[1./4]]
            # Weight the SDV of every integration point in element by the intergration point weight
            SDV_Element_Weighted = np.dot(SDV_field,intergration_point_Weight)
            Element_Weight = self.get_element_weight()
            SDV_Average = np.sum(SDV_Element_Weighted*Element_Weight)
        
        # Write the solution dependent variable output field data into file
        SDV_file_folder = os.path.join(extracted_data_folder, odb_name + '_' + self.output_option)
        if os.path.exists(SDV_file_folder) == False:
            os.mkdir(SDV_file_folder)
        # Write the data for each frame to the document
        np.savetxt(SDV_file_folder + '\%s_%s.dat' %(self.output_option, frame_count), 
            SDV_field, fmt='%-4.5e', delimiter=',')
        return SDV_Average

    def get_stress_value(self):
        keyword_map = {'11':1, '22':2, '33':3, '12':4, '13':5, '23':6}
        keyword = None
        output_index = self.output_index
        if self.output_index in keyword_map.keys():
            logwrite('Output stress component: %s' %self.output_index)
            keyword = keyword_map[self.output_index]
            output_index = None
        Stress_Tensor_at_integ_nodes = Get_field_output(self.frame, self.output_option, output_index, field_vector_length=self.element_output_length).get_values_array()
        if keyword != None:
            stress_field = Stress_Tensor_at_integ_nodes[:, keyword].reshape(len(Stress_Tensor_at_integ_nodes)/self.nodes_per_element, self.nodes_per_element)
        else:
            stress_field = Stress_Tensor_at_integ_nodes[:,1].reshape(len(Stress_Tensor_at_integ_nodes)/self.nodes_per_element, self.nodes_per_element)
        if self.method == 'gauss_integration':
            Stress_integ, Stress_total = self.gauss_integ(stress_field)
            IVOL_Total = self.get_integration_point_volume()
            Stress_Average = Stress_total / IVOL_Total
            print('Total integrated %s: %f' %(self.output_option, Stress_total))
            print('Total integration point volume: %f' %IVOL_Total)
        # Write the user defined output field data into a txt file
        Stress_file_folder = os.path.join(extracted_data_folder, odb_name + '_' + self.output_option)
        if os.path.exists(Stress_file_folder) == False:
            os.mkdir(Stress_file_folder)
        # Write the data for each frame to the document
        np.savetxt(Stress_file_folder + '\%s_%s.dat' %(self.output_option+self.output_index, frame_count), 
            stress_field, fmt='%-4.5e', delimiter=',')
        return Stress_Average
    
    def get_strain_value(self):
        keyword_map = {'11':1, '22':2, '33':3, '12':4, '13':5, '23':6}
        keyword = None
        output_index = self.output_index
        if self.output_index in keyword_map.keys():
            logwrite('Output stress component: %s' %self.output_index)
            keyword = keyword_map[self.output_index]
            output_index = None
        Strain_Tensor_at_integ_nodes = Get_field_output(self.frame, self.output_option, output_index, field_vector_length=self.element_output_length).get_values_array()
        if keyword != None:
            strain_field = Strain_Tensor_at_integ_nodes[:, keyword].reshape(len(Strain_Tensor_at_integ_nodes)/self.nodes_per_element, self.nodes_per_element)
        else:
            strain_field = Strain_Tensor_at_integ_nodes[:,1].reshape(len(Strain_Tensor_at_integ_nodes)/self.nodes_per_element, self.nodes_per_element)
        if self.method == 'gauss_integration':
            Strain_integ, Strain_total = self.gauss_integ(strain_field)
            IVOL_Total = self.get_integration_point_volume()
            Strain_Average = Strain_total / IVOL_Total
            print('Total integrated %s: %f' %(self.output_option, Strain_total))
            print('Total integration point volume: %f' %IVOL_Total)
        # Write the user defined output field data into a txt file
        Strain_file_folder = os.path.join(extracted_data_folder, odb_name + '_' + self.output_option)
        if os.path.exists(Strain_file_folder) == False:
            os.mkdir(Strain_file_folder)
        # Write the data for each frame to the document
        np.savetxt(Strain_file_folder + '\%s_%s.dat' %(self.output_option+self.output_index, frame_count), 
            strain_field, fmt='%-4.5e', delimiter=',')
        return Strain_Average
    
    def get_displacement_value(self):
        field_array = Get_field_output(self.frame, self.output_option, field_vector_length=self.node_output_length).get_values_array()
        Displacement = field_array[:,1:]
        [ rows_num, columns_num] = np.shape(Displacement)
        Displacement_Total_Weighted = np.sum(Displacement[:,1])/float(rows_num)
        # Write the displacement data into a txt file
        Displacement_file_folder = os.path.join(extracted_data_folder,'Displacement')
        if os.path.exists(Displacement_file_folder) == False:
            os.mkdir(Displacement_file_folder)
        with open(Displacement_file_folder+'\%s_%s.dat'%(self.output_option, frame_count),'w') as file:
            for row in field_array:
                file.write(str(row)+'\n')
        return Displacement_Total_Weighted
    
    def get_nt_value(self):
        # TODO: Replace with the elemental integration method
        Nodal_Temperature_at_nodes = Get_field_output(self.frame, self.output_option, field_vector_length=self.node_output_length).get_values_array()
        # Write the nodal temperature field data into file
        NT_file_folder = os.path.join(extracted_data_folder, odb_name + '_' + self.output_option)
        if os.path.exists(NT_file_folder) == False:
            os.mkdir(NT_file_folder)
        np.savetxt(NT_file_folder + '\%s_%s.dat'%(self.output_option, frame_count), Nodal_Temperature_at_nodes, fmt='%-4.5e', delimiter=',')
        
        return 0
    
    def get_et_value(self):
        Elemental_Temperature = Get_field_output(self.frame, self.output_option, field_vector_length=self.node_output_length).get_values_array()
        Elemental_Temperature_Average = np.mean(Elemental_Temperature[:,1])
        TEMP_file_folder = os.path.join(extracted_data_folder, odb_name + '_' + self.output_option)
        if os.path.exists(TEMP_file_folder) == False:
            os.mkdir(TEMP_file_folder)
        np.savetxt(TEMP_file_folder + '\%s_%s.dat'%(self.output_option, frame_count), Elemental_Temperature, fmt='%-4.5e', delimiter=',')
        return Elemental_Temperature_Average
    def get_temp_gradient_value(self):
        Temp_grad = Get_field_output(self.frame, self.output_option, field_vector_length=self.element_output_length).get_values_array()
        # Calculate the temperature gradient of every element by the gauss integration method
        if self.method == 'gauss_integration':
            Temp_grad_integ, Temp_grad_total = self.gauss_integ(Temp_grad[:,1:])
            EVOL_Total = self.get_element_volume()
            Temp_grad_Average = Temp_grad_total / EVOL_Total
        elif self.method == 'volume_average':
            [rows_num, columns_num] = np.shape(Temp_grad)
            Temp_grad_integ = np.zeros([rows_num,1])
            intergration_point_Weight = [[1./4],[1./4],[1./4],[1./4]]
            # Weight the temperature gradient of every integration point in element by the intergration point weight
            Temp_grad_integ = np.dot(Temp_grad[:,1:],intergration_point_Weight)
            Element_Weight = self.get_element_weight()
            Temp_grad_total = np.sum(Temp_grad_integ*Element_Weight)
        # Save the temperature gradient data into a .dat file
        Temp_grad_file_folder = os.path.join(extracted_data_folder,'Temperature_Gradient_' + odb_name)
        if os.path.exists(Temp_grad_file_folder) == False:
            os.mkdir(Temp_grad_file_folder)
        np.savetxt(Temp_grad_file_folder + '\%s_%s.dat'%(self.output_option, frame_count), Temp_grad, fmt='%-4.5e', delimiter=',')
        return Temp_grad_Average

    def get_hflux_value(self):
        Heat_flux = Get_field_output(self.frame, self.output_option, field_vector_length=self.element_output_length).get_values_array()
        # Calculate the heat flux of every element by the gauss integration method
        # print(Heat_flux.shape)
        # print(Heat_flux[:8,:])
        flux_field = Heat_flux[:,1:].reshape(len(Heat_flux)/self.nodes_per_element, self.nodes_per_element)
        # print(flux_field.shape)
        # print(flux_field[0,:])
        if self.method == 'gauss_integration':
            Heat_flux_integ, Heat_flux_total = self.gauss_integ(flux_field)
            IVOL_Total = self.get_integration_point_volume()
            Total_HFL_Average = Heat_flux_total / IVOL_Total
        elif self.method == 'volume_average':
            [rows_num, columns_num] = np.shape(Heat_flux)
            Heat_flux_integ = np.zeros([rows_num,1])
            intergration_point_Weight = [[1./4],[1./4],[1./4],[1./4]]
            # Weight the heat flux of every integration point in element by the intergration point weight
            Heat_flux_integ = np.dot(Heat_flux,intergration_point_Weight)
            Element_Weight = self.get_element_weight()
            Heat_flux_total = np.sum(Heat_flux_integ*Element_Weight)
        # Save the heat flux data into a .dat file
        HFL_file_folder = os.path.join(extracted_data_folder,'Heat_Flux_' + odb_name)
        if os.path.exists(HFL_file_folder) == False:
            os.mkdir(HFL_file_folder)
        np.savetxt(HFL_file_folder + '\%s_%s.dat'%(self.output_option, frame_count), Heat_flux, fmt='%-4.5e', delimiter=',')
        
        return Total_HFL_Average
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

def Output_filed_result(frame,output_option, output_index, dimension):
    if output_option == 'EVOL':
        EVOL_Total = output_field_access(frame,output_option, dimension, output_index).get_element_volume()
        return EVOL_Total
    elif output_option == 'IVOL':
        IVOL_Total = output_field_access(frame,output_option, dimension, output_index).get_integration_point_volume()
        return IVOL_Total
    elif 'UVARM' in output_option:
        UVARM_Total_Weighted = output_field_access(frame,output_option, dimension, output_index).get_uvarm_value()
        return UVARM_Total_Weighted
    elif 'SDV' in output_option:
        SDV_Total_Weighted = output_field_access(frame,output_option, dimension, output_index).get_sdv_value()
        return SDV_Total_Weighted
    elif output_option in ['S','S11','S22','S33','S12','S13','S23','SEQV','Mises']:
        Stress_Total_Weighted = output_field_access(frame,output_option, dimension, output_index).get_stress_value()
        return Stress_Total_Weighted
    elif output_option in ['U','U1','U2','U3']:
        Displacement_Total_Weighted = output_field_access(frame,output_option, dimension, output_index).get_displacement_value()
        return Displacement_Total_Weighted
    elif output_option in ['LE', 'E']:
        Linear_Strain_Total_Weighted = output_field_access(frame,output_option, dimension, output_index).get_strain_value()
        return Linear_Strain_Total_Weighted
    elif 'NT' in output_option:
        Nodal_Temperature_Total_Weighted = output_field_access(frame, output_option, dimension, output_index).get_nt_value()
        return Nodal_Temperature_Total_Weighted
    elif 'TEMP' in output_option:
        Element_Temperature_Total_Weighted = output_field_access(frame,output_option, dimension, output_index).get_et_value()
        return Element_Temperature_Total_Weighted
    elif 'HFL' in output_option:
        Heat_Flux_Total_Weighted = output_field_access(frame,output_option, dimension, output_index).get_hflux_value()
        return Heat_Flux_Total_Weighted
    elif 'GRADT' in output_option:
        Temp_grad_Total_Weighted = output_field_access(frame,output_option, dimension, output_index).get_temp_gradient_value()
        return Temp_grad_Total_Weighted


def extract_odb(odb_dir, odb_file_name, output_option, output_index, dimension, 
                specified_step, specified_frame,
                xmin, xmax, ymin, ymax):
    global extracted_data_folder
    global odb_name
    os.chdir(odb_dir)
    # Create the folder to store the extracted data
    extracted_data_folder = os.path.join(odb_dir,'Extracted Results')
    odb_name = odb_file_name.split('.')[0]  # Get the name of the ODB file without the extension
    if os.path.exists(extracted_data_folder) == False:
        os.mkdir(extracted_data_folder)
    # *********************************************************************************************************************
    logwrite('*'*70 + '\n' +
             '*' + '{0:^68}'.format('Extracting Data from ODB') + '*' + '\n' +
             '*'*70)
    files_list = os.listdir(odb_dir)
    Remove_file(odb_dir, 'lck')
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
        if output_option in ['S', 'U', 'LE', 'E']:
            logwrite('{0:<25s}'.format('Output option') + ':' + output_option + '\n' + 
                     '{0:<25s}'.format('Output index') + ':' + output_index + '\n' + '-'*60)
            pass
        else:
            output_option = output_option + output_index
            output_index = None
            logwrite('{0:<25s}'.format('Output field parameter') + ':' + output_option + '\n' +
                    '-'*60)

    # Open the ODB file
    try:
        odb_file_path = os.path.join(odb_dir,odb_file_name)
        logwrite('Odb file path: ' + odb_file_path)
        odb = openOdb(odb_file_path, readOnly=True)
        logwrite("Open the ODB file successfully")
    except:      
        try:
            # The ODB file may be created by the old version of Abaqus, so we need to upgrade the ODB file
            upgrade_odb_file_name = 'upgraded_' + odb_file_name
            upgrade_odb_path = os.path.join(odb_dir,upgrade_odb_file_name)
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
    # Data_plot = np.zeros((frame_total,2))
    Data_plot = []  # Use a list to store the data
    for i, name in enumerate(steps_name):
        logwrite('{0:=^70}'.format('Step: %s'%name))
        step = odb.steps[name]
        keywords = step.frames[0].fieldOutputs.keys()
        logwrite('Available field outputs: %s' % keywords)
        # ----------------------------------------------------------------------------------------------------------------
        # Extract the data in the every frame
        # ----------------------------------------------------------------------------------------------------------------
        if specified_frame == 0:
            # TODO: Use the append method to store the data in the 'Data_plot' array
            logwrite('Extract the data from every frame...')
            for frame in step.frames:
                current_time = frame.frameValue + step_time
                # Extract the data corresponding to the output option
                logwrite('{0:-^70}'.format('Frame: %d'%frame.frameId))
                Weighted_Variable = Output_filed_result(frame, output_option, output_index, dimension)
                # prettyprint(frame.frameValue)
                # prettyprint(frame.frameId)
                prettyprint(frame.description + '\n')
                logwrite("Progress: Frame:%d\t%.2f s" %(frame_count,current_time))
                Data_plot.append(np.array([current_time, Weighted_Variable]))
                
                logwrite('The weighted %s is %.8f' %(output_option, Weighted_Variable))
                frame_count += 1
                # TODO: Replot the figure every time a frame is processed
            step_time = current_time
        # ----------------------------------------------------------------------------------------------------------------
        # Extract the data in the last frame
        # ----------------------------------------------------------------------------------------------------------------
        elif specified_frame == -1:
            Data_plot = None
            logwrite('Extract the data from the last frame...')
            Weighted_Variable = Output_filed_result(step.frames[-1],output_option, output_index, dimension)
            if Weighted_Variable != None:
                logwrite('The weighted %s is %.8f' %(output_option, Weighted_Variable))
        # ----------------------------------------------------------------------------------------------------------------
        # Extract the data in the specified frame
        # ----------------------------------------------------------------------------------------------------------------
        else:
            Data_plot = None
            logwrite('{0:-^70}'.format('Frame: %d'%step.frames[specified_frame-1].frameId))
            logwrite('Extract the data from the ' + str(specified_frame) + 'th frame...')
            Weighted_Variable = Output_filed_result(step.frames[specified_frame-1], output_option, output_index, dimension)
            if Weighted_Variable != None:
                logwrite('The weighted %s is %.8f' %(output_option, Weighted_Variable))
    # *********************************************************************************************************************
    # Plot the curve of the output field parameter w.r.t. time
    # *********************************************************************************************************************
    Data_plot = np.array(Data_plot)  # Convert the list to a numpy array
    if Data_plot.any():
        plot_x_y_curve(xmin, xmax, ymin, ymax, data_array=Data_plot, 
                    x_label='Time', y_label=output_option, title='%s-Time Curve'%output_option, 
                    file_name=odb_file_name+'-%s-Time Curve'%output_option)
        np.savetxt(os.path.join(extracted_data_folder,odb_file_name+'%s-Time.dat'%output_option),Data_plot,fmt='%-4.5e',delimiter=',')

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
    [odb_dir, odb_file]= os.path.split(odb_file_path)
    # Check if the log file exists
    # logfile_path = os.path.join(extracted_data_folder,'extract_odb.log')
    # if os.path.exists(logfile_path):
    #     os.remove(logfile_path)
    extract_odb(odb_dir, odb_file, output_option, output_index, dimension, 
                specified_step, specified_frame,
                xmin, xmax, ymin, ymax)
  
