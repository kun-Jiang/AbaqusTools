import os
import numpy as np
from utilities.logConfig import logconfig, logging

class Num_Offset:
    def __init__(self,input_file_path:str, Offset_Magnitude:int, layer_info:dict) -> None:
        # *******************************************************
        print("{0:=^80}".format("EleNum_Offset Process"))
        # *******************************************************
        #               Extract the inp file name
        # *******************************************************
        # Get the working directory: x:xx/xx/xx.inp
        self.input_file_path = input_file_path
        # Get the magnitude of the offset, default is 10000
        self.Offset_Magnitude = Offset_Magnitude
        # Remove the space before the string
        # E.g. Directory = ' D:/Desktop' --> 'D:/Desktop' 
        self.input_file_path = self.input_file_path.strip()
        # Get the layer information from the dictionary
        # E.g. layer_info = {'Mech':1, 'Chem':2}
        self.layer_info = layer_info
        # ***************************************************************************
        # Split the file path into directory and file name
        # E.g. Offset_origin_file_path = 'D:/Desktop/xx.inp'
        #      Directory = 'D:/Desktop', InpFile = 'xx.inp'
        [self.Origin_Directory,self.InpFile] = os.path.split(self.input_file_path) 
        # Create the folder to store the offseted inp file
        self.EleNum_Offset_folder = os.path.join(self.Origin_Directory,'EleNum_Offset')
        if os.path.exists(self.EleNum_Offset_folder) == False:
            os.mkdir(self.EleNum_Offset_folder)
        # ***************************************************************************
        #               Set the log file
        # ***************************************************************************
        logconfig('NumOffset.log',self.EleNum_Offset_folder)
        # ***************************************************************************
        #              Automatically calculate the offset magnitude
        # ***************************************************************************
        # If the offset magnitude is given as 0, then the offset magnitude will be
        # the maximum element number in the inp file.
        connectivity = np.loadtxt(self.input_file_path,delimiter=',',dtype=int)
        elements_num = connectivity[:,0]
        if self.Offset_Magnitude == 0:
            self.Offset_Magnitude = np.max(elements_num)
        elif self.Offset_Magnitude < np.max(elements_num):
            self.Offset_Magnitude = np.max(elements_num)
            logging.warning(
                "The offset magnitude is smaller than the maximum element"
                "number in the inp file, the offset magnitude will be set to the maximum element number.")
        # ***************************************************************************
        #               Adding the offset to the element number
        # ***************************************************************************
        # origin_file_folder : the folder of the origin inp file
        # Offset_Magnitude          : the magnitude of the offset
        # EleNum_Offset_folder      : the folder to store the offseted inp file
        # InpFile                   : the name of the origin inp file
        logging.info("{0:=^80}".format("Start EleNum_Offset Process"))
        logging.info("The input file is: %s"%self.input_file_path)
        logging.info("The output directory is: %s"%self.EleNum_Offset_folder)
        self.EleNum_Offset()
        # *******************************************************
        #               Finished
        # *******************************************************
        logging.info("{0:=^80}".format("Finished") + '\n')
    
    def EleNum_Offset(self):
        os.chdir(self.EleNum_Offset_folder)
        # ***************************************************************************
        # Extracting the inp file name
        InpFile_Name = self.InpFile.split('.inp')[0]
        # E.g. 'Eleset_Grain_broundraies' --> 'Eleset', 'Grain', 'broundraies'
        # prefix = 'Eleset'
        InpFile_prefix   = InpFile_Name.split('_')[0]
        # suffix = 'Grain', 'broundraies'
        InpFile_suffix_temp   = InpFile_Name.split('_')[1:]
        # There are usually more than one suffix in the element set name, so 
        # we need to add '_' to connect them.
        # E.g. 'Eleset_Grain_broundraies'
        InpFile_suffix = '_'.join(InpFile_suffix_temp)
        # ***************************************************************************
        # Offset_mode: The parameter to define the mode of identify
        #    Offset_mode = 1 : the file just contains elements number
        #                      (so every number in line should be added an offset)
        #    Offset_mode = 2 : the file contains elements number and nodes of element
        #                      (so only the first number in line should be added an offset)
        # --------------------------------------------------------------------------
        #              Judge the kind of element set
        # --------------------------------------------------------------------------
        # If there is string 'Element' in the name of file, indicating that
        # the file contains elements define information, so node number
        # in line don't need to be offseted.
        Offset_mode = 1
        if 'element' in InpFile_prefix.lower():
            Offset_mode = 2
        logging.info("{0:<40}:{1:^15}".format('Input file prefix',InpFile_prefix))
        logging.info("{0:<40}:{1:^15}".format('Offset mode',Offset_mode))
        logging.info("{0:<40}:{1:^15}".format('Element number offset magnitude',self.Offset_Magnitude))
        # ===========================================================================
        #
        #                      Element number offset
        #
        # ===========================================================================  
        logging.info("{0:-^80}".format("Output"))
        logging.info("{0:^40}{1:^15}".format("Output file name","Multiplier"))
        with open(self.input_file_path,'r') as Inp_Origin_File:
            Inp_Origin_lines = Inp_Origin_File.readlines()
        # for layer in layer_info:
        for key, value in self.layer_info.items():
            # Get the layer name and layer multiplier
            layer_name = key
            layer_multiplier = int(value)
            EleNumOffset     = self.Offset_Magnitude*layer_multiplier
            if Offset_mode == 1:
                layer_file_name  = '%s_%s_%s'%(InpFile_prefix,InpFile_suffix,layer_name)
            elif Offset_mode == 2:
                layer_file_name = '%s_%s'%(InpFile_prefix,layer_name)
            logging.info("{0:^40}{1:^15}".format(layer_file_name,layer_multiplier))
            # with open('%s.inp'%layer_file_name,'w') as Inp_Offset_File:
            if Offset_mode == 1:
                # Due to the last line of the file isn't always contain same number of 
                # objects as other lines, so the numpy is not suitable for manipulating the file.
                # -------------------------------------------------------
                # Considering that the file just contains elements number in element set,
                # so every elemnts in line should be offseted.
                with open("%s.inp"%layer_file_name,'w') as InpFile:
                    for Real_line in Inp_Origin_lines:
                        Real_line = Real_line.strip('\n')
                        offset_line = ''
                        # Spliting the single element label from the line
                        Eles = Real_line.split(',')
                        EleNum_list = []
                        for EleNum_temp in Eles:
                            if EleNum_temp == '':
                                continue
                            # Adding the offset on every element labels
                            EleNum = int(EleNum_temp) + EleNumOffset
                            EleNum_list.append(str(EleNum))
                        offset_line = ','.join(EleNum_list)
                        InpFile.write(offset_line + '\n')
            elif Offset_mode == 2:
                # Operating offset through numpy
                output_file = '%s.inp'%layer_file_name
                Element_Connectivity = np.loadtxt(self.input_file_path, delimiter=',', dtype=int)
                # Adding the offset on the element label
                Element_Connectivity[:,0] += EleNumOffset
                np.savetxt(output_file, Element_Connectivity, fmt='%d', delimiter=',')
                # Operating offset through string manipulation which seems to be faster, but I don't know why.
                # with open("%s.inp"%layer_file_name,'w') as InpFile:
                #     for Real_line in Inp_Origin_lines:
                #         Real_line = Real_line.strip('\n')
                #         offset_line = ''
                #         # Spliting the single element label from the line
                #         EleNum_temp = Real_line.split(',')
                #         # Solely adding the offset on the element label
                #         EleNum = int(EleNum_temp[0]) + EleNumOffset
                #         # Merge the element label and node label
                #         offset_line = str(EleNum) + ',' + ','.join(EleNum_temp[1:])
                #         InpFile.write(offset_line + '\n')    