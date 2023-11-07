import os
import argparse

# *

def EleNum_Offset(Offset_origin_file_folder, Offset_Magnitude, EleNum_Offset_folder, InpFile):
    os.chdir(EleNum_Offset_folder)
    # ***************************************************************************
    # Extracting the inp file name
    InpFile_Name     = InpFile.split('.inp')[0]
    print('Input file name   : %s.inp'%InpFile_Name)
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
    print('Input file prefix : %s\tOffset mode : %d'%(InpFile_prefix,Offset_mode))
    print('Element number offset magnitude : %d'%Offset_Magnitude)
    # ***************************************************************************
    # Reading the layer information from the layer_info.txt file
    layer_info = []
    layer_info_file = open(os.path.join(Offset_origin_file_folder,'layer_info.txt'),'r')
    layer_info_file_lines = layer_info_file.readlines()
    for line in layer_info_file_lines:
        # Extracting the layer name and layer multiplier
        # E.g. line = 'Mech,1' --> layer_name = 'Mech', layer_multiplier = 1
        line_split = line.split(',')
        layer_name = line_split[0]
        layer_multiplier = int(line_split[1].split('\n')[0])
        layer_info.append([layer_name,layer_multiplier])
    # ===========================================================================
    #
    #                      Element number offset
    #
    # ===========================================================================  
    print("{0:-^70}".format("Output"))
    print("{0:^40}{1:^15}".format("Output file name","Multiplier"))
    Offset_origin_file_path = os.path.join(Offset_origin_file_folder,InpFile)
    with open(Offset_origin_file_path,'r') as Inp_Origin_File:
        Inp_Origin_lines = Inp_Origin_File.readlines()
    for layer in layer_info:
        # Get the layer name and layer multiplier
        layer_name       = layer[0]
        layer_multiplier = layer[1]
        EleNumOffset     = Offset_Magnitude*layer_multiplier
        if Offset_mode == 1:
            layer_file_name  = '%s_%s_%s'%(InpFile_prefix,InpFile_suffix,layer_name)
        elif Offset_mode == 2:
            layer_file_name = '%s_%s'%(InpFile_prefix,layer_name)
        print("{0:^40}{1:^15}".format(layer_file_name,layer_multiplier))
        # with open('%s.inp'%layer_file_name,'w') as Inp_Offset_File:
        if Offset_mode == 1:
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
            with open("%s.inp"%layer_file_name,'w') as InpFile:
                for Real_line in Inp_Origin_lines:
                    Real_line = Real_line.strip('\n')
                    offset_line = ''
                    # Spliting the single element label from the line
                    EleNum_temp = Real_line.split(',')
                    # Solely adding the offset on the element label
                    EleNum = int(EleNum_temp[0]) + EleNumOffset
                    # Merge the element label and node label
                    offset_line = str(EleNum) + ',' + ','.join(EleNum_temp[1:])
                    InpFile.write(offset_line + '\n')    
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The script of adding offset to the element number in the element set file.')
    parser.add_argument('-w', '--work_directory', type=str, metavar='', required=True, help='working directory')
    parser.add_argument('-m', '--Magnitude', type=float, metavar='', required=True, help='Magnitude of the offset')
    args = parser.parse_args()
    # *******************************************************
    print("{0:=^80}".format("EleNum_Offset Process"))
    # *******************************************************
    #               Extract the inp file name
    # *******************************************************
    # Get the working directory: x:xx/xx/xx.inp
    Eleset_origin_file_path = args.work_directory
    # Get the magnitude of the offset, default is 10000
    Offset_Magnitude = int(args.Magnitude)
    # Remove the space before the string
    # E.g. Directory = ' D:/Desktop' --> 'D:/Desktop' 
    Eleset_origin_file_path = Eleset_origin_file_path.strip()
    print('Input file path   : %s'%Eleset_origin_file_path)
    # ***************************************************************************
    # Split the file path into directory and file name
    # E.g. Offset_origin_file_path = 'D:/Desktop/xx.inp'
    #      Directory = 'D:/Desktop', InpFile = 'xx.inp'
    [Offset_origin_file_folder,InpFile] = os.path.split(Eleset_origin_file_path) 
    # Create the folder to store the offseted inp file
    EleNum_Offset_folder = os.path.join(Offset_origin_file_folder,'EleNum_Offset')
    if os.path.exists(EleNum_Offset_folder) == False:
        os.mkdir(EleNum_Offset_folder)
    # ***************************************************************************
    #               Adding the offset to the element number
    # ***************************************************************************
    # Offset_origin_file_folder : the folder of the origin inp file
    # Offset_Magnitude          : the magnitude of the offset
    # EleNum_Offset_folder      : the folder to store the offseted inp file
    # InpFile                   : the name of the origin inp file
    EleNum_Offset(Offset_origin_file_folder, Offset_Magnitude, EleNum_Offset_folder, InpFile)
    print("{0:=^80}".format("Finished") + '\n')
    