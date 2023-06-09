import os
import shutil
import argparse


def Inp_Split(Inp_Split_folder,Inp_origin_file_path,Directory):
    os.chdir(Inp_Split_folder)
    InpFile_Name     = InpFile.split('.inp')[0]
    print('Input file name : %s.inp'%InpFile_Name)
    # *******************************************************
    InpFile_origin      = open(Inp_origin_file_path,'r')
    InpFile_origin_lines    = InpFile_origin.readlines()
    InpFile_origin.close()
    # Rewrite the Inp file
    InpFile_New_name = "%s-Split.inp"%InpFile_Name
    InpFile_New_path = os.path.join(Inp_Split_folder,InpFile_New_name)
    InpFile_New = open(InpFile_New_path,'w')
    print('Output file name: %s-Split.inp'%InpFile_Name)
    Node_num_Total = []
    Ele_num_Total = []
    Ele_num_Grain_boundaries = []
    Node_num_Grain_boundaries = []
    # *******************************************************
    print("{:-^80}".format('Running log'))
    iter_num = 0
    Inp_set_count = 0
    while len(InpFile_origin_lines) > 0:
        iter_num += 1
        line = InpFile_origin_lines[0]
        InpFile_origin_lines.remove(line)
        InpFile_New.write(line)
        line = line.lower()
        if '**' in line:
            # If there is a '**' in line, indicating that this line is a comment
            continue
        # *******************************************************
        if '*node' in line:
            if 'output' in line:
                continue
            # Writing the Node.inp file
            Inp_set_count += 1
            NodeFile    = open("Node_%s.inp"%Inp_set_count,'w')
            temp_lines = list(InpFile_origin_lines)
            for subline in temp_lines:
                if '*' in subline:
                    # If there is a '*' in subline, this section is finished
                    InpFile_New.write("*include,input=Model\\Node_%s.inp\n"%Inp_set_count)
                    break
                else:
                    # Writing the node and coordinates
                    NodeFile.write(subline)
                    InpFile_origin_lines.remove(subline)
                    # Reading all node number
                    # The list Ele_num_Total will be used in the later
                    Node_num = int(subline.split(',')[0])
                    Node_num_Total.append(Node_num)       
                            
            NodeFile.close()
            continue
        # *******************************************************
        if '*element' in line:
            if 'output' in line:
                continue
            # There must be a elset after the elements define (For example: *Element, type=***, elset=***)
            # Writing the ***.inp file which define the element(*** is the name of elset)
            if 'elset' in line:
                elset_Name_temp  = line.split('elset=')[1]
                elset_Name       = elset_Name_temp.split('\n')[0]
            else:
                # If there not defines a elset, using the name of element type temporarily
                elset_Name_temp  = line.split('type=')[1]
                elset_Name       = elset_Name_temp.split('\n')[0]
                print('There is not a name for the element set, using the name of element type temporarily\n%s'%elset_Name)
            elementFile_path = "Element_%s.inp"%elset_Name
            elementFile      = open(elementFile_path,'w')
            temp_lines       = list(InpFile_origin_lines)
            for subline in temp_lines:
                if '*' in subline:
                    # If there is a '*' in subline, this section is finished
                    InpFile_New.write("*include,input=Model\\Element_%s.inp\n"%elset_Name)
                    break
                else:
                    # Writing the element number and nodes
                    elementFile.write(subline)
                    InpFile_origin_lines.remove(subline)
                    # Reading all element number
                    # The list Ele_num_Total will be used in the later
                    Ele_num = int(subline.split(',')[0])
                    Ele_num_Total.append(Ele_num)
                    
            elementFile.close()
        # *******************************************************
        if '*elset' in line:
            # Writing the ***.inp file (*** is the name of elset)
            elset_Name_temp = line.split('elset=')[1]
            if 'generate' in elset_Name_temp:
                # *Elset, elset=***, generate
                # This situation is ignored
                continue
            elif ',' in elset_Name_temp:
                # *Elset, elset=***, instance=***\n
                elset_Name = elset_Name_temp.split(',')[0]
            elif '"' in elset_Name_temp:
                # *Nset, nset="***"\n
                elset_Name = elset_Name_temp.split('"')[1]
            else:
                # *Elset, elset=***\n
                elset_Name = elset_Name_temp.split('\n')[0]
            elsetFile_path = "Ele_%s.inp"%elset_Name
            elsetFile = open(elsetFile_path,'w')
            temp_lines = list(InpFile_origin_lines)
            for subline in temp_lines:
                if '*' in subline:
                    # If there is a '*' in subline, this section is finished
                    InpFile_New.write("*include,input=Model\\Ele_%s.inp\n"%elset_Name)
                    if 'grain_boundaries' in elset_Name:
                        # If there is a element set named grain_boundaries, creat a set named grain_interior
                        InpFile_New.write("*Elset, elset=Ele_grain_interior\n*include,input=Model\\Ele_grain_interior.inp\n")
                    break
                else:
                    # Writing the element number into the new .inp file
                    elsetFile.write(subline)
                    InpFile_origin_lines.remove(subline)
                    # Reading the element number corresponding to the grain boundaries
                    # The list Ele_num_Grain_boundaries will be used in the later
                    if 'grain_boundaries' in elset_Name:
                        subline_split = subline.split(',')
                        for Ele_num in subline_split:
                            # Remove the '\n' in the end of the element number
                            if '\n' in Ele_num:
                                Ele_num = Ele_num.split('\n')[0]
                            Ele_num_Grain_boundaries.append(int(Ele_num))
            elsetFile.close()
        # *******************************************************
        if '*nset' in line:
            # Writing the ***.inp file (*** is the name of nset)
            nset_Name_temp = line.split('nset=')[1]
            if 'generate' in nset_Name_temp:
                # *Nset, nset=***, generate
                # This situation is ignored
                continue
            elif ',' in nset_Name_temp:
                # *Nset, nset=***, instance=***\n
                nset_Name = nset_Name_temp.split(',')[0]
            elif '"' in nset_Name_temp:
                # *Nset, nset="***"\n
                nset_Name = nset_Name_temp.split('"')[1]
            else:
                # *Nset, nset=***\n
                nset_Name = nset_Name_temp.split('\n')[0]
            nsetFile_path = "Nset_%s.inp"%nset_Name
            nsetFile = open(nsetFile_path,'w')
            temp_lines = list(InpFile_origin_lines)
            for subline in temp_lines:
                if '*' in subline:
                    # If there is a '*' in subline, this section is finished
                    InpFile_New.write("*include,input=Model\\Nset_%s.inp\n"%nset_Name)
                    if 'grain_boundaries' in nset_Name:
                        # If there is a element set named Grain_boundaries, creat a set named Grain_interior
                        InpFile_New.write("*Nset, nset=Nset_grain_interior\n*include,input=Model\\Nset_grain_interior.inp\n")
                    break
                else:
                    # Writing the node number
                    nsetFile.write(subline)
                    InpFile_origin_lines.remove(subline)
                    # Reading the node number corresponding to the grain boundaries
                    # The list Node_num_Grain_boundaries will be used in the later
                    if 'grain_boundaries' in nset_Name:
                        subline_split = subline.split(',')
                        for Node_num in subline_split:
                            # Remove the '\n' in the end of the element number
                            if '\n' in Node_num:
                                Node_num = Node_num.split('\n')[0]
                            Node_num_Grain_boundaries.append(int(Node_num))
            nsetFile.close()
            
    # *******************************************************
    # Create the new .inp file contains the element number of grain interior
    # *******************************************************
    # Remove the repeated element number
    Ele_num_Total = set(Ele_num_Total)
    Ele_num_Grain_boundaries = set(Ele_num_Grain_boundaries)
    Ele_num_Grain_interior = Ele_num_Total - Ele_num_Grain_boundaries
    # print(Ele_num_Grain_boundaries)
    # os._exit(0)
    elset_Name = 'grain_interior'
    elsetFile_path = "Ele_%s.inp"%elset_Name
    elsetFile = open(elsetFile_path,'w')
    for Ele_num in Ele_num_Grain_interior:
        elsetFile.write("%d\n"%Ele_num)
    elsetFile.close()
    # *******************************************************
    # Create the new .inp file contains the node number of grain interior
    # *******************************************************
    # Remove the repeated node number
    Node_num_Total = set(Node_num_Total)
    Node_num_Grain_boundaries = set(Node_num_Grain_boundaries)
    Node_num_Grain_interior = Node_num_Total - Node_num_Grain_boundaries
    # print(Node_num_Grain_boundaries)
    # os._exit(0)
    nset_Name = 'grain_interior'
    nsetFile_path = "Nset_%s.inp"%nset_Name
    nsetFile = open(nsetFile_path,'w')
    for Node_num in Node_num_Grain_interior:
        nsetFile.write("%d\n"%Node_num)
    nsetFile.close()
    # *******************************************************
    InpFile_New.close()
    # *******************************************************
    # Move the new.inp file to the working directory
    os.chdir(Directory)
    InpFile_split_old = os.path.join(Directory,InpFile_New_name)
    if os.path.exists(InpFile_split_old):
        os.remove(InpFile_split_old)
    shutil.move(InpFile_New_path,Directory)
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The script of spliting the ABAQUS inp file')
    parser.add_argument('-w', '--work_directory', type=str, metavar='', required=True, help='working directory')
    args = parser.parse_args()
    # *******************************************************
    print("{0:=^80}".format("Inp_Split Process"))
    # *******************************************************
    #               Extract the inp file name
    # *******************************************************
    # Get the working directory: x:xx/xx/xx.inp
    Inp_origin_file_path = args.work_directory
    # Remove the space before the string
    # E.g. Directory = ' D:/Desktop' --> 'D:/Desktop' 
    Inp_origin_file_path = Inp_origin_file_path.strip()
    # Split the file path into directory and file name
    # E.g. Inp_origin_file_path = 'D:/Desktop/xx.inp'
    #      Directory = 'D:/Desktop', InpFile = 'xx.inp'
    [Directory,InpFile] = os.path.split(Inp_origin_file_path)
    # Create the folder to store the split inp file
    Inp_Split_folder = os.path.join(Directory,'Model')
    if os.path.exists(Inp_Split_folder) == False:
        os.mkdir(Inp_Split_folder)
    # *******************************************************
    #               Split the inp file
    # *******************************************************
    Inp_Split(Inp_Split_folder,Inp_origin_file_path,Directory)
    print("{0:=^80}".format("Finished") + '\n')
    