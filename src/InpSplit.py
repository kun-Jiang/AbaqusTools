import os
import math
import logging
from utilities.logConfig import logconfig

class Inpsplit:
    
    def __init__(self,input_file_path:str):
        # *******************************************************
        #               Extract the inp file name
        # *******************************************************
        self.input_file_path = input_file_path
        # Remove the space before the string
        # E.g. Directory = ' D:/Desktop' --> 'D:/Desktop' 
        self.input_file_path = self.input_file_path.strip()
        # Split the file path into directory and file name
        # E.g. input_file_path = 'D:/Desktop/xx.inp'
        #      Directory = 'D:/Desktop', InpFile = 'xx.inp'
        [self.Directory,self.InpFile] = os.path.split(self.input_file_path)
        # Create the new .inp file
        self.FileName = os.path.splitext(self.InpFile)[0]
        self.Split_inp_file = os.path.join(self.Directory,self.FileName+'_split.inp')
        # Create the folder to store the split inp file
        Inp_Split_folder = os.path.join(self.Directory,'Model')
        os.makedirs(Inp_Split_folder, exist_ok=True)
        os.chdir(Inp_Split_folder)
        # *******************************************************
        #               Set the log file
        # *******************************************************
        logconfig('InpSplit.log',self.Directory)
        # *******************************************************
        #               Start split the inp file
        # *******************************************************
        logging.info("{0:=^80}".format("Start inp_Split Process"))
        logging.info("The input file is: %s"%self.input_file_path)
        logging.info("The output file is: %s"%self.Split_inp_file)
        self.InpSplit()
        # *******************************************************
        #               Finished
        # *******************************************************
        logging.info("{0:=^80}".format("Finished") + '\n')

    def InpSplit(self):
        with open(self.input_file_path,'r') as InpFile:
            InpFile_origin_lines = InpFile.readlines()
        with open(self.Split_inp_file,'w') as InpFile_split:
            # **************************************************************************************************
            #                             loop over the original inp file
            # **************************************************************************************************
            InpFile_lines_iter = iter(InpFile_origin_lines)
            self.Node_count = 0
            self.Element_count = 0
            line = next(InpFile_lines_iter)
            while True:
                # line = next(InpFile_lines_iter, None)
                if line == None:
                    # If the end of the file is reached, break the loop
                    break
                line,judge = self.KeywordJudge(line,InpFile_lines_iter,InpFile_split)
                if line == None:
                    break
                # If the line doesn't contain the specific keyword, write the line to the new inp file
                # and iterate to the next line.
                if judge == False:
                    InpFile_split.write(line)
                    line = next(InpFile_lines_iter, None)

    def KeywordJudge(self,line_start:str,InpFile_lines_iter:iter,InpFile_split:str)->str:
        judge = False
        line_temp = line_start.strip()
        keyword = line_temp.split(',')[0].lower()
        if line_start != None and keyword == '*node':
            # If the line contains '*node' and not contains 'output', 
            # create a new node file and write the node information to the new file
            self.Node_count += 1
            # Write the line which contains the keyword '*node' to the new inp file
            InpFile_split.write(line_start)
            line = self.NodeAccess(self.Node_count,InpFile_lines_iter)
            InpFile_split.write('*include, input=Model\\Node_%s.inp\n'%self.Node_count)
            judge = True
            
        elif line_start != None and keyword == '*element':
            # If the line contains '*element' and not contains 'output',
            # create a new element file and write the element information to the new file
            self.Element_count += 1
            # Write the line which contains the keyword '*element' to the new inp file
            InpFile_split.write(line_start)
            line = self.ElementAccess(self.Element_count,InpFile_lines_iter)
            InpFile_split.write('*include, input=Model\\Element_%s.inp\n'%self.Element_count)
            judge = True
        elif line_start != None and keyword == '*elset':
            # If the line contains 'elset', write the line to the new inp file
            # InpFile_split.write(line)
            line, Elset_Name = self.ElsetAccess(InpFile_lines_iter, line_start, InpFile_split)
            if Elset_Name != None:
                InpFile_split.write('*include, input=Model\\Ele_%s.inp\n'%Elset_Name)
            judge = True
        elif line_start != None and keyword == '*nset':
            # If the line contains 'nset', write the line to the new inp file
            # InpFile_split.write(line)
            line, Nset_Name = self.NsetAccess(InpFile_lines_iter, line_start, InpFile_split)
            if Nset_Name != None:
                InpFile_split.write('*include, input=Model\\Nset_%s.inp\n'%Nset_Name)
            judge = True
        # If the input line does not contain the specified keyword (judge = False), 
        # just return the input line and don't iterate the next line, 
        # else return the iterated line (judge = True)
        if judge == False:
            line = line_start
        return line,judge
        
    def NodeAccess(self,Node_count:int,InpFile_lines_iter:iter)->str:
        """Extract the node information from the original inp file and write it to the new node file

        Args:
            Node_count (int): the number of the node file
            InpFile_lines_iter (iter): an iterator of the original inp file

        Returns:
            line: the current line of the iterator
        """        
        with open("Node_%s.inp"%Node_count,'w') as NodeFile:
            while True:
                # Iterate each line of the original inp file to extract the 
                # node information until the next '*' in the line
                line = next(InpFile_lines_iter, None)
                if line == None or '*' in line:
                    # line == None: The iterator reaches the end of the file
                    # '*' in line : The line contains '*', return the line and break the loop
                    break
                # Write the node information to the new node file
                NodeFile.write(line)
            return line
    def ElementAccess(self,Element_count:int,InpFile_lines_iter:iter)->str:
        """Extract the element information from the original inp file and write it to the new element file

        Args:
            Element_count (int): the number of the element file
            InpFile_lines_iter (iter): an iterator of the original inp file

        Returns:
            line: the current line of the iterator
        """        
        with open("Element_%s.inp"%Element_count,'w') as ElementFile:
            while True:
                # Iterate each line of the original inp file to extract the 
                # element information until the next '*' in the line
                line = next(InpFile_lines_iter, None)
                if line == None or '*' in line:
                    # line == None: The iterator reaches the end of the file
                    # '*' in line : The line contains '*', return the line and break the loop
                    # return line
                    break
                # Write the element information to the new element file
                ElementFile.write(line)
            return line
        
    def ElsetAccess(self,InpFile_lines_iter:iter, line_start:str,InpFile_split:str)->str:
        """Extract the elset information from the original inp file and write it to the new elset file

        Args:
            InpFile_lines_iter (iter): an iterator of the original inp file
            line_start (str): the start line when the function is called

        Returns:
            line (str): the current line of the iterator
            Elset_Name (str): _description_
        """        
        # Remove the '\n' at the end of the string
        line_temp = line_start.strip()
        line_split = line_temp.split(',')
        try:
            if '*elset' in line_split[0].lower():
                Elset_Name = line_split[1].split('=')[1]
                if 'generate' in line_start.lower():
                    # Remove the 'generate' keyword from the line and write the line
                    # which contains the keyword 'elset' to the new inp file
                    line_temp = ','.join(line_split[:-1]) + '\n'
                    InpFile_split.write(line_temp)
                    # Extract the start, end, and step information from the line
                    line = next(InpFile_lines_iter, None)
                    line_split = line.split(',')
                    start = int(line_split[0])
                    end = int(line_split[1])
                    step = int(line_split[2])
                    ele_sequence = self.Sequence_Generate(start, end, step)
                    with open('Ele_' + Elset_Name + '.inp','w') as ElsetFile:
                        self.DataLineWrite(ele_sequence, ElsetFile)
                    line = next(InpFile_lines_iter, None)
                    return line, Elset_Name
            else:
                line = line_start
                return line, None
        except:
            logging.error('There is an error when extracting the elset name\n' +
                          'Error in line: ' + line_start,
                          exc_info=True)
        InpFile_split.write(line_start)
        with open('Ele_' + Elset_Name + '.inp','w') as ElsetFile:
            while True:
                # Iterate each line of the original inp file to extract the
                # nset information until the next '*' in the line
                line = next(InpFile_lines_iter, None)
                if line == None:
                    # The iterator reaches the end of the file
                    break
                if '*' in line:
                    # If the line contains '*', return the line and Elset_Name and break the loop
                    return line, Elset_Name
                # Write the nset information to the new nset file
                ElsetFile.write(line)
    
    def NsetAccess(self,InpFile_lines_iter:iter, line_start:str,InpFile_split:str)->str:
        # Remove the '\n' at the end of the string
        line_temp = line_start.strip()
        line_split = line_temp.split(',')
        try:
            if '*nset' in line_split[0].lower():
                Nset_Name = line_split[1].split('=')[1]
                if 'generate' in line_start.lower():
                    # Remove the 'generate' keyword from the line and write the line
                    # which contains the keyword 'nset' to the new inp file
                    line_temp = ','.join(line_split[:-1]) + '\n'
                    InpFile_split.write(line_temp)
                    # Extract the start, end, and step information from the line
                    line = next(InpFile_lines_iter, None)
                    line_split = line.split(',')
                    start = int(line_split[0])
                    end = int(line_split[1])
                    step = int(line_split[2])
                    node_sequence = self.Sequence_Generate(start, end, step)
                    with open('Nset_' + Nset_Name + '.inp','w') as NsetFile:
                        self.DataLineWrite(node_sequence, NsetFile)
                    line = next(InpFile_lines_iter, None)
                    return line, Nset_Name
            else:
                line = line_start
                return line, None
        except:
            logging.error('There is an error when extracting the nset name\n' +
                          'Error in line: ' + line_start,
                          exc_info=True)
        InpFile_split.write(line_start)
        with open('Nset_' + Nset_Name + '.inp','w') as NsetFile:
            while True:
                # Iterate each line of the original inp file to extract the
                # nset information until the next '*' in the line
                line = next(InpFile_lines_iter, None)
                if line == None:
                    # The iterator reaches the end of the file
                    break
                if '*' in line:
                    # If the line contains '*', return the line and Nset_Name and break the loop
                    return line, Nset_Name
                # Write the nset information to the new nset file
                NsetFile.write(line)
    
    def Sequence_Generate(self,start:int,end:int,inc:int)->list:
        """This function is used to generate a sequence of numbers

        Args:
            start (int): the start number of the sequence
            end (int): the end number of the sequence
            inc (int): the increment of the sequence

        Returns:
            Sequence (list): A list contains the sequence of numbers
        """    
        Sequence = []
        for i in range(start,end+1,inc):
            Sequence.append(str(i))
        return Sequence
    
    def DataLineWrite(self,vector:list,file,max_num:int = 16):
        # Writing the Element/Node number into the ***.inp file
        # 
        """Writing the Element/Node number into the ***.inp file. Following the 
        convention of ABAQUS, there are no more than 16 element numbers in a line
            
        Args:
            vector (list): A vector contains the element/node number
            file (file): The file object of the ***.inp file
            max_num (int) (optional): The maximum number of element/node numbers 
                            allowed in a single line. Defaults to 16 if not provided.
        """    
        # Convert elements of vector from int to str
        vector = list(map(str,vector))
        Num = len(vector)
        lines_num = math.ceil(Num / max_num)
        if lines_num == 1:
            line = ','.join(vector)
            file.write(line + '\n')
        else:
            for j in range(lines_num):
                line = ', '.join(vector[(j*16):(j+1)*max_num])
                file.write(line + '\n')

                
    # def GrainGB_Extract(Inp_Split_folder:str,Ele_num_Grain_boundaries:list,Node_num_Grain_boundaries:list):
    #     file_list = os.listdir(Inp_Split_folder)
    #     # Find out all input file that contains elements' number of grain
    #     grain_file_list = []
    #     for file_name in file_list:
    #         file_name_temp = file_name.lower()
    #         if 'ele' in file_name_temp and 'grain' in file_name_temp and 'gb' not in file_name_temp:
    #             grain_file_list.append(file_name)

    #     for file in grain_file_list:
    #         file_name = file.split('.')[0]
    #         file_path = os.path.join(Inp_Split_folder,file)
    #         Ele_num_grain = []
    #         with open(file_path,'r') as grain_file:
    #             lines = grain_file.readlines()
    #             for line in lines:
    #                 line = line.strip()
    #                 line_split = line.split(',')
    #                 for ele_num in line_split:
    #                     Ele_num_grain.append(int(ele_num))
    #         # Compare the element number of grain with the element number of grain boundaries, and find
    #         # out the intersection of them.
    #         # Ele_num_grain_GB: the element number of grain boundary part in a grain
    #         Ele_num_grain_GB = list(set(Ele_num_grain) & set(Ele_num_Grain_boundaries))
    #         # Writing the list into file
    #         Grain_GB_file_path = os.path.join(Inp_Split_folder,file_name+'_GB.inp')
    #         with open(Grain_GB_file_path,'w') as Grain_GB_file:
    #             DataLineWrite(Ele_num_grain_GB,Grain_GB_file)
        

                
        # # *******************************************************
        # # Create the new .inp file contains the element number of grain interior
        # # *******************************************************
        # if Ele_num_Grain_boundaries != None:
        #     # Remove the repeated element number
        #     Ele_num_Total = set(Ele_num_Total)
        #     Ele_num_Grain_boundaries = set(Ele_num_Grain_boundaries)
        #     Ele_num_Grain_interior = list(Ele_num_Total - Ele_num_Grain_boundaries)
        #     elset_Name = 'gi'
        #     elsetFile_path = "Ele_%s.inp"%elset_Name
        #     elsetFile = open(elsetFile_path,'w')
        #     DataLineWrite(Ele_num_Grain_interior,elsetFile)
        # # *******************************************************
        # # Create the new .inp file contains the node number of grain interior
        # # *******************************************************
        # if Node_num_Grain_boundaries != None:
        #     # Remove the repeated node number
        #     Node_num_Total = set(Node_num_Total)
        #     Node_num_Grain_boundaries = set(Node_num_Grain_boundaries)
        #     Node_num_Grain_interior = Node_num_Total - Node_num_Grain_boundaries
        #     nset_Name = 'gi'
        #     nsetFile_path = "Nset_%s.inp"%nset_Name
        #     nsetFile = open(nsetFile_path,'w')
        #     DataLineWrite(Node_num_Grain_interior,nsetFile)
        # # *******************************************************
        # GrainGB_Extract(Inp_Split_folder,Ele_num_Grain_boundaries,Node_num_Grain_boundaries)
        # # *******************************************************
        # InpFile_New.close()
        # # *******************************************************
        # # Move the new.inp file to the working directory
        # os.chdir(Directory)
        # InpFile_split_old = os.path.join(Directory,InpFile_New_name)
        # if os.path.exists(InpFile_split_old):
        #     os.remove(InpFile_split_old)
        # shutil.move(InpFile_New_path,Directory)
        # print('-'*80)
        # return Elset_name_list