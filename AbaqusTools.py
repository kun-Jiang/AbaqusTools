import threading
import tkinter as tk
from tkinter import ttk
import subprocess
from tkinter import filedialog
import tkinter.messagebox
import os

class AbaqusGUI(tk.Tk):
    def __init__(self):
        # 设置窗口标题和窗口大小
        super().__init__()
        self.title("ABAQUS GUI")
        self.geometry("520x400")
        self.create_widgets()
        
    def create_widgets(self):
        # 创建界面切换按钮
        self.style = ttk.Style()
        self.style.configure('TNotebook.Tab', padding=(30, 8))
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True)

        # 创建CMD选项卡
        self.tab_cmd = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_cmd, text='CMD')
        create_cmd_gui(self.tab_cmd)

        # 创建InpSplit选项卡
        self.tab_InpSplit = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_InpSplit,text='InpSplit')
        create_InpSplit_gui(self.tab_InpSplit)        
        
        # 创建Element_Num_Offset选项卡
        self.tab_Ele_num_offset = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_Ele_num_offset, text='Offset')
        create_Ele_num_offset_gui(self.tab_Ele_num_offset)
        
        # Create a tabwidget for extracting data from odb file
        self.tab_odb_extract = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_odb_extract, text='ODB Extract')
        create_odb_extract_gui(self.tab_odb_extract)
        
class create_cmd_gui():
    def __init__(self, parent):
        self.tab_cmd = parent
        # Get all the files in the folder
        self.CMD_working_directory_label = ttk.Label(self.tab_cmd, text="Working directory:")
        self.CMD_working_directory_entry = ttk.Entry(self.tab_cmd, width=40)
        self.button_browse_working_directory = ttk.Button(self.tab_cmd, text="Browse", 
                                            command=self.browse_working_directory)
        
        # Create a combobox for displaying inp files
        self.CMD_inp_file_lable = ttk.Label(self.tab_cmd, text="Select inp file:")
        self.CMD_inp_combobox = ttk.Combobox(self.tab_cmd, width=37)

        # Create a combobox for displaying for files
        self.CMD_for_file_label = ttk.Label(self.tab_cmd, text="Select for file:")
        self.CMD_for_combobox = ttk.Combobox(self.tab_cmd, width=37)
        self.button_update_combobox = ttk.Button(self.tab_cmd, text="Update",
                                            command=self.update_combobox)
         
        # Create a button to run abaqus, and run it in a new thread by the function Run_in_new_thread
        self.CMD_button_run_abaqus = ttk.Button(self.tab_cmd, text="Call ABAQUS", 
                                            command=lambda: Run_in_new_thread(self.run_abaqus))
        # Create a spinbox, default value is 2, indicate 2 cores
        self.cpu_cores_label = ttk.Label(self.tab_cmd, text="CPU cores:")
        self.cpu_cores_value = tk.IntVar(value=2)
        self.cpu_cores_spinbox = tk.Spinbox(self.tab_cmd, from_=1, to=1000, width=5, textvariable=self.cpu_cores_value)
        # Create a check button, default value is 1
        self.check_Int_value = tk.IntVar(value=1)
        self.check_Int_button = ttk.Checkbutton(self.tab_cmd, text="Display running status in CMD window",variable=self.check_Int_value)
        # Create a check button, default value is 1
        self.check_ask_value = tk.IntVar(value=1)
        self.check_ask_button = ttk.Checkbutton(self.tab_cmd, text="Overwrite old files",variable=self.check_ask_value)
        # *********************************************************************** #
        # Layout the GUI
        # *********************************************************************** #
        # Get working directory
        self.CMD_working_directory_label.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.CMD_working_directory_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W)
        self.button_browse_working_directory.grid(row=0, column=2, padx=5, pady=5)
        # Get inp file
        self.CMD_inp_file_lable.grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.CMD_inp_combobox.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W)
        # Get for file
        self.CMD_for_file_label.grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.CMD_for_combobox.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W)
        self.button_update_combobox.grid(row=2, column=2, padx=5, pady=5)
        # Set CPU cores
        self.cpu_cores_label.grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
        self.cpu_cores_spinbox.grid(row=3, column=1, padx=5, pady=5, sticky=tk.W)
        # Set check button
        self.check_Int_button.grid(row=4, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W )
        self.check_ask_button.grid(row=5, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W )
        # Run button
        self.CMD_button_run_abaqus.grid(row=6, column=1, padx=5, pady=5)
        
    def browse_working_directory(self):
        # 打开文件对话框
        file_path = filedialog.askdirectory()
        
        # 将选择的文件路径显示在输入框中
        self.CMD_working_directory_entry.delete(0, tk.END)
        self.CMD_working_directory_entry.insert(0, file_path)
        
    def update_combobox(self):
        All_files = os.listdir(self.CMD_working_directory_entry.get())

        # Filter out only the .inp files
        CMD_inp_files = [f for f in All_files if f.endswith(".inp")]        
        # Update the combobox with the inp files
        self.CMD_inp_combobox["values"] = CMD_inp_files   
        # Filter out only the .for files
        CMD_for_files = [f for f in All_files if f.endswith(".for")]
        # Update the combobox with the for files
        self.CMD_for_combobox["values"] = CMD_for_files
        
    def run_abaqus(self):
        # Get working directory from input box
        CMD_working_directory = self.CMD_working_directory_entry.get()
        # Get inp file name from combobox
        CMD_inp_file_name = self.CMD_inp_combobox.get()
        # Get for file name from combobox
        CMD_for_file_name = self.CMD_for_combobox.get()
        scratch_path    = os.path.join(CMD_working_directory, 'Temp files')
        cpu_cores       = self.cpu_cores_value.get()
        check_Int       = self.check_Int_value.get()
        check_ask       = self.check_ask_value.get()

        os.chdir(CMD_working_directory)
        abaqus_bat_path = 'D:/ABAQUS/Commands/abaqus.bat'
        abaqus_job      = 'job='+CMD_inp_file_name
        abaqus_user     = 'user='+CMD_for_file_name
        abaqus_scratch_path = 'scratch='+scratch_path
        abaqus_cpu_cores = 'cpus='+str(cpu_cores)
        if os.path.exists(scratch_path) == False:
            try:
                os.makedirs(scratch_path)
            except:
                tk.messagebox.showerror("Error", scratch_path + '\n'
                    "You do not have the permission to create the scratch folder! \
                    Please check the path and create the folder manually. If these don\'t work,\
                    modifying the AbaqusTools.py file to set a appropriate scratch path.")
        # Construct the command line
        abaqus_command_line = [abaqus_bat_path, abaqus_job, abaqus_user, 
                               abaqus_cpu_cores, abaqus_scratch_path]
        if check_Int == 1:
            # Display the progress of the calculation
            abaqus_Int = 'int'
            abaqus_command_line.append(abaqus_Int)
        if check_ask == 1:
            # Overwrite the existing file
            abaqus_ask = 'ask=off'
            abaqus_command_line.append(abaqus_ask)
        # print(abaqus_command_line)
        try:
            subprocess.call(abaqus_command_line)
            # subprocess.call(['D:/ABAQUS/Commands/abaqus.bat', 'job='+input_file, 'user='+for_file,
            #                  'scratch="D:/ABAQUS/temp/CMD_Temp"', 'int', 'ask=off'])
        except FileNotFoundError:
            tk.messagebox.showerror("Error", "Execute command line failed. Please check the path of the abaqus.bat file and modify it in AbaqusTools.py")
            print(abaqus_command_line)
class create_InpSplit_gui():
    def __init__(self, parent):
        self.tab_InpSplit = parent
        # 创建ttk风格的界面元素
        self.label_InpSplit_file = ttk.Label(self.tab_InpSplit, text='Inp file path:')
        self.entry_InpSplit = ttk.Entry(self.tab_InpSplit,width=40)
        self.button_browse_InpSplit = ttk.Button(self.tab_InpSplit, text="Select", 
                                                command=self.browse_InpSplit_file)
        
        self.button_run_InpSplit = ttk.Button(self.tab_InpSplit,text='Execute split',
                                                command=lambda: Run_in_new_thread(self.run_InpSplit))
        # 界面元素布局
        self.label_InpSplit_file.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.entry_InpSplit.grid(row=0, column=1, padx=5, pady=5)
        self.button_browse_InpSplit.grid(row=0, column=2, padx=5, pady=5)

        self.button_run_InpSplit.grid(row=2, column=1, padx=5, pady=5)
        
    def run_InpSplit(self):
        # 从输入框中获取文件路径
        input_file_path = self.entry_InpSplit.get()
        InpSplit_py_file_path = os.path.join(working_directory,'src\InpSplit.py')
        # 执行InpSplit.py
        try:
            subprocess.call(['python', InpSplit_py_file_path, '-w '+input_file_path])
        except FileNotFoundError:
            tk.messagebox.showerror("错误", "找不到Inp文件或Python执行文件！")
                    
    def browse_InpSplit_file(self):
        # 打开文件对话框
        file_path = filedialog.askopenfilename(filetypes=[("Inp files", "*.inp"), ("All files", "*.*")])

        # 将选择的文件路径显示在输入框中
        self.entry_InpSplit.delete(0, tk.END)
        self.entry_InpSplit.insert(0, file_path)

class create_Ele_num_offset_gui():
    def __init__(self, parent):
        self.tab_Ele_num_offset = parent
        # *********************************************************************** #
        #                             Get input file                              #
        # *********************************************************************** #
        self.Offset_origin_file_label = ttk.Label(self.tab_Ele_num_offset, text='Inp file path:')
        self.Offset_origin_file_entry = ttk.Entry(self.tab_Ele_num_offset,width=40)
        self.Button_browse_offset = ttk.Button(self.tab_Ele_num_offset, text="Select", 
                                                 command=self.browse_Offset_file)
        # *********************************************************************** #
        #                       Creating the text input box                       #
        # *********************************************************************** #
        self.Offset_layer_name = ttk.Label(self.tab_Ele_num_offset, text="Layer")
        self.Offset_file_name = ttk.Label(self.tab_Ele_num_offset, text="Layer name")
        self.Offset_multiplier = ttk.Label(self.tab_Ele_num_offset, text="Multiplier")
        # Creating the first layer input
        self.First_layer_label = ttk.Label(self.tab_Ele_num_offset, text="First layer")
        self.First_layer_input_entry = ttk.Entry(self.tab_Ele_num_offset, width=40)
        self.First_layer_input_entry.insert(tk.END, "Mech")
        self.First_layer_multiplier_entry = ttk.Entry(self.tab_Ele_num_offset, width=10, justify='center')
        self.First_layer_multiplier_entry.insert(tk.END, "0")
        # Creating the second layer input
        self.Second_layer_label = ttk.Label(self.tab_Ele_num_offset, text="Second layer")
        self.Second_layer_input_entry = ttk.Entry(self.tab_Ele_num_offset, width=40)
        self.Second_layer_input_entry.insert(tk.END, "Chem")
        self.Second_layer_multiplier_entry = ttk.Entry(self.tab_Ele_num_offset, width=10, justify='center')
        self.Second_layer_multiplier_entry.insert(tk.END, "1")
        # Creating the third layer input
        self.Third_layer_label = ttk.Label(self.tab_Ele_num_offset, text="Third layer")
        self.Third_layer_input_entry = ttk.Entry(self.tab_Ele_num_offset, width=40)
        self.Third_layer_input_entry.insert(tk.END, "Phase")
        self.Third_layer_multiplier_entry = ttk.Entry(self.tab_Ele_num_offset, width=10, justify='center')
        self.Third_layer_multiplier_entry.insert(tk.END, "2")
        # Creating the fourth layer input
        self.Fourth_layer_label = ttk.Label(self.tab_Ele_num_offset, text="Fourth layer")
        self.Fourth_layer_input_entry = ttk.Entry(self.tab_Ele_num_offset, width=40)
        self.Fourth_layer_input_entry.insert(tk.END, "Dummy")
        self.Fourth_layer_multiplier_entry = ttk.Entry(self.tab_Ele_num_offset, width=10, justify='center')
        self.Fourth_layer_multiplier_entry.insert(tk.END, "3")
        # Creatinng the magnitude of offset
        self.Offest_Magnitude_label = ttk.Label(self.tab_Ele_num_offset, text="Offset magnitude")
        self.Offset_Magnitude_entry = ttk.Entry(self.tab_Ele_num_offset, width=10, justify='center')
        # The default value of offset magnitude is 10000
        self.Offset_Magnitude_entry.insert(tk.END, "10000")
        # Creating the execute button
        self.button_run_EleNum_Offset = ttk.Button(self.tab_Ele_num_offset,text='Run',
                                              command=lambda: Run_in_new_thread(self.run_Ele_num_offset))
        # Text of tips
        self.label_Ele_num_offset = ttk.Label(self.tab_Ele_num_offset, text='Tips:\n'
            +'The first layer elements use similar element number with the input file, so its\n'
            +'multiplier is equal to ZERO. And the offset is equal to the number of elements')
        # *********************************************************************** #
        #                        Layout the text input box                        #
        # *********************************************************************** #
        # Get input file
        self.Offset_origin_file_label.grid(row=0, column=0, padx=5, pady=5, sticky=tk.N)
        self.Offset_origin_file_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.N)
        self.Button_browse_offset.grid(row=0, column=2, padx=5, pady=5, sticky=tk.N)
        # Table head
        self.Offset_layer_name.grid(row=1, column=0, padx=5, pady=5, sticky=tk.N)
        self.Offset_file_name.grid(row=1, column=1, padx=5, pady=5, sticky=tk.N)
        self.Offset_multiplier.grid(row=1, column=2, padx=5, pady=5, sticky=tk.N)
        # First layer
        self.First_layer_label.grid(row=2, column=0, padx=5, pady=5, sticky=tk.N)
        self.First_layer_input_entry.grid(row=2, column=1, padx=5, pady=5)
        self.First_layer_multiplier_entry.grid(row=2, column=2, padx=5, pady=5)
        # Second layer
        self.Second_layer_label.grid(row=3, column=0, padx=5, pady=5, sticky=tk.N)
        self.Second_layer_input_entry.grid(row=3, column=1, padx=5, pady=5)
        self.Second_layer_multiplier_entry.grid(row=3, column=2, padx=5, pady=5)
        # Third layer
        self.Third_layer_label.grid(row=4, column=0, padx=5, pady=5, sticky=tk.N)
        self.Third_layer_input_entry.grid(row=4, column=1, padx=5, pady=5)
        self.Third_layer_multiplier_entry.grid(row=4, column=2, padx=5, pady=5)
        # Fourth layer
        self.Fourth_layer_label.grid(row=5, column=0, padx=5, pady=5, sticky=tk.N)
        self.Fourth_layer_input_entry.grid(row=5, column=1, padx=5, pady=5)   
        self.Fourth_layer_multiplier_entry.grid(row=5, column=2, padx=5, pady=5)
        # Magnitude of offset
        self.Offest_Magnitude_label.grid(row=6, column=0, padx=5, pady=5, sticky=tk.W)
        self.Offset_Magnitude_entry.grid(row=6, column=1, padx=5, pady=5, sticky=tk.W)
        # Execute button
        self.button_run_EleNum_Offset.grid(row=7, column=1, padx=5, pady=5, sticky=tk.N)
        # Noting text
        self.label_Ele_num_offset.grid(row=8, column=0, columnspan=3, padx=5, pady=5, sticky=tk.W)
    
    def browse_Offset_file(self):
        # 打开文件对话框
        file_path = filedialog.askopenfilenames(filetypes=[("Inp files", "*.inp"), ("All files", "*.*")])
        # 将选择的文件路径显示在输入框中
        self.Offset_origin_file_entry.delete(0, tk.END)
        self.Offset_origin_file_entry.insert(0, list(file_path))

    def run_Ele_num_offset(self):
        # Get file path from input box
        Offset_files_path = self.Offset_origin_file_entry.get()
        Offset_files_path_list = Offset_files_path.split(' ')
        print(Offset_files_path_list)
        # Get the magnitude of offset
        Offset_Magnitude = self.Offset_Magnitude_entry.get()
        if Offset_files_path == None:
            tk.messagebox.showerror("错误", "找不到Inp文件！")
            return
        elif len(Offset_files_path_list) == 1:
            Offset_file_path = Offset_files_path_list[0]
            [Offset_file_folder, Offset_file]= os.path.split(Offset_file_path)        
            # *******************************************************************************************
            # Writing layer information into a txt file, which will be read by EleNum_Offset.py
            layer_info = [[self.First_layer_input_entry.get(),  self.First_layer_multiplier_entry.get()],
                        [self.Second_layer_input_entry.get(), self.Second_layer_multiplier_entry.get()],
                        [self.Third_layer_input_entry.get(),  self.Third_layer_multiplier_entry.get()],
                        [self.Fourth_layer_input_entry.get(), self.Fourth_layer_multiplier_entry.get()]]
            layer_info_file_path = os.path.join(Offset_file_folder,'layer_info.txt')
            with open(layer_info_file_path, 'w') as layer_info_file:
                for layer in layer_info:
                    # layer[0] is the layer number, layer[1] is the multiplier
                    layer_info_file.write(str(layer[0]) + ',' + str(layer[1]) + '\n')
            # *******************************************************************************************
            EleNum_Offset_py_file_path = os.path.join(working_directory,'src\EleNum_Offset.py')
            # Execute EleNum_Offset.py
            try:
                subprocess.call(['python', EleNum_Offset_py_file_path, '-w '+Offset_file_path, '-m '+Offset_Magnitude])
            except FileNotFoundError:
                tk.messagebox.showerror("错误", "找不到Inp文件或Python执行文件！")
        elif len(Offset_files_path_list) > 1:
            print('There are more than one inp files')
            for Offset_file_path in Offset_files_path_list:
                [Offset_file_folder, Offset_file]= os.path.split(Offset_file_path)
                # *******************************************************************************************
                # Writing layer information into a txt file, which will be read by EleNum_Offset.py
                layer_info = [[self.First_layer_input_entry.get(),  self.First_layer_multiplier_entry.get()],
                            [self.Second_layer_input_entry.get(), self.Second_layer_multiplier_entry.get()],
                            [self.Third_layer_input_entry.get(),  self.Third_layer_multiplier_entry.get()],
                            [self.Fourth_layer_input_entry.get(), self.Fourth_layer_multiplier_entry.get()]]
                layer_info_file_path = os.path.join(Offset_file_folder,'layer_info.txt')
                with open(layer_info_file_path, 'w') as layer_info_file:
                    for layer in layer_info:
                        # layer[0] is the layer number, layer[1] is the multiplier
                        layer_info_file.write(str(layer[0]) + ',' + str(layer[1]) + '\n')
                # *******************************************************************************************
                EleNum_Offset_py_file_path = os.path.join(working_directory,'src\EleNum_Offset.py')
                # Execute EleNum_Offset.py
                try:
                    subprocess.call(['python', EleNum_Offset_py_file_path, '-w '+Offset_file_path, '-m '+Offset_Magnitude])
                except FileNotFoundError:
                    tk.messagebox.showerror("错误", "找不到Inp文件或Python执行文件！")                

class create_odb_extract_gui():
    def __init__(self, parent):
        self.tab_odb_extract = parent
        # Get all the files in the folder
        self.odb_extract_file_label = ttk.Label(self.tab_odb_extract, text="ODB file path")
        self.odb_extract_file_entry = ttk.Entry(self.tab_odb_extract, width=40)
        self.button_select_odb_file = ttk.Button(self.tab_odb_extract, text="Select", 
                                            command=self.browse_odb_file)        
        # Create a combobox for displaying available output options
        self.Output_options_label = ttk.Label(self.tab_odb_extract, text="Output options:")
        self.Output_options_combobox = ttk.Combobox(self.tab_odb_extract, width=10)
        self.Output_options_combobox['values'] = ('S','U','LE','EVOL','UVARM','SDV')
        self.Output_suboptions_label = ttk.Label(self.tab_odb_extract, text="Suboptions:")
        self.Output_suboptions_entry = ttk.Entry(self.tab_odb_extract, width=10)
        # Get necessary information for odb extract
        #---------------------------------Dimension---------------------------------#
        self.Dimension_label = ttk.Label(self.tab_odb_extract, text="Dimension:")
        self.Dimension_entry = ttk.Entry(self.tab_odb_extract, width=10)
        self.Dimension_entry.insert(tk.END, "2")
        #-----------------------------------step------------------------------------#
        self.Step_label = ttk.Label(self.tab_odb_extract, text="Step:")
        self.Step_entry = ttk.Entry(self.tab_odb_extract, width=10)
        self.Step_entry.insert(tk.END, "-1")
        self.Step_description_label = ttk.Label(self.tab_odb_extract, text="0: All step; -1: Last step")
        #-----------------------------------frame-----------------------------------#
        self.Frame_label = ttk.Label(self.tab_odb_extract, text="Frame:")
        self.Frame_entry = ttk.Entry(self.tab_odb_extract, width=10)
        self.Frame_entry.insert(tk.END, "0")
        self.Frame_description_label = ttk.Label(self.tab_odb_extract, text="0: All frame; -1: Last frame")
        #--------------------------------Dividing line------------------------------#
        self.Dividing_line_label = ttk.Label(self.tab_odb_extract, text="-----------------------------------------------------------------------------------")
        #-----------------------Parameters for figure output------------------------#
        self.Figure_x_axis_min_label = ttk.Label(self.tab_odb_extract, text="X-min:")
        self.Figure_x_axis_min_entry = ttk.Entry(self.tab_odb_extract, width=10)
        self.Figure_x_axis_min_entry.insert(tk.END, "0")
        self.Figure_x_axis_max_label = ttk.Label(self.tab_odb_extract, text="X-max:")
        self.Figure_x_axis_max_entry = ttk.Entry(self.tab_odb_extract, width=10)
        self.Figure_x_axis_max_entry.insert(tk.END, "0")
        self.Figure_y_axis_min_label = ttk.Label(self.tab_odb_extract, text="Y-min:")
        self.Figure_y_axis_min_entry = ttk.Entry(self.tab_odb_extract, width=10)
        self.Figure_y_axis_min_entry.insert(tk.END, "0")
        self.Figure_y_axis_max_label = ttk.Label(self.tab_odb_extract, text="Y-max:")
        self.Figure_y_axis_max_entry = ttk.Entry(self.tab_odb_extract, width=10)
        self.Figure_y_axis_max_entry.insert(tk.END, "0")
        self.Figure_description_label = ttk.Label(self.tab_odb_extract, text="Don't input anything if you want to use the default value")
        # Creating the execute button
        self.Odb_extract_button = ttk.Button(self.tab_odb_extract,text='Execute',
                                              command=lambda: Run_in_new_thread(self.execute_odb_extract))
        # *********************************************************************** #
        # Layout the GUI
        # *********************************************************************** #
        # Get odb file
        self.odb_extract_file_label.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.odb_extract_file_entry.grid(row=0, column=1, columnspan=3, padx=5, pady=5, sticky=tk.W)
        self.button_select_odb_file.grid(row=0, column=4, padx=5, pady=5, sticky=tk.W)
        # Get output options
        self.Output_options_label.grid(row=1, column=0, padx=5, pady=5,   sticky=tk.W)
        self.Output_options_combobox.grid(row=1, column=1, padx=5, pady=5,   sticky=tk.W)
        self.Output_suboptions_label.grid(row=1, column=2, padx=5, pady=5,   sticky=tk.W)
        self.Output_suboptions_entry.grid(row=1, column=3, padx=5, pady=5,   sticky=tk.W)
        # Get necessary information for odb extract
        # Dimension
        self.Dimension_label.grid(row=2, column=0, padx=5, pady=5,   sticky=tk.W)
        self.Dimension_entry.grid(row=2, column=1, padx=5, pady=5,   sticky=tk.W)
        # Step
        self.Step_label.grid(row=3, column=0, padx=5, pady=5,   sticky=tk.W)
        self.Step_entry.grid(row=3, column=1, padx=5, pady=5,   sticky=tk.W)
        self.Step_description_label.grid(row=3, column=2, columnspan=4, padx=5, pady=5,  sticky=tk.W)
        # Frame
        self.Frame_label.grid(row=4, column=0, padx=5, pady=5,   sticky=tk.W)
        self.Frame_entry.grid(row=4, column=1, padx=5, pady=5,   sticky=tk.W)
        self.Frame_description_label.grid(row=4, column=2, columnspan=4, padx=5, pady=5,  sticky=tk.W)
        # Dividing line
        self.Dividing_line_label.grid(row=5, column=0, columnspan=4, padx=5, pady=5,  sticky=tk.W)
        # Parameters for figure output
        self.Figure_x_axis_min_label.grid(row=6, column=0, padx=5, pady=5,   sticky=tk.W)
        self.Figure_x_axis_min_entry.grid(row=6, column=1, padx=5, pady=5,   sticky=tk.W)
        self.Figure_x_axis_max_label.grid(row=6, column=2, padx=5, pady=5,   sticky=tk.W)
        self.Figure_x_axis_max_entry.grid(row=6, column=3, padx=5, pady=5,   sticky=tk.W)
        self.Figure_y_axis_min_label.grid(row=7, column=0, padx=5, pady=5,   sticky=tk.W)
        self.Figure_y_axis_min_entry.grid(row=7, column=1, padx=5, pady=5,   sticky=tk.W)
        self.Figure_y_axis_max_label.grid(row=7, column=2, padx=5, pady=5,   sticky=tk.W)
        self.Figure_y_axis_max_entry.grid(row=7, column=3, padx=5, pady=5,   sticky=tk.W)
        self.Figure_description_label.grid(row=8, column=0, columnspan=4, padx=5, pady=5,  sticky=tk.W)
        # Execution button
        self.Odb_extract_button.grid(row=9, column=2, padx=5, pady=5,   sticky=tk.W)
        # # Create a combobox for displaying inp files
        # self.odb_file_label = ttk.Label(self.tab_odb_extract, text="Select odb file:")
        # self.odb_file_combobox = ttk.Combobox(self.tab_odb_extract, width=37)
        
    def browse_odb_file(self):
        # 打开文件对话框
        file_path = filedialog.askopenfilename(filetypes=[("ODB files", "*.odb"), ("All files", "*.*")])
        # 将选择的文件路径显示在输入框中
        self.odb_extract_file_entry.delete(0, tk.END)
        self.odb_extract_file_entry.insert(0, file_path)
            
    def execute_odb_extract(self):
        odb_file_path = self.odb_extract_file_entry.get()
        output_option = self.Output_options_combobox.get()
        output_suboption = self.Output_suboptions_entry.get()
        dimension = self.Dimension_entry.get()
        specified_step = self.Step_entry.get()
        specified_frame = self.Frame_entry.get()
        x_axis_min = self.Figure_x_axis_min_entry.get()
        x_axis_max = self.Figure_x_axis_max_entry.get()
        y_axis_min = self.Figure_y_axis_min_entry.get()
        y_axis_max = self.Figure_y_axis_max_entry.get()
        if odb_file_path == None:
            tk.messagebox.showerror("Error", "Please select an file")
            return
        else:
            odb_extract_py_file_path = os.path.join(working_directory,'src\Function_ODB_Extract.py')
            script_path = odb_extract_py_file_path
            abaqus_command_line = ['D:/ABAQUS/Commands/abaqus.bat','python', script_path,
                                   '-f',  odb_file_path,
                                   '-o',  output_option,
                                   '-so', output_suboption,
                                   '-d',  dimension,
                                   '-s',  specified_step,
                                   '-r',  specified_frame,
                                   '-x1', x_axis_min,
                                   '-x2', x_axis_max,
                                   '-y1', y_axis_min,
                                   '-y2', y_axis_max]
            
            try:
                subprocess.call(abaqus_command_line)
                # subprocess.call(['D:/ABAQUS/Commands/abaqus.bat', 'job='+input_file, 'user='+for_file,
                #                  'scratch="D:/ABAQUS/temp/CMD_Temp"', 'int', 'ask=off'])
            except FileNotFoundError:
                tk.messagebox.showerror("Error", "Execute command line failed. Please check the path of the abaqus.bat file and modify it in AbaqusTools.py")
                print(abaqus_command_line)
        
        pass
    # def update_combobox(self):
        # Display all field output variables in the combobox
        # Get the keywoards of the field output variables
        # self.current_frame.fieldOutputs.keys()
        
        # All_files = os.listdir(self.odb_extract_working_directory_entry.get())

        # # Filter out only the .inp files
        # CMD_inp_files = [f for f in All_files if f.endswith(".inp")]        
        # # Update the combobox with the inp files
        # self.CMD_inp_combobox["values"] = CMD_inp_files   
        # # Filter out only the .for files
        # CMD_for_files = [f for f in All_files if f.endswith(".for")]
        # # Update the combobox with the for files
        # self.CMD_for_combobox["values"] = CMD_for_files
        
class Run_in_new_thread():
    def __init__(self, func):
        self.func = func
        # Create a thread to run the function
        threading.Thread(target=func).start()

                
if __name__ == '__main__':
    working_directory = os.getcwd()
    os.chdir(working_directory)
    app = AbaqusGUI()
    app.mainloop()
