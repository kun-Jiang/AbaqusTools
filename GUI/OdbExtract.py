import os
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import subprocess
from utilities.RunInNewThread import Run_in_new_thread
from utilities.GetEnv import Get_ABAQUS_env


class create_odb_extract_gui():
    def __init__(self, parent, working_directory):
        self.tab_odb_extract = parent
        self.working_directory = working_directory
        # Get all the files in the folder
        self.odb_extract_file_label = ttk.Label(self.tab_odb_extract, text="ODB file path")
        self.odb_extract_file_entry = ttk.Entry(self.tab_odb_extract, width=40)
        self.button_select_odb_file = ttk.Button(self.tab_odb_extract, text="Select", 
                                            command=self.browse_odb_file)        
        # Create a combobox for displaying available output options
        self.Output_options_label = ttk.Label(self.tab_odb_extract, text="Output options:")
        self.Output_options_combobox = ttk.Combobox(self.tab_odb_extract, width=10)
        self.Output_options_combobox['values'] = ('S','U','E','LE','EVOL','UVARM','SDV')
        self.Output_suboptions_label = ttk.Label(self.tab_odb_extract, text="Suboptions:")
        self.Output_suboptions_entry = ttk.Entry(self.tab_odb_extract, width=10)
        self.Output_suboptions_entry.insert(tk.END, '1')
        self.Output_options_combobox.current(5)
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
        ABAQUS_env, ABAQUS_Execute = Get_ABAQUS_env()
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
            odb_extract_py_file_path = os.path.join(self.working_directory,'src\Function_ODB_Extract.py')
            script_path = odb_extract_py_file_path
            abaqus_command_line = [ABAQUS_Execute,'python', script_path,
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