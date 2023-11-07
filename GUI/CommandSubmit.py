import os
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import subprocess
from Src import RunInNewThread

# from AbaqusTools.Src import Run_in_new_thread



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
                                            command=lambda: RunInNewThread(run_abaqus()))
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

