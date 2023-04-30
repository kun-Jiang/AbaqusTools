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

        # 创建界面切换按钮
        self.style = ttk.Style()
        self.style.configure('TNotebook.Tab', padding=(30, 8))
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True)

        # 创建CMD选项卡
        self.tab_cmd = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_cmd, text='CMD')
        self.create_cmd_gui()

        # 创建InpSplit选项卡
        self.tab_InpSplit = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_InpSplit,text='InpSplit')
        self.create_InpSplit_gui()
        
        # 创建Element_Num_Offset选项卡
        self.tab_Ele_num_offset = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_Ele_num_offset, text='Offset')
        self.create_Ele_num_offset_gui()

    def create_Ele_num_offset_gui(self):
        # *********************************************************************** #
        #                             Get input file                              #
        # *********************************************************************** #
        self.Offset_origin_file_label = ttk.Label(self.tab_Ele_num_offset, text='Inp文件路径:')
        self.Offset_origin_file_entry = ttk.Entry(self.tab_Ele_num_offset,width=40)
        self.Button_browse_offset = ttk.Button(self.tab_Ele_num_offset, text="选择文件", 
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
        self.Second_layer_input_entry.insert(tk.END, "Mass")
        self.Second_layer_multiplier_entry = ttk.Entry(self.tab_Ele_num_offset, width=10, justify='center')
        self.Second_layer_multiplier_entry.insert(tk.END, "1")
        # Creating the third layer input
        self.Third_layer_label = ttk.Label(self.tab_Ele_num_offset, text="Third layer")
        self.Third_layer_input_entry = ttk.Entry(self.tab_Ele_num_offset, width=40)
        self.Third_layer_input_entry.insert(tk.END, "Frac")
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
                                              command=lambda: self.Run_in_new_thread(self.run_Ele_num_offset))
        # Text of tips
        self.label_Ele_num_offset = ttk.Label(self.tab_Ele_num_offset, text='Tips:\n'
            +'The first layer elements use similar element number with the input file, so its\n'
            +'multiplier is equal to ZERO. And the offset is equal to the number of elements')
        # *********************************************************************** #
        #                        Layout the text input box                        #
        # *********************************************************************** #
        # Get input file
        self.Offset_origin_file_label.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.Offset_origin_file_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W)
        self.Button_browse_offset.grid(row=0, column=2, padx=5, pady=5, sticky=tk.W)
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
             
    def create_cmd_gui(self):
        # Get all the files in the folder
        self.CMD_working_directory_label = ttk.Label(self.tab_cmd, text="Working directory:")
        self.CMD_working_directory_entry = ttk.Entry(self.tab_cmd, width=40)
        self.button_browse_working_directory = ttk.Button(self.tab_cmd, text="Browse", 
                                            command=self.browse_working_directory)
        
        # Create a combobox for displaying inp files
        self.CMD_inp_file_lable = ttk.Label(self.tab_cmd, text="Select inp file:")
        self.CMD_inp_combobox = ttk.Combobox(self.tab_cmd, width=37)
        self.button_update_inp_combobox = ttk.Button(self.tab_cmd, text="Update",
                                            command=self.update_inp_combobox)

        # Create a combobox for displaying for files
        self.CMD_for_file_label = ttk.Label(self.tab_cmd, text="Select for file:")
        self.CMD_for_combobox = ttk.Combobox(self.tab_cmd, width=37)
        self.button_update_for_combobox = ttk.Button(self.tab_cmd, text="Update",
                                            command=self.update_for_combobox)
         
        # Create a button to run abaqus, and run it in a new thread by the function Run_in_new_thread
        self.CMD_button_run_abaqus = ttk.Button(self.tab_cmd, text="执行有限元计算", 
                                            command=lambda: self.Run_in_new_thread(self.run_abaqus))
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
        self.button_update_inp_combobox.grid(row=1, column=2, padx=5, pady=5)
        # Get for file
        self.CMD_for_file_label.grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.CMD_for_combobox.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W)
        self.button_update_for_combobox.grid(row=2, column=2, padx=5, pady=5)
        # Set CPU cores
        self.cpu_cores_label.grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
        self.cpu_cores_spinbox.grid(row=3, column=1, padx=5, pady=5, sticky=tk.W)
        # Set check button
        self.check_Int_button.grid(row=4, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W )
        self.check_ask_button.grid(row=5, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W )
        # Run button
        self.CMD_button_run_abaqus.grid(row=6, column=1, padx=5, pady=5)
        
        
    def create_InpSplit_gui(self):
        # 创建ttk风格的界面元素
        self.label_InpSplit_file = ttk.Label(self.tab_InpSplit, text='Inp文件路径:')
        self.entry_InpSplit = ttk.Entry(self.tab_InpSplit,width=40)
        self.button_browse_InpSplit = ttk.Button(self.tab_InpSplit, text="选择文件", 
                                                command=self.browse_InpSplit_file)
        
        self.button_run_InpSplit = ttk.Button(self.tab_InpSplit,text='进行Inp文件切分',
                                                command=lambda: self.Run_in_new_thread(self.run_InpSplit))
        # 界面元素布局
        self.label_InpSplit_file.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.entry_InpSplit.grid(row=0, column=1, padx=5, pady=5)
        self.button_browse_InpSplit.grid(row=0, column=2, padx=5, pady=5)

        self.button_run_InpSplit.grid(row=2, column=1, padx=5, pady=5)
        
    def run_InpSplit(self):
        # 从输入框中获取文件路径
        input_file_path = self.entry_InpSplit.get()
        InpSplit_py_file_path = os.path.join(working_directory,'Source code\InpSplit.py')
        # 执行InpSplit.py
        try:
            subprocess.call(['python', InpSplit_py_file_path, '-w '+input_file_path])
        except FileNotFoundError:
            tk.messagebox.showerror("错误", "找不到Inp文件或Python执行文件！")
            
    def run_Ele_num_offset(self):
        # Get file path from input box
        Offset_file_path = self.Offset_origin_file_entry.get()
        # Get the magnitude of offset
        Offset_Magnitude = self.Offset_Magnitude_entry.get()
        if Offset_file_path == None:
            tk.messagebox.showerror("错误", "找不到Inp文件！")
            return
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
        EleNum_Offset_py_file_path = os.path.join(working_directory,'Source code\EleNum_Offset.py')
        # Execute EleNum_Offset.py
        try:
            subprocess.call(['python', EleNum_Offset_py_file_path, '-w '+Offset_file_path, '-m '+Offset_Magnitude])
        except FileNotFoundError:
            tk.messagebox.showerror("错误", "找不到Inp文件或Python执行文件！")

    def update_inp_combobox(self):
        All_inp_files = os.listdir(self.CMD_working_directory_entry.get())

        # Filter out only the .inp files
        inp_files = [f for f in All_inp_files if f.endswith(".inp")]        
        # Update the combobox with the inp files
        self.CMD_inp_combobox["values"] = inp_files   
    
    def update_for_combobox(self):
        All_for_files = os.listdir(self.CMD_working_directory_entry.get())
        # Filter out only the .for files
        CMD_for_files = [f for f in All_for_files if f.endswith(".for")]
        # Update the combobox with the for files
        self.CMD_for_combobox["values"] = CMD_for_files

    def browse_working_directory(self):
        # 打开文件对话框
        file_path = filedialog.askdirectory()
        
        # 将选择的文件路径显示在输入框中
        self.CMD_working_directory_entry.delete(0, tk.END)
        self.CMD_working_directory_entry.insert(0, file_path)
        
    def browse_inp_file(self):
        # 打开文件对话框
        file_path = filedialog.askopenfilename()

        # 将选择的文件路径显示在输入框中
        self.CMD_inp_file_entry.delete(0, tk.END)
        self.CMD_inp_file_entry.insert(0, file_path)

    def browse_for_file(self):
        # 打开文件对话框
        file_path = filedialog.askopenfilename()

        # 将选择的文件路径显示在输入框中
        self.CMD_for_file_entry.delete(0, tk.END)
        self.CMD_for_file_entry.insert(0, file_path)
        
    def browse_InpSplit_file(self):
        # 打开文件对话框
        file_path = filedialog.askopenfilename()

        # 将选择的文件路径显示在输入框中
        self.entry_InpSplit.delete(0, tk.END)
        self.entry_InpSplit.insert(0, file_path)
        
    def browse_Offset_file(self):
        # 打开文件对话框
        file_path = filedialog.askopenfilename()

        # 将选择的文件路径显示在输入框中
        self.Offset_origin_file_entry.delete(0, tk.END)
        self.Offset_origin_file_entry.insert(0, file_path)

        
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
        abaqus_bat_path = 'D:/ABAQUS/2021/Commands/abaqus.bat'
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
            tk.messagebox.showerror("错误", "找不到ABAQUS执行文件！")
    
    def Run_in_new_thread(self, func):
        # Create a thread to run the function
        threading.Thread(target=func).start()
        
if __name__ == '__main__':
    working_directory = os.getcwd()
    os.chdir(working_directory)
    app = AbaqusGUI()
    app.mainloop()
