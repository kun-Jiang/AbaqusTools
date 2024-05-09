import os
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import subprocess
from utilities.RunInNewThread import Run_in_new_thread
from src.NumOffset import Num_Offset


class create_offset_gui():
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
                layer_info_file.write('Offset magnitude' + ',' + str(Offset_Magnitude) + '\n')
            # *******************************************************************************************
            Num_Offset(Offset_file_path, int(Offset_Magnitude))
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
                Num_Offset(Offset_file_path, int(Offset_Magnitude))         
