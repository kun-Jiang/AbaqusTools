import os
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import subprocess
from utilities.RunInNewThread import Run_in_new_thread
from src.InpSplit import Inpsplit


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
        # self.button_run_InpSplit = ttk.Button(self.tab_InpSplit,text='Execute split',
        #                                         command=self.run_InpSplit)
        # 界面元素布局
        self.label_InpSplit_file.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.entry_InpSplit.grid(row=0, column=1, padx=5, pady=5)
        self.button_browse_InpSplit.grid(row=0, column=2, padx=5, pady=5)

        self.button_run_InpSplit.grid(row=2, column=1, padx=5, pady=5)
        
    def run_InpSplit(self):
        # 从输入框中获取文件路径
        input_file_path = self.entry_InpSplit.get()
        # InpSplit_py_file_path = os.path.join(working_directory,'src\InpSplit.py')
        # 执行InpSplit.py
        print(input_file_path)
        Inpsplit(input_file_path)
        # try:
        #     subprocess.call(['python', InpSplit_py_file_path, '-w '+input_file_path])
        # except FileNotFoundError:
        #     tk.messagebox.showerror("错误", "找不到Inp文件或Python执行文件！")
                    
    def browse_InpSplit_file(self):
        # 打开文件对话框
        file_path = filedialog.askopenfilename(filetypes=[("Inp files", "*.inp"), ("All files", "*.*")])
        # 将选择的文件路径显示在输入框中
        self.entry_InpSplit.delete(0, tk.END)
        self.entry_InpSplit.insert(0, file_path)