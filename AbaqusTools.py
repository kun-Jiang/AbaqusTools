import sys
import os
# Obtain the main folder of the package
main_folder = os.path.dirname(os.path.abspath(__file__))
# Add the main folder to the path
sys.path.insert(0, main_folder)

import threading
import tkinter as tk
from tkinter import ttk
import subprocess
from tkinter import filedialog
import tkinter.messagebox

from gui import CommandSubmit
from gui.InpFileSplit import create_InpSplit_gui
from gui.ElementNumOffset import create_offset_gui
from gui.OdbExtract import create_odb_extract_gui

class AbaqusGUI(tk.Tk):
    def __init__(self):
        # 设置窗口标题和窗口大小
        super().__init__()
        self.title("ABAQUS GUI")
        self.geometry("520x400")
        self.style = ttk.Style()
        self.style.configure('TNotebook.Tab', padding=(20, 5))
        self.create_widgets()
    def create_widgets(self):
        # 创建界面切换按钮

        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True)
        # 创建CMD选项卡
        self.tab_cmd = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_cmd, text='CMD')
        JobSubmit = CommandSubmit.create_cmd_gui(self.tab_cmd)
        # Create InpSplit widget
        self.tab_inp_split = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_inp_split, text='InpSplit')
        InpSplit = create_InpSplit_gui(self.tab_inp_split)
        # Create Element_Num_Offset widget
        self.tab_offset = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_offset, text='Offset')
        Element_Num_Offset = create_offset_gui(self.tab_offset)
        # Create odb_extract_gui widget
        self.tab_OdbExtract = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_OdbExtract, text='OdbExtract')
        Odb_Extract = create_odb_extract_gui(self.tab_OdbExtract, working_directory)
    

if __name__ == '__main__':
    working_directory = main_folder
    os.chdir(working_directory)
    app = AbaqusGUI()
    app.mainloop()
    