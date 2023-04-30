import sys
import subprocess
import os
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QMessageBox

class AbaqusGUI(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        # 设置窗口标题和窗口大小
        self.setWindowTitle("ABAQUS GUI")
        self.setGeometry(0, 0, 500, 300)

        # 创建界面切换按钮
        self.notebook = QtWidgets.QTabWidget(self)
        self.notebook.setGeometry(0, 0, 500, 300)

        # 创建CMD选项卡
        self.tab_cmd = QtWidgets.QWidget(self)
        self.notebook.addTab(self.tab_cmd, "CMD")
        self.create_cmd_gui()

        # 创建InpSplit选项卡
        self.tab_InpSplit = QtWidgets.QWidget(self)
        self.notebook.addTab(self.tab_InpSplit, "InpSplit")
        self.create_InpSplit_gui()

        # 创建Element_Num_Offset选项卡
        self.tab_Ele_num_offset = QtWidgets.QWidget(self)
        self.notebook.addTab(self.tab_Ele_num_offset, "Offset")
        self.create_Ele_num_offset_gui()

    def create_cmd_gui(self):
        pass

    def create_InpSplit_gui(self):
        pass

    def browse_Offset_file(self):
        pass

    def create_Ele_num_offset_gui(self):
        # *********************************************************************** #
        #                             Get input file                              #
        # *********************************************************************** #
        self.Offset_origin_file_label = QtWidgets.QLabel(self.tab_Ele_num_offset)
        self.Offset_origin_file_label.setText("Inp文件路径:")
        self.Offset_origin_file_label.setGeometry(10, 10, 80, 20)

        self.Offset_origin_file_entry = QtWidgets.QLineEdit(self.tab_Ele_num_offset)
        self.Offset_origin_file_entry.setGeometry(100, 10, 300, 20)

        self.Button_browse_offset = QtWidgets.QPushButton(self.tab_Ele_num_offset)
        self.Button_browse_offset.setText("选择文件")
        self.Button_browse_offset.setGeometry(410, 10, 80, 20)
        self.Button_browse_offset.clicked.connect(self.browse_Offset_file)

        # *********************************************************************** #
        #                       Creating the text input box                       #
        # *********************************************************************** #
        self.Offset_layer_name = QtWidgets.QLabel(self.tab_Ele_num_offset)
        self.Offset_layer_name.setText("Layer")
        self.Offset_layer_name.setGeometry(10, 50, 40, 20)

        self.Offset_file_name = QtWidgets.QLabel(self.tab_Ele_num_offset)
        self.Offset_file_name.setText("File name")
        self.Offset_file_name.setGeometry(90, 50, 60, 20)

        self.Offset_multiplier = QtWidgets.QLabel(self.tab_Ele_num_offset)
        self.Offset_multiplier.setText("Multiplier")
        self.Offset_multiplier.setGeometry(200, 50, 60, 20)

        # Creating the first layer input
        self.First_layer_label = QtWidgets.QLabel(self.tab_Ele_num_offset)
        self.First_layer_label.setText("First layer")
        self.First_layer_label.setGeometry(10, 80, 60, 20)

        self.First_layer_input_entry = QtWidgets.QLineEdit(self.tab_Ele_num_offset)
        self.First_layer_input_entry.setGeometry(90, 80, 300, 20)
        self.First_layer_input_entry.insert("Mech")

        self.First_layer_multiplier_entry = QtWidgets.QLineEdit(self.tab_Ele_num_offset)
        self.First_layer_multiplier_entry.setGeometry(200, 80, 100, 20)
        self.First_layer_multiplier_entry.insert("0")

        # Creating the second layer input
        self.Second_layer_label = QtWidgets
        
def InpSplit(self):
    # 获取InpSplit界面的输入值
    origin_file = self.InpSplit_origin_file_entry.get()
    file_name = self.InpSplit_file_name_entry.get()
    output_directory = self.InpSplit_output_directory_entry.get()
    layer_name = self.InpSplit_layer_name_entry.get()
    start_element = self.InpSplit_start_element_entry.get()
    end_element = self.InpSplit_end_element_entry.get()

    # 检查输入值是否有效
    if not os.path.exists(origin_file):
        self.display_message('请输入有效的文件路径')
        return
    if not os.path.exists(output_directory):
        self.display_message('请输入有效的输出路径')
        return
    if not file_name:
        self.display_message('请输入文件名')
        return
    if not layer_name:
        self.display_message('请输入层名称')
        return
    if not start_element or not end_element:
        self.display_message('请输入元素范围')
        return

    # 将输入值传递给InpSplit程序并运行
    args = ['InpSplit.exe', origin_file, output_directory, file_name, layer_name, start_element, end_element]
    subprocess.call(args)
    self.display_message('InpSplit 运行完毕')

def browse_file(self, entry):
    # 打开文件选择对话框并将选中的文件路径填充到对应的输入框中
    file_path = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(0, file_path)

def browse_directory(self, entry):
    # 打开文件夹选择对话框并将选中的文件夹路径填充到对应的输入框中
    dir_path = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, dir_path)

def display_message(self, message):
    # 在消息框中显示指定的消息
    messagebox.showinfo('消息', message)
    
if __name__ == '__main__':
    app = QApplication(sys.argv)
    app = AbaqusGUI()
    app.mainloop()
