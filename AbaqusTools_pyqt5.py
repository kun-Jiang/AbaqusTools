import sys
import subprocess
import os
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QMessageBox, QHBoxLayout, QFormLayout

class AbaqusGUI(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        # Setting the title and size of the window
        self.setWindowTitle("ABAQUS GUI")
        self.setGeometry(500, 500, 620, 400)

        # Creating the tab widget
        self.MainWindow = QtWidgets.QTabWidget(self)
        self.MainWindow.setGeometry(0, 0, 620, 400)
        self.MainWindow.setStyleSheet("font-size: 20px")
        # Creating tab Ele_num_offset
        create_Ele_num_offset_gui(self.MainWindow)

class create_Ele_num_offset_gui(QtWidgets.QWidget):
    def __init__(self, parent):
        # Inherit from QWidget
        super().__init__(parent)
        formLayout = QFormLayout()
        # Creating tab Ele_num_offset, and add it to MainWindow
        self.tab_Ele_num_offset = QtWidgets.QWidget(parent)
        parent.addTab(self.tab_Ele_num_offset, "Offset")
        # Setting the font size of tab Ele_num_offset
        # self.tab_Ele_num_offset.setStyleSheet("font-size: 22px")
        
        self.Button_browse_offset = QtWidgets.QPushButton(self.tab_Ele_num_offset)
        self.Button_browse_offset.setText("选择文件")
        self.Button_browse_offset.clicked.connect(self.browse_file_button)
        self.Offset_origin_file_entry = QtWidgets.QLineEdit(self.tab_Ele_num_offset)
        # self.Offset_origin_file_entry.setGeometry(100, 10, 200, 40)
        line = QHBoxLayout()
        line.addWidget(self.Offset_origin_file_entry)
        line.addWidget(self.Button_browse_offset)
        # self.setLayout(line)
        
        # formLayout.addRow(self.tab_Ele_num_offset)
        formLayout.addRow(line)
        self.setLayout(formLayout)
        
        QtWidgets.QGridLayout(self.tab_Ele_num_offset)
    def browse_file_button(self):
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(self.tab_Ele_num_offset, '选择文件', '', 'Inp Files (*.inp)')
        if file_name:
            self.Offset_origin_file_entry.setText(file_name)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    # app.setStyleSheet("""
    # QPushButton {
    #     border-style: outset;
    #     border-width: 1px;
    #     border-radius: 6px;
    #     border-color: rgb(208, 208, 208);
    #     font: bold 14px;
    #     min-width: 10em;
    #     padding: 6px;
    # }
    # QPushButton:hover {
    #     background-color: rgb(233, 233, 233);
    # }
    # QPushButton:pressed {
    #     background-color: rgb(219, 210, 210);
    #     border-style: inset;
    # }""")
    gui = AbaqusGUI()
    gui.show()
    sys.exit(app.exec_())
