# -*- coding: utf-8 -*-
from odbAccess import *
from textRepr import prettyPrint as prettyprint
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import threading
import argparse

def logwrite(message):
    log_file_path = os.path.join(extract_folder, 'extract_odb.log')
    with open(log_file_path,'a') as log_file:
        log_file.write(message+'\n')
        print(message)
        log_file.close()
        
def Remove_file(folder, extension):
    for file in os.listdir(folder):
        if file.endswith(extension):
            logwrite('Remove file: %s'%file)
            os.remove(os.path.join(folder, file))

class plot_x_y_curve():
    def __init__(self,xmin, xmax, ymin, ymax,data_array,x_label,y_label,title,file_name):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.data_x = data_array[:,0]
        self.data_y = data_array[:,1]
        self.label_x = x_label
        self.label_y = y_label
        self.title_name = title
        self.export_file_name = file_name
        self.plot()
    def plot(self):
        # Plot the curve
        logwrite('Plot the curve')
        fig=plt.figure(num=1,figsize=(6,6))
        plt.plot(self.data_x,self.data_y,'r-',linewidth=1)
        plt.title(self.title_name,fontsize=16, weight='bold', fontname='sans-serif')
        plt.xlabel(self.label_x,fontsize=16, weight='bold', fontname='sans-serif')
        plt.ylabel(self.label_y,fontsize=16, weight='bold', fontname='sans-serif')
        if (self.xmax - self.xmax) != 0:
            plt.xlim(self.xmin, self.xmax)  # 设置 x 轴范围
        if (self.ymax - self.ymin) != 0:
            plt.ylim(self.ymin, self.ymax)  # 设置 y 轴范围 
        # Save the figure in each sub folder
        fig_path = os.path.join(extract_folder,'%s.png'%self.export_file_name)
        plt.savefig(fig_path,dpi=300,bbox_inches='tight')
        plt.close(fig)
        plt.show()
        plt.show(block=True)

def odb_access(odb_dir, odb_file, output_option, output_index, dimension, 
                specified_step, specified_frame):
    os.chdir(odb_dir)
    odb_name = odb_file.split('.')[0]
    # *********************************************************************************************************************
    logwrite('*'*70 + '\n' +
             '*' + '{0:^68}'.format('Extracting Data from ODB') + '*' + '\n' +
             '*'*70)
    # Remove the lock file to prevent access denial.
    Remove_file(odb_dir, 'lck')
    # *********************************************************************************************************************
    # Display the information of ODB file and output field parameter
    logwrite('{0:<25s}'.format('ODB file name') + ':' + odb_file)
    if output_index != '0':
        output_option = output_option + output_index
    logwrite('{0:<25s}'.format('Output field parameter') + ':' + output_option + '\n' +
             '-'*60)
    # Open the ODB file
    try:
        odb_file_path = os.path.join(odb_dir,odb_file)
        logwrite('Odb file path: ' + odb_file_path)
        odb = openOdb(odb_file_path, readOnly=True)
        logwrite("Open the ODB file successfully")
    except:      
        try:
            # The ODB file may be created by the old version of Abaqus, so we need to upgrade the ODB file
            upgrade_odb_file = 'upgraded_' + odb_file
            upgrade_odb_path = os.path.join(odb_dir, upgrade_odb_file)
            logwrite("Fail to open the ODB file, try to upgrade the ODB file" + '\n' +
                     'New odb file path: ' + upgrade_odb_path)
            if os.path.exists(upgrade_odb_path):
                logwrite('The upgraded Odb file already exists.')
            else:
                upgradeOdb(existingOdbPath=odb_file_path, upgradedOdbPath=upgrade_odb_path)
            logwrite("Upgrade the ODB file successfully")
            odb_file_path = upgrade_odb_path
            odb = openOdb(odb_file_path, readOnly=True)
        except:
            logwrite('Fail to open the Odb file, please check it.')
    # *********************************************************************************************************************
    # If the specified_step = 0, extract the data in every step
    #                       = -1, extract the data in the last step
    #                       = n, extract the data in the nth step
    # We can obtain the step name by odb.steps.keys()
    # This code is to get string of step name
    if specified_step == 0:
        steps_name = odb.steps.keys()
        step_num = len(steps_name)
        logwrite('Extract the data from every step...')
    else:
        step_num = 1
        # Add the specified step name into the list
        if specified_step == -1:
            steps_name = [odb.steps.keys()[-1]]
            logwrite('Extract the data from the last setp...')
        else:
            steps_name = [odb.steps.keys()[specified_step-1]]
            logwrite('Extract the data from the ' + str(specified_frame) + 'th step...')
            
            
            
    
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The script of spliting the ABAQUS inp file')
    parser.add_argument('-f', '--odb_file', type=str, metavar='', required=True, help='odb file path')
    parser.add_argument('-o', '--output', type=str, metavar='', required=True, help='output option')
    parser.add_argument('-so', '--suboption', type=str, metavar='', required=False, help='suboption')
    parser.add_argument('-d', '--dimension', type=int, metavar='', required=True, help='dimension of the model')
    parser.add_argument('-s', '--step', type=int, metavar='', required=False, help='specified step')
    parser.add_argument('-r', '--frame', type=int, metavar='', required=False, help='specified frame')
    parser.add_argument('-x1', '--xmin', type=float, metavar='', required=False, help='minimum value of x axis')
    parser.add_argument('-x2', '--xmax', type=float, metavar='', required=False, help='maximum value of x axis')
    parser.add_argument('-y1', '--ymin', type=float, metavar='', required=False, help='minimum value of y axis')
    parser.add_argument('-y2', '--ymax', type=float, metavar='', required=False, help='maximum value of y axis')
    args = parser.parse_args()
    # *********************************************************************************************************************
    odb_file_path = args.odb_file
    output_option = args.output
    output_index = args.suboption
    dimension = args.dimension
    specified_step = args.step
    specified_frame = args.frame
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    [odb_dir, odb_file]= os.path.split(odb_file_path)
    # Check if the log file exists
    logfile_path = os.path.join(odb_dir,'extract_odb.log')
    if os.path.exists(logfile_path):
        os.remove(logfile_path)
    global root_dir
    global extract_folder
    root_dir = odb_dir
    # Create the folder to store the extracted data
    extract_folder = os.path.join(root_dir, 'Extracted Results')
    if not os.path.exists(extract_folder):
        os.mkdir(extract_folder)
    # Acess the .odb database
    integrated_list = odb_access(odb_dir, odb_file, output_option, output_index, dimension,
               specified_step, specified_frame)
    # Generate a plot showing the integral of the variable over time.
    if integrated_list is not None:
        plot_x_y_curve(xmin, xmax, ymin, ymax)