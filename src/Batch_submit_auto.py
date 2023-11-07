import subprocess
import os
import psutil


# *************************************************************************************************
# Initialize parameters
# *************************************************************************************************
# Parameters corresponding to the Abaqus command line
# Default regard there is only one inp file and one for file in each sub folder, and the program
# will read them automatically.
# ------------------------------------------------------------------------------------------------
cpu_cores = 2
# The path of the abaqus.bat file, you must modify it to your own path
command_line_bat_path = 'abaqus'
command_line_cpu_cores = 'cpus='+str(cpu_cores)
command_line_int = 'int'
command_line_ask = 'ask=off'
# Parameters corresponding to the job sequence
# ------------------------------------------------------------------------------------------------
# If you want to submit the job in a specific order or just submit some of the jobs, you can modify
# the sub_folders list manually. Otherwise, the program will submit all the jobs in the primary working directory.
# E.g. sub_folders = ['B2','B3','B4','B5'], these string are the name of the sub folders.
sub_folders = []
# *************************************************************************************************
# Get the primary working directory
# *************************************************************************************************
global working_directory
working_directory = os.getcwd()
# Create the scratch folder
scratch_path = os.path.join(working_directory, 'Temp')
if os.path.exists(scratch_path) == False:
    try:
        os.makedirs(scratch_path)
    except:
        print("Error", scratch_path + '\n'
            "You do not have the permission to create the scratch folder! \
            Please check the path and create the folder manually. If these don\'t work,\
            modifying the python file to set a appropriate scratch path.")
        os._exit(0)
# *************************************************************************************************
# Get sub folders in the primary working directory, each sub folder contains a inp file and a for file for submitting
if sub_folders == []:
    sub_folders = os.listdir(working_directory)
    sub_folders = [folder for folder in sub_folders if os.path.isdir(folder)]
# *************************************************************************************************
# Create the initial CMD window
cmd = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
# *************************************************************************************************
# Collect the inp file information and for file information in each sub folder
inp_file_list = []
for_file = []
Job_label = []
sub_folder_temp = list(sub_folders)
for sub_folder in sub_folder_temp:
    sub_folder_path = os.path.join(working_directory, sub_folder)
    # Get the inp file and for file in each sub folder
    inp_file = [file for file in os.listdir(sub_folder_path) if file.endswith('.inp')]
    for_file_name = [file for file in os.listdir(sub_folder_path) if file.endswith('.for')]
    if len(inp_file) == 1 and len(for_file_name) == 1:
        inp_file_list.append(inp_file[0])
        for_file.append(for_file_name[0])
        Job_label.append(sub_folder + '\\' + inp_file[0])
    else:
        sub_folders.remove(sub_folder)
# *************************************************************************************************
# Batch submit the job
# *************************************************************************************************
Job_status   = ['Waiting'] * len(sub_folders)
Job_time     = ['0'] * len(sub_folders)
Job_Inc      = ['0'] * len(sub_folders)
Job_sequence = 0
Job_progress = [0]*len(sub_folders)
processes    = [[]]*len(sub_folders)
# ------------------------------------------------------------------------------------------------
# The program will obtain the CPU usage at intervals, if the CPU usage is less than the warning value,
# it will continue to submit the job, otherwise it will wait for the CPU usage to drop below the warning value.
interval_time = 10
cpu_waring_value = 90
while True:
    # Get the CPU usage at intervals
    cpu_percent = psutil.cpu_percent(interval=interval_time)
    # Optimize the intervel_time
    if cpu_percent < cpu_waring_value:
        # get the CPU usage per 2 seconds
        interval_time = 10
    else:
        # get the CPU usage per one minute
        interval_time = 60
    # Clear the CMD window
    os.system('cls' if os.name == 'nt' else 'clear')
    # print("\033[F\033[K"*20, end='\r')
    print(working_directory)
    if cpu_percent <= cpu_waring_value and Job_sequence < len(sub_folders):

        # *************************************************************************************************
        # Set the working directory
        os.chdir(os.path.join(working_directory, sub_folders[Job_sequence]))
        # Remove the lock file
        all_file = os.listdir(os.getcwd())
        lock_file = [file for file in all_file if file.endswith('.lck')]
        if lock_file != []:
            os.remove(lock_file[0])
        # Remove the log file
        log_file = [file for file in all_file if file.endswith('.log')]
        if log_file != []:
            os.remove(log_file[0])
        # *************************************************************************************************
        # Construct the command line
        # *************************************************************************************************
        command_line_job = 'job='+inp_file_list[Job_sequence]
        command_line_user = 'user='+for_file[Job_sequence]
        command_line_scratch_path = 'scratch='+scratch_path
        abaqus_command_line = ('start /wait cmd /C ' + command_line_bat_path + ' ' + command_line_job + ' ' + command_line_user + ' ' +  
                                command_line_cpu_cores + ' ' + command_line_scratch_path + ' ' + command_line_ask)
        # *************************************************************************************************
        try:
            job = subprocess.Popen(abaqus_command_line,creationflags=subprocess.CREATE_NEW_CONSOLE, shell=True)
            processes[Job_sequence] = job
            Job_sequence += 1
        except FileNotFoundError:
            print('start /wait cmd /C ' + '\n' +
                    command_line_bat_path + ' ' + '\n' +
                    '\t' + command_line_job + ' ' + '\n' +
                    '\t' + command_line_user + ' ' + '\n' +
                    '\t' + command_line_cpu_cores + ' ' + '\n' +
                    '\t' + command_line_scratch_path + ' ' + '\n' +
                    # '\t' + command_line_int + ' ' + '\n' +
                    '\t' + command_line_ask)
            Job_status[Job_sequence] = 'The Abaqus command line is fault!'
            
    # *************************************************************************************************
    # Check the status of jobs
    # *************************************************************************************************
    for process_index,process in enumerate(processes):
        if process != []:
            if process.poll() == None:
                Job_status[process_index] = 'Running'
            else:
                sub_folder = sub_folders[process_index]
                sub_folder_path = os.path.join(working_directory, sub_folder)
                inp_file = inp_file_list[process_index]
                inp_file_name = inp_file.split('.')[0]
                log_file = os.path.join(sub_folder_path, inp_file_name + '.log')
                try:
                    with open(log_file, 'r') as f:
                        status = f.readlines()[-1]
                        status = status.strip('\n')
                        if 'error' in status.lower():
                            Job_status[process_index] = status
                        elif 'completed' in status.lower():
                            Job_status[process_index] = status
                except:
                    pass

                status_file = os.path.join(sub_folder_path, inp_file_name + '.sta')
                try:
                    with open(status_file, 'r') as f:
                        status = f.readlines()[-1]
                        if 'analysis' in status.lower():
                            continue
                        status = status.strip('\n')
                        status = status.split(' ')
                        status = list(filter(lambda x: x != '', status))
                        Job_time[process_index] = str(status[-2])
                        Job_Inc[process_index] = str(status[-1])
                except:
                    Job_time[process_index] = 'Fault'
                    Job_Inc[process_index] = 'Fault'
            # elif process.poll() == 0:
            #     Job_status[process_index] = 'Finished'
            # elif process.poll() > 0:
            #     Job_status[process_index] = 'Fault: ' + str(process.poll())
            # elif process.poll() < 0:
            #     Job_status[process_index] = 'Killed: ' + str(process.poll())
    # ------------------------------------------------------------------------------------------------
    # Display the job sequence
    print('*'*60 + '\n' +
            '*' + '{0:^58}'.format('Abaqus batch submisson') + '*' + '\n' +
            '*'*60)
    # ------------------------------------------------------------------------------------------------
    # Display the job status
    print('{0:^50}'.format('Job name') + '|' + '{0:^30}'.format('Job status') + '|' + '{0:^8}'.format('Time') + '|' + '{0:^8}'.format('Inc')
          + '\n' + '-'*96)
    for Job_index,Job in enumerate(Job_status):
        print('{0:^50}'.format(Job_label[Job_index]) + '|' + '{0:^30}'.format(Job) + '|' + '{0:^8}'.format(Job_time[Job_index]) + '|' + '{0:^8}'.format(Job_Inc[Job_index]) )
        if Job != 'Waiting' and Job != 'Running':
            Job_progress[Job_index] = 1
    print('-'*50 + '\n' +
          'Don\'t close the window! The job is running ... ' + '\n' +
          'CPU usage: ' + str(cpu_percent) + '%' + '\t' + 'Limit: ' + str(cpu_waring_value) + '%' + '\t' + 'Refresh time: ' + str(interval_time) + 's' )
    if sum(Job_progress) == len(sub_folders):
        print('All jobs have been completed!')
        break