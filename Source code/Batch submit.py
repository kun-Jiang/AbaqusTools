import subprocess
import os
# *************************************************************************************************
# Set the Abaqus command line parameters
# *************************************************************************************************
inp_file_name = 'Model.inp'
for_file_name = 'subroutine.for'
abaqus_job    = 'job='+inp_file_name
abaqus_user   = 'user='+for_file_name
cpu_cores = 2
abaqus_bat_path = 'D:/ABAQUS/2021/Commands/abaqus.bat'
abaqus_cpu_cores = 'cpus='+str(cpu_cores)
abaqus_Int = 'int'
abaqus_ask = 'ask=off'
# *************************************************************************************************
# Create the working directory
# *************************************************************************************************
dir = []
Job_label = []
working_directory = 'D:\Desktop\Diffusion\Simple_Model'
# The temporary folder for the scratch files
scratch_path = os.path.join(working_directory, 'Temp')
sub_job_folder = ['DPf1e-6','DPf1e-8','DPf1e-10','DPf1e-12','kr1e-2','kr1e-4','kr1e-6','kr1e-8']
for folder in sub_job_folder:
    dir.append(os.path.join(working_directory, folder))
    # This label is just for the display of the job status, it's not actually the name of the job
    Job_label.append(folder + '\\'+ '%s'%inp_file_name)

# *************************************************************************************************
# Create the initial CMD window
# 
cmd = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
# *************************************************************************************************
# Create the scratch folder
# *************************************************************************************************
if os.path.exists(scratch_path) == False:
    try:
        os.makedirs(scratch_path)
    except:
        print("Error", scratch_path + '\n'
            "You do not have the permission to create the scratch folder! \
            Please check the path and create the folder manually. If these don\'t work,\
            modifying the AbaqusTools.py file to set a appropriate scratch path.")
# *************************************************************************************************
# Batch submit the job
# *************************************************************************************************
Job_status = ['Waiting']*len(dir)
for i in range(len(dir)):
    # *************************************************************************************************
    # Display the job sequence
    # *************************************************************************************************
    print('*'*60 + '\n' +
            '*' + '{0:^58}'.format('Abaqus batch submisson') + '*' + '\n' +
            '*' + '{0:^58}'.format('Job sequence:%d of %d' %((i+1),len(dir))) + '*' + '\n' +
            '*'*60)
    # *************************************************************************************************
    # Display the job status
    # *************************************************************************************************
    Job_status[i] = 'Running'
    print('{0:^25}'.format('Job name') + '{0:^25}'.format('Job status') + '\n' + '-'*50)
    for Job_index,Job in enumerate(Job_status):
        print('{0:^25}'.format(Job_label[Job_index] + ' :') + '{0:^25}'.format(Job))
    # *************************************************************************************************
    # Set the working directory
    os.chdir(dir[i])
    # Remove the lock file
    all_file = os.listdir(dir[i])
    lock_file = [file for file in all_file if file.endswith('.lck')]
    if lock_file == []:
        pass
    else:
        os.remove(lock_file[0])
        print('The lock file has been removed!')
    # *************************************************************************************************
    # Construct the command line
    # *************************************************************************************************
    abaqus_scratch_path = 'scratch='+scratch_path
    abaqus_command_line = ('start /wait cmd /C ' + abaqus_bat_path + ' ' + abaqus_job + ' ' + abaqus_user + ' ' +  
                            abaqus_cpu_cores + ' ' + abaqus_scratch_path + ' ' + abaqus_Int + ' ' + abaqus_ask)
    # Display the command line
    print('{0:-^50}'.format('The ABAQUS command line'))
    print('start /wait cmd /C ' + '\n' +
            abaqus_bat_path + ' ' + '\n' +
            '\t' + abaqus_job + ' ' + '\n' +
            '\t' + abaqus_user + ' ' + '\n' +
            '\t' + abaqus_cpu_cores + ' ' + '\n' +
            '\t' + abaqus_scratch_path + ' ' + '\n' +
            '\t' + abaqus_Int + ' ' + '\n' +
            '\t' + abaqus_ask)
    print('{0:-^50}'.format('Submit the job'))
    # *************************************************************************************************
    try:
        job = subprocess.Popen(abaqus_command_line,creationflags=subprocess.CREATE_NEW_CONSOLE, shell=True)
        print('Don\'t close the window! The job is running ... ')
        # Wait for the job to finish
        job.wait()
        if job == 0:
            Job_status[i] = 'Success'
        else:
            Job_status[i] = 'Fail'
    except FileNotFoundError:
        print('The Abaqus command line is fault!')
        os._exit(0)
    os._exit(0)