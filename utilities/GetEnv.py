import os

def Get_ABAQUS_env():
    env_dist = os.environ
    PATH_address = env_dist['PATH']
    # print(PATH_address)
    PATH_address = PATH_address.split(';')
    for address in PATH_address:
        if ('aba' or 'abaqus') in address.lower():
            if 'command' in address.lower():
                ABAQUS_env = address
            for root, dirs, files in os.walk(address):
                try:
                    for dir in dirs:
                        if 'command' in dir.lower():
                            ABAQUS_env = address
                            break
                except:
                    ABAQUS_env = address
            break
    files = os.listdir(ABAQUS_env)
    for file in files:
        if 'abaqus' in file.lower():
            ABAQUS_Executable = os.path.join(ABAQUS_env, file)
    return ABAQUS_env, ABAQUS_Executable


if __name__ == '__main__':
    # This is for testing
    test = Get_ABAQUS_env()
    print(test)
