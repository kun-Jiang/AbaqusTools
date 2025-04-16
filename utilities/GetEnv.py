import os
import shutil

def Get_ABAQUS_env():
    # Search for the abaqus.bat executable
    abaqus_executable = shutil.which("abaqus.bat")
    if abaqus_executable:
        abaqus_env = os.path.dirname(abaqus_executable)
        return abaqus_env, abaqus_executable
    else:
        raise FileNotFoundError("Abaqus executable (abaqus.bat) not found in PATH.")


if __name__ == '__main__':
    # This is for testing
    test = Get_ABAQUS_env()
    print(test)
