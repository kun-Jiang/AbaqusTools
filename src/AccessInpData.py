"""
Access specific information from the .inp file.

functions:
    AccessNodes(Inp_file_path:str)->dict
    AccessElements(Inp_file_path:str)->dict
    AccessSets(Inp_file_path:str)->dict
"""

__author__ = 'Kun Jiang'
__email__  = '2305050680jk@gmail.com'
__date__   = '2023-09-03'
__url__ = 'https://github.com/kun-Jiang/AbaqusTools'

def AccessNodes(inp_file_path:str)->dict:
    """
    Access the nodes' information from the .inp file.

    Args:
        Inp_file_path (str): the path of the .inp file.
    """
    print('Accessing the nodes information from the .inp file...')
    return
    coords = dict()
    with open(inp_file_path,'r') as Inpfile:
        Inplines = Inpfile.readlines()
        for line in Inplines:
            # Remove whitespace characters from both ends of the string
            line.strip()
            if line.startswith('*NODE'):
                break
    
    
    
    
    
    return coords
    
def AccessElements(inp_file_path:str)->dict:
    """Access the elements' information from the .inp file.

    Args:
        Inp_file_path (str): the path of the .inp file.
    """
    print('Accessing the elements information from the .inp file...')
    return

def AccessSets(inp_file_path:str)->dict:
    """Access the sets' information from the .inp file.

    Args:
        Inp_file_path (str): the path of the .inp file.
    """
    print('Accessing the sets information from the .inp file...')
    return