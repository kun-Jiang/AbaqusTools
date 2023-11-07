import threading



class Run_in_new_thread():
    def __init__(self, func):
        self.func = func
        # Create a thread to run the function
        threading.Thread(target=func).start()