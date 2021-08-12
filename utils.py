import os

def create_if_not_exist(path):
    if os.path.exists(path):
        return
    os.makedirs(path)
