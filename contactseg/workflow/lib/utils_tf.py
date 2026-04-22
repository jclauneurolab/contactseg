import os
from appdirs import AppDirs

def get_download_dir():
    if "TEMPLATEFLOW_HOME" in os.environ.keys():
        download_dir = os.environ["TEMPLATEFLOW_HOME"]
    else:
        # Create local download dir for atlases if it doesn't exist
        # This keeps it separate from nnUNet and autoafids
        dirs = AppDirs("templateflow", "khanlab")
        download_dir = dirs.user_cache_dir
    return download_dir