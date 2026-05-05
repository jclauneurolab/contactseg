import os
from appdirs import AppDirs


def get_cache_dir():
    if "CONTACTSEG_CACHE_DIR" in os.environ.keys():
        cache_dir = os.environ["CONTACTSEG_CACHE_DIR"]
    else:
        # Create local download dir for atlases if it doesn't exist
        # This keeps it separate from nnUNet and autoafids
        dirs = AppDirs("contactseg", "khanlab")
        cache_dir = dirs.user_cache_dir
    return cache_dir
