import os
import stat
import shutil

# code taken from:
# http://stackoverflow.com/questions/1889597/deleting-directory-in-python


def _remove_readonly(fn, path_, excinfo):
    # Handle read-only files and directories
    if fn is os.rmdir:
        os.chmod(path_, stat.S_IWRITE)
        os.rmdir(path_)
    elif fn is os.remove:
        os.lchmod(path_, stat.S_IWRITE)
        os.remove(path_)


def force_remove_file_or_symlink(path_):
    try:
        os.remove(path_)
    except OSError:
        os.lchmod(path_, stat.S_IWRITE)
        os.remove(path_)


# Code from shutil.rmtree()
def is_regular_dir(path_):
    try:
        mode = os.lstat(path_).st_mode
    except os.error:
        mode = 0
    return stat.S_ISDIR(mode)


def clear_dir(path_):
    if not os.path.isdir(path_):
        print('Path: {} does not exist, creating new'.format(path_))
        os.makedirs(path_)
    elif is_regular_dir(path_):
        # Given path is a directory, clear its content
        for name in os.listdir(path_):
            fullpath = os.path.join(path_, name)
            if is_regular_dir(fullpath):
                shutil.rmtree(fullpath, onerror=_remove_readonly)
            else:
                force_remove_file_or_symlink(fullpath)
    else:
        # Given path is a file or a symlink.
        # Raise an exception here to avoid accidentally clearing the content
        # of a symbolic linked directory.
        raise OSError("Cannot call clear_dir() on a symbolic link")