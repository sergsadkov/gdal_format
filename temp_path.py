import os
import shutil
import time


__temp_dir = os.path.join(os.environ['TMP'], 'image_processor')
__shp_extensions = ['shp', 'dbf', 'shx', 'prj', 'sbn', 'sbx', 'cpg']


def splitPath(path):
    if os.path.isdir(path):
        return path, '', ''
    folder, file = os.path.split(path)
    name, ext = os.path.splitext(file)
    return folder, name, ext[1:]


def fullPath(folder, file, ext=''):
    return f'{folder}\\{file}{("",".")[bool(ext)]}{ext.lstrip(".")}'


def sureDir(folder):
    if os.path.isdir(folder) and not os.path.exists(folder):
        os.makedirs()


def copyFile(path_in, path_out, overwrite=False):
    if (not os.path.exists(path_out)) or overwrite:
        if os.path.exists(path_in):
            try:
                shutil.copyfile(path_in, path_out)
            except Exception as e:
                pass


def deleteFile(file):
    if os.path.exists(file):
        try:
            os.remove(file)
        except Exception as e:
            pass


def copySHP(file_in, file_out, overwrite=False, ext_list=__shp_extensions):
    folder_in, name_in, ext_in = splitPath(file_in)
    folder_out, name_out, ext_out = splitPath(file_out)
    for ext in ext_list:
        copyFile(fullPath(folder_in, name_in, ext),
                 fullPath(folder_out, name_out, ext), overwrite=overwrite)


def deleteSHP(file, ext_list=__shp_extensions):
    folder, name, ext = splitPath(file)
    for ext in ext_list:
        deleteFile(fullPath(folder, name, ext))


def listFiles(folder, extensions=None, target_path=None, miss_path=None):

    files = []

    if extensions is not None:
        if isinstance(extensions, (tuple, list)):
            extensions = list(extensions)
        else:
            extensions = [extensions]
        exts = ['.' + str(ext).lower().lstrip('.') for ext in extensions]

    for corner, _folders, _files in os.walk(folder):
        if miss_path:
            if re.search(miss_path, corner):
                continue
        for file in _files:
            if miss_path:
                if re.search(miss_path, file):
                    continue
            if extensions is not None:
                if all(not file.lower().endswith(ext) for ext in exts):
                    continue
            if target_path:
                if not re.search(target_path, file):
                    continue
            files.append(corner + '\\' + file)

    return files


def tempPath(name='tmp', ext='', temp_dir=__temp_dir):
    i = 0
    temp_folder = f'{temp_dir}\\{i}'
    while os.path.exists(temp_folder):
        i += 1
        temp_folder = f'{temp_dir}\\{i}'
    if not name:
        name = str(i)
    os.makedirs(temp_folder)
    return fullPath(temp_folder, name, ext)


# ----------------------------------------------------------------------------


sureDir(__temp_dir)


# Delete all files with last change over 1 day ago from __temp_dir
for corner, folders, names in os.walk(__temp_dir):

    for name in names:
        file = fullPath(corner, name)
        if (time.time()-os.path.getmtime(file))/86400 > 1:
            deleteFile(file)

    if len(os.listdir(corner))==0:
        if corner != __temp_dir:
            os.rmdir(corner)
