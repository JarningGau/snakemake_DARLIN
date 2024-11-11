import os
import fnmatch
import shutil

##########
# you can run this script after you downloaded new data:
# 
# the folder directory structure should  be as follows:
# 20240101_DARLIN_Mon_Gr_XXX
#     |
#     |___ Data
#         |___Gr_A1
#             |___Gr_A1_R1.fq.gz
#             |___Gr_A1_R2.fq.gz
# or
# 20240101_DARLIN_Mon_Gr_XXX
#     |
#     |___ raw_fastq
#         |___Gr_A1
#             |___Gr_A1_R1.fq.gz
#             |___Gr_A1_R2.fq.gz

# run rename script in the folder 20240101_DARLIN_Mon_Gr_XXX:
# cd 20240101_DARLIN_Mon_Gr_XXX then run this rename script
##########


def check_cwd(cwd):
    items = os.listdir(cwd)

    if not 'Data' in items or 'raw_fastq' in items:
        new_folder = os.path.join(cwd, 'raw_fastq')
        os.makedirs(new_folder, exist_ok=True)
        files_items = [item for item in items if ( item.endswith('.fq.gz') or item.endswith('.fastq.gz') and 'md5' not in item and 'MD5' not in item)]
        for item in files_items:
            if os.path.isfile(os.path.join(cwd, item)):
                shutil.move(os.path.join(cwd, item), os.path.join(new_folder, item))

# get files
def find_files(directory, pattern):
    matched_files = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, pattern):
            matched_files.append(os.path.join(root, filename))
    matched_files = [item for item in matched_files if ('checkpoint' not in item and 'md5' not in item and 'MD5' not in item)]
    return matched_files

def process_format(pattern, remove_str, replace_str, process_folders):
    files = find_files(cwd, pattern)

    if len(files):
        for item in files:
            file_dir = os.path.dirname(item)
            file_basename = os.path.basename(item)
            filter_name = file_basename.replace(remove_str, replace_str)
            shutil.move(item, os.path.join(file_dir,filter_name))

        if process_folders:
            folders = set(os.path.dirname(item) for item in files)
            for item in folders:
                filter_name = item.replace('-', '_')
                shutil.move(item, filter_name)

    # rename folder: Data -> raw_astq
    if os.path.exists(os.path.join(cwd, 'Data')):
        shutil.move(os.path.join(cwd, 'Data'), os.path.join(cwd, 'raw_fastq'))

def trim_path(full_path, target_dir='raw_fastq'):
    parts = full_path.split(os.sep)
    if target_dir in parts:
        index = parts.index(target_dir)
        trimmed_path = os.sep.join(parts[:index + 2])
        return trimmed_path
    else:
        return full_path

def folders_pairing(cwd):
    raw_fastq_path = cwd + '/raw_fastq'
    files = find_files(raw_fastq_path, '*.fastq.gz')

    folders = files
    folders = [item.replace('_R1.fastq.gz', '') for item in folders]
    folders = set([item.replace('_R2.fastq.gz', '') for item in folders])
    # print(raw_fastq_path)
    for item in folders:
        # print('item',item)
        trimmed_path = trim_path(item, target_dir='raw_fastq')
        # target_folder = os.path.join(raw_fastq_path, item)
        if not os.path.exists(trimmed_path):
            
            # print(trimmed_path)
            os.makedirs(trimmed_path)
            sources_file_R1 = item +'_R1.fastq.gz'
            sources_file_R2 = item +'_R2.fastq.gz'
            # print(sources_file_R2)
            shutil.move(sources_file_R1, trimmed_path)
            shutil.move(sources_file_R2, trimmed_path)


def main(cwd):
    patterns = ['*.fq.gz', '*.fastq.gz', '*.fastq.gz', '*.fastq.gz', '*.fastq.gz']
    remove_strs = ['fq.gz', '-R1', '-R2', '_L001', '_001']
    replace_strs = ['fastq.gz', '_R1', '_R2', '', '']


    for remove_str in remove_strs:
        process_folders = False
        index = remove_strs.index(remove_str)
        pattern = patterns[index]
        replace_str = replace_strs[index]
        if index == len(patterns)-1:
            process_folders = True
        process_format(pattern, remove_str, replace_str, process_folders)


if __name__ == '__main__':
    cwd = os.getcwd()
    check_cwd(cwd)
    main(cwd)
    folders_pairing(cwd)