def check_file_for_newline(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if '\n' in line:
                return True
    return False

# 使用示例
file_path = 'D:/研/KDD16_HAM-master/KDD16_HAM-master/karatecrn.txt'
contains_newline = check_file_for_newline(file_path)

if contains_newline:
    print("文件包含 '\\n' 字符.")
else:
    print("文件不包含 '\\n' 字符.")

# def remove_newline_from_file(file_path):
#     with open(file_path, 'r') as file:
#         lines = file.readlines()
#
#     # 移除每一行中的 '\n' 字符
#     new_lines = [line.replace('\n', '') for line in lines]
#
#     with open(file_path, 'w') as file:
#         file.writelines(new_lines)
#
# # 使用示例
# file_path = 'D:/研/KDD16_HAM-master/KDD16_HAM-master/footballcrn.txt'
# remove_newline_from_file(file_path)