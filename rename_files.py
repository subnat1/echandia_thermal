import os

def rename_files(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith('.jpg'):
            old_name = filename
            # first part of the filename is either Cell Body, Cell Center-section
            split1 = filename.split("_")
            # print(split1[1])
            # print("First split done")
            # Split the second part at "." and take the first part and if it is a number convert that into a 5 digit no str.
            split2 = split1[1].split(".")
            # print("second split done")
            if split2[0].isdigit():
                five_digt = f"{split2[0].zfill(5)}"
                # print(five_digt)
                new_name = f"{split1[0]}_{five_digt}.{split2[2]}"
                # print(new_name)
                curr_path = os.getcwd()
                os.rename(f"{folder_path}/{old_name}", f"{folder_path}/{new_name}")

            else:
                continue

            # Club the first part from 1, the five digit no. and the remaining part of the filename and create this into a new name
            # Rename the oldname into newname
            

# Usage example:
folder_path = "output_fdm_laplace_20241027_141636/jpegs"
subfolders  = sorted(next(os.walk(folder_path))[1])
for x in subfolders:
    print(f"Renaming {x}")
    try:
        rename_files(f"{folder_path}/{x}")
    except:
        continue

# rename_files_1('output_fdm_laplace_20241027_141636/jpegs/6.5C')


