from collections import defaultdict
import regex
import sys
import gzip
import os
import json
BLUE, RED, WHITE, YELLOW, MAGENTA, GREEN, END = '\33[94m', '\033[91m', \
                                                '\33[97m', '\33[93m', \
                                                '\033[1;35m', '\033[1;32m', \
                                                '\033[0m'

def clear():
    os.system('clear')

def heading(direction='R1'):
    clear()

    sys.stdout.write(GREEN + f'''
                                   8
                        .,,aadd88P=8=Y88bbaa,,.
                  .,ad88888P:a8P:d888b:Y8a:Y88888ba,.
              ,ad888888P:a8888:a8888888a:8888a:Y888888ba,
           ,a8888888:d8888888:d888888888b:8888888b:8888888a,
        ,a88888888:d88888888:d88888888888b:88888888b:88888888a,
      ,d88888888:d888888888:d8888888888888b:888888888b:88888888b,
    ,d88888888:d8888888888I:888888888888888:I8888888888b:88888888b,
  ,d888888888:d88888888888:88888888888888888:88888888888b:888888888b,
 d8888888888:I888888888888:88888888888888888:888888888888I:8888888888b
d8P"'   `"Y8:8P"'     `"Y8:8P"'    8    `"Y8:8P"'     `"Y8:8P"'   `"Y8b
"           "             "        8        "             "           "
                                   8
                                   8  ''' + RED + '''    ![Virus analysis]! v1.0''' + GREEN + '''
                                   8
        [''' + WHITE + 'C' + GREEN + '''] Check virus sequence   8     ''' + WHITE + 'by:' + MAGENTA + ''' Chanju Jung (''' + YELLOW + 'Jung handsome' + GREEN + ''')
        [''' + WHITE + 'R' + GREEN + '''] Set R direction:''' + WHITE + f''' {direction}''' + GREEN + '''    8           BAELAB
        [''' + WHITE + 'F' + GREEN + '''] Set input file folder  8             
        [''' + WHITE + 'S' + GREEN + '''] Set save folder        8
        [''' + WHITE + 'P' + GREEN + '''] Parameters check       8
        [''' + WHITE + 'A' + GREEN + '''] Execute program        8
        [''' + WHITE + 'Q' + GREEN + '''] Quit                   8
\n''' + END)
    print("Select an option from menu:\n")
    
########

def virus_classify(virus_seq_dict, seq):
    """"
    This function classifies the virus sequence in the input sequence
    output: virus name, if no virus is found, return None
    """
    for virus in virus_seq_dict:
        virus_match = regex.search(rf"({virus_seq_dict[virus].upper()}){{s<=1}}", seq)
        if virus_match:
            return virus
    return None

def reverse_complement(seq):
    """
    This function returns the reverse complement of the input sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq = seq.upper()
    return "".join(complement[base] for base in reversed(seq))






def check_virus_sequence():
    try:
        with open('virus_seq.json', 'r') as file:
            virus_data = json.load(file)
            for virus, sequence in virus_data.items():
                print(f"{virus}: {sequence}")
    except FileNotFoundError:
        print("virus_seq.json file not found.")
    except json.JSONDecodeError:
        print("Error reading the virus_seq.json file.")

def set_file_directory(direction):
    date = input("Enter the NGS date (YYYYMMDD): ")
    base_directory = "/home/baelab/MGEL/"
    directory = base_directory + date

    if not os.path.exists(directory):
        print(f"{directory} does not exists. Check the NGS date.")
        return
    index_range = input("Enter the index range (e.g., 25-30): ")
    if not index_range:
        alert_val_error()
        return

    input_list = []
    for i in index_range.split(","):
        if "-" in i:
            input_list.extend(range(int(i.split("-")[0]), int(i.split("-")[1]) + 1))
        else:
            input_list.append(int(i))
            

    files = [f"{directory}/{i}_S{i}_L001_R{direction[-1]}_001.fastq.gz" for i in input_list]
    return files


def set_save_folder():
    save_folder = input(f"Enter the save folder under current directory: Default is {os.path.join(os.getcwd(), 'virus_analysis')}\n \
    Press enter to use default.")
    if not save_folder:
        save_folder = os.path.join(os.getcwd(), 'virus_analysis')
    else:
        save_folder = os.path.join(os.getcwd(), save_folder)
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
        print(f"Directory {save_folder} created.")
    print(f"Save folder set to: {save_folder}")
    return save_folder


def execute_program(files = None, save_folder = None, direction = 'R1', qc_length = 20, reference_seq = None, virus_seq_dict = None, virus_count = None, virus_data=None):
    '''
    This function executes the program
    It returns the virus count dictionary and number of lines failed QC
    '''
    result_dict = defaultdict(list)

    if not files:
        print("Input file directory not set.")
        return
    if not save_folder:
        print("Save folder not set.")
        return

    for file in files:
        if direction == 'R1':
            virus_count = {key: 0 for key in virus_data["R1_virus_seq_dict"]}
        if direction == 'R2':
            virus_count = {key: 0 for key in virus_data["R2_virus_seq_dict"]}
        virus_count["No_virus"] = 0
        lines = read_fastq(file)
        result_list = count_sequence(lines, direction=direction, virus_seq_dict=virus_seq_dict, virus_count=virus_count, reference_seq=reference_seq)  # [sequence_count, qc_failed]
        result_dict[file] = result_list
    return result_dict
    
def save_as_txt(save_folder, result_dict):
    if not save_folder:
        print("Save folder not set.")
        return
    if not result_dict:
        print("No result to save.")
        return

    print(result_dict.keys())

    file_name = "virus_analysis.txt"
    with open(os.path.join(save_folder, file_name), 'a') as f:
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(\
            "WT reads (no virus)",
            "Delta_only",
            "Omicron_only",
            "Delta-Delta",
            "Delta-Omicron",
            "Omicron-Omicron",
            "Omicron-Delta",
            "NGS_date",
            "Index_num"
            ))
        f.write("\n")
        for file, result_list in result_dict.items():
            class_count = {virus: count for virus, count in result_list[2].items()}
            f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format( 
                class_count['No_virus'],
                class_count['Delta_only'],
                class_count['Omicron_only'],
                class_count['Delta-Delta'],
                class_count['Delta-Omicron'],
                class_count['Omicron-Omicron'],
                class_count['Omicron-Delta'],
                file.split('/')[-2],
                file.split('/')[-1].split('_')[0]
                    ))
            f.write("\n")
    print("Files saved successfully.")
    

def read_fastq(file_path):
    """
    This function reads the fastq file and returns the sequence count
    """
    # read the fastq file
    with gzip.open(file_path, 'rt') as f:
        lines = f.readlines()
        f.close()
    return lines


def count_sequence(lines, direction='R1', reference_seq=None, virus_seq_dict=None, qc_length=20, virus_count=None):
    """
    This function counts the sequence in the fastq file
    """

    sequence_count = defaultdict(int)
    if direction == 'R1':
        indicator = reference_seq[20:20+qc_length]
    elif direction == 'R2':
        indicator = reverse_complement(reference_seq[-20-qc_length:-20])

    for i in range(0, len(lines), 4):
        seq = lines[i+1].strip()
        sequence_count[seq] += 1

    # QC the sequence using the indicator

    qc_failed = 0


    for seq in sequence_count:
        if indicator not in seq:
            qc_failed += sequence_count[seq]
            continue
        virus = virus_classify(virus_seq_dict, seq)
        if virus:
            virus_count[virus] += sequence_count[seq]
        else:
            virus_count["No_virus"] += sequence_count[seq]
    return [sequence_count, qc_failed, virus_count]

heading()


# print(f"Indicator: {indicator}")
# print(f"total reads: {sum(sequence_count.values())}")
# print(f"lines with no QC indicator: {qc_failed}")
# print(virus_count)


def alert_val_error():
    print(RED + "Invalid value. Please enter a valid value.")



def main():
    direction = 'R1'
    qc_length = 20
    save_folder = None
    files = None


    files = [] 
    clear()
    heading(direction=direction)
    try:
        while True:
            choice = input(BLUE + "Enter your choice: ").upper()
            if choice == 'C':
                clear()
                heading(direction=direction)
                check_virus_sequence()

            elif choice == 'R':
                clear()
                heading(direction=direction)
                print("Enter the direction (R1 or R2): ")
                choice = input(GREEN + '''[''' + WHITE + '1' + GREEN + '''] R1 ''' +"\n"
                                + GREEN + '''[''' + WHITE + '2' + GREEN + '''] R2 ''' + "\n") .upper()
                if choice == '1':
                    direction = 'R1'
                elif choice == '2':
                    direction = 'R2'
                else:
                    alert_val_error()
                heading(direction=direction)

            elif choice == 'F':
                clear()
                heading(direction=direction)
                files = set_file_directory(direction=direction)

            elif choice == 'S':
                clear()
                heading(direction=direction)
                save_folder = set_save_folder()
            
            elif choice == 'P':
                clear()
                heading(direction=direction)
                print(f"Direction: {direction}")
                print(f"Input files: {files}")
                print(f"Save folder: {save_folder}")
            
            elif choice == 'A':
                clear()
                heading(direction=direction)
                virus_data = json.load(open('virus_seq.json', 'r'))
                if direction == 'R1':
                    virus_count = {key: 0 for key in virus_data["R1_virus_seq_dict"]}
                    reference_seq = virus_data["R1_REFERENCE_SEQ"].upper()
                    virus_seq_dict = virus_data["R1_virus_seq_dict"]
                elif direction == 'R2':
                    virus_count = {key: 0 for key in virus_data["R2_virus_seq_dict"]}
                    reference_seq = virus_data["R2_REFERENCE_SEQ"].upper()
                    virus_seq_dict = virus_data["R2_virus_seq_dict"]
                virus_count["No_virus"] = 0
                print("Executing the program.. ")


                

                result_dict = execute_program(files=files, save_folder=save_folder, direction=direction, reference_seq=reference_seq, virus_seq_dict=virus_seq_dict, qc_length=qc_length, virus_count=virus_count, virus_data=virus_data)
                save_as_txt(save_folder, result_dict)
                print("Program executed successfully.")
                sys.exit(0)

            elif choice == "Q":
                clear()
                heading(direction=direction)
                print("Exiting the program")
                sys.exit(0)
            else:
                alert_val_error()

    except KeyboardInterrupt:
        print("Exiting the program")
        sys.exit(0)
                
            
if __name__ == '__main__':
    main()
