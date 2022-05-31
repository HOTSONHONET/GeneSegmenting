from pprint import pprint
import os
import argparse
from Bio import SeqIO
from glob import glob
import json
from datetime import datetime


class Shatz:
    def readGBfile(self, file_name: str):
        to_strings = []
        for rec in SeqIO.parse(file_name, "gb"):
            if rec.features:
                for feature in rec.features:
                    if feature.type == "CDS":
                        # print(feature.location)
                        # print(type(feature.location)
                        to_strings.append(feature.location.__str__())

        list_of_pairs = []
        for line in to_strings:
            if line.startswith("join"):
                comp1, comp2 = line.split(",")
                # splitting the comp1
                comp1_temp1, comp1_temp2 = comp1.split(":")
                comp1_start_value = int(comp1_temp1[1+comp1_temp1.find("["):])
                compl_end_value = int(comp1_temp2[:comp1_temp2.find("]")])
                # Splitting the comp2
                comp2_temp1, comp2_temp2 = comp2.split(":")
                comp2_start_value = int(comp2_temp1[2:])
                comp2_end_value = int(comp2_temp2[:comp2_temp2.find("]")])

                pair = (comp2_start_value, comp2_end_value)
                list_of_pairs.append(pair)
            else:
                pair = tuple(map(int, line.split("]")[
                             0].split("[")[1].split(":")))
                list_of_pairs.append(pair)

        # pprint(list_of_pairs)
        list_of_pairs = sorted(list_of_pairs)
        s, l = list_of_pairs[0]
        mergeList = []
        for i in range(1, len(list_of_pairs)):
            tmp_s, tmp_e = list_of_pairs[i]
            if l >= tmp_s:
                l = tmp_e
            else:
                mergeList.append((s, l))
                s = tmp_s
                l = tmp_e
        mergeList.append((s, l))

        # pprint(mergeList)
        return mergeList

    def readSeq(self, file_name: str, sequences):
        main = ''
        with open(file_name, 'r') as reader:
            line = reader.readline()
            while line != "":
                line = reader.readline()
                main = main + line[:-1]
            reader.close()

        sequence_len = len(sequences)

        matrix = [0]*(sequence_len)
        # print("Enter your choice: 0 = include the sequence between indices, 1 = exclude the sequence between indices")
        # choice = int(input())
        choice = 1
        if (choice == 0):
            for i in range(sequence_len):
                # adding subsequences from sequence to rows in matrix
                matrix[i] = list(
                    main[int(sequences[i][0])-1:int(sequences[i][1])])
            print("\n\nThe Sequences are:")
            str1 = ""
            for i in range(sequence_len):
                print(str1.join((matrix[i])))
        elif (choice == 1):
            result = main
            for i in range(sequence_len-1, -1, -1):
                result = result[0:int(sequences[i][0])-1:] + \
                    "  "+result[int(sequences[i][1])::]
            # print("\n\nThe sequence is:")
            # print(result.split("  "))
            f_res = result.split("  ")
            to_return = {}
            for idx, seq in enumerate(f_res):
                to_return[f"seq_{idx + 1}"] = seq
            return to_return

    def __call__(self, gb_file_dir_path, seq_file_dir_path):
        gb_files = sorted(glob(f"{gb_file_dir_path}/*gb"))
        seq_files = sorted(glob(f"{seq_file_dir_path}/*txt"))
        json_op = {}
        for idx, (gb_file_path, seq_file_path) in enumerate(list(zip(gb_files, seq_files))):
            sequences = self.readGBfile(gb_file_path)
            dict_val = self.readSeq(seq_file_path, sequences)
            json_op[f"file_{idx + 1}"] = dict_val
        return json_op


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gb_dir_path", type=str, required=True,
                        help="path to the directory containing the .gb formatted files")
    parser.add_argument("--txt_dir_path", type=str, required=True,
                        help="path to the directory containing the .txt formatted files")
    parser.add_argument("--name", type=str,
                        default=f"shatz_results")

    args = parser.parse_args()

    # Running the Class
    shatabdi = Shatz()
    json_object = shatabdi(args.gb_dir_path, args.txt_dir_path)
    if not os.path.exists("./results"):
        os.mkdir("./results")

    date = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")
    with open(f"./results/{args.name}_{date}.json", "w") as jfile:
        json.dump(json_object, jfile, indent=4)

    print("[INFO] Successfully completed. You can check out the result in ./results directory...")


if __name__ == "__main__":
    main()
