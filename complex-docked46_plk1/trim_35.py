import sys
import numpy as np
import re

input_filename = sys.argv[1]
output_filename = sys.argv[2]

def extract_text(input_filename, start_string, end_string):
    with open(input_filename, "r") as f:
        text = f.read()
    return re.findall(start_string + "(.*?)" + end_string, text, re.DOTALL)

if __name__ == "__main__":
    filename = input_filename
    start_string = "MODEL 4"
    end_string = "MODEL 5"
    text = extract_text(input_filename, start_string, end_string)
#with open(output_filename, "w") as f:
    #f.write(text)
    np.savetxt(output_filename, text, "%s")
#    print(text)
