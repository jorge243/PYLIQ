from PYLIQ import analyze
import os

input_file = 'PYLIQ Input.xlsx'
file_dir = os.path.dirname(__file__)

analyze(input_file, file_dir)
