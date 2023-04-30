import sys
import pandas as pd

# Check if the required modules are installed
required_modules = ['data.table', 'dplyr', 'ggplot2', 'tidyr']
for module in required_modules:
    try:
        __import__(module)
    except ImportError:
        print(f"Module {module} not found")
        sys.exit(1)

input_file = None
interval = "all"
sample = "all"
original_file_name = None

# Parse command line arguments
for i in range(1, len(sys.argv), 2):
    if sys.argv[i] == "--input":
        input_file = sys.argv[i+1]
    elif sys.argv[i] == "--interval":
        interval = sys.argv[i+1]
    elif sys.argv[i] == "--sample":
        sample = sys.argv[i+1]
    elif sys.argv[i] == "--original_file_name":
        original_file_name = sys.argv[i+1]

# Read input file using pandas
if input_file:
    data = pd.read_csv(input_file, sep='\t', header=0)
else:
    print('VEP annotated input file not specified')
    sys.exit(1)

annotsv = data[~data.duplicated(subset='AnnotSV_ID')][['AnnotSV_ID', 'SV_type', 'GnomAD_pLI', 'B_loss_source', 'B_gain_source']]
annotsv['novel_annotsv'] = annotsv.apply(lambda row: True if pd.isna(row['GnomAD_pLI']) and row['B_loss_source']=='' and row['B_gain_source']=='' else False, axis=1)
annotsv = annotsv.groupby(['novel_annotsv', 'SV_type']).size().reset_index(name='Count')
annotsv = annotsv[['SV_type', 'Count', 'novel_annotsv']]
annotsv['interval'] = interval
annotsv['sample'] = sample
annotsv['original_file_name'] = original_file_name

# Write the output to a file
output_file = f"{interval}.{sample}.annotsv.counted"
annotsv.to_csv(output_file, sep='\t', index=False)
