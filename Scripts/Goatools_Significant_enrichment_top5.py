'''
@author: hayoung kim
'''

import pandas as pd
import os
import glob
type_list = ['ALL','CIR','COM','COP','COR','FB2U','FBU','FFL']

for type in type_list:
    module_list = os.listdir('/Users/hayoungkim/Desktop/Masters_SEM2/Design_Project/goatools-main/Atha_goatools_outputs/SCHype'+type)
    for module in module_list:
        if '.txt' in module:
            data = pd.read_csv('/Users/hayoungkim/Desktop/Masters_SEM2/Design_Project/goatools-main/Atha_goatools_outputs/SCHype'+type+'/'+module, sep="\t", header=0)
            data['type_module'] = type+'_module'+module[:-4]
            print(data.shape[0])
            new_data = data[data['enrichment'].str.contains('e', na=False)]
            if new_data.shape[0] >= 5:
                top5_data = new_data.nlargest(5, 'p_fdr_bh', keep='all')
            else:
                top5_data = new_data
            top5_data.to_csv('/Users/hayoungkim/Desktop/Masters_SEM2/Design_Project/goatools-main/Atha_goatools_outputs/GO_csv/' + 'Atha_' + 'SCHype' + type + '_' + module[:-4] +'_top5.csv', index = True)
os.chdir('/Users/hayoungkim/Desktop/Masters_SEM2/Design_Project/goatools-main/Atha_goatools_outputs/GO_csv')
extension = 'csv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

#combine all files in the list
combined_csv = pd.concat([pd.read_csv(f, index_col = 'type_module', usecols = ["type_module","# GO","NS","enrichment","name","ratio_in_study","ratio_in_pop","p_uncorrected","depth","study_count","p_fdr_bh","study_items"]) for f in all_filenames])
#export to csv
combined_csv.to_csv( "/Users/hayoungkim/Desktop/Masters_SEM2/Design_Project/goatools-main/Atha_goatools_outputs/GOATOOLS_summary.csv",  encoding='utf-8-sig')

