import pandas as pd
import os
import glob
import re
from pathlib import Path
from openpyxl import load_workbook
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

#Define the directory containing the files to reorganize
os.chdir('/Users/eber373/OneDrive - PNNL/Documents/Data/Test Data Aggregation/')
cwd=os.getcwd()

# Glob the files into a list
FileList = glob.glob('*.csv')

def iter_dir(FileList, out_dir=cwd):
    ''' Iterate over the files in a directory to reorder them into a new output folder'''
    
    #create out directory
    directory = "reorg_output"
    path = os.path.join(out_dir, directory)
    print(path)
    if Path(path).exists()==False:
        os.mkdir(path)
        print("Directory '% s' created" % directory)

    print("Status: Organizing")
    print(FileList)
    for file in FileList:
        #get name of working file
        match = re.match("^[^_]*_[^_]*_", file)
        prefix = match.group(0)

        #get name of existing files
        outmatch = prefix+"reorg.xlsx"
        outpath = os.path.join(out_dir, directory , outmatch)
        print(outpath)
        #If the outfile doesnt exist create a new file to populate
        if Path(outpath).exists()==False:
            
            #import current file into pandas df
            whole_data = pd.read_csv(file)

            #drop NaN from working df
            clean_data = whole_data.dropna(how='any', subset=['Calculated m/z', 'm/z Error (ppm)'])

            # data frame for each sheet
            mz = clean_data[['Molecular Formula','Ion Type',
                                'm/z']]
            Calibrated_mz = clean_data[['Molecular Formula','Ion Type',
                                'Calibrated m/z']]	
            Calculated_mz = clean_data[['Molecular Formula','Ion Type',
                                'Calculated m/z']]	
            PeakHeight = clean_data[['Molecular Formula','Ion Type',
                                'Peak Height']]	
            ResolvingPower = clean_data[['Molecular Formula','Ion Type',
                                'Resolving Power']]	
            SN = clean_data[['Molecular Formula','Ion Type',
                                'S/N']]	
            mzErrorPPM = clean_data[['Molecular Formula','Ion Type',
                                'm/z Error (ppm)']]
            mzErrorScore = clean_data[['Molecular Formula','Ion Type',
                                'm/z Error Score']]	
            IsotopologueSimilarity = clean_data[['Molecular Formula','Ion Type',
                                'Isotopologue Similarity']]	
            ConfidenceScore = clean_data[['Molecular Formula','Ion Type',
                                'Confidence Score']]

            # supress SettingWithCopyWarning
            pd.options.mode.chained_assignment = None  # default='warn'

            #rename the field of interest to include the sample it comes from
            mz.rename(columns={'m/z': file[:-3]+'m/z'}, inplace=True)
            Calibrated_mz.rename(columns={'Calibrated m/z': file[:-3]+'Calibrated m/z'}, inplace=True)
            Calculated_mz.rename(columns={'Calculated m/z': file[:-3]+'Calculated m/z'}, inplace=True)
            PeakHeight.rename(columns={'Peak Height': file[:-3]+'Peak Height'}, inplace=True)
            ResolvingPower.rename(columns={'Resolving Power': file[:-3]+'Resolving Power'}, inplace=True)
            SN.rename(columns={'S/N': file[:-3]+'S/N'}, inplace=True)
            mzErrorPPM.rename(columns={'m/z Error (ppm)': file[:-3]+'m/z Error (ppm)'}, inplace=True)
            mzErrorScore.rename(columns={'m/z Error Score': file[:-3]+'m/z Error Score'}, inplace=True)
            IsotopologueSimilarity.rename(columns={'Isotopologue Similarity': file[:-3]+'Isotopologue Similarity'}, inplace=True)
            ConfidenceScore.rename(columns={'Confidence Score': file[:-3]+'Confidence Score'}, inplace=True)

            # Create a Pandas Excel writer using XlsxWriter as the engine.
            writer = pd.ExcelWriter(path+'/'+prefix+'reorg.xlsx', mode='w', engine='openpyxl')

            # Write each dataframe to a different worksheet.
            clean_data.to_excel(writer, sheet_name='All', index=False)
            mz.to_excel(writer,sheet_name='mz', index=False)
            Calibrated_mz.to_excel(writer,sheet_name='Calibrated mz', index=False)	
            Calculated_mz.to_excel(writer,sheet_name='Calculated mz', index=False)	
            PeakHeight.to_excel(writer,sheet_name='Peak Height', index=False)	
            ResolvingPower.to_excel(writer,sheet_name='Resolving Power', index=False)	
            SN.to_excel(writer,sheet_name='SN', index=False)	
            mzErrorPPM.to_excel(writer,sheet_name='mz Error ppm', index=False)
            mzErrorScore.to_excel(writer,sheet_name='mz Error Score', index=False)	
            IsotopologueSimilarity.to_excel(writer,sheet_name='Isotopologue Similarity', index=False)
            ConfidenceScore.to_excel(writer,sheet_name='Confidence Score', index=False)

            # Close the Pandas Excel writer and output the Excel file.
            writer.save()

        #If the outfile exists and the prefix matches go ahead and append to that file
        else:
            #open new data file as pandas df
            new_data = pd.read_csv(file, index_col=0)

            #drop NaN from newly imported data
            clean_data = new_data.dropna(how='any', subset=['Calculated m/z', 'm/z Error (ppm)'])

            #open existing excel file as a dict with sheets to append to as pandas dataframes
            sheetAll = pd.read_excel(outpath, sheet_name="All", index_col=0)
            sheetMZ = pd.read_excel(outpath, sheet_name="mz",index_col=0)
            sheetCalibrated_mz = pd.read_excel(outpath, sheet_name="Calibrated mz",index_col=0)
            sheetCalculated_mz = pd.read_excel(outpath, sheet_name="Calculated mz",index_col=0)
            sheetPeakHeight = pd.read_excel(outpath, sheet_name="Peak Height",index_col=0)
            sheetResolvingPower = pd.read_excel(outpath, sheet_name="Resolving Power",index_col=0)
            sheetSN = pd.read_excel(outpath, sheet_name="SN",index_col=0)
            sheetmzErrorPPM = pd.read_excel(outpath, sheet_name="mz Error ppm",index_col=0)
            sheetmzErrorScore = pd.read_excel(outpath, sheet_name="mz Error Score",index_col=0)
            sheetIsotopologueSimilarity = pd.read_excel(outpath, sheet_name="Isotopologue Similarity",index_col=0)
            sheetConfidenceScore = pd.read_excel(outpath, sheet_name="Confidence Score",index_col=0)

            #update All
            updateAll = pd.concat((sheetAll,clean_data), sort=False)
            # update MZ
            updateMZ = pd.merge(left=sheetMZ, right=clean_data[['Molecular Formula','Ion Type','m/z']], how="outer", on=["Molecular Formula",'Ion Type'])
            updateMZ.rename(columns={'m/z': file[:-3]+'m/z'}, inplace=True)
            #update Calibrated m/z
            updateCalibrated_mz = pd.merge(left=sheetCalibrated_mz, right=clean_data[['Molecular Formula','Ion Type','Calibrated m/z']], how="outer", on=["Molecular Formula",'Ion Type'])
            updateCalibrated_mz.rename(columns={'Calibrated m/z': file[:-3]+'Calibrated m/z'}, inplace=True) 
            #update Calculated m/z    
            updateCalculated_mz = pd.merge(left=sheetCalculated_mz, right=clean_data[['Molecular Formula','Ion Type','Calculated m/z']], how="outer", on=["Molecular Formula",'Ion Type'])
            updateCalculated_mz.rename(columns={'Calculated m/z': file[:-3]+'Calculated m/z'}, inplace=True)
            #update Peak Height  
            updatePeakHeight = pd.merge(left=sheetPeakHeight, right=clean_data[['Molecular Formula','Ion Type','Peak Height']], how="outer", on=["Molecular Formula",'Ion Type'])
            updatePeakHeight.rename(columns={'Peak Height': file[:-3]+'Peak Height'}, inplace=True)
            #update Resolving Power         
            updateResolvingPower = pd.merge(left=sheetResolvingPower, right=clean_data[['Molecular Formula','Ion Type','Resolving Power']], how="outer", on=["Molecular Formula",'Ion Type'])
            updateResolvingPower.rename(columns={'Resolving Power': file[:-3]+'Resolving Power'}, inplace=True)
            #update S/N  
            updateSN = pd.merge(left=sheetSN, right=clean_data[['Molecular Formula','Ion Type','S/N']], how="outer", on=["Molecular Formula",'Ion Type'])
            updateSN.rename(columns={'S/N': file[:-3]+'S/N'}, inplace=True)
            #update m/z Error (ppm)  
            updatemzErrorPPM = pd.merge(left=sheetmzErrorPPM, right=clean_data[['Molecular Formula','Ion Type','m/z Error (ppm)']], how="outer", on=["Molecular Formula",'Ion Type'])
            updatemzErrorPPM.rename(columns={'m/z Error (ppm)': file[:-3]+'m/z Error (ppm)'}, inplace=True)
            #update m/z Error Score
            updatemzErrorScore = pd.merge(left=sheetmzErrorScore, right=clean_data[['Molecular Formula','Ion Type','m/z Error Score']], how="outer", on=["Molecular Formula",'Ion Type'])
            updatemzErrorScore.rename(columns={'m/z Error Score': file[:-3]+'/z Error Score'}, inplace=True)
            #update Isotopologue Similarity
            updateIsotopologueSimilarity = pd.merge(left=sheetIsotopologueSimilarity, right=clean_data[['Molecular Formula','Ion Type','Isotopologue Similarity']], how="outer", on=["Molecular Formula",'Ion Type'])
            updateIsotopologueSimilarity.rename(columns={'Isotopologue Similarity': file[:-3]+'Isotopologue Similarity'}, inplace=True)
            #update Confidence Score
            updateConfidenceScore = pd.merge(left=sheetConfidenceScore, right=clean_data[['Molecular Formula','Ion Type','Confidence Score']], how="outer", on=["Molecular Formula",'Ion Type'])
            updateConfidenceScore.rename(columns={'Confidence Score': file[:-3]+'Confidence Score'}, inplace=True)

            #read the existing sheets so that openpyxl won't create a new one later
            #write out
            with pd.ExcelWriter(outpath,engine = "openpyxl",  mode='a') as writer:
                workBook = writer.book
                try:
                    workBook.remove(workBook['All'])
                    workBook.remove(workBook['mz'])
                    workBook.remove(workBook['Calibrated mz'])
                    workBook.remove(workBook['Calculated mz'])
                    workBook.remove(workBook['Peak Height'])
                    workBook.remove(workBook['Resolving Power'])
                    workBook.remove(workBook['SN'])
                    workBook.remove(workBook['mz Error ppm'])
                    workBook.remove(workBook['mz Error Score'])
                    workBook.remove(workBook['Isotopologue Similarity'])
                    workBook.remove(workBook['Confidence Score'])
                except:
                    print("worksheet doeS_N't exist")
                finally:
                    updateAll.to_excel(writer, sheet_name='All')
                    updateMZ.to_excel(writer, sheet_name='mz')
                    updateCalibrated_mz.to_excel(writer, sheet_name='Calibrated mz')
                    updateCalculated_mz.to_excel(writer, sheet_name='Calculated mz')
                    updatePeakHeight.to_excel(writer, sheet_name='Peak Height')
                    updateResolvingPower.to_excel(writer, sheet_name='Resolving Power')
                    updateSN.to_excel(writer, sheet_name='SN')
                    updatemzErrorPPM.to_excel(writer, sheet_name='mz Error ppm')
                    updatemzErrorScore.to_excel(writer, sheet_name='mz Error Score')
                    updateIsotopologueSimilarity.to_excel(writer, sheet_name='Isotopologue Similarity')
                    updateConfidenceScore.to_excel(writer, sheet_name='Confidence Score')
                writer.save()

    return (print("All done!"))

iter_dir(FileList)
