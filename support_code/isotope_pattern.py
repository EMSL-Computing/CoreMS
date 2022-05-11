
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from support_code import AtomsDescription_standardized as atom

#from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
#from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas


#Create plot of (a) nmass error vs. ratio, (b) r-squared cutoff, (c) slope cutoff based on results. 
def isotopehunter_qc_plots(results,pattern,file_name,correlation):
    mzdif=[]
    r_dif=[]
    qc=[]

    elementpatterns=pd.DataFrame(atom.get_elementpattern())

    #Determine 2 most abundant isotopologues for correlation analysis.  
    isotope1=pattern.sort_values(by='ratio',ascending=False).isotope[0]
    isotope2=pattern.sort_values(by='ratio',ascending=False).isotope[1]

    for result in results[0]:
        mzdif.append(result[isotope2]['mz']-result[isotope1]['mz'])
        r_dif.append(result[isotope2]['intense']/result[isotope1]['intense'])
        qc.append(result['qc'])

    plot_df=pd.DataFrame({'mzdif':mzdif,'rdif':r_dif,'qc':qc})
    clean_results=results[1]
    final_results=results[2]

    fig, (ax1,ax2,ax3) = plt.subplots(3,1)
    fig.set_size_inches(8,11)

    allcolor='lightgray'
    validcolor='black'
    maxcolor='red'
    ax1.scatter('mzdif','rdif',color=allcolor,data=plot_df)
    ax1.scatter('mzdif','rdif',color=validcolor,data=plot_df[plot_df.qc=='qc'])
    ax1.scatter('mzdif','rdif',color=maxcolor,data=plot_df[plot_df.qc=='max'])
    ax1.scatter('mdiff','ratio',color='blue',data=elementpatterns)
    for i, txt in enumerate(elementpatterns.element):
        ax1.annotate(txt, (elementpatterns.mdiff.iloc[i], elementpatterns.ratio.iloc[i]),color='blue')
    ax1.set(xlabel='m/z difference (Da)',ylabel='Ratio',title=file_name + ': '+ isotope2 + ' / ' + isotope1)
    ax1.set(xlim=[min(mzdif),max(mzdif)], ylim=[min(r_dif),max(r_dif)])
    ax1.legend(['all','valid','max','true'],bbox_to_anchor=(1.05, 1.0), loc='upper left')

    ax2.hist('corr',bins=np.arange(0,1,0.05),data=pd.DataFrame(results[0]),color=allcolor)
    ax2.hist('corr',bins=np.arange(0,1,0.05),data=pd.DataFrame(clean_results),color=validcolor)
    ax2.hist('corr',bins=np.arange(0,1,0.05),data=pd.DataFrame(final_results),color=maxcolor)
    ax2.legend(['all','valid','max'],bbox_to_anchor=(1.05, 1.0), loc='upper left')

    ax2.set(xlabel='R-squared',ylabel='Frequency')
    ax2.axvline(x=correlation, color='gray', linestyle='dashed')

    ax3.hist('slope',bins=20,data=pd.DataFrame(clean_results),color=validcolor)
    ax3.hist('slope',bins=20,data=pd.DataFrame(final_results),color=maxcolor)
    ax3.set(xlabel='Normalized Slope',ylabel='Frequency')
    ax3.legend(['valid','max'],bbox_to_anchor=(1.05, 1.0), loc='upper left')

def metal_assignment(lcms,result,pattern,file_name):

    mass_spectrum=lcms.get_average_mass_spectrum_by_scanlist(result['scan'])

    # mass_spectrum.filter_by_max_resolving_power(15, 2)
    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

    ms_df=mass_spectrum.to_dataframe()

    ms_df_range=ms_df[(abs(ms_df['m/z']-result['mass'])<6)]
    mf=ms_df_range[ms_df_range['m/z'].round(3)==result['mass']]['Molecular Formula'].iloc[0]
    error=ms_df_range[ms_df_range['m/z'].round(3)==result['mass']]['m/z Error (ppm)'].iloc[0]
    score=ms_df_range[ms_df_range['m/z'].round(3)==result['mass']]['Confidence Score'].iloc[0]

    result['Molecular Formula']=mf
    result['Error']=error
    result['Score']=score
    result['file_name']=file_name


def metal_chromatogram(ax,lcms,result,pattern,timerange,file_name):
    
    #Determine most abundant isotopologues for correlation analysis.  
    isotope1=pattern.sort_values(by='ratio',ascending=False).isotope[0]
    
    plotmasses=pattern.mdiff+result[isotope1]['mz']-pattern[pattern.isotope==isotope1].mdiff.iloc[0]

    eics=lcms.get_eics(target_mzs=plotmasses.to_list(),tic_data={},peak_detection=False,smooth=False)

    for i in plotmasses.index:
        ax.plot(eics[0][plotmasses[i]].time, eics[0][plotmasses[i]].eic/pattern.ratio[i])
    ax.set(xlabel='Retention time (min)',ylabel='Scaled Intensity',xlim=timerange)
    ax.set_title(file_name)
    ax.legend(pattern.isotope.to_list(),bbox_to_anchor=(1.05, 1.0), loc='upper left',frameon=False)
    ax.text(1.05,0.7,'EIC:',transform=ax.transAxes)
    ax.text(1.05,0.6,plotmasses.round(4).to_list(),transform=ax.transAxes)
    ax.text(1.05,0.5,'R-Squared:'+ result['corr'].round(3).astype(str),transform=ax.transAxes)
    ax.text(1.05,0.4,'Slope:'+ result['slope'].round(2).astype(str),transform=ax.transAxes)
    ax.text(1.05,0.3,'d m/z:'+ result['dmz'].round(4).astype(str),transform=ax.transAxes)
    ax.axvline(x=result['time_peak'],color='gray',linestyle='dashed')

#This works on assigned spectra. 
def MS_pattern_plot(ax,lcms,result,pattern):

    mass_spectrum=lcms.get_average_mass_spectrum_by_scanlist(result['scan'])

    ms_df=mass_spectrum.to_dataframe()

    ms_df_range=ms_df[(abs(ms_df['m/z']-result['mass'])<6)]

    ax.stem('m/z','Peak Height',data=ms_df_range,  markerfmt=' ', basefmt=' ')
    
    theor_mz=pattern.mdiff+result['mass']
    theor_int=pattern.ratio*result['abundance']
    ax.stem(theor_mz,theor_int, basefmt=' ',linefmt='gray')

    for isotope in pattern.isotope[pattern.requirement=='Y']:
        ax.stem('mz','intense',data=result[isotope],  markerfmt=' ', basefmt=' ',linefmt='red')

    ax.legend(('all peaks','theoretical','pattern'),bbox_to_anchor=(1.05, 1.0), loc='upper left',frameon=False)
    ax.set(xlabel='m/z',ylabel='Intensity')
    ax.set_title('Time: '+result['time_peak'].round(2).astype(str)+' min')
    if('Score' in result.keys()):
        ax.text(1.05,0.7,result['Molecular Formula'],transform=ax.transAxes)
        ax.text(1.05,0.6,'Error (ppm) = '+round(result['Error'],3).astype(str),transform=ax.transAxes)
        ax.text(1.05,0.5,'Score = '+round(result['Score'],3).astype(str),transform=ax.transAxes)
    ax.axhline(y=0.0, color='black')

def apoplot(ax,lcms,result,pattern,timerange,file_name,metalcharge):
            
    #Determine most abundant isotopologues for correlation analysis.  
    isotope1=pattern.sort_values(by='ratio',ascending=False).isotope[0]
    
    apodiff=pattern.mass[0]-atom.atoms['1H'].atomic_mass*metalcharge
        
    plotmass=result[isotope1]['mz']-apodiff
    #print(plotmass)
    eics=lcms.get_eics(target_mzs=[plotmass],tic_data={},peak_detection=False,smooth=False)

    ax.plot(eics[0][plotmass].time, eics[0][plotmass].eic)
    ax.set(xlabel='Retention time (min)',ylabel='Intensity',xlim=timerange)
    ax.set_title(file_name)
    #ax.legend(pattern.isotope.to_list(),bbox_to_anchor=(1.05, 1.0), loc='upper left',frameon=False)
    ax.text(1.05,0.7,'Apo EIC:',transform=ax.transAxes)
    ax.text(1.05,0.6,[plotmass.round(4)],transform=ax.transAxes)
    #ax.text(1.05,0.5,'R-Squared:'+ result['corr'].round(3).astype(str),transform=ax.transAxes)
    #ax.text(1.05,0.4,'Slope:'+ result['slope'].round(2).astype(str),transform=ax.transAxes)
    #ax.text(1.05,0.3,'d m/z:'+ result['dmz'].round(4).astype(str),transform=ax.transAxes)
    #ax.axvline(x=result['time'],color='gray',linestyle='dashed')
