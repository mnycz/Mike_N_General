#!/usr/bin/env python
# coding: utf-8

# exec(open('Waveform_Ring_Analysis.py').read())


import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gc
import matplotlib as mt
import timeit
from datetime import datetime
import time
import cProfile
from scipy.signal import chirp, find_peaks, peak_widths, argrelextrema
from scipy.ndimage import label, generate_binary_structure, iterate_structure,center_of_mass
#from scipy.integrate import simps
#from numpy import trapz
#import matplotlib.gridspec as gridspec
#import matplotlib.ticker as ticker

print('Hello')
startTime = datetime.now()
print(startTime)

pd.set_option('mode.chained_assignment',None)
#file = uproot.open("/lustre19/expphy/volatile/halla/triton/mnycz/SOLID/Cerenkov_160.root")
#file = uproot.open("/lustre19/expphy/volatile/halla/triton/mnycz/SOLID/Cerenkov_280.root")
file = uproot.open("/lustre19/expphy/volatile/halla/triton/mnycz/SOLID/Cerenkov_302.root")

EvTree = file['EvTree']

#%%time
#branches = EvTree.arrays(namedecode='utf-8')


# # Just read in portion of the rootfile 


#get_ipython().run_cell_magic('time', '', "branches = EvTree.arrays(namedecode='utf-8',entrystop=20000)\n#branches = EvTree.arrays(namedecode='utf-8',entrystop=50000)")
branches = EvTree.arrays(namedecode='utf-8',entrystop=100000)



Branch_Names=list(branches.keys())



sorted_Names = sorted(Branch_Names)



def Branch_To_Array(x,Generic_DF):
    Generic_DF.append(branches[x])

sorted_Names[40]


Undefined_Entry = -1


# In[10]:


Cer_String = ['Cer']
Calo_String = ['C1','C2','C3','C4','C5','C6','C7','C8','C9']
S2_String = ['S2']
#S2_String = ['SKIP']



def Cer_Calo_S2_Split(x,DF_S,DF_Det,DF_Ele):
    DF_S=(x.split("_"))
    #print(len(DF_S))
    #if (len(DF_S)==3 and DF_S[0]!=('C1'or 'C2'or 'C3' or 'C4')  ):
    if (len(DF_S)==3 and (DF_S[0][0:3] in Cer_String)  ):
        DF_Det.append(DF_S[0])
        DF_Ele.append(DF_S[2])
    elif (len(DF_S)==4 and (DF_S[0][0:3] in Cer_String)  ):
        DF_Det.append(DF_S[0])
        DF_Ele.append(DF_S[2] + DF_S[3])
   
    elif (len(DF_S)==3 and (DF_S[0][0:3] in S2_String)  ):
        DF_Det.append(DF_S[0])
        DF_Ele.append(DF_S[2])
    elif (len(DF_S)==4 and (DF_S[0][0:3] in S2_String)  ):
        DF_Det.append(DF_S[0])
        DF_Ele.append(DF_S[2] + DF_S[3])

    elif ((DF_S[0] in Calo_String) and len(DF_S)==2 ):
        DF_Det.append(DF_S[0])
        DF_Ele.append(DF_S[1])
    
    elif ((DF_S[0] in Calo_String) and len(DF_S)==3 ):
        #print(DF_S[1])
        if ((DF_S[1]) in ['1','2','3','4','5']):         
            DF_Det.append(DF_S[0])
            DF_Ele.append(DF_S[2])
        else:
            DF_Det.append(DF_S[0])
            DF_Ele.append(DF_S[1]+DF_S[2])
    elif ((DF_S[0] in Calo_String) and len(DF_S)==4 ):
        DF_Det.append(DF_S[0])
        DF_Ele.append(DF_S[2] + DF_S[3])
    else:
        DF_Det.append("Skip")
        DF_Ele.append("Skip")




test_s=[]
test_Det=[]
test_Ele=[]
Cer_Calo_S2_Split(sorted_Names[40],test_s,test_Det,test_Ele)

Praw = []
Ptime= []
Ptime_sum = []
Peak= []
Peak_sum = []
Ped =[]
Ped_err =[]
for i in range(1,len(sorted_Names)):
    temp_raw = []
    temp_time_sum=[]
    temp_peak_sum=[]
    NPtime = []
    NPeak=[]
    NRaw=[]
    NInt=[]
    NPed=[]
    NPed_err=[]
    temp_string= []
    temp_Detector = []
    temp_Element = []
    Cer_Calo_S2_Split(sorted_Names[i],temp_string,temp_Detector,temp_Element)
    if (temp_Element[0] == 'raw'):
        Branch_To_Array(sorted_Names[i],NRaw)
        NRaw=pd.DataFrame(NRaw[0])
        for k in range(len(NRaw.columns)):
            NRaw = NRaw.rename(columns={list(NRaw)[k]:"%i"%k})
        NRaw.fillna(Undefined_Entry,inplace=True)
        Praw.append(NRaw)
        
    elif (temp_Element[0] == 'Ptime'):
        Branch_To_Array(sorted_Names[i],NPtime)
        temp_time_sum= NPtime[0].flatten()
        temp_time_sum=pd.DataFrame(temp_time_sum,columns=[sorted_Names[i]])

        NPtime=pd.DataFrame(NPtime[0])
        if (len(NPtime.columns)==0):
            NPtime['0'] = Undefined_Entry
        for k in range(len(NPtime.columns)):
            NPtime = NPtime.rename(columns={list(NPtime)[k]:sorted_Names[i]+"_%i"%k})
            NPtime.fillna(Undefined_Entry,inplace=True)
        Ptime.append(NPtime)
        Ptime_sum.append(temp_time_sum)
    elif (temp_Element[0] == 'Ppeak'):
        #print("hello!",i)
        Branch_To_Array(sorted_Names[i],NPeak)
        temp_peak_sum= NPeak[0].flatten()
        temp_peak_sum=pd.DataFrame(temp_peak_sum,columns=[sorted_Names[i]])

        NPeak=pd.DataFrame(NPeak[0])
        if (len(NPeak.columns)==0):
            NPeak['0'] = Undefined_Entry
        for k in range(len(NPeak.columns)):
            NPeak = NPeak.rename(columns={list(NPeak)[k]:sorted_Names[i]+"_%i"%k})
            NPeak.fillna(Undefined_Entry,inplace=True)
        Peak.append(NPeak)
        Peak_sum.append(temp_peak_sum)

    elif (temp_Element[0]=='pedmean'):
        Branch_To_Array(sorted_Names[i],NPed)
        NPed = pd.DataFrame(NPed[0],columns=[sorted_Names[i]])
        Ped.append(NPed)
        
    elif (temp_Element[0]=='pederr'):
        Branch_To_Array(sorted_Names[i],NPed_err)
        NPed_err = pd.DataFrame(NPed_err[0],columns=[sorted_Names[i]])
        Ped_err.append(NPed_err)

    elif (temp_Element[0]=='Skip'):
        continue 


Praw = Praw[0:104]
Ped = Ped[0:104]
Ped_err = Ped_err[0:104]
Peak = Peak[0:104]
Ptime = Ptime[0:104]
Ptime_sum = Ptime_sum[0:104]

Time_Names =[]
for i in range(len(Ptime)):
    Time_Names .append(list(Ptime[i].columns))
Peak_Names =[]
for i in range(len(Peak)):
    Peak_Names.append(list(Peak[i].columns))


Combined_Time_Data=[]
for i in range(len(Ptime)):
    Combined_Time_Data.append(pd.concat([Ptime[i]],axis=1))

Combined_Cer_Data=[]
for i in range(0,len(Praw)):
    #Combined_Cer_Data.append(pd.concat([Praw[i],Peak[i][Peak_Names[i][0]],Ptime[i],Ped[i]],axis=1))
    Combined_Cer_Data.append(pd.concat([Praw[i],Peak[i][Peak_Names[i][0]],Ptime[i],Ped[i],Ped_err[i]],axis=1))
    
    tmp_det = Combined_Cer_Data[i].columns.tolist()[65]
    tmp_det_1=tmp_det.split('_')
    if (len(tmp_det_1)==4 and tmp_det_1[0][0:3] in Cer_String):
        tmp_combined = tmp_det_1[0][0:3]
        tmp_number= int(tmp_det_1[0][3:])
        ID_number = int(tmp_det_1[1][0:])
    elif (len(tmp_det_1)==3 and (tmp_det_1[0][0:2] in S2_String)):
        tmp_combined =tmp_det_1[0][0:1]
        tmp_number= int(tmp_det_1[0][1:])
        ID_number = 1
    elif (len(tmp_det_1)==3 and (tmp_det_1[0][0:2] in Calo_String)):
        tmp_combined =tmp_det_1[0][0:1]
        tmp_number= int(tmp_det_1[0][1:])
        ID_number = 1
    elif (len(tmp_det_1)==4 and (tmp_det_1[0][0:2] in Calo_String)):
        tmp_combined =tmp_det_1[0][0:1]
        tmp_number= int(tmp_det_1[0][1:])
        ID_number = int(tmp_det_1[1][0:])

    elif (len(tmp_det_1)==3 and tmp_det_1[0][0:3] in Cer_String):
        tmp_combined = tmp_det_1[0][0:3]
        tmp_number= int(tmp_det_1[0][3:])
        ID_number = int(tmp_det_1[1][0:])
    else:
        tmp_combined = tmp_det_1[0][0:1]
        tmp_number=int(tmp_det_1[0][1:])
        ID_number=int(tmp_det_1[1][0:])

    Combined_Cer_Data[i]['detector']= tmp_combined
    Combined_Cer_Data[i]['PMT'] = tmp_combined+'_'+str(tmp_number)+'_'+ str(ID_number)
    Combined_Cer_Data[i]['PMT_ID'] = tmp_number
    Combined_Cer_Data[i]['Pixel_ID'] = ID_number

Max_Val = [[],]*len(Praw)

for i in range(len(Praw)):
    Max_Val[i] = Praw[i].max(axis=1)
    Max_Val[i] = pd.DataFrame(Max_Val[i],columns=['Max_Peak'])
    Max_Val[i]['Max_Peak_Pos'] = Praw[i].idxmax(axis=1)
    Max_Val[i]['Max_Peak_Pos'] = Max_Val[i]['Max_Peak_Pos'].astype(int)

Filter_Ped =[[],]*len(Ped)
Ped_Error = [[],]*len(Ped_err)
for i in range(len(Praw)):
    Praw[i] =pd.concat([Praw[i],Max_Val[i]],axis=1)
    Ped[i] = pd.concat([Ped[i],Max_Val[i]],axis=1)
    Ped_err[i] = pd.concat([Ped_err[i],Max_Val[i]],axis=1)
    Combined_Cer_Data[i] = pd.concat([Combined_Cer_Data[i],Max_Val[i]],axis=1)
    Praw[i].query('Max_Peak!=-1',inplace=True)
    Filter_Ped[i] = Ped[i][Ped[i]['Max_Peak']!=-1]
    Ped_Error[i] = Ped_err[i][Ped_err[i]['Max_Peak']!=-1]
    Combined_Cer_Data[i] = Combined_Cer_Data[i][Combined_Cer_Data[i]['Max_Peak']!=-1]

for i in range(len(Praw)):
    Praw[i].drop(['Max_Peak','Max_Peak_Pos'],axis=1,inplace=True)
    Filter_Ped[i].drop(['Max_Peak','Max_Peak_Pos'],axis=1,inplace=True)
    Ped_Error[i].drop(['Max_Peak','Max_Peak_Pos'],axis=1,inplace=True)
    Combined_Cer_Data[i].drop(['Max_Peak','Max_Peak_Pos'],axis=1,inplace=True)
    Combined_Cer_Data[i].reset_index(inplace=True)
def Pedestal_Sub(x,y):
     return(x.sub(y.squeeze(),axis=0))


Pedestal_Subtracted = [[],]*len(Praw)
for i in range(len(Praw)):
    Pedestal_Subtracted[i] = Pedestal_Sub(Praw[i],Filter_Ped[i])
    Pedestal_Subtracted[i].reset_index(inplace=True)
    Ped_Error[i].reset_index(inplace=True)

def Peak_Finder(x,y):
    #return (find_peaks(x,height=20,prominence=20)[0])
    #return (find_peaks(x,height=50,prominence=30)[0])
    return (find_peaks(x,height=2*y,prominence=10)[0])

def Full_Width(x,y):
    return(peak_widths(x,y,rel_height=0.95)[1:] )

Peaks = [[],]*len(Pedestal_Subtracted)
for j in range(len(Pedestal_Subtracted)):
    temp=[]
    for i in range(len(Pedestal_Subtracted[j])):
        temp.append(Peak_Finder(Pedestal_Subtracted[j].iloc[i][1:65],Ped_Error[j].iloc[i][1]))
    Peaks[j]=temp


Integration = [[],]*len(Pedestal_Subtracted)
for j in range(len(Pedestal_Subtracted)):
    temp_int=[]
    for i in range(len(Pedestal_Subtracted[j])):
        temp_int.append(Full_Width(Pedestal_Subtracted[j].iloc[i][1:65],Peaks[j][i]))
    Integration[j]=temp_int



for i in range(len(Peaks)):
    Peaks[i] = pd.DataFrame(Peaks[i])
    Integration[i] = pd.DataFrame(Integration[i])

for i in range(len(Integration)):
    tmp_Int_name = Integration[i].columns.tolist()
    Integration[i].rename(columns={tmp_Int_name[0]:'Int_1'},inplace=True)
    Integration[i].rename(columns={tmp_Int_name[1]:'Int_2'},inplace=True)
    Integration[i].rename(columns={tmp_Int_name[2]:'Int_3'},inplace=True)

Peak_Names = ['peak_1','peak_2','peak_3','peak_4','peak_5','peak_6','peak_7','peak_8','peak_9','peak_10','peak_11']

for i in range(len(Peaks)):
    tmp_name = Peaks[i].columns.tolist()

    if (len(tmp_name) == 1):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)

    elif (len(tmp_name) == 2):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)

    elif (len(tmp_name) == 3):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)
        Peaks[i].rename(columns={tmp_name[2]:Peak_Names[2]},inplace=True)

    elif (len(tmp_name) == 4):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)
        Peaks[i].rename(columns={tmp_name[2]:Peak_Names[2]},inplace=True)
        Peaks[i].rename(columns={tmp_name[3]:Peak_Names[3]},inplace=True)

    elif (len(tmp_name) == 5):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)
        Peaks[i].rename(columns={tmp_name[2]:Peak_Names[2]},inplace=True)
        Peaks[i].rename(columns={tmp_name[3]:Peak_Names[3]},inplace=True)
        Peaks[i].rename(columns={tmp_name[4]:Peak_Names[4]},inplace=True)

    elif (len(tmp_name) == 6):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)
        Peaks[i].rename(columns={tmp_name[2]:Peak_Names[2]},inplace=True)
        Peaks[i].rename(columns={tmp_name[3]:Peak_Names[3]},inplace=True)
        Peaks[i].rename(columns={tmp_name[4]:Peak_Names[4]},inplace=True)
        Peaks[i].rename(columns={tmp_name[5]:Peak_Names[5]},inplace=True)


    elif (len(tmp_name) == 7):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)
        Peaks[i].rename(columns={tmp_name[2]:Peak_Names[2]},inplace=True)
        Peaks[i].rename(columns={tmp_name[3]:Peak_Names[3]},inplace=True)
        Peaks[i].rename(columns={tmp_name[4]:Peak_Names[4]},inplace=True)
        Peaks[i].rename(columns={tmp_name[5]:Peak_Names[5]},inplace=True)
        Peaks[i].rename(columns={tmp_name[6]:Peak_Names[6]},inplace=True)

    elif (len(tmp_name) == 8):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)
        Peaks[i].rename(columns={tmp_name[2]:Peak_Names[2]},inplace=True)
        Peaks[i].rename(columns={tmp_name[3]:Peak_Names[3]},inplace=True)
        Peaks[i].rename(columns={tmp_name[4]:Peak_Names[4]},inplace=True)
        Peaks[i].rename(columns={tmp_name[5]:Peak_Names[5]},inplace=True)
        Peaks[i].rename(columns={tmp_name[6]:Peak_Names[6]},inplace=True)
        Peaks[i].rename(columns={tmp_name[7]:Peak_Names[7]},inplace=True)

    elif (len(tmp_name) == 9):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)
        Peaks[i].rename(columns={tmp_name[2]:Peak_Names[2]},inplace=True)
        Peaks[i].rename(columns={tmp_name[3]:Peak_Names[3]},inplace=True)
        Peaks[i].rename(columns={tmp_name[4]:Peak_Names[4]},inplace=True)
        Peaks[i].rename(columns={tmp_name[5]:Peak_Names[5]},inplace=True)
        Peaks[i].rename(columns={tmp_name[6]:Peak_Names[6]},inplace=True)
        Peaks[i].rename(columns={tmp_name[7]:Peak_Names[7]},inplace=True)
        Peaks[i].rename(columns={tmp_name[8]:Peak_Names[8]},inplace=True)

    elif (len(tmp_name) == 10):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)
        Peaks[i].rename(columns={tmp_name[2]:Peak_Names[2]},inplace=True)
        Peaks[i].rename(columns={tmp_name[3]:Peak_Names[3]},inplace=True)
        Peaks[i].rename(columns={tmp_name[4]:Peak_Names[4]},inplace=True)
        Peaks[i].rename(columns={tmp_name[5]:Peak_Names[5]},inplace=True)
        Peaks[i].rename(columns={tmp_name[6]:Peak_Names[6]},inplace=True)
        Peaks[i].rename(columns={tmp_name[7]:Peak_Names[7]},inplace=True)
        Peaks[i].rename(columns={tmp_name[8]:Peak_Names[8]},inplace=True)
        Peaks[i].rename(columns={tmp_name[9]:Peak_Names[9]},inplace=True)

    elif (len(tmp_name) == 11):
        Peaks[i].rename(columns={tmp_name[0]:Peak_Names[0]},inplace=True)
        Peaks[i].rename(columns={tmp_name[1]:Peak_Names[1]},inplace=True)
        Peaks[i].rename(columns={tmp_name[2]:Peak_Names[2]},inplace=True)
        Peaks[i].rename(columns={tmp_name[3]:Peak_Names[3]},inplace=True)
        Peaks[i].rename(columns={tmp_name[4]:Peak_Names[4]},inplace=True)
        Peaks[i].rename(columns={tmp_name[5]:Peak_Names[5]},inplace=True)
        Peaks[i].rename(columns={tmp_name[6]:Peak_Names[6]},inplace=True)
        Peaks[i].rename(columns={tmp_name[7]:Peak_Names[7]},inplace=True)
        Peaks[i].rename(columns={tmp_name[8]:Peak_Names[8]},inplace=True)
        Peaks[i].rename(columns={tmp_name[9]:Peak_Names[9]},inplace=True)
        Peaks[i].rename(columns={tmp_name[10]:Peak_Names[10]},inplace=True)

# Add event index to Peaks and Integration info.
Pedestal_Sub_Combined = [[],]*len(Pedestal_Subtracted)
for i in range(len(Pedestal_Subtracted)):
    Pedestal_Sub_Combined[i] = pd.concat([Pedestal_Subtracted[i],Combined_Cer_Data[i]['PMT'],Combined_Cer_Data[i]['detector'],Combined_Cer_Data[i]['PMT_ID'],Combined_Cer_Data[i]['Pixel_ID'],Peaks[i],Integration[i]],axis=1)
    #Peaks[i] = pd.concat([Pedestal_Subtracted[i]['index'],Peaks[i]])
    #Integration[i] = pd.concat([Pedestal_Subtracted[i]['index'],Integration[i]])

#del Combined_Cer_Data
#gc.collect()

Combined_Detectors = pd.concat(Pedestal_Sub_Combined,axis=0)
Combined_Detectors.replace(np.nan,-1,inplace=True)

Grouped_Events = Combined_Detectors.groupby('index')
Combined_group_index = list(Grouped_Events.groups.keys())

Grouped_By_Event = [[],]*len(Combined_group_index)
for i in range(len(Combined_group_index)):
    Grouped_By_Event[i] = Grouped_Events.get_group(Combined_group_index[i])


Grouped_By_Detector = []
for i in range(len(Grouped_By_Event)):
    Grouped_By_Detector.append(Grouped_By_Event[i].groupby('detector'))

Sub_Detectors = [[],]*len(Grouped_By_Detector)
Sub_Detectors_Comp = [[],]*len(Grouped_By_Detector)

for i in range(len(Grouped_By_Detector)):
    tmp=[[],]*3
    tmp_Keys = list(Grouped_By_Detector[i].groups.keys())
    for j in range(len(tmp_Keys)):
        tmp[j]=(Grouped_By_Detector[i].get_group(tmp_Keys[j]))
    Sub_Detectors[i]=(tmp)
    Sub_Detectors_Comp[i] = (tmp)

Detectors_To_Del = []
Detectors_To_Del_Comp =[]
for i in range(len(Sub_Detectors)):
    if (len(Sub_Detectors[i][0])==0):
        Detectors_To_Del.append(i)
        Detectors_To_Del_Comp.append(i)
    else:
        temp_sub_Calo = Sub_Detectors[i][0][Sub_Detectors[i][0]['peak_1']>27]
        Sub_Detectors[i][0] = temp_sub_Calo[temp_sub_Calo['peak_1']<32]
        if (len(Sub_Detectors[i][0])==0):
            Detectors_To_Del.append(i)
    if (len(Sub_Detectors[i][1])==0):
        Detectors_To_Del.append(i)
        #Detectors_To_Del_Comp.append(i)

    else:
        temp_sub_Cer_1 = Sub_Detectors[i][1][Sub_Detectors[i][1]['peak_1']>6]
        temp_sub_Cer_2 = temp_sub_Cer_1[temp_sub_Cer_1['peak_1']<30]

        temp_sub_Cer_3 = Sub_Detectors[i][1][Sub_Detectors[i][1]['peak_2']>6]
        temp_sub_Cer_4 = temp_sub_Cer_3[temp_sub_Cer_3['peak_2']<30]

        temp_sub_Cer_5 = Sub_Detectors[i][1][Sub_Detectors[i][1]['peak_3']>6]
        temp_sub_Cer_6 = temp_sub_Cer_5[temp_sub_Cer_5['peak_3']<30]


        temp_sub_Cer_7 = Sub_Detectors[i][1][Sub_Detectors[i][1]['peak_4']>6]
        temp_sub_Cer_8 = temp_sub_Cer_7[temp_sub_Cer_7['peak_4']<30]

        temp_sub_Cer_9 = Sub_Detectors[i][1][Sub_Detectors[i][1]['peak_5']>6]
        temp_sub_Cer_10 = temp_sub_Cer_9[temp_sub_Cer_9['peak_5']<30]

        if (len(temp_sub_Cer_2)==0 and len(temp_sub_Cer_4)==0 and len(temp_sub_Cer_6)==0  and len(temp_sub_Cer_8)==0 and len(temp_sub_Cer_10)==0 ):
            Detectors_To_Del.append(i)
    Sub_Detectors[i][1] = pd.concat([temp_sub_Cer_2,temp_sub_Cer_4,temp_sub_Cer_6],axis=0)

Sorted_Detectors_To_Del = sorted(set(Detectors_To_Del))
for i in sorted(Sorted_Detectors_To_Del, reverse=True):
    del Sub_Detectors[i]

for i in sorted(Detectors_To_Del_Comp, reverse=True):
    del Sub_Detectors_Comp[i]


for i in range(len(Sub_Detectors)):
    Sub_Detectors[i][0]['Num_Pixels'] = Sub_Detectors[i][0]['Pixel_ID'].count()
    Sub_Detectors[i][0]['Num_PMTs'] = Sub_Detectors[i][0]['PMT_ID'].nunique()
for i in range(len(Sub_Detectors)):
    Sub_Detectors[i][1]['Num_Pixels'] = Sub_Detectors[i][1][Sub_Detectors[i][1]['Pixel_ID']!=5].count()['Pixel_ID']
    Sub_Detectors[i][1]['Num_PMTs'] = Sub_Detectors[i][1][Sub_Detectors[i][1]['Pixel_ID']!=5]['PMT_ID'].nunique()

Calorimeter_Signals_1 = []
for i in range(len(Sub_Detectors)):
    Calorimeter_Signals_1.append(Sub_Detectors[i][0])
Calorimeter_Signals= pd.concat(Calorimeter_Signals_1)
Calorimeter_Signals.reset_index(inplace=True,drop=True)

Group_Calorimeter_PMT = Calorimeter_Signals.groupby('PMT')
Calo_Keys = list(Group_Calorimeter_PMT.groups.keys())

Calo_PMT = []
for i in range(len(Calo_Keys)):
    Calo_PMT.append(Group_Calorimeter_PMT.get_group(Calo_Keys[i]))
Waveform_Calo_PMT = Calo_PMT.copy()
#Copy Raw Waveform
for i in range(len(Waveform_Calo_PMT)):
    Waveform_Calo_PMT[i] = Waveform_Calo_PMT[i].loc[:,'0':'63']
######

Good_Cerenkov_Signals_1=[]
for i in range(len(Sub_Detectors)):
    Good_Cerenkov_Signals_1.append(Sub_Detectors[i][1])

Good_Cerenkov_Signals= pd.concat(Good_Cerenkov_Signals_1)
#Good_Cerenkov_Signals.reset_index(inplace=True)
Good_Cerenkov_Signals.reset_index(inplace=True,drop=True)

Group_Cerenkov_PMT = Good_Cerenkov_Signals.groupby('PMT')
Cer_Keys = list(Group_Cerenkov_PMT.groups.keys())

Cer_PMT = []
for i in range(len(Cer_Keys)):
    Cer_PMT.append(Group_Cerenkov_PMT.get_group(Cer_Keys[i]))

Waveform_Cer_PMT = Cer_PMT.copy()

# Copy only Raw Waveform data ---- Use to sum signals
for i in range(len(Waveform_Cer_PMT)):
    Waveform_Cer_PMT[i] = Waveform_Cer_PMT[i].loc[:,'0':'63'] 


# Need to separate the Integration Information based on the correct peak

def Integration_1(w,x,y,z,m,n,q,a,b):
    if (len(w)==1):
        if (x>=(q-4) and x<=(q+4)):
            return(w[0],int(a[0]),int(b[0]),int(x))
        else:
            return(0,0,0,0)
    elif (len(w)==2):
        #if (x>=7 and x<=30):
        if (x>=(q-4) and x<=(q+4)):
            return(w[0],int(a[0]),int(b[0]),int(x))
        elif (y>=(q-4) and y<=(q+4)):
            return(w[1],int(a[1]),int(b[1]),int(y))
        elif (z>=(q-4) and z<=(q+4)):
            return(w[2],int(a[2]),int(b[2]),int(z))
        elif (m>=(q-4) and m<=(q+4)):
            return(w[3],int(a[3]),int(b[3]),int(m))
        elif (n>=(q-4) and n<=(q+4)):
            return(w[4],int(a[4]),int(b[4]),int(n))
        else:
            return(0,0,0,0)
    elif (len(w)==3):
        if (x>=(q-4) and x<=(q+4)):
            return(w[0],int(a[0]),int(b[0]),int(x))
        elif (y>=(q-4) and y<=(q+4)):
            return(w[1],int(a[1]),int(b[1]),int(y))
        elif (z>=(q-4) and z<=(q+4)):
            return(w[2],int(a[2]),int(b[2]),int(z))
        elif (m>=(q-4) and m<=(q+4)):
            return(w[3],int(a[3]),int(b[3]),int(m))
        elif (n>=(q-4) and n<=(q+4)):
            return(w[4],int(a[4]),int(b[4]),int(n))
        else:
            return(0,0,0,0)
    
    elif (len(w)==4):
        if (x>=(q-4) and x<=(q+4)):
            return(w[0],int(a[0]),int(b[0]),int(x))
        elif (y>=(q-4) and y<=(q+4)):
            return(w[1],int(a[1]),int(b[1]),int(y))
        elif (z>=(q-4) and z<=(q+4)):
            return(w[2],int(a[2]),int(b[2]),int(z))
        elif (m>=(q-4) and m<=(q+4)):
            return(w[3],int(a[3]),int(b[3]),int(m))
        elif (n>=(q-4) and n<=(q+4)):
            return(w[4],int(a[4]),int(b[4]),int(n))
        else:
            return(0,0,0,0)
    else:
        return(0,0,0,0)

#Int_tmp_1 = [[],]*len(Cer_PMT)
#Int_tmp_2 = [[],]*len(Cer_PMT)
#Int_tmp_3 = [[],]*len(Cer_PMT)
############
Int_tmp_1 = []
Int_tmp_2 = []
Int_tmp_3 = []
Int_tmp_4 =[]
for i in range(len(Cer_PMT)):
    #print(i)
    Mean_Peak = Cer_PMT[i][Cer_PMT[i]['peak_1']!=-1].mean()['peak_1']
    Mode_tmp = Cer_PMT[i][Cer_PMT[i]['peak_1']!=-1]
    Mode_Peak = Mode_tmp['peak_1'].mode()[0]
    Avg_Value = (Mean_Peak+Mode_Peak)/2
    tmp_1,tmp_2,tmp_3,tmp_4=(np.vectorize(Integration_1)(Cer_PMT[i]['Int_1'],Cer_PMT[i]['peak_1'],Cer_PMT[i]['peak_2'],Cer_PMT[i]['peak_3'],Cer_PMT[i]['peak_4'],Cer_PMT[i]['peak_5'],Avg_Value,Cer_PMT[i]['Int_2'],Cer_PMT[i]['Int_3']))
    #Int_tmp_1[i],Int_tmp_2[i],Int_tmp_3[i]=(np.vectorize(Integration_1)(Cer_PMT[i]['Int_1'],Cer_PMT[i]['peak_1'],Cer_PMT[i]['peak_2'],Cer_PMT[i]['peak_3'],Cer_PMT[i]['Int_2'],Cer_PMT[i]['Int_3']))
    #print(tmp_1,tmp_2,tmp_3)
    Int_tmp_1.append(tmp_1)
    Int_tmp_2.append(tmp_2)
    Int_tmp_3.append(tmp_3)
    Int_tmp_4.append(tmp_4)

Int_1 = [[],]*len(Cer_PMT) 
Int_2 = [[],]*len(Cer_PMT) 
Int_3 = [[],]*len(Cer_PMT) 
Int_4 = [[],]*len(Cer_PMT) 

for i in range(len(Cer_PMT)):
    Int_1[i] = pd.DataFrame(Int_tmp_1[i],columns=['Int_1'])
    Int_2[i] = pd.DataFrame(Int_tmp_2[i],columns=['Int_2'])
    Int_3[i] = pd.DataFrame(Int_tmp_3[i],columns=['Int_3'])
    Int_4[i] = pd.DataFrame(Int_tmp_4[i],columns=['Peak'])


Combined_Integration = [[],]*len(Cer_PMT)
for i in range(len(Cer_PMT)):
    Combined_Integration[i] = pd.concat([Int_1[i],Int_2[i],Int_3[i],Int_4[i]],axis=1)
################
Int_tmp_5 = []
Int_tmp_6 = []
Int_tmp_7 = []
Int_tmp_8 =[]
for i in range(len(Calo_PMT)):
    #print(i)
    Mean_Peak = Calo_PMT[i][Calo_PMT[i]['peak_1']!=-1].mean()['peak_1']
    Mode_tmp = Calo_PMT[i][Calo_PMT[i]['peak_1']!=-1]
    Mode_Peak = Mode_tmp['peak_1'].mode()[0]
    Avg_Value = (Mean_Peak+Mode_Peak)/2
    tmp_5,tmp_6,tmp_7,tmp_8=(np.vectorize(Integration_1)(Calo_PMT[i]['Int_1'],Calo_PMT[i]['peak_1'],Calo_PMT[i]['peak_2'],Calo_PMT[i]['peak_3'],Calo_PMT[i]['peak_4'],Calo_PMT[i]['peak_5'],Avg_Value,Calo_PMT[i]['Int_2'],Calo_PMT[i]['Int_3']))

    Int_tmp_5.append(tmp_5)
    Int_tmp_6.append(tmp_6)
    Int_tmp_7.append(tmp_7)
    Int_tmp_8.append(tmp_8)

Int_5 = [[],]*len(Calo_PMT) 
Int_6 = [[],]*len(Calo_PMT) 
Int_7 = [[],]*len(Calo_PMT) 
Int_8 = [[],]*len(Calo_PMT) 

for i in range(len(Calo_PMT)):
    Int_5[i] = pd.DataFrame(Int_tmp_5[i],columns=['Int_1'])
    Int_6[i] = pd.DataFrame(Int_tmp_6[i],columns=['Int_2'])
    Int_7[i] = pd.DataFrame(Int_tmp_7[i],columns=['Int_3'])
    Int_8[i] = pd.DataFrame(Int_tmp_8[i],columns=['Peak'])


Combined_Calo_Integration = [[],]*len(Calo_PMT)
for i in range(len(Calo_PMT)):
    Combined_Calo_Integration[i] = pd.concat([Int_5[i],Int_6[i],Int_7[i],Int_8[i]],axis=1)


###############
def Signal_Sum(x,y,z):
    return(x[y:z].sum())
Summation = [[],]*len(Waveform_Cer_PMT)
for i in range(len(Waveform_Cer_PMT)):
    temp_int=[]
    for j in range(len(Waveform_Cer_PMT[i])):
        #Summation.append(np.vectorize(Signal_Sum)(Waveform_Cer_PMT[i].iloc[j],Combined_Integration[i]['Int_2'].iloc[j],Combined_Integration[i]['Int_3'].iloc[j]))
        temp_int.append(Waveform_Cer_PMT[i].iloc[j][Combined_Integration[i]['Int_2'].iloc[j]:Combined_Integration[i]['Int_3'].iloc[j]].sum())
    Summation[i] = temp_int


for i in range(len(Summation)):
    Cer_PMT[i].reset_index(inplace=True,drop=True)
    Summation[i] = pd.DataFrame(Summation[i],columns=['sum'])
    Summation[i]= pd.concat([Summation[i],Cer_PMT[i]['index'],Cer_PMT[i]['Num_PMTs'],Cer_PMT[i]['Num_Pixels'],Cer_PMT[i]['PMT'],Cer_PMT[i]['PMT_ID'],Cer_PMT[i]['Pixel_ID'],Cer_PMT[i]['detector'],Combined_Integration[i]['Peak']],axis=1)
####
#Sum the Calo Signals
Calo_Summation = [[],]*len(Waveform_Calo_PMT)
for i in range(len(Waveform_Calo_PMT)):
    temp_int2=[]
    for j in range(len(Waveform_Calo_PMT[i])):
        temp_int2.append(Waveform_Calo_PMT[i].iloc[j][Combined_Calo_Integration[i]['Int_2'].iloc[j]:Combined_Calo_Integration[i]['Int_3'].iloc[j]].sum())
    Calo_Summation[i] = temp_int2
for i in range(len(Calo_Summation)):
    Calo_PMT[i].reset_index(inplace=True,drop=True)
    Calo_Summation[i] = pd.DataFrame(Calo_Summation[i],columns=['sum'])
    Calo_Summation[i]= pd.concat([Calo_Summation[i],Calo_PMT[i]['index'],Calo_PMT[i]['Num_PMTs'],Calo_PMT[i]['Num_Pixels'],Calo_PMT[i]['PMT'],Calo_PMT[i]['PMT_ID'],Calo_PMT[i]['Pixel_ID'],Calo_PMT[i]['detector'],Combined_Calo_Integration[i]['Peak']],axis=1)

#####
# Group the events together ----- This is redundant now, but want to include Calo Signal in the Neural Network. Have to make same cuts on same events
tmp_summation_1 = pd.concat(Summation,axis=0)
tmp_summation_2 = pd.concat(Calo_Summation,axis=0)
Combined_Signals = pd.concat([tmp_summation_1,tmp_summation_2],axis=0)
Combined_Signals_Group = Combined_Signals.groupby('index')
Combined_Signals_Keys = list(Combined_Signals_Group.groups.keys())



Combined_Summation=[[],]*len(Combined_Signals_Keys)
for i in range(len(Combined_Signals_Keys)):
    Combined_Summation[i] = Combined_Signals_Group.get_group(Combined_Signals_Keys[i])

Combined_Events_To_Del =[]
for i in range(len(Combined_Summation)):
    mode_CS = Combined_Summation[i]['Num_PMTs'].mode()[0]
    #if (mode_CS >12 and mode_CS <24):
    if (mode_CS <3 and mode_CS >5):
        Combined_Events_To_Del.append(i)

Ring_Events = Combined_Summation.copy()
for i in sorted(Combined_Events_To_Del, reverse=True):
    #del Combined_Summation[i]
    del Ring_Events[i]
Event_Ring_Sum =[[],]*len(Ring_Events)
Calo_Ring_Sum = [[],]*len(Ring_Events)
for i in range(len(Ring_Events)):
    Event_Ring_Sum[i]=Ring_Events[i][Ring_Events[i]['detector']!='C']
    Calo_Ring_Sum[i] = Ring_Events[i][Ring_Events[i]['detector']!='Cer']


Event_Wise_Sum =[[],]*len(Combined_Summation)
Calo_Event_Wise_Sum = [[],]*len(Combined_Summation)
for i in range(len(Combined_Summation)):
    Event_Wise_Sum[i]=Combined_Summation[i][Combined_Summation[i]['detector']!='C']
    Calo_Event_Wise_Sum[i] = Combined_Summation[i][Combined_Summation[i]['detector']!='Cer']


for i in range(len(Event_Ring_Sum)):
    Event_Ring_Sum[i]['Map_ID'] = (Event_Ring_Sum[i]['PMT_ID'].astype(str)+ Event_Ring_Sum[i]['Pixel_ID'].astype(str)).astype(int)
    Calo_Ring_Sum[i]['Map_ID'] = (Calo_Ring_Sum[i]['PMT_ID'].astype(str)+ Calo_Ring_Sum[i]['Pixel_ID'].astype(str)).astype(int)

for i in range(len(Event_Wise_Sum)):
    Event_Wise_Sum[i]['Map_ID'] = (Event_Wise_Sum[i]['PMT_ID'].astype(str)+ Event_Wise_Sum[i]['Pixel_ID'].astype(str)).astype(int)
    Calo_Event_Wise_Sum[i]['Map_ID'] = (Calo_Event_Wise_Sum[i]['PMT_ID'].astype(str)+ Calo_Event_Wise_Sum[i]['Pixel_ID'].astype(str)).astype(int)

Array_1 = [141, 142, 131, 132, 121, 122, 111, 112]
Array_2 = [143, 144, 133, 134, 123, 124, 113, 114]
Array_3 = [241, 242, 231, 232, 221, 222, 211, 212]
Array_4 = [243, 244, 233, 234, 223, 224, 213, 214]
Array_5 = [341, 342, 331, 332, 321, 322, 311, 312]
Array_6 = [343, 344, 333, 334, 323, 324, 313, 314]
Array_7 = [441, 442, 431, 432, 421, 422, 411, 412] 
Array_8 = [443, 444, 433, 434, 423, 424, 413, 414] 

Array = [Array_1,Array_2,Array_3,Array_4,Array_5,Array_6,Array_7,Array_8]
HitMap_Array=np.asarray(Array)
###
Calo_Array_1 =[41,41,61,64,71,74]
Calo_Array_2 =[41,41,62,63,72,73]
Calo_Array_3 =[51,54,91,92,82,83]
Calo_Array_4 =[52,53,94,93,81,84]
Calo_Array_5 =[11,11,21,21,31,31]

Calo_Array = [Calo_Array_1,Calo_Array_2,Calo_Array_3,Calo_Array_4,Calo_Array_5]
Calo_HitMap_Array = np.asarray(Calo_Array)

####

Hit_Map_Ring = [[],]*len(Event_Ring_Sum)
Hit_Map_Ring_Time = [[],]*len(Event_Ring_Sum)
for i in range(len(Event_Ring_Sum)):
    tmp_Array=HitMap_Array.copy()
    tmp_Array.fill(0)
    tmp_Array_2=HitMap_Array.copy()
    tmp_Array_2.fill(0)
    for j in range(len(Event_Ring_Sum[i])):
        tmp = np.asarray(np.where(HitMap_Array == Event_Ring_Sum[i]['Map_ID'].iloc[j]))
        if ((tmp.size)!=0):
            tmp=tmp.flatten()
            tmp_Array[tmp[0]][tmp[1]]=Event_Ring_Sum[i]['sum'].iloc[j]
            tmp_Array_2[tmp[0]][tmp[1]]=Event_Ring_Sum[i]['Peak'].iloc[j]
        
    Hit_Map_Ring[i]=tmp_Array
    Hit_Map_Ring_Time[i]=tmp_Array_2
Calo_Hit_Map_Ring = [[],]*len(Event_Ring_Sum)
Calo_Hit_Map_Ring_Time = [[],]*len(Event_Ring_Sum)
for i in range(len(Calo_Ring_Sum)):
    tmp_Array_3=Calo_HitMap_Array.copy()
    tmp_Array_3.fill(0)
    tmp_Array_4=Calo_HitMap_Array.copy()
    tmp_Array_4.fill(0)
    for j in range(len(Calo_Ring_Sum[i])):
        tmp_2 = np.asarray(np.where(Calo_HitMap_Array == Calo_Ring_Sum[i]['Map_ID'].iloc[j]))
        if ((tmp_2.size)!=0):
            tmp_2=tmp_2.flatten()
            tmp_Array_3[tmp_2[0]][tmp_2[1]]=Calo_Ring_Sum[i]['sum'].iloc[j]
            tmp_Array_4[tmp_2[0]][tmp_2[1]]=Calo_Ring_Sum[i]['Peak'].iloc[j]

    Calo_Hit_Map_Ring[i]=tmp_Array_3
    Calo_Hit_Map_Ring_Time[i]=tmp_Array_4
Calo_Ring_Hit_Map = np.stack(Calo_Hit_Map_Ring)
Calo_Ring_Hit_Time_Map = np.stack(Calo_Hit_Map_Ring_Time)
Combined_Ring_Hit_Map = np.stack(Hit_Map_Ring)
Combined_Ring_Time_Hit_Map = np.stack(Hit_Map_Ring_Time)
#####
Hit_Map = [[],]*len(Event_Wise_Sum)
Hit_Map_2 = [[],]*len(Event_Wise_Sum)
for i in range(len(Event_Wise_Sum)):
     tmp_Array=HitMap_Array.copy()
     tmp_Array.fill(0)
     tmp_Array_2=HitMap_Array.copy()
     tmp_Array_2.fill(0)
     for j in range(len(Event_Wise_Sum[i])):
             tmp = np.asarray(np.where(HitMap_Array == Event_Wise_Sum[i]['Map_ID'].iloc[j]))
             if ((tmp.size)!=0):
                     tmp=tmp.flatten()
                     tmp_Array[tmp[0]][tmp[1]]=Event_Wise_Sum[i]['sum'].iloc[j]
                     tmp_Array_2[tmp[0]][tmp[1]]=Event_Wise_Sum[i]['Peak'].iloc[j]
     Hit_Map[i]=tmp_Array
     Hit_Map_2[i]=tmp_Array_2
###
Calo_Hit_Map = [[],]*len(Calo_Event_Wise_Sum)
Calo_Hit_Map_2 = [[],]*len(Calo_Event_Wise_Sum)
for i in range(len(Calo_Event_Wise_Sum)):
     tmp_Array_3=Calo_HitMap_Array.copy()
     tmp_Array_3.fill(0)
     tmp_Array_4=Calo_HitMap_Array.copy()
     tmp_Array_4.fill(0)
     for j in range(len(Calo_Event_Wise_Sum[i])):
             tmp_2 = np.asarray(np.where(Calo_HitMap_Array == Calo_Event_Wise_Sum[i]['Map_ID'].iloc[j]))
             if ((tmp_2.size)!=0):
                     tmp_2=tmp_2.flatten()
                     tmp_Array_3[tmp_2[0]][tmp_2[1]]=Calo_Event_Wise_Sum[i]['sum'].iloc[j]
                     tmp_Array_4[tmp_2[0]][tmp_2[1]]=Calo_Event_Wise_Sum[i]['Peak'].iloc[j]
     Calo_Hit_Map[i]=tmp_Array_3
     Calo_Hit_Map_2[i]=tmp_Array_4







# Combine to output to read in for ML
Good_Hit_Map= np.stack(Hit_Map_Ring)
Good_Hit_Time_Map= np.stack(Hit_Map_Ring_Time)
#Good_Calo_Hit_Map= np.stack(Good_Events_Calo)
#Good_Calo_Hit_Time_Map= np.stack(Good_Events_Calo_Time)


np.save('/lustre19/expphy/volatile/halla/triton/mnycz/SOLID/ML/Total_Hit_Map_302.npy',Good_Hit_Map)
np.save('/lustre19/expphy/volatile/halla/triton/mnycz/SOLID/ML/Total_Hit_Map_302_Time.npy',Good_Hit_Time_Map)



print(datetime.now() - startTime)


