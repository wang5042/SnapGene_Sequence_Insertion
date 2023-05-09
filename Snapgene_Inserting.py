from snapgene_reader import snapgene_file_to_seqrecord, snapgene_file_to_gbk
from Bio import SeqIO
import os
import sys
sys.path.append("D:/PycharmProjects/SnapgeneProcessing/snapgenefile-genbankfile-Trans")
import Script


filepath = "D:/PycharmProjects/SnapgeneProcessing/"
Insert_seq = "D:/PycharmProjects/SnapgeneProcessing/Pre-editing/Insert_Sequence/UAS-Pxyl.dna"
'''
The following above is the parameter you need to change,while the below is the program
'''


Insertion_Start_Site= 528 #Starting site of backbone for Insert_seq
Insertion_End_Site = 941 #
insert_seqrecord = Script.snapgene_file_to_dict(Insert_seq)
#print(insert_seqrecord)
#print(len(insert_seqrecord['seq']))
insert_seqrecord_feature = []
for feature in insert_seqrecord['features']:
    #print(feature)
    feature['start']=feature['start']+Insertion_Start_Site
    feature['end']=feature['end']+Insertion_Start_Site-1
    #print(feature)
    insert_seqrecord_feature.append(feature)
del insert_seqrecord['features']
for root, dirs, files in os.walk(filepath+'Pre-editing/Backbone_Sequence/'):
    print(files)
for file in files:
    bone_seqrecord = Script.snapgene_file_to_dict(filepath+'Pre-editing/Backbone_Sequence/'+file)
#print(bone_seqrecord)
#print(len(bone_seqrecord['seq']))
#print(type(bone_seqrecord['features']))
    feature_0 = []
    for feature in bone_seqrecord['features']:
        #print(type(feature))
        #print(feature)
        #print(type(feature['start']))
        if feature['start']>Insertion_End_Site:
            feature['start'] = feature['start']-(Insertion_End_Site-Insertion_Start_Site+1-len(insert_seqrecord['seq'])) #+57
            feature['end'] = feature['end']-(Insertion_End_Site-Insertion_Start_Site+1-len(insert_seqrecord['seq']))#+56
            #print(feature)
            feature_0.append(feature)
        if feature['end']<Insertion_Start_Site:
            feature_0.append(feature)
#print(bone_seqrecord['features'])
#print(feature_0)
    del bone_seqrecord['features']
    bone_seqrecord['features'] = feature_0+insert_seqrecord_feature
    #print(type(bone_seqrecord['features']))
    #print(bone_seqrecord['features'])

    insertion_site= Insertion_Start_Site-1
    new_sequence = bone_seqrecord['seq'][:insertion_site] + insert_seqrecord['seq'] + bone_seqrecord["seq"][Insertion_End_Site+1:]
    del bone_seqrecord['seq']
    bone_seqrecord['seq']=new_sequence
    Script.snapgene_file_to_gbk(bone_seqrecord,filepath+"Already-Editing/"+file[0:12]+"-UASpxyl-"+file[18:-4]+".gbk")
#print()

