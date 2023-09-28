import numpy as np
from Feature.read_fasta_sequences import read_nucleotide_sequences
from Feature.Kmer import Kmer
from Feature.RCKmer import RCKmer
from Feature.PseEIIP import PseEIIP
import joblib
import pandas as pd
from sklearn.decomposition import PCA
                

file = 'Inpactor2_library.fasta'
fastas = read_nucleotide_sequences(file)
total_result = []
for sample in fastas:
  result = []
  sample = [sample]
  ## ------------------------------------------feature---------------------------------------------
  fea_stage1_RCKmer = np.array(RCKmer(sample, k=3, upto=False, normalize=True))
  fea_stage1_Kmer = np.array(Kmer(sample, k=3, type="DNA", upto=False, normalize=True))
  fea_stage1_PseEIIP = np.array(PseEIIP(sample))
  fea_stage1 = np.hstack((fea_stage1_RCKmer,fea_stage1_Kmer,fea_stage1_PseEIIP))
  # print(fea_stage1.shape)  # (1, 160)

  feature = np.array(Kmer(sample, k=3, type="DNA", upto=False, normalize=True))
  # print(feature.shape)   # (1, 64)

  ## ------------------------------------------model---------------------------------------------
  model_stage1 = joblib.load('./Model/model1_lgb')

  pre_stage1 = model_stage1.predict(fea_stage1)

  if pre_stage1 >= 0.5:
      result.append('LTR-RT')
      print('LTR-RT')
      model_sup = joblib.load('./Model/Sup_Kmer(64)_knn')
      pre_sup = model_sup.predict(feature)
  
      if pre_sup >= 0.5:
          print('Superfamily: RLC')
          result.append('Superfamily: RLC')
          Model = joblib.load('./Model/Kmer(64)_1_knn')
          pre_1 = Model.predict(feature)
          if pre_1 >= 0.5:
              print('Lineage: ALE/RETROFIT')
              result.append('Lineage: ALE/RETROFIT')              
          else:
              model_2 = joblib.load('./Model/Kmer(64)_2_knn')
              pre_2 = model_2.predict(feature)
              if pre_2 >= 0.5:
                  print('Lineage: TORK/TAR')
                  result.append('Lineage: TAR/TORK')
              else:
                  model_3 = joblib.load('./Model/Kmer(64)_3_knn')
                  pre_3 = model_3.predict(feature)
                  if pre_3 >= 0.5:
                      model_4 = joblib.load('./Model/Kmer(64)_4_knn')
                      pre_4 = model_4.predict(feature)
                      if pre_4 >= 0.5:
                          print('Lineage: BIANCA')
                          result.append('Lineage: BIANCA')
                      else:
                          print('Lineage: SIRE')
                          result.append('Lineage: SIRE')
                  else:
                      model_5 = joblib.load('./Model/Kmer(64)_5_svm')
                      pre_5 = model_5.predict(feature)
                      if pre_5 >= 0.5:
                          print('Lineage: IVANA-ORYCO')
                          result.append('Lineage: IVANA-ORYCO')                          
                      else:
                          model_6 = joblib.load('./Model/Kmer(64)_6_rf_SMOTH')
                          pre_6 = model_6.predict(feature)
                          if pre_6 >= 0.5:
                              print('Lineage: ANGELA')
                              result.append('Lineage: ANGELA')                                                        
                          else:
                              model_7 = joblib.load('./Model/Kmer(64)_7_knn')
                              pre_7 = model_7.predict(feature)
                              if pre_7 >= 0.5:
                                  print('Lineage: IKEROS')
                                  result.append('Lineage: IKEROS')                                                                                          
                              else:
                                  print('Lineage: IVANA-ORYCO')
                                  result.append('Lineage: IVANA-ORYCO')
  
  
      else:
          print('Superfamily: RLG')
          result.append('Superfamily: RLG')          
          model_8 = joblib.load('./Model/Kmer(64)_8_cab')
          pre_8 = model_8.predict(feature)
          if pre_8 >= 0.5:
              print('Lineage: TAT')
              result.append('Lineage: TAT')                        
          else:
              model_9 = joblib.load('./Model/Kmer(64)_9_knn')
              pre_9 = model_9.predict(feature)
              if pre_9 >= 0.5:
                  print('Lineage: TEKAY/DEL')
                  result.append('Lineage: TEKAY/DEL')
              else:
                  Model0 = joblib.load('./Model/Kmer(64)_10_knn')
                  pre_10 = Model0.predict(feature)
                  if pre_10 >= 0.5:
                      Model1 = joblib.load('./Model/Kmer(64)_11_cab')
                      pre_11 = Model1.predict(feature)
                      if pre_11 >= 0.5:
                          print('Lineage: ATHILA')
                          result.append('Lineage: ATHILA')                                                                    
                      else:
                          print('Lineage: CRM')
                          result.append('Lineage: CRM')                                                                                              
                  else:
                      Model2 = joblib.load('./Model/Kmer(64)_12_cab_SMOTH')
                      pre_12 = Model2.predict(feature)
                      if pre_12 >= 0.5:
                          print('Lineage: GALADRIEL')
                          result.append('Lineage: GALADRIEL')                          
                      else:
                          print('Lineage: REINA')
                          result.append('Lineage: REINA')                                                    

  else:
      print('nLTR-RT')
      result.append('nLTR-RT')
  print(result)      
  total_result.append(result)
  
# print(total_result)


with open(r'result_class1.txt', 'w') as f:
    for i in total_result:
        f.write(str(i) + '\n')







