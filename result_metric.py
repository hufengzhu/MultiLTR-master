from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import classification_report


def metrics(Y_validation,predictions):
    print('Accuracy:', accuracy_score(Y_validation, predictions))
    print('F1 score:', f1_score(Y_validation, predictions,average='weighted'))
    print('Recall:', recall_score(Y_validation, predictions,average='weighted'))
    print('Precision:', precision_score(Y_validation, predictions, average='weighted'))
    
    print('\n clasification report:\n', classification_report(Y_validation, predictions))
    print('\n confusion matrix:\n',confusion_matrix(Y_validation, predictions))
    # To create confusion matrix
    confusion_matrix(Y_validation, predictions)


lineages_names_dic = {'nLTR-RT':0,
                      'LTR/RLC/ALE/RETROFIT':1,
                      'LTR/RLC/ANGELA':2,
                      'LTR/RLC/BIANCA':3,
                      'LTR/RLC/IKEROS':4,
                      'LTR/RLC/IVANA/ORYCO':5,
                      'LTR/RLC/SIRE':6,
                      'LTR/RLC/TAR/TORK':7,
                      'LTR/RLG/ATHILA':8,
                      'LTR/RLG/CRM':9,
                      'LTR/RLG/TEKAY/DEL':10,
                      'LTR/RLG/GALADRIEL':11,
                      'LTR/RLG/REINA':12,
                      'LTR/RLG/TAT':13}


true_path = "Inpactor2_library.fasta"
predic_path = "result_class.txt"



from Bio import SeqIO

true_list = []
for seq_record in SeqIO.parse("Inpactor2_library.fasta", "fasta"):
    na = seq_record.id
    # se = seq_record.seq
    b = na.split('_')
    true_ = b[3].split('#')[1].replace("-", "/")
    # print(lineages_names_dic.get(true_))
    true_list.append(lineages_names_dic.get(true_))


predict_list = []
with open(predic_path, 'r') as f:
    samples = f.readlines()
    for sample in samples:
        splited_list = sample[1:-2].split(",")
        if splited_list[0].split("'")[1] == 'nLTR-RT':
            predic_ = splited_list[0].split("'")[1]
        else:
            predic_ = splited_list[0].split("'")[1].split("-")[0] +"/"+splited_list[1].split("'")[1].split()[-1]+"/"+splited_list[2].split("'")[1].split()[-1].replace("-", "/")
        predict_list.append(lineages_names_dic.get(predic_))



metrics(true_list, predict_list)

