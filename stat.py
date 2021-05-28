import vcf
import sys

vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))
total = 0
tp = 0
fp = 0
p = 0
n = 0
for record in vcf_reader:
    total += 1
    x = int(record.INFO['HG2count'])
    y = 0
    for sample in record.samples:
        y = int(sample['NIA'])
    if x == 0:
        n += 1
        if y > 0:
            fp += 1
    else:
        p += 1
        if y > 0:
            tp += 1
fn = p - tp
tn = n - fp
print("Tot: " + str(total) + " Pos: " + str(p) + " Neg: " + str(n) + " TP: " + str(tp) + " FP: " + str(fp) + " FN: " + str(fn) + " TN: " + str(tn))
recall = float(tp)/float(tp+fn)
precision = float(tp)/(float)(tp+fp)
accuracy = float(tp+fn)/float(tp+fp+fn+tn)
fdr = float(fp)/(float)(tp+fp)
print("Recall: ", round(recall, 2))
print("Precision: ", round(precision, 2))
print("Accuracy: ", round(accuracy, 2))
print("False discovery rate: ", round(fdr, 2))
