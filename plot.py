import vcf
import sys
import matplotlib.pyplot as plt
import pandas as pd

vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))
total = 0
tp = 0
fp = 0
p = 0
n = 0
NIA = []
NIR = []
TP_ID = []

NIA_1 = []
NIR_1 = []
FP_ID = []

f = open("result.tsv", "w")
line = "ID" +"\t" +"TP/FP" + "\t" + "alt support" + "\t" + "ref support"
f.write(line)
f.write("\n")
tpid = 0
fpid = 0
#TF = sys.argv[2]
for record in vcf_reader:
    total += 1
    x = int(record.INFO['HG2count'])
    y = 0
    for sample in record.samples:
        y = int(sample['NIA'])
        z = int(sample['NIR'])
    if x == 0:
        n += 1
        if y > 0:
            f.write(record.ID + "\tFP\t" + str(y) + "\t" + str(z))
            f.write("\n")
            if y > 50000 or z > 50000:
                print(record.ID, y, z)
            #if TF == "TP":
            fpid += 1
            FP_ID.append(fpid) 
            NIA_1.append(y)
            NIR_1.append(z)
            #color.append("FP")
            fp += 1
    else:
        p += 1
        if y > 0:
            f.write(record.ID + "\tTP\t" + str(y) + "\t" + str(z))
            f.write("\n")
            if y > 50000 or z > 50000:
                print(record.ID, y, z)
            #if TF == "FP":
            tpid += 1
            TP_ID.append(tpid)
            NIA.append(y)
            NIR.append(z)
            #color.append("TP")
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
plot1 = plt.figure(1)
plt.scatter(NIA, NIR)
plt.title("alt support vs ref support for TPs")
plt.xlabel("alt support")
plt.ylabel("ref support")
plt.savefig('NIA_vs_NIR.TP.png')
plot1 = plt.figure(2)
plt.scatter(NIA_1, NIR_1)
plt.title("alt support vs ref support for FPs")
plt.xlabel("alt support")
plt.ylabel("ref support")
plt.savefig('NIA_vs_NIR.FP.png')
'''
plot3 = plt.figure(3)
plt.plot(NIA,'g*', NIR, 'ro')
plt.savefig('NIA_vs_NIR.TP.line.png')
plot4 = plt.figure(4)
label=['NIA', 'NIR']
plt.plot(NIA_1,'g*', NIR_1, 'ro', label=label)
plt.legend()
plt.savefig('NIA_vs_NIR.FPP.line.png')
'''
plot3 = plt.figure(3)
#plt.xlim(-50, 1000)
plt.ylim(-50, 1000)
plt.title("alt support vs ref support for TPs")
plt.plot(TP_ID, NIA, marker='*', color = 'red', label='NIA')
plt.plot(TP_ID, NIR, marker='o', color = 'green',label='NIR')
plt.legend()
plt.xlabel("# TPs")
plt.ylabel("count")
plt.savefig('NIA_vs_NIR.TP.line.zoom.png')

plot4 = plt.figure(4)
#plt.xlim(-50, 1000)
plt.title("alt support vs ref support for FPs")
plt.ylim(-50, 1000)
plt.plot(FP_ID, NIA_1, marker='*', color = 'red', label='NIA')
plt.plot(FP_ID, NIR_1, marker='o', color = 'green', label='NIR')
plt.legend()
plt.xlabel("# FPs")
plt.ylabel("count")
plt.savefig('NIA_vs_NIR.FP.line.zoom.png')

plot5 = plt.figure(5)
plt.title("alt support vs ref support for TPs")
plt.plot(TP_ID, NIA, label='alt support')
plt.plot(TP_ID, NIR, label='ref support')
plt.legend()
plt.savefig('NIA_vs_NIR.TP.line1.png')
plot6 = plt.figure(6)
plt.title("alt support vs ref support for FPs")
plt.plot(FP_ID, NIA_1, label='alt support')
plt.plot(FP_ID, NIR_1, label='ref support')
plt.legend()
plt.savefig('NIA_vs_NIR.FP.line1.png')

plot7 = plt.figure(7)
plt.title("alt support vs ref support for TPs")
label=['alt support: green, ref support: red']
plt.plot(NIA,'g*', NIR, 'ro', label=label)
plt.xlabel("# TPs")
plt.ylabel("count")
plt.legend()
plt.savefig('NIA_vs_NIR.TP.point.png')

plot8 = plt.figure(8)
plt.title("alt support vs ref support for FPs")
label=['alt support: green, ref support: red']
plt.plot(NIA_1,'g*', NIR_1, 'ro', label=label)
plt.xlabel("# FPs")
plt.ylabel("count")
plt.legend()
plt.savefig('NIA_vs_NIR.FP.point.png')


plot9 = plt.figure(9)
plt.title("alt support vs ref support for TPs(<1000)")
plt.ylim(-50, 1000)
label=['alt support: green, ref support: red']
plt.plot(NIA,'g*', NIR, 'ro', label=label)
plt.xlabel("# TPs")
plt.ylabel("count")
plt.legend()
plt.savefig('NIA_vs_NIR.TP.point.zoom.png')

plot10 = plt.figure(10)
plt.title("alt support vs ref support for FPs(<1000)")
plt.ylim(-50, 1000)
label=['alt support: green, ref support: red']
plt.plot(NIA_1,'g*', NIR_1, 'ro', label=label)
plt.xlabel("# FPs")
plt.ylabel("count")
plt.legend()
plt.savefig('NIA_vs_NIR.FP.point.zoom.png')
#plt.xlim(-2, 1000)
#plt.ylim(-2, 5000)
#plt.xscale('log')
#plt.yscale('log')

plt.subplot(1, 2, 1)
plt.scatter(NIA, NIR)
plt.title("alt vs ref support for TPs")
plt.xlabel("alt support")
plt.ylabel("ref support")
plt.xlim(-50, 1000)
plt.ylim(-50, 1000)
#plt.xscale('log')
#plt.yscale('log')

#plt.savefig('NIA_vs_NIR.TP.png')
plt.subplot(1, 2, 2)
plt.scatter(NIA_1, NIR_1)
plt.title("alt vs ref support for FPs")
plt.xlabel("alt support")
#plt.ylabel("ref support")
plt.xlim(-50, 1000)
plt.ylim(-50, 1000)
plt.savefig('NIA_vs_NIR.png')

f.close()
