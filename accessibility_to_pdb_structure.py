import sys

fhi = open(sys.argv[1])
cluster_choice = sys.argv[2]

print("attribute: percentAcc")
print("match mode: 1-to-1")
print("recipient: residues")

cnt_f = 1
cnt_r = 147

for line in fhi:
    split = line.split()
    idx = int(split[0])
    val = float(split[1])
    cluster = split[2]
    if idx < -76 or idx > 74:
        continue
    if cluster != cluster_choice:
        continue
    print("\t/J:%s\t%s" % (cnt_f, val))
    print("\t/I:%s\t%s" % (cnt_r, val))
    cnt_f += 1
    cnt_r -= 1

fhi.close()
