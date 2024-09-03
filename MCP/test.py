chi2 = []
p = []
p1 = []
p2 = []
p3 = []

with open ( "chi2.txt", "r", encoding="utf-8" ) as f:
    for line in f:
        liste = line.split(" ")
        try:
            chi2.append(float(liste[0]))
        except:
            continue
        p.append(float(liste[1]))
        p1.append(float(liste[2]))
        p2.append(float(liste[3]))
        p3.append(float(liste[4]))


# sort chi2 and sort p with this sorting
chi2_sorted = sorted(chi2)
p_sorted = [x for _,x in sorted(zip(chi2,p))]

# sort p1, p2, p3 with this sorting
p1_sorted = [x for _,x in sorted(zip(chi2,p1))]
p2_sorted = [x for _,x in sorted(zip(chi2,p2))]
p3_sorted = [x for _,x in sorted(zip(chi2,p3))]

for i in range(10):
    print(chi2_sorted[i], p_sorted[i], p1_sorted[i], p2_sorted[i], p3_sorted[i])