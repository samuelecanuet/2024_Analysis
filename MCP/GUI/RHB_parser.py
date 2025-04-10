lines=[]
with open("../fit_params_2025.txt", 'r') as file:
    for line in file : 
        print(line)
        try:
            lines.append([float(i) for i in line.split(' ') if i != "\n"])
        except:
            break

bin = 400
#### Spacial params
print("FCondition;X_IN;(X_norm_log<[0])&&(X_norm_log>[1])")
print("[0];0.4")
print("[1];-0.4")

print("FCondition;Y_IN;(Y_norm_log<[0])&&(Y_norm_log>[1])")
print("[0];0.4")
print("[1];-0.4")

print("FCondition;SpaceCondition;X_IN && Y_IN")


####CalParam
print("#")
print("# X correction")
print("#")
#0X
print(f"FParamCalc;a;[0]+[1]*Y_norm_log+[2]*Y_norm_log*Y_norm_log+[3]*Y_norm_log*Y_norm_log*Y_norm_log")
print(f"[0];{lines[0][0]}")
print(f"[1];{lines[0][1]}")
print(f"[2];{lines[0][2]}")
print(f"[3];{lines[0][3]}")
print("#")
#1X
print(f"FParamCalc;b;[0]*X_norm_log+[1]*X_norm_log*Y_norm_log+[2]*X_norm_log*Y_norm_log*Y_norm_log+[3]*X_norm_log*Y_norm_log*Y_norm_log*Y_norm_log")
print(f"[0];{lines[0][4]}")
print(f"[1];{lines[0][5]}")
print(f"[2];{lines[0][6]}")
print(f"[3];{lines[0][7]}")
print("#")
#2X
print(f"FParamCalc;c;[0]*X_norm_log*X_norm_log+[1]*X_norm_log*X_norm_log*Y_norm_log+[2]*X_norm_log*X_norm_log*Y_norm_log*Y_norm_log+[3]*X_norm_log*X_norm_log*Y_norm_log*Y_norm_log*Y_norm_log")
print(f"[0];{lines[0][8]}")
print(f"[1];{lines[0][9]}")
print(f"[2];{lines[0][10]}")
print(f"[3];{lines[0][11]}")
print("#")
#3X
print(f"FParamCalc;d;[0]*X_norm_log*X_norm_log*X_norm_log+[1]*X_norm_log*X_norm_log*X_norm_log*Y_norm_log+[2]*X_norm_log*X_norm_log*X_norm_log*Y_norm_log*Y_norm_log+[3]*X_norm_log*X_norm_log*X_norm_log*Y_norm_log*Y_norm_log*Y_norm_log")
print(f"[0];{lines[0][12]}")
print(f"[1];{lines[0][13]}")
print(f"[2];{lines[0][14]}")
print(f"[3];{lines[0][15]}")
print("#")
print("#")
print("# Y correction")
print("#")
#0X
print(f"FParamCalc;h;[0]+[1]*Y_norm_log+[2]*Y_norm_log*Y_norm_log+[3]*Y_norm_log*Y_norm_log*Y_norm_log")
print(f"[0];{lines[1][0]}")
print(f"[1];{lines[1][1]}")
print(f"[2];{lines[1][2]}")
print(f"[3];{lines[1][3]}")
print("#")
#1X
print(f"FParamCalc;k;[0]*X_norm_log+[1]*X_norm_log*Y_norm_log+[2]*X_norm_log*Y_norm_log*Y_norm_log+[3]*X_norm_log*Y_norm_log*Y_norm_log*Y_norm_log")
print(f"[0];{lines[1][4]}")
print(f"[1];{lines[1][5]}")
print(f"[2];{lines[1][6]}")
print(f"[3];{lines[1][7]}")
print("#")
#2X
print(f"FParamCalc;l;[0]*X_norm_log*X_norm_log+[1]*X_norm_log*X_norm_log*Y_norm_log+[2]*X_norm_log*X_norm_log*Y_norm_log*Y_norm_log+[3]*X_norm_log*X_norm_log*Y_norm_log*Y_norm_log*Y_norm_log")
print(f"[0];{lines[1][8]}")
print(f"[1];{lines[1][9]}")
print(f"[2];{lines[1][10]}")
print(f"[3];{lines[1][11]}")
print("#")
#3X
print(f"FParamCalc;m;[0]*X_norm_log*X_norm_log*X_norm_log+[1]*X_norm_log*X_norm_log*X_norm_log*Y_norm_log+[2]*X_norm_log*X_norm_log*X_norm_log*Y_norm_log*Y_norm_log+[3]*X_norm_log*X_norm_log*X_norm_log*Y_norm_log*Y_norm_log*Y_norm_log")
print(f"[0];{lines[1][12]}")
print(f"[1];{lines[1][13]}")
print(f"[2];{lines[1][14]}")
print(f"[3];{lines[1][15]}")
print("FParamCalc;X_rec;a+b+c+d")
print("FParamCalc;Y_rec;h+k+l+m")
print("#")
print("FH2F;anode_image_rec;anode_image_rec")
print("Condition:SpaceCondition")
print(f"X_rec;{bin};-8;8")
print(f"Y_rec;{bin};-8;8")
print("#")
print("FH1F;proj_x_rec;proj_x_rec")
print(f"Condition:SpaceCondition")
print(f"X_rec;{bin};-8;8")
print("#")
print("FH1F;proj_y_rec;proj_y_rec")
print(f"Condition:SpaceCondition")
print(f"Y_rec;{bin};-8;8")
