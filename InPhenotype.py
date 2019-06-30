import os
import random
import sys

import numpy as np
from scipy import stats

from IP_Classes import Forest
from IP_Classes import Ind
from IP_Classes import SNP
from IP_Classes import Tree

global Total_Leaves_Num
global Ind_Threshold
global Imprv_Threshold
global Category1
global Category2
Tuple_Error_Threshold = 5


def InPhenotype(Train_Path, Test_Path, Only_Prediction, Gene_Num, SNP_Num, Ind_Threshold_Given, Imprv_Threshold_Given, Tree_Num, Output_Path, Repr_SNPs_Num, Is_Weight, Final_Gene_Num = 1):

    global Total_Leaves_Num
    global Ind_Threshold
    global Imprv_Threshold
    global Category1
    global Category2

    random.seed(4321)

    Gene2Values = {}
    Gene2Forest = {}
    Perm_Gene2Values = {}
    Permuted_Gene2Forest = {}
    Gene2Pvalue = {}
    Gene2Scores_Lists = {}

    Forests = []
    Permuted_Forests = []
    List_Of_SNPs=[]
    Ind_Array=[]
    Perm_Ind_Array = []
    Total_Trees_List = []
    Total_Perm_Trees_List = []
    Pred_Inds = []

    Temp_SNP=SNP()
    Temp_Ind= Ind()

    Valid = Check_Args(Train_Path, Test_Path, Only_Prediction, Gene_Num, SNP_Num, Ind_Threshold_Given, Imprv_Threshold_Given, Output_Path, Final_Gene_Num)
    if(Valid == 0):
        exit()

    Imprv_Threshold = 1- Imprv_Threshold_Given
    Ind_Threshold = Ind_Threshold_Given

    Train_Data = open(Train_Path, 'r').readlines()
    Test_Data = open(Test_Path, 'r').readlines()
    Test_Output_File = open(Output_Path + "IP_Test_Results.txt", 'w')

    Genes_List = [row.rstrip().split("\t")[0] for row in Train_Data[(1 + SNP_Num):(1 + SNP_Num + Gene_Num)]]
    Names = ["Ind" + (str)(i+1) for i in range(len(Train_Data[0].rstrip().split("\t")[1:]))]
    for row in Train_Data[(1 + SNP_Num):(1 + SNP_Num + Gene_Num)]:
        Data = row.rstrip().split("\t")
        Vals_Across_Inds = []
        Gene = Data[0]
        for val in Data[1:]:
            Vals_Across_Inds.append((float)(val))
        Gene2Values[Gene] = Vals_Across_Inds

    Phenotype = Convert2Float(Train_Data[0].rstrip().split("\t")[1:])

    for i in range(1,(1 + SNP_Num)):
        Curr_SNP_Data = Train_Data[i].rstrip().split("\t")
        Is_Valid, Opt_Cat1, Opt_Cat2 = Check_SNP(Curr_SNP_Data, 1, Ind_Threshold)
        if(Is_Valid == 0):
            print("Error - There are more than two SNP values")
            exit()
        elif(Is_Valid == 1):
            continue
        Category1 = Opt_Cat1
        Category2 = Opt_Cat2
        Temp_SNP.Name = Curr_SNP_Data[0]
        Temp_SNP.Chromosome = Curr_SNP_Data[1]
        Temp_SNP.Position = Curr_SNP_Data[2]
        for j in range(1,len(Curr_SNP_Data)):
            Temp_SNP.IndMap[Names[j-1]] = Curr_SNP_Data[j]
        List_Of_SNPs.append(Temp_SNP)
        Temp_SNP = SNP()

    for i in range(len(Names)):
        GE_Vals = {}
        Temp_Ind.Name = Names[i]
        Temp_Ind.Phenotype = Phenotype[i]
        for key in Gene2Values:
            GE_Vals[key] = Gene2Values[key][i]
        Temp_Ind.GEmap = GE_Vals
        Ind_Array.append(Temp_Ind)
        Temp_Ind = Ind()

    Repr_SNPs_Num = min(len(List_Of_SNPs), Repr_SNPs_Num)
    SNPs_Indices = range(0, len(List_Of_SNPs))

    for j in range(Tree_Num):

        Chosen_Samples = []
        Tree_Inds = []

        Curr_Tree = Tree()
        Curr_Tree.Errors = []
        Curr_Tree.Total_Error = 0

        Chosen_Gene = Genes_List[random.randint(0,len(Genes_List)-1)]

        for k in range(len(Names)):
            Chosen_Sample_Idx = random.randint(0,len(Names) - 1)
            Chosen_Sample = Names[Chosen_Sample_Idx]
            Chosen_Samples.append(Chosen_Sample)
        for samp in Chosen_Samples:
            for Curr_Ind in Ind_Array:
                if(Curr_Ind.Name == samp):
                    Tree_Inds.append(Curr_Ind)

        Curr_Tree.Unused_Samples = [sample for sample in Ind_Array if sample.Name not in Chosen_Samples]

        random.shuffle(SNPs_Indices)
        SNP_Reprs = [List_Of_SNPs[SNPs_Indices[i]] for i in range(Repr_SNPs_Num)]

        Gene_Across_Inds,Phen_Across_Inds = Get_Gene_BioData(Tree_Inds, Chosen_Gene)

        if(len(set(Gene_Across_Inds))) < Ind_Threshold:
                continue

        Root_SSR = Get_Reg_Measure(Gene_Across_Inds, Phen_Across_Inds, "SSR")
        Root_SST = Get_Reg_Measure(Gene_Across_Inds, Phen_Across_Inds, "SST")
        Curr_Tree.SSR = Root_SSR
        Curr_Tree.SST = Root_SST
        Curr_Tree.Primary_Params = Get_Reg_Measure(Gene_Across_Inds, Phen_Across_Inds)

        Root = Create_Branch(None, SNP_Reprs, Tree_Inds, (float)(Root_SSR)/Root_SST, Chosen_Gene)

        Curr_Tree.Root = Root
        Curr_Tree.Gene = Chosen_Gene
        Curr_Tree.Uniqe_Inds = len(set(Tree_Inds))
        Curr_Tree.Leaves_Num = Get_Leaves_Num(Curr_Tree.Root, True)

        Total_Trees_List.append(Curr_Tree)

    for tree in Total_Trees_List:
        Curr_Gene = tree.Gene
        if(Curr_Gene in Gene2Forest):
            Curr_Forest = Gene2Forest[Curr_Gene]
            Curr_Forest.Tree_List.append(tree)
        else:
            New_Forest = Forest()
            New_Forest.Gene = Curr_Gene
            New_Forest.Tree_List = [tree]
            Gene2Forest[Curr_Gene] = New_Forest
    for key in Gene2Forest.keys():
        Forests.append(Gene2Forest[key])

    for i in range(len(Ind_Array)):

        for tree in Total_Trees_List:
            if(Ind_Array[i] not in tree.Unused_Samples):
                continue
            Gene_Name = tree.Gene
            if(tree.Leaves_Num == 1):
                Prediction = tree.Primary_Params[0] * Ind_Array[i].GEmap[Gene_Name] + tree.Primary_Params[1]
                tree.Errors.append((Prediction, Ind_Array[i].Phenotype))
            else:
                Predict(tree.Root, Ind_Array[i], Gene_Name, Ind_Array[i].Phenotype)

    for tree in Total_Trees_List:
        if(tree.Leaves_Num == 1):
            Errors = tree.Errors
            Total_Error = 0
            for tup in Errors:
                Total_Error += (tup[0] - tup[1]) ** 2
            Total_Error /= len(Errors)
            tree.Total_Error = Total_Error ** 0.5
        else:
            Calc_Regerssion_Error(tree.Root)

    for tree in Total_Trees_List:
        if(tree.Leaves_Num == 1):
            tree.Score = 1.0/(tree.Total_Error)
        else:
            Inds_Tuples = (Agg_Tree_PR_Tuples(tree.Root, []))
            Error = 0
            for tup in Inds_Tuples:
                Error += (tup[0] - tup[1]) ** 2
            Error = (float)(Error)/(len(Inds_Tuples))
            Error = Error ** 0.5
            tree.Score = 1.0/Error

    Sorted_Trees = sorted(Total_Trees_List, key = lambda x:x.Score)[::-1]
    Print_Trees_To_File(Sorted_Trees, Output_Path, False)

    for row in Train_Data[(1 + SNP_Num):(1 + SNP_Num + Gene_Num)]:
        Data = row.rstrip().split("\t")
        Vals_Across_Inds = []
        Gene = Data[0]
        for val in Data[1:]:
            Vals_Across_Inds.append((float)(val))
        random.shuffle(Vals_Across_Inds)
        Perm_Gene2Values[Gene] = Vals_Across_Inds

    Phenotype_Permuted = Convert2Float(Train_Data[0].rstrip().split("\t")[1:])
    random.shuffle(Phenotype_Permuted)

    for i in range(len(Names)):
        Temp_Ind.Name = Names[i]
        Temp_Ind.Phenotype = Phenotype_Permuted[i]
        GE_Vals = {}
        for key in Perm_Gene2Values:
            GE_Vals[key] = Perm_Gene2Values[key][i]
        Temp_Ind.GEmap = GE_Vals
        Perm_Ind_Array.append(Temp_Ind)
        Temp_Ind = Ind()

    for j in range(Tree_Num):

        Chosen_Samples = []
        Tree_Inds = []

        Total_Leaves_Num = 0

        Curr_Tree = Tree()
        Curr_Tree.Errors = []
        Curr_Tree.Total_Error = 0

        Chosen_Gene = Genes_List[random.randint(0, len(Genes_List) - 1)]

        for k in range(len(Names)):
            Chosen_Sample_Idx = random.randint(0,len(Names) -1)
            Chosen_Sample = Names[Chosen_Sample_Idx]
            Chosen_Samples.append(Chosen_Sample)
        for samp in Chosen_Samples:
            for Curr_Ind in Perm_Ind_Array:
                if(Curr_Ind.Name == samp):
                    Tree_Inds.append(Curr_Ind)

        Curr_Tree.Unused_Samples = [sample for sample in Perm_Ind_Array if sample.Name not in Chosen_Samples]

        random.shuffle(SNPs_Indices)
        SNP_Reprs = [List_Of_SNPs[SNPs_Indices[i]] for i in range(Repr_SNPs_Num)]

        Gene_Across_Inds,Phen_Across_Inds = Get_Gene_BioData(Tree_Inds, Chosen_Gene)

        if(len(set(Gene_Across_Inds))) < Ind_Threshold:
                continue
        Root_SSR = Get_Reg_Measure(Gene_Across_Inds, Phen_Across_Inds, "SSR")
        Root_SST = Get_Reg_Measure(Gene_Across_Inds, Phen_Across_Inds, "SST")
        Curr_Tree.SSR = Root_SSR
        Curr_Tree.SST = Root_SST
        Curr_Tree.Primary_Params = Get_Reg_Measure(Gene_Across_Inds, Phen_Across_Inds)

        Root = Create_Branch(None, SNP_Reprs, Tree_Inds, (float)(Root_SSR)/Root_SST, Chosen_Gene)

        Curr_Tree.Root = Root
        Curr_Tree.Gene = Chosen_Gene
        Curr_Tree.Uniqe_Inds = len(set(Tree_Inds))
        Curr_Tree.Leaves_Num = Get_Leaves_Num(Curr_Tree.Root, True)

        Total_Perm_Trees_List.append(Curr_Tree)

    for tree in Total_Perm_Trees_List:
        Curr_Gene = tree.Gene
        if(Curr_Gene in Permuted_Gene2Forest):
            Curr_Forest = Permuted_Gene2Forest[Curr_Gene]
            Curr_Forest.Tree_List.append(tree)
        else:
            New_Forest = Forest()
            New_Forest.Gene = Curr_Gene
            New_Forest.Tree_List = [tree]
            Permuted_Gene2Forest[Curr_Gene] = New_Forest
    for key in Permuted_Gene2Forest.keys():
        Permuted_Forests.append(Permuted_Gene2Forest[key])

    for i in range(len(Perm_Ind_Array)):

        for tree in Total_Perm_Trees_List:
            if(Perm_Ind_Array[i] not in tree.Unused_Samples):
                continue
            Gene_Name = tree.Gene
            if(tree.Leaves_Num == 1):
                Prediction = tree.Primary_Params[0] * Perm_Ind_Array[i].GEmap[Gene_Name] + tree.Primary_Params[1]
                tree.Errors.append((Prediction, Perm_Ind_Array[i].Phenotype))
            else:
                Predict(tree.Root, Perm_Ind_Array[i], Gene_Name, Perm_Ind_Array[i].Phenotype)

    for tree in Total_Perm_Trees_List:
        if(tree.Leaves_Num == 1):
            Errors = tree.Errors
            Total_Error = 0
            for tup in Errors:
                Total_Error += (tup[0] - tup[1]) ** 2
            Total_Error /= len(Errors)
            tree.Total_Error = Total_Error ** 0.5
        else:
            Calc_Regerssion_Error(tree.Root)

    for tree in Total_Perm_Trees_List:
        if(tree.Leaves_Num == 1):
            tree.Score = 1.0/(tree.Total_Error)
        else:
            Inds_Tuples = (Agg_Tree_PR_Tuples(tree.Root, []))
            Error = 0
            for tup in Inds_Tuples:
                Error += (tup[0] - tup[1]) ** 2
            Error = (float)(Error)/(len(Inds_Tuples))
            Error = Error ** 0.5
            tree.Score = 1.0/Error

    Sorted_Perm_Trees = sorted(Total_Perm_Trees_List, key = lambda x:x.Score)[::-1]
    Print_Trees_To_File(Sorted_Perm_Trees, Output_Path, True)

    Permuted_Gene2Scores = Extract_Gene2Score(Total_Perm_Trees_List)
    Real_Gene2Scores = Extract_Gene2Score(Total_Trees_List)
    Relevant_Genes = Unify(Permuted_Gene2Scores.keys(), Real_Gene2Scores.keys())
    for rg in Relevant_Genes:
        try:
            Curr_Real_Scores = Real_Gene2Scores[rg]
        except  KeyError:
            Curr_Real_Scores = []
        try:
            Curr_Permuted_Scores = Permuted_Gene2Scores[rg]
        except  KeyError:
            Curr_Permuted_Scores = []
        Gene2Scores_Lists[rg] = (Curr_Real_Scores, Curr_Permuted_Scores)
    for key in Gene2Scores_Lists.keys():
        S_Lists = Gene2Scores_Lists[key]
        if(len(S_Lists[0]) == 0 or len(S_Lists[1]) == 0):
            Gene2Pvalue[key] = (float)("nan")
            continue
        temp, Gene2Pvalue[key] = stats.mannwhitneyu(S_Lists[0], S_Lists[1], use_continuity = False, alternative = "greater")

    Gene_Scores_Output_File = open(Output_Path + "Gene_Scores.txt", 'w')
    for key in Gene2Pvalue.keys():
        Gene_Scores_Output_File.write(key + "\t" + (str)(Gene2Pvalue[key]) + "\n")
    Gene_Scores_Output_File.close()

    for forest in Forests:
        Gene = forest.Gene
        forest.Score = (float)(Gene2Pvalue[Gene])

    Forests = sorted(Forests, key= lambda x:x.Score)

    Gene2Values = {}
    SNP2Values = {}

    if(Only_Prediction):
        Feature_Start_Idx = 0
        Test_Output_File.write("Individual\tPrediction\n")
    else:
        Test_Output_File.write("Individual\tPhenotype\tPrediction\tError\n")
        Known_Phenotypes = Test_Data[0].rstrip().split("\t")[1:]
        Feature_Start_Idx = 1
        Error_Tuples = []

    for row in Test_Data[Feature_Start_Idx:(SNP_Num + Feature_Start_Idx)]:
        Features = row.rstrip().split("\t")
        Vals_Across_Inds = []
        Snip = Features[0]
        for val in Features[1:]:
            Vals_Across_Inds.append((int)(val))
        SNP2Values[Snip] = Vals_Across_Inds

    for row in Test_Data[(SNP_Num + Feature_Start_Idx): (SNP_Num + Feature_Start_Idx + Gene_Num)]:
        Features = row.rstrip().split("\t")
        Vals_Across_Inds = []
        Gene = Features[0]
        for val in Features[1:]:
            Vals_Across_Inds.append((float)(val))
        Gene2Values[Gene] = Vals_Across_Inds

    Pred_Inds_Names = ["Ind" + (str)(i+1) for i in range(len(Test_Data[0].rstrip().split("\t")[1:]))]

    for i in range(len(Pred_Inds_Names)):
        Curr_Ind = Ind()
        GE_Vals = {}
        SNP_Vals = {}
        Curr_Ind.Name = Pred_Inds_Names[i]
        if(Only_Prediction):
            Curr_Ind.Phenotype = None
        else:
            Curr_Ind.Phenotype = Known_Phenotypes[i]
        for key in Gene2Values:
            GE_Vals[key] = Gene2Values[key][i]
        Curr_Ind.GEmap = GE_Vals
        for key in SNP2Values:
            SNP_Vals[key] = SNP2Values[key][i]
        Curr_Ind.SNPmap = SNP_Vals
        Pred_Inds.append(Curr_Ind)

    for Curr_Ind in Pred_Inds:

        Test_Output_File.write(Curr_Ind.Name + "\t")

        Total_Weight = 0.0
        Total_Prediction = 0.0
        for i in range(Final_Gene_Num):
            Prediction_Entries = []
            for tree in sorted(Forests[i].Tree_List, key=lambda x:x.Score)[::-1]:
                Gene_Name = tree.Gene
                Prediction_Values = Predict_By_Tree(tree.Root, Curr_Ind, Gene_Name, tree.Primary_Params, True, tree.Total_Error)
                if(Prediction_Values[0] == -1):
                    continue
                if(Prediction_Values[1] == 0):
                    Prediction_Entries.append([Prediction_Values[0], tree.Score, 0.0, tree.Leaves_Num])
                else:
                    Prediction_Entries.append([Prediction_Values[0], tree.Score, (float)(1/Prediction_Values[1]), tree.Leaves_Num])
            if(Is_Weight == "1"):
                for ent in Prediction_Entries:
                    Total_Weight += ent[2]
                    Total_Prediction += ent[0] * ent[2]
            else:
                for ent in Prediction_Entries:
                    Total_Weight += 1
                    Total_Prediction += ent[0]
        try:
            Final_Prediction = (float)(Total_Prediction)/(Total_Weight)
        except ZeroDivisionError:
            Final_Prediction = ((float)("inf"))
        if(Only_Prediction):
            Test_Output_File.write((str)(Final_Prediction) + "\n")
        else:
            Test_Output_File.write((str)(Curr_Ind.Phenotype) + "\t" + (str)(Final_Prediction) + "\t")
            Test_Output_File.write((str)(abs((float)(Curr_Ind.Phenotype) - Final_Prediction)) + "\n")
            Error_Tuples.append(((float)(Curr_Ind.Phenotype), Final_Prediction))

    if(not Only_Prediction):
        Total_Error = 0.0
        for tup in Error_Tuples:
            Total_Error += ((tup[0] - tup[1]) ** 2)
        Total_Error = ((float)(Total_Error)/len(Error_Tuples)) ** 0.5
        Test_Output_File.write("Test RMSE is: \t" + (str)(Total_Error))

    Test_Output_File.close()

    return


def Get_Gene_BioData(Ind_Array, Gene_Name):

    x = []
    y = []

    for i in range(len(Ind_Array)):
        x.append(Ind_Array[i].GEmap[Gene_Name])
        y.append(Ind_Array[i].Phenotype)

    return (x,y)


def Get_Reg_Measure(x, y, SST_Or_SSR = None):

    Tuple_x_y = []
    x_new = []
    y_new = []

    for i in range(len(x)):
        Tuple_x_y.append((x[i], y[i]))
    Tuple_x_y = set(Tuple_x_y)
    for tup in Tuple_x_y:
        x_new.append(tup[0])
        y_new.append(tup[1])
    x_new = np.array(x_new)
    y_new = np.array(y_new)

    if(len(x_new) != len(y_new)):
        print("Data problem - number of gene entries different than number of phenotype entries")
        exit()

    y_Avg = (float)(np.sum(y_new))/len(y_new)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_new,y_new)
    Predicted_y = intercept + slope * x_new
    Pred_Error = y_new - Predicted_y
    Avg_Error = y_new - y_Avg
    SSE = np.sum(Pred_Error**2)
    SST = np.sum(Avg_Error**2)
    SSR = SST - SSE

    if(SST_Or_SSR == None):
        return [slope,intercept]
    if(SST_Or_SSR == "SST"):
        return SST
    elif(SST_Or_SSR == "SSR"):
        return SSR


def Create_Branch(Parent, List_Of_SNPs, Inds, Parent_Measure, Gene_Name):

    global Ind_Threshold
    global Imprv_Threshold
    global OT_Global
    global Category1
    global Category2

    Chosen_SNP = SNP()
    Results = []

    Selected_SNPs = list(set(List_Of_SNPs))
    SNP_subset_num = (int)(2 * round(len(List_Of_SNPs)**0.5))
    for i in range(len(List_Of_SNPs) - SNP_subset_num):
        Selected_SNPs.pop(random.randint(0,len(Selected_SNPs)-1))

    for i in range(len(Selected_SNPs)):

        Curr_Ind_Left = []
        Curr_Ind_Right = []
        Reg1_x=[]
        Reg2_x=[]
        Reg1_y=[]
        Reg2_y=[]

        Curr_SNP = Selected_SNPs[i]

        for j in range(len(Inds)):
            Name = Inds[j].Name
            Origin = Curr_SNP.IndMap[Name]
            if(Origin == Category1):
                Reg1_x.append(Inds[j].GEmap[Gene_Name])
                Reg1_y.append(Inds[j].Phenotype)
                Curr_Ind_Left.append(Inds[j])
            elif(Origin == Category2):
                Reg2_x.append(Inds[j].GEmap[Gene_Name])
                Reg2_y.append(Inds[j].Phenotype)
                Curr_Ind_Right.append(Inds[j])

        if(len(set(Reg1_x)) < Ind_Threshold):
            R_sq1 = (-1 * (float)("inf"))
        else:
            R_sq1 = (float)(Get_Reg_Measure(Reg1_x, Reg1_y, "SSR")) / Get_Reg_Measure(Reg1_x, Reg1_y, "SST")

        if(len(set(Reg2_x)) < Ind_Threshold):
            R_sq2 = (-1 * (float)("inf"))
        else:
            R_sq2 = (float)(Get_Reg_Measure(Reg2_x, Reg2_y, "SSR")) / Get_Reg_Measure(Reg2_x, Reg2_y, "SST")

        Score = (float)((((R_sq1 + R_sq2))/2) * Imprv_Threshold)
        if(Score > Parent_Measure):
            Results.append((i, Score))

    if(len(Results) == 0):
        return None

    Results = sorted(Results, key = lambda x:x[1])[::-1]
    Best_SNP = Selected_SNPs[Results[0][0]]

    Curr_Ind_Left = []
    Curr_Ind_Right = []
    Reg1_x = []
    Reg2_x = []
    Reg1_y = []
    Reg2_y = []
    for j in range(len(Inds)):
        Name = Inds[j].Name
        Origin = Best_SNP.IndMap[Name]
        if(Origin == Category1):
            Reg1_x.append(Inds[j].GEmap[Gene_Name])
            Reg1_y.append(Inds[j].Phenotype)
            Curr_Ind_Left.append(Inds[j])
        elif(Origin == Category2):
            Reg2_x.append(Inds[j].GEmap[Gene_Name])
            Reg2_y.append(Inds[j].Phenotype)
            Curr_Ind_Right.append(Inds[j])

    Chosen_SNP.Name = Best_SNP.Name
    Chosen_SNP.IndMap = Best_SNP.IndMap
    Chosen_SNP.Position = Best_SNP.Position
    Chosen_SNP.Chromosome = Best_SNP.Chromosome
    Chosen_SNP.Left_Sample_Size = len(set(Reg1_x))
    Chosen_SNP.Right_Sample_Size = len(set(Reg2_x))
    Chosen_SNP.Left_Inds_List = [sample.Name for sample in Curr_Ind_Left]
    Chosen_SNP.Right_Inds_List = [sample.Name for sample in Curr_Ind_Right]
    Chosen_SNP.Left_Rsq = (float)(Get_Reg_Measure(Reg1_x, Reg1_y, "SSR"))/(Get_Reg_Measure(Reg1_x, Reg1_y, "SST"))
    Chosen_SNP.right_Rsq = (float)(Get_Reg_Measure(Reg2_x, Reg2_y, "SSR"))/(Get_Reg_Measure(Reg2_x, Reg2_y, "SST"))
    Chosen_SNP.Parent = Parent

    Left_Inds = Curr_Ind_Left
    Right_Inds = Curr_Ind_Right

    Best_SNP_Left_x, Best_SNP_Left_y, Best_SNP_Right_x, Best_SNP_Right_y = Reg1_x, Reg1_y, Reg2_x, Reg2_y

    Chosen_SNP.Left_Branch = Create_Branch(Chosen_SNP, List_Of_SNPs, Left_Inds, Chosen_SNP.Left_Rsq, Gene_Name)
    Chosen_SNP.Right_Branch = Create_Branch(Chosen_SNP, List_Of_SNPs, Right_Inds, Chosen_SNP.right_Rsq, Gene_Name)

    Chosen_SNP.Left_Params = Get_Reg_Measure(Best_SNP_Left_x, Best_SNP_Left_y)
    Chosen_SNP.Right_Params = Get_Reg_Measure(Best_SNP_Right_x, Best_SNP_Right_y)

    return Chosen_SNP


def Get_Leaves_Num(Current_SNP, Is_Root_Level):

    global Total_Leaves_Num

    if(Is_Root_Level == True):
        Total_Leaves_Num = 0
        if(Current_SNP == None):
            return 1

    if(Current_SNP == None):
        return -1

    retval = Get_Leaves_Num(Current_SNP.Left_Branch, False)
    if(retval == -1):
        Total_Leaves_Num += 1
    retval = Get_Leaves_Num(Current_SNP.Right_Branch, False)
    if(retval == -1):
        Total_Leaves_Num += 1
    if(Is_Root_Level == True):
        return Total_Leaves_Num
    return None


def Predict_By_Tree(SNP_Node, Sample, Gene_Name, Primary_Params, Initial_Call_Flag, Primary_Total_Error):

    global Tuple_Error_Threshold
    global Category1
    global Category2

    if(Initial_Call_Flag == True):

        if(SNP_Node == None):
            Prediction = Primary_Params[0] * Sample.GEmap[Gene_Name] + Primary_Params[1]
            return (Prediction, Primary_Total_Error)

    if(SNP_Node == None):
        return None

    Origin = Sample.SNPmap[SNP_Node.Name]

    if(Origin == Category1):
        retval = Predict_By_Tree(SNP_Node.Left_Branch, Sample, Gene_Name, None, False, None)
        if(retval == None):
            Params = SNP_Node.Left_params
            Prediction = Params[0] * Sample.GEmap[Gene_Name] + Params[1]
            if(len(SNP_Node.Left_PR_Tuple) < Tuple_Error_Threshold):
                Score = 0
            else:
                Score = SNP_Node.left_Total_Error
            return (Prediction, Score)
        else:
            return retval
    elif(Origin == Category2):
        retval = Predict_By_Tree(SNP_Node.Right_Branch, Sample, Gene_Name, None, False, None)
        if(retval == None):
            Params = SNP_Node.Right_params
            Prediction = Params[0] * Sample.GEmap[Gene_Name] + Params[1]
            if(len(SNP_Node.Right_PR_Tuple) < Tuple_Error_Threshold):
                Score = 0
            else:
                Score = SNP_Node.Right_Total_Error
            return (Prediction, Score)
        else:
            return retval
    else:
        return (-1, None)


def Print_Trees_To_File(Curr_Trees, Output_Path, Is_Perm):

    Count = 1

    if(Is_Perm):
        Full = "_Permuted"
    else:
        Full = ""
    Trees_File = open(Output_Path + "/Trees" + Full + ".txt", 'w')
    Trees_File.write("Tree Num\tScore\tGene\tLeaves Num\t#Unused Samples\n")

    for tree in Curr_Trees:
        Trees_File.write("Tree " + (str)(Count) + "\t")
        Trees_File.write((str)(tree.Score) + "\t")
        Trees_File.write(tree.Gene + "\t")
        Trees_File.write((str)(tree.Leaves_Num) + "\t")
        Trees_File.write((str)(len(tree.Unused_Samples)))
        Print_Tree(tree.Root, Trees_File)
        Trees_File.write("\n")
        Count += 1


def Print_Tree(Current_SNP, Trees_File):

    if(Current_SNP == None):
        return

    Trees_File.write("\t" + Current_SNP.Name)
    Print_Tree(Current_SNP.Left_Branch, Trees_File)
    Print_Tree(Current_SNP.Right_Branch, Trees_File)


def Predict(Node, Sample, Gene_Name, Real_Phen_Value):

    global Category1
    global Category2

    if(Node == None):
        return None

    Origin = Node.IndMap[Sample.Name]
    if(Origin == Category1):
        retval = Predict(Node.Left_Branch, Sample, Gene_Name, Real_Phen_Value)
        if(retval == None):
            Params = Node.Left_Params
            Prediction = (Params[0] * Sample.GEmap[Gene_Name]) + Params[1]
            Node.Left_PR_Tuple.append((Prediction, Real_Phen_Value))
            return 0
        else:
            return retval

    elif(Origin == Category2):
        retval = Predict(Node.Right_Branch, Sample, Gene_Name, Real_Phen_Value)
        if(retval == None):
            Params = Node.Right_Params
            Prediction = Params[0]*Sample.GEmap[Gene_Name] + Params[1]
            Node.Right_PR_Tuple.append((Prediction, Real_Phen_Value))
            return 0
        else:
            return retval
    else:
        return -1


def Calc_Regerssion_Error(SNP_Node):

    if(SNP_Node == None):
        return None

    retval = Calc_Regerssion_Error(SNP_Node.Left_Branch)
    if(retval == None):
        Errors = SNP_Node.Left_PR_Tuple
        Total_Error = 0
        if(len(Errors) == 0):
            SNP_Node.Left_Total_Error = 0.0
        else:
            for tup in Errors:
                Total_Error += (tup[0] - tup[1]) ** 2
            Total_Error /= len(Errors)
            SNP_Node.Left_Total_Error = Total_Error ** 0.5
    retval = Calc_Regerssion_Error(SNP_Node.Right_Branch)
    if(retval == None):
        Errors = SNP_Node.Right_PR_Tuple
        Total_Error = 0
        if(len(Errors) == 0):
            SNP_Node.Right_Total_Error = 0.0
        else:
            for tup in Errors:
                Total_Error += (tup[0] - tup[1]) ** 2
            Total_Error /= len(Errors)
            SNP_Node.Right_Total_Error = Total_Error ** 0.5
    return 0


def Agg_Tree_PR_Tuples(SNP_Node, Total_Tuples):

    if(SNP_Node == None):
        return None

    retval = Agg_Tree_PR_Tuples(SNP_Node.Left_Branch, Total_Tuples)
    if(retval == None):
        Total_Tuples = Total_Tuples + SNP_Node.Left_PR_Tuple
    else:
        Total_Tuples = Total_Tuples + retval
    retval = Agg_Tree_PR_Tuples(SNP_Node.Right_Branch, Total_Tuples)
    if(retval == None):
        Total_Tuples = Total_Tuples + SNP_Node.Right_PR_Tuple
    else:
        Total_Tuples = retval
    return Total_Tuples


def Get_SNPs_From_Tree(SNP_Node):

    SNPs_List = []
    if(SNP_Node == None):
        return []

    SNPs_List.extend(Get_SNPs_From_Tree(SNP_Node.Left_Branch))
    SNPs_List.extend(Get_SNPs_From_Tree(SNP_Node.Right_Branch))
    SNPs_List.extend([SNP_Node.Name])
    return SNPs_List

def Check_SNP(SNP_Data, Start_Idx, Ind_Threshold):
    Category1_Count = 0
    Category2_Count = 0
    Unique_Categories = []
    for entry in SNP_Data[Start_Idx:]:
        if(entry not in Unique_Categories):
            if(len(Unique_Categories)) > 2:
                return (0, None, None)
            Unique_Categories.append(entry)
    for entry in SNP_Data[Start_Idx:]:
        if(entry == Unique_Categories[0]):
            Category1_Count += 1
        elif(entry == Unique_Categories[1]):
            Category2_Count += 1
    if(Category1_Count < Ind_Threshold or Category2_Count < Ind_Threshold):
        return (1, None, None)
    return (2, Unique_Categories[0], Unique_Categories[1])


def Convert2Float(Data_List):
    for i in range(len(Data_List)):
        Data_List[i] = (float)(Data_List[i])
    return Data_List


def Extract_Gene2Score(Tree_List):

    Gene2Score = {}

    for tree in Tree_List:
        Gene = tree.Gene
        Score = tree.Score
        if(Gene not in Gene2Score):
            Gene2Score[Gene] = [Score]
        else:
            Scores_List = Gene2Score[Gene]
            Scores_List.append(Score)
            Gene2Score[Gene] = Scores_List

    return Gene2Score


def Unify(List1, List2):
    Final_List = list(set(List1) | set(List2))
    return Final_List


def Check_Args(Train_Path, Test_Path, Only_Prediction, Gene_Num, SNP_Num, Ind_Threshold_Given, Imprv_Threshold_Given, Output_Path, Final_Gene_Num):
    try:
        Curr_File = open(Train_Path, 'r').readlines()
    except IOError:
        print("Train file not found")
        return 0
    if(len(Curr_File) != SNP_Num + Gene_Num + 1):
        print("Error - Gene and SNP numbers not suitable with train file")
        return 0

    try:
        Curr_File = open(Test_Path, 'r').readlines()
    except IOError:
        print("Test file not found")
        return 0
    if(Only_Prediction):
        if(len(Curr_File) != SNP_Num + Gene_Num):
            print("Error - test file not suitable for prediction mode")
            return 0
    else:
        if(len(Curr_File) != SNP_Num + Gene_Num + 1):
            print("Error - test file not suitable for error mode")
            return 0
    Ind_Num = len(open(Train_Path).readlines()[1].rstrip().split("\t")[1:])
    if(Ind_Threshold_Given > Ind_Num/2):
        print "Error - individual threshold too high"
        return 0
    if(Imprv_Threshold_Given > 1 or Imprv_Threshold_Given < 0):
        print "Error - improvement threshold should be between 0 and 1"
        return 0
    if(Final_Gene_Num > Gene_Num):
        print "Error - number of selected forests cannot be larger than number of genes"
        return 0

    if(not os.path.exists(Output_Path)):
        print "Error - output path not found"
        return 0

    return 1

InPhenotype(sys.argv[1], sys.argv[2], (bool)((int)(sys.argv[3])), (int)(sys.argv[4]), (int)(sys.argv[5]), (int)(sys.argv[6]),
            (float)(sys.argv[7]), (int)(sys.argv[8]), sys.argv[9], (int)(sys.argv[10]), int(sys.argv[12]), (int)(sys.argv[11]))


