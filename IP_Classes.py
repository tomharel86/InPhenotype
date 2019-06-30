class Ind:
    def __init__(self):
        self.Name=None
        self.Phenotype=0
        self.GEmap={}
        self.SNPmap = {}

    def __repr__(self):
        return "Ind({},{})".format(self.Name,(str)(self.Pexp))


class SNP:

    def __init__(self):
        self.Name=None
        self.Chromosome=0
        self.Position=0
        self.IndMap={}
        self.Left_Branch = None
        self.Right_Branch = None
        self.Parent = None
        self.Left_Sample_Size = 0
        self.Right_Sample_Size = 0
        self.Left_Params = []
        self.Left_Inds_List =[]
        self.Right_Inds_List = []
        self.Right_Params = []
        self.Left_PR_Tuple = []
        self.Right_PR_Tuple = []
        self.Left_Total_Error = 0
        self.Right_Total_Error = 0


    def __repr__(self):
        return "SNP({},{})".format(self.Name,(str)(self.Chromosome))


class Tree:

    def __init__(self):
        self.Root = None
        self.Score = 0.0          #Root_SSE / SSE_Leaves = best score is the highest
        self.Unused_Samples = []
        self.Errors = []

    def __repr__(self):
        return "Tree({},{})".format(self.Root,(str)(self.Score))



class Forest:

    def __init__(self):
        self.Tree_List = []
        self.Score = 0.0
        self.Gene = None

    def __repr__(self):
        return "Forest({},{},{})".format((str)(self.Score),(str)(self.Gene),(str)(len(self.Tree_List)))



