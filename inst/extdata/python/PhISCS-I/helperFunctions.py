import numpy as np
import sys
import os

def getEntryStringValue(keyValueString, entryID):
    keyValuePairs = keyValueString.split(";")
    for keyValuePair in keyValuePairs:
        keyString   = keyValuePair.strip().split("=")[0]
        valueString = keyValuePair.strip().split("=")[1]
        if keyString == str(entryID):
            return valueString
    print("ERROR in getEntryStringValue. There does not exit entry with ID " + entryID + " in string" + keyValueString + ". EXITING !!!")
    sys.exit(2)


def getEntryStrValue(keyValueString, entryID):
    return getEntryStringValue(keyValueString, entryID)


def isFloat(value):
    try:
        float(value)
        return True
    except ValueError:
            return False


def nearestInt(x):
    return int(x+0.5)


def floatToStr(floatNumber, numDecimals = 2):
    if isFloat(floatNumber):
        return ("{:." + str(numDecimals) + "f}").format(floatNumber)
    else:
        return floatNumber


def floatToString(floatNumber, numDecimals = 2):
    return floatToStr(floatNumber, numDecimals)


def simulateMutantReferenceCounts(coverage, VAF, coverageStdev = 0.01):
    totalCount     = coverage
    mutantCount    = np.random.binomial(totalCount, VAF/2)
    referenceCount = totalCount - mutantCount
    return [mutantCount, referenceCount]





class Mutation:

    def getINFOEntryStringValue(self, entryID):
        INFO_columns = self.INFO.split(";")
        for column in INFO_columns:
            if column.strip().split("=")[0] == str(entryID):
                return column.strip().split("=")[1]
        print("There does not exit " + entryID + " in " + self.ID + ". EXITING !!!")
        sys.exit(2)


    def getINFOEntryStrValue(self, entryID):
        return self.getINFOEntryStringValue(entryID)


    def __init__(self, ID = "NA", chromosome="none", position="0", mutReads={}, refReads={}, INFO=""):
        self.ID        = ID
        self.chromosome = chromosome
        self.position    = position
        self.INFO    = INFO
        assert set(mutReads.keys()) == set(refReads.keys()), "In Mutation constructor, set of samples different for variant and reference reads"
        self.mutReads = {}
        for sampleID in mutReads.keys():
            self.mutReads[sampleID] = mutReads[sampleID]
        self.refReads = {}
        for sampleID in refReads.keys():
            self.refReads[sampleID] = refReads[sampleID]
        if len(self.refReads.keys()) > 0:
            ERROR_MESSAGE = "In Mutation constructor, sampleIDs field in INFO does not agree with sampleIDs present in mutReads"
            assert set(self.getINFOEntryStringValue("sampleIDs").split(",")) == set(refReads.keys()), ERROR_MESSAGE
        else:
            ERROR_MESSAGE = "If no mut/ref reads INFO is provided then sampleIDs should be empty or absent in INFO field"
            assert ("sampleIDs" not in INFO) or len(self.getINFOEntryStringValue("sampleIDs")) == 0, ERROR_MESSAGE


    def setMutReadsInSample(self, newReadCount, sampleID):
        assert sampleID in self.mutReads.keys(), "ERROR in Mutation.py in setMutReadsInSample. Sample ID not present in hash."
        self.mutReads[sampleID] = newReadCount

    def setRefReadsInSample(self, newReadCount, sampleID):
        assert sampleID in self.refReads.keys(), "ERROR in Mutation.py in setRefReadsInSample. Sample ID not present in hash."
        self.refReads[sampleID] = newReadCount        

    def getVAF(self, sampleID=""):
        if sampleID == "":
            assert len(self.refReads.keys()) == 1, "ERROR. sampleID must be provided in getVAF of Mutation class if num samples is not 1"
            return self.getVAF(self.mutReads.keys()[0])
        assert sampleID in self.mutReads.keys(), "ERROR. Incorrect sampleID " + sampleID + " in class Mutation, function getVAF"
        if self.mutReads[sampleID] + self.refReads[sampleID] == 0:
            return "NA"
        return float(2.0*self.mutReads[sampleID])/(self.mutReads[sampleID] + self.refReads[sampleID])


    def getChromosome(self):
        return self.chromosome


    def getPosition(self):
        return self.position

        
    def getSampleIDs(self):
        # we expect to have something like sampleIDs=S1,S2,S3;
        return self.getINFOEntryStringValue("sampleIDs").strip().split(",")


    def updateINFOEntryValue(self, entryID, newValue):
        INFO_pairs = self.INFO.strip(";").split(";")
        updatedINFO = ""
        IDfound = False
        for pair in INFO_pairs: # pair is (ID, value) pair in form "ID=value"
            pairID = pair.split("=")[0]
            if pairID == entryID:
                updatedINFO += pairID + "=" + str(newValue) + ";"
                IDfound = True
            else:
                updatedINFO += pair + ";"
        assert IDfound == True, "ERROR in function updateINFOEntryValue in Mutation class. INFO does not contain entry with ID " + entryID
        self.INFO = updatedINFO


    def addINFOEntry(self, entryID, value):
        INFO_pairs = self.INFO.strip(";").split(";")
        IDfound = False
        for pair in INFO_pairs: # pair is (ID, value) pair in form "ID=value"
            pairID = pair.split("=")[0]
            if pairID == entryID:
                IDfound = True
        assert IDfound == False, "ERROR in function addINFOEntry in Mutation class. INFO already contains entry with ID " + entryID
        self.INFO = self.INFO.rstrip(";") + ";" + entryID + "=" + value + ";"
                 

    def updateID(self, newID):
        self.updateINFOEntryValue("ID", newID)
    

    def toString(self):
        strRepresentation = ""
        strRepresentation += str(self.ID) + "\t"
        strRepresentation += str(self.chromosome) + "\t"
        strRepresentation += str(self.position) + "\t"
        strRepresentation += ';'.join([str(self.mutReads[sampleID]) for sampleID in self.getSampleIDs()]) + "\t"
        strRepresentation += ';'.join([str(self.refReads[sampleID]) for sampleID in self.getSampleIDs()]) + "\t"
        strRepresentation += self.INFO
        return strRepresentation


    def __lt__(self, other):
        if int(self.chromosome) < int(other.chromosome):
            return True
        if int(self.chromosome) == int(other.chromosome):
            return self.position < other.position
        return False


    def getID(self):
        return self.ID
    
    def getGeneID(self):
        if (self.INFO.startswith("geneID=") or ";geneID=" in self.INFO) and (self.INFO.startswith("ID=") or ";ID=" in self.INFO):
            print("ERROR in function getGeneID in Mutation class")
            print("Mutation " + self.toString() + " has both geneID and ID in its INFO.")
            print("This causes issues in plotting functions (ID used in phiscs, geneID in BSCITE")
            sys.exit(2)
        if self.INFO.startswith("geneID") or ";geneID=" in self.INFO:
            return self.getINFOEntryStringValue("geneID")
        elif self.INFO.startswith("ID=") or ";ID=" in self.INFO:
            return self.getINFOEntryStringValue("ID")
        else:
            print("ERROR in function getGeneID in Mutation class. INFO entry does not contain geneID nor ID")
            print(self.toString())
            sys.exit(2)
        

    def getRefReads(self, sampleID=""):
        if sampleID == "":
            ERROR_MESSAGE =  "ERROR. In getRefReads() function in Mutation class, sampleID argument must be provided, "
            ERROR_MESSAGE += "unless the number of samples is equal to 1."
            assert len(self.refReads.keys()) == 1, ERROR_MESSAGE
            return self.getRefReads(self.refReads.keys()[0])
        return self.refReads[sampleID]


    def getMutReads(self, sampleID=""):
        if sampleID == "":
            ERROR_MESSAGE =  "ERROR. In getMutReads() function in Mutation class, sampleID argument must be provided, "
            ERROR_MESSAGE += "unless the number of samples is equal to 1."
            assert len(self.mutReads.keys()) == 1, ERROR_MESSAGE
            return self.getMutReads(self.mutReads.keys()[0])
        return self.mutReads[sampleID]

    
    def getTotalReads(self, sampleID=""):
        return self.getMutReads(sampleID) + self.getRefReads(sampleID)


    def getTrueVAF(self, sampleID=""):
        if sampleID == "":
            ERROR_MESSAGE =  "ERROR. In getTrueVAF() function in Mutation class, sampleID argument must be provided, "
            ERROR_MESSAGE += "unless the number of samples is equal to 1."
            assert len(self.mutReads.keys()) == 1, ERROR_MESSAGE
            return self.getTrueVAF(self.mutReads.keys()[0])
        sampleIDs = self.getSampleIDs()
        assert sampleID in sampleIDs, "ERROR in function getTrueVAF() in Mutation class. sampleID " + sampleID + " not found!"
        trueVAFs = [float(x) for x in self.getINFOEntryStringValue("trueVAF").split(",")]
        assert len(trueVAFs) == len(sampleIDs), "ERROR in function getTrueVAF() in Mutation class. Lengths of trueVAFs and sampleIDs unequal"
        for i in range(len(sampleIDs)):
            if sampleIDs[i] == sampleID:
                return trueVAFs[i]
        assert True, "ERROR. Last assert in getTrueVAF() in Mutation class failed."
        


    def updateMutRefReads(self, newCoverage, sampleID=""):
        if sampleID == "":
            ERROR_MESSAGE =  "ERROR. In updateMutRefReads() function in Mutation class, sampleID argument must be provided, "
            ERROR_MESSAGE += "unless the number of samples is equal to 1."
            assert len(self.mutReads.keys()) == 1, ERROR_MESSAGE
            return self.updateMutRefReads(newCoverage, self.mutReads.keys()[0])

        VAF = self.getTrueVAF(sampleID)
        #newTotalReads    = int(np.random.normal(newCoverage, newCoverageStdev) + 0.5)  # number of reads spanning mutation locus
        newTotalReads   = newCoverage
        newMutReads     = np.random.binomial(newTotalReads, VAF/2)
        newRefReads    = newTotalReads - newMutReads
        self.mutReads[sampleID]    = newMutReads
        self.refReads[sampleID]    = newRefReads
        return True


    def updateMutRefReadsAllSamples(self, newCoverage):
        for sampleID in self.mutReads.keys():
            self.updateMutRefReads(newCoverage, sampleID)


    def updateTrueVAFinSample(self, sampleID, newVAF):
        assert sampleID in self.getSampleIDs(), "ERROR in function updateTrueVAFinSample. No sample with ID " + sampleID
        updateSampleVAFindex = self.getSampleIDs().index(sampleID)
        currentValues = self.getINFOEntryStringValue("trueVAF").strip().split(",")
        newValues     = []
        for i in range(len(currentValues)): #i is sample index
            if i == updateSampleVAFindex:
                newValues.append(str(newVAF))
            else:
                newValues.append(currentValues[i])
        self.updateINFOEntryValue("trueVAF", ",".join(newValues))
            
    

    def addSample(self, sampleID, mutReads, refReads, trueVAF = ""):
        ERROR_MESSAGE =  "\nEROR in addSample in Mutation.py. SampleID " + sampleID + " is already present."
        ERROR_MESSAGE += "\nMutation ID is " + self.getID()
        assert sampleID not in self.getSampleIDs(), ERROR_MESSAGE
        if trueVAF != "":
            currentTrueVAFstring = self.getINFOEntryStringValue("trueVAF")
            self.updateINFOEntryValue("trueVAF", currentTrueVAFstring.rstrip(";") + "," + floatToString(trueVAF, 5))
        currentSampleIDsString = self.getINFOEntryStringValue("sampleIDs")
        self.updateINFOEntryValue("sampleIDs", currentSampleIDsString.rstrip(";") + "," + sampleID)
        self.mutReads[sampleID] = mutReads
        self.refReads[sampleID] = refReads


    def getAverageCoverageInAllSamples(self):
        assert set(self.mutReads.keys()) == set(self.getSampleIDs()), "ERROR in getAverageCoverageInAllSamples. Unequal sets of sample IDs."
        assert set(self.mutReads.keys()) == set(self.refReads.keys()), "ERROR in getAverageCoverageInAllSamples. Unequal sets of sample IDs."
        totalCoverage = 0
        for sampleID in self.getSampleIDs():
            totalCoverage += self.mutReads[sampleID]
            totalCoverage += self.refReads[sampleID]
        return int(totalCoverage/len(self.getSampleIDs()))


    def reorderSampleIDs(self, desiredOrderOfSamples): # this function was used in parsing data from Andrew's paper
        sampleIDs = self.getINFOEntryStringValue("sampleIDs").strip().split(",")
        assert set(sampleIDs) == set(desiredOrderOfSamples), "ERROR in reorderSampleIDs. Different sets of sampleIDs."

        newINFOvalue = ""
        processedEntries = []

        for entryID in [x.strip().split("=")[0] for x in self.INFO.rstrip(";").split(";")]:
            assert entryID not in processedEntries, "ERROR. Entry " + entryID + " repeated in reorderSampleIDs. Mut ID is " + self.getID()
            entryValueInSample = {}
            valuesInSamples = self.getINFOEntryStringValue(entryID).strip().split(",")

            if entryID in ["refNuc", "altNuc", "geneID"]:
                assert len(valuesInSamples) == 1, "ERROR in reorderSampleIDs. Reference or altered nucleotide takes only one value."
                newINFOvalue += entryID + "=" + valuesInSamples[0] + ";"
            else:
                assert len(valuesInSamples) == len(sampleIDs), "ERROR in reorderSampleIDs. Different lengths of sample IDs."
                for i in range(len(sampleIDs)):
                    entryValueInSample[sampleIDs[i]] = valuesInSamples[i]
                newINFOvalue += entryID
                newINFOvalue += "="
                newINFOvalue += ",".join([entryValueInSample[sampleID] for sampleID in desiredOrderOfSamples])
                newINFOvalue += ";"
                processedEntries.append(entryID)

        self.INFO = newINFOvalue


def strToMutation(inputString):
    stringColumns     = inputString.strip().split()
    ID        = stringColumns[0]
    chromosome     = stringColumns[1]
    position     = int(stringColumns[2])
    assert len(stringColumns) > 5, "ERROR in function strToMutation. INFO column empty. sampleIDs= required in this column"
    INFO = stringColumns[5]

    sampleIDs = getEntryStringValue(INFO, "sampleIDs").strip().split(",")
    assert len(sampleIDs) > 0, "sampleIDs field absent in " + inputString + " and therefore can not convert this string to Mutation."    
    
    mutReads = {} # column[3]
    refReads = {} # column[4]
    ERROR_MESSAGE = "ERROR! Number of sampleIDs and number of mutReads (thirdColumn), as well as number of refReads (column 4) must be equal."
    ERROR_MESSAGE += "\nSome of these numbers are not equal and therefore string " + inputString + " can not be converted to Mutation."
    assert len(stringColumns[3].strip().split(";")) == len(sampleIDs), ERROR_MESSAGE
    assert len(stringColumns[4].strip().split(";")) == len(sampleIDs), ERROR_MESSAGE

    for i in range(len(sampleIDs)):
        sampleID = sampleIDs[i]
        mutReads[sampleID] = int(stringColumns[3].split(";")[i])
        refReads[sampleID] = int(stringColumns[4].split(";")[i])

    return Mutation(ID, chromosome, position, mutReads, refReads, INFO)


def readMutationsFromBulkFile(pathBulkFile):
    assert os.path.exists(pathBulkFile), "ERROR in function readMutationsFromBulkFile!!! There does not exist bulk file " + pathBulkFile
    bulkFile = open(pathBulkFile, "r")
    bulkFile.readline()
    bulkMutations = []
    for line in bulkFile:
        bulkMutations.append(strToMutation(line))
    bulkFile.close()
    return bulkMutations


def validateSampleIDsConcordance(bulkMutations):
    sampleIDs = bulkMutations[0].getSampleIDs()
    numSamples = len(sampleIDs)
    assert numSamples > 0, "ERROR in validateSampleIDsConcordance. Number of samples zero!"
    for i in range(1, len(bulkMutations)):
        currentSampleIDs = bulkMutations[i].getSampleIDs()
        assert len(currentSampleIDs) == len(sampleIDs), "ERROR in function validateSampleIDsConcordance for mutations " + bulkMutations[0].getID() + " " + bulkMutations[i].getID()
        for sampleIndex in range(numSamples):
            assert sampleIDs[sampleIndex] == currentSampleIDs[sampleIndex], "ERROR in function validateSampleIDsConcordance. Sample IDs are not concordant"



def generateArtificialNullMutation(existingBulkMutation):
    sampleIDs = existingBulkMutation.getSampleIDs()
    refReads = {}
    mutReads = {}
    for sampleID in sampleIDs:
        mutReads[sampleID] = existingBulkMutation.getTotalReads(sampleID)
        refReads[sampleID] = 0
    INFO = ""
    INFO += "trueVAF=" + ",".join(['1.0' for s in sampleIDs]) + ";"
    INFO += "sampleIDs=" + ",".join(sampleIDs) + ";" 
    return Mutation("NULL Mutation abcdef", "NULLchromosome", 0, mutReads, refReads, INFO)


def draw_tree(filename, addBulk, bulkfile):
    import pandas as pd
    import pygraphviz as pyg

    graph = pyg.AGraph(strict=False, directed=True)
    font_name = 'Avenir'

    class Node:
        def __init__(self, name, parent):
            self.name = name
            self.parent = parent
            self.children = []
            if parent:
                parent.children.append(self)

    def print_tree(node):
        graph.add_node(node.name, label=node.name, fontname=font_name, color='black', penwidth=3.5)
        for child in node.children:
            graph.add_edge(node.name, child.name)
            print_tree(child)

    def contains(col1, col2):
        for i in range(len(col1)):
            if not col1[i] >= col2[i]:
                return False
        return True

    def write_tree(matrix, names):
        i = 0
        while i < matrix.shape[1]:
            j = i + 1
            while j < matrix.shape[1]:
                if np.array_equal(matrix[:,i], matrix[:,j]):
                    matrix = np.delete(matrix, j, 1)
                    x = names.pop(j)
                    names[i] += '<br/><br/>' + x
                    j -= 1
                j += 1
            names[i] = '<'+names[i]+'>'
            i += 1

        rows = len(matrix)
        cols = len(matrix[0])
        dimensions = np.sum(matrix, axis=0)
        # ordered indeces
        indeces = np.argsort(dimensions)
        dimensions = np.sort(dimensions)
        mutations_name = []
        for i in range(cols):
            mutations_name.append(names[indeces[i]])

        root = Node(mutations_name[-1], None)
        mut_nod = {}
        mut_nod[mutations_name[cols-1]] = root

        i = cols - 2
        while i >=0:
            if dimensions[i] == 0:
                break
            attached = False
            for j in range(i+1, cols):
                if contains(matrix[:, indeces[j]], matrix[:, indeces[i]]):
                    node = Node(mutations_name[i], mut_nod[mutations_name[j]])
                    mut_nod[mutations_name[i]] = node
                    attached = True
                    break
            if not attached:
                node = Node(mutations_name[i], root)
                mut_nod[mutations_name[i]] = node
            i -=1
        print_tree(root)

    if addBulk:
        vafs = {}
        bulkMutations = readMutationsFromBulkFile(bulkfile)
        sampleIDs = bulkMutations[0].getSampleIDs()
        for mut in bulkMutations:
            temp_vaf = []
            for sample in sampleIDs:
                temp_vaf.append('<font color="blue">' + str(mut.getVAF(sampleID=sample)) + '</font>')
            vafs[mut.getID()] = '{} ({})'.format(mut.getID(), ','.join(temp_vaf))

    inp = np.genfromtxt(filename, skip_header=1, delimiter='\t')
    with open(filename, 'r') as fin:
        if addBulk:
            mutation_names = [vafs[x] for x in fin.readline().strip().split('\t')[1:]]
        else:
            mutation_names = fin.readline().strip().split('\t')[1:]
    sol_matrix = np.delete(inp, 0, 1)
    write_tree(sol_matrix, mutation_names)
    graph.layout(prog='dot')
    outputpath = filename[:-len('.CFMatrix')]
    graph.draw('{}.png'.format(outputpath))



def draw_farid(filename, addBulk, bulkfile):
    add_cells=False

    import pandas as pd
    import pygraphviz as pyg
    import networkx as nx
    from networkx.drawing.nx_agraph import graphviz_layout, to_agraph

    def contains(col1, col2):
        for i in range(len(col1)):
            if not col1[i] >= col2[i]:
                return False
        return True

    df = pd.read_csv(filename, sep='\t').set_index('cellID/mutID')
    splitter_mut = '\n'
    matrix = df.values
    names_mut = list(df.columns)

    i = 0
    while i < matrix.shape[1]:
        j = i + 1
        while j < matrix.shape[1]:
            if np.array_equal(matrix[:,i], matrix[:,j]):
                matrix = np.delete(matrix, j, 1)
                x = names_mut.pop(j)
                names_mut[i] += splitter_mut + x
                j -= 1
            j += 1
        i += 1

    rows = matrix.shape[0]
    cols = matrix.shape[1]
    dimensions = np.sum(matrix, axis=0)
    indices = np.argsort(dimensions)
    dimensions = np.sort(dimensions)
    names_mut = [names_mut[indices[i]] for i in range(cols)]

    G = nx.DiGraph()
    G.add_node(cols)
    G.add_node(cols-1)
    G.add_edge(cols, cols-1, label=names_mut[cols-1])
    node_mud = {}
    node_mud[names_mut[cols-1]] = cols-1

    i = cols - 2
    while i >= 0:
        if dimensions[i] == 0:
            break
        attached = False
        for j in range(i+1, cols):
            if contains(matrix[:, indices[j]], matrix[:, indices[i]]):
                G.add_node(i)
                G.add_edge(node_mud[names_mut[j]], i, label=names_mut[i])
                node_mud[names_mut[i]] = i
                attached = True
                break
        if not attached:
            G.add_node(i)
            G.add_edge(cols, i, label=names_mut[i])
            node_mud[names_mut[i]] = i
        i -=1

    clusters = {}
    for node in G:
        if node == cols:
            G.node[node]['label'] = '<<b>germ<br/>cells</b>>'
            G.node[node]['fontname'] = 'Helvetica'
            G.node[node]['width'] = 0.4
            G.node[node]['style'] = 'filled'
            G.node[node]['penwidth'] = 3
            G.node[node]['fillcolor'] = 'gray60'
            continue
        untilnow_mut = []
        sp = nx.shortest_path(G, cols, node)
        for i in range(len(sp)-1):
            untilnow_mut += G.get_edge_data(sp[i], sp[i+1])['label'].split(splitter_mut)
        untilnow_cell = df.loc[(df[untilnow_mut] == 1).all(axis=1) & \
                               (df[[x for x in df.columns if x not in untilnow_mut]] == 0).all(axis=1)].index
        if len(untilnow_cell) > 0:
            clusters[node] = '\n'.join(untilnow_cell)
        else:
            clusters[node] = '-'
        
        if add_cells:
            G.node[node]['label'] = clusters[node]
        else:
            G.node[node]['label'] = ''
            G.node[node]['shape'] = 'circle'
        G.node[node]['fontname'] = 'Helvetica'
        G.node[node]['width'] = 0.4
        G.node[node]['style'] = 'filled'
        G.node[node]['penwidth'] = 2
        G.node[node]['fillcolor'] = 'gray90'
    i = 1
    for k, v in clusters.items():
        if v == '-':
            clusters[k] = i*'-'
            i += 1

    header = ''
    if addBulk:
        vafs = {}
        bulkMutations = readMutationsFromBulkFile(bulkfile)
        sampleIDs = bulkMutations[0].getSampleIDs()
        for mut in bulkMutations:
            temp_vaf = []
            for sample in sampleIDs:
                temp_vaf.append(str(mut.getVAF(sampleID=sample)))
            vafs[mut.getID()] = '<font color="blue">'+','.join(temp_true)+'</font>'        
        for edge in G.edges():
            temp = []
            for mut in G.get_edge_data(edge[0],edge[1])['label'].split(splitter_mut):
                mut = '<u>' + mut + '</u>' + ': ' + vafs_true[mut] + '; ' + vafs_noisy[mut]
                temp.append(mut)
            temp = '<' + '<br/>'.join(temp) + '>'
            G.get_edge_data(edge[0],edge[1])['label'] = temp

        for mut in bulkMutations:
            try:
                isatype = mut.getINFOEntryStringValue('ISAVtype')
                header += mut.getID() + ': ' + isatype + '<br/>'
            except:
                pass
    
    temp = df.columns[(df==0).all(axis=0)]
    if len(temp) > 0:
        header += 'Became Germline: ' + ','.join(temp) + '<br/>'
    
    with open(filename[:-len('.CFMatrix')]+'.log') as fin:
        i = 0
        for line in fin:
            i += 1
            if i > 10 and i < 18:
                header += line.rstrip() + '<br/>'
    

    H = nx.relabel_nodes(G, clusters)
    html = '''<{}>'''.format(header)
    H.graph['graph'] = {'label':html, 'labelloc':'t', 'resolution':300, 'fontname':'Helvetica', 'fontsize':8}
    H.graph['node'] = {'fontname':'Helvetica', 'fontsize':8}
    H.graph['edge'] = {'fontname':'Helvetica', 'fontsize':8}
    
    mygraph = to_agraph(H)
    mygraph.layout(prog='dot')
    outputpath = filename[:-len('.CFMatrix')]
    mygraph.draw('{}.png'.format(outputpath))
