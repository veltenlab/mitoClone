#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==============================================================================
# Written by : Salem Malikic
# Modified by: Farid Rashidi
# Last Update: Apr 25, 2019
# ==============================================================================


from gurobipy import *
from helperFunctions import *
from datetime import datetime
import argparse
import errno
import pandas as pd
from faridFunctions import *


# COMMAND LINE ARGUMENTS PARSING
parser = argparse.ArgumentParser(description='PhISCS-I by Gurobi solver', add_help=True)
# Required arguments:
parser.add_argument('-SCFile', '--SCFile', required=True, type=str,
                    help='Path to single cell data matrix file')
parser.add_argument('-fn', '--fnProbability', required=True, type=float,
                    help='Probablity of false negative')
parser.add_argument('-fp', '--fpProbability', required=True, type=float,
                    help='Probablity of false positive')

# Optional:
parser.add_argument('-o', '--outDir', default='.', type=str,
                    help='Output directory')
parser.add_argument('-kmax', '--maxMutationsToEliminate', default=0, type=int,
                    help='Max number of mutations to be eliminated [default value is 0]')
parser.add_argument('-bulkFile', '--bulkFile', default=None, type=str,
                    help='Path to bulk data file')
parser.add_argument('-delta', '--delta', default=0.20, type=float,
                    help='Delta parameter accounting for VAF variance [default value is 0.20]')
parser.add_argument('-time', '--time', type=int, default=86400,
                    help='Max time (in seconds) allowed for the computation [default value is 24 hours]')
parser.add_argument('--drawTree', action='store_true',
                    help='Draw output tree by Graphviz')
parser.add_argument('--drawFarid', action='store_true',
                    help='Draw output tree by Graphviz')
# https://stackoverflow.com/questions/15753701/argparse-option-for-passing-a-list-as-option
# parser.add_argument('--candidateISAV',
#                     help='',
#                     type=str)

parser.add_argument('-w', '--colEliminationWeight', default=0, type=float,
                    help='Weight of column elimination [default value is 0]')
parser.add_argument('-threads', '--threads', default=1, type=int,
                    help='Number of threads [default value is 1]')


args = parser.parse_args()


# assert os.path.exists(args.outDir) == False, "ERROR!!! There already exists file or folder with name " + args.outDir + ". Exiting."
try:
    os.makedirs(args.outDir)
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(args.outDir):
        pass
    else:
        raise

filename = os.path.splitext(os.path.basename(args.SCFile))[0]
outfile = os.path.join(args.outDir, filename)
gurobi_log = "{}.gurobi".format(outfile)
start_model = datetime.now()
verbose = False
tree = False

# ======= Reading SC data input (cell x mut matrix with headers)
df = pd.read_csv(args.SCFile, sep='\t', index_col=0)
df = df.replace('?', 3)
df = df.astype(int)
mutIDs = df.columns
cellIDs = df.index
I = df.values
# SCFile  = open(args.SCFile, "r")
# mutIDs  = SCFile.readline().rstrip().split()[1:]
# cellIDs = []
# I = []
# for line in SCFile:
#     lineColumns = line.strip().split()
#     cellID = lineColumns[0]
#     cellIDs.append(cellID)
#     I.append([int(x) for x in lineColumns[1:]])
# SCFile.close()


numCells     = len(cellIDs)
numMutations    = len(mutIDs)
assert numMutations == len(set(mutIDs)), "ERROR!!! Some of the mutation IDs appear multiple times. Mutation IDs are " + str(mutIDs)
assert numCells == len(set(cellIDs)), "ERROR!!! Some of the cell IDs appear multiple times. Cell IDs are " + str(cellIDs)


beta  = args.fnProbability
alpha = args.fpProbability


isTrueVAF = False # trueVAF was used for the internal purposes during the method development
usingBulk = False
delta = None
bulkMutations = None
if args.bulkFile:
    delta = args.delta
    usingBulk = True
    bulkMutations = readMutationsFromBulkFile(args.bulkFile)
    assert len(bulkMutations) == numMutations, "ERROR!!! Single-cell and bulk data do not have the same number of mutations"
    for i in range(numMutations):
        assert bulkMutations[i].ID == mutIDs[i], "Mutations must be sorted in the same order in single cell and bulk data"


# =========== VARIABLES
model = Model('PhISCS_ILP')
# model.Params.LogFile = gurobi_log
model.Params.LogFile = ''
model.Params.Threads = args.threads
model.setParam('TimeLimit', args.time)

print('Generating variables...')


# ===== Matrix Y is matrix of corrected (i.e. true) genotypes w.r.t. input SC matrix I
Y = {}
for c in range(numCells):
    for m in range(numMutations):
            Y[c, m] = model.addVar(vtype=GRB.BINARY, name='Y({0},{1})'.format(c, m))



# ===== Variables B control the existence of conflict between columns
B = {}
for p in range(numMutations+1):
    for q in range(numMutations+1):
        B[p, q, 1, 1] = model.addVar(vtype=GRB.BINARY, obj=0,
                                     name='B[{0},{1},1,1]'.format(p, q))
        B[p, q, 1, 0] = model.addVar(vtype=GRB.BINARY, obj=0,
                                     name='B[{0},{1},1,0]'.format(p, q))
        B[p, q, 0, 1] = model.addVar(vtype=GRB.BINARY, obj=0,
                                     name='B[{0},{1},0,1]'.format(p, q))


# ==== Variable K[j] is set to 1 if and only if mutation j is among eliminated mutations
K = {}
# candidateISAV = [int(item.replace('mut','')) for item in args.candidateISAV.split(',')]
# f_i = args.SCFile
# f_o = '/data/frashidi/PhISCS/_result/TP_FP/noBulk_k_0/' + filename + '.CFMatrix'
# candidateISAV = give_me_muts_to_filter(f_i, f_o, args.maxMutationsToEliminate)
for m in range(numMutations+1):
    '''
    if m in list(set(range(numMutations)) - set(candidateISAV)):
        K[m] = model.addVar(vtype=GRB.BINARY, name='K[{0}]'.format(m), lb=0, ub=0)
    else:
        K[m] = model.addVar(vtype=GRB.BINARY, name='K[{0}]'.format(m))
    '''
    K[m] = model.addVar(vtype=GRB.BINARY, name='K[{0}]'.format(m))
model.addConstr(K[numMutations] == 0) # null mutation can not be eliminated 


# ==== A[p,q] = 1 if p is ancestor of q
A = {}
if usingBulk:
    for p in range(numMutations + 1): # mutation with index numMutation is null mutation
        for q in range(numMutations + 1):
            A[p,q] = model.addVar(vtype=GRB.BINARY, obj=0, name='A[{0},{1}]'.format(p,q))    


model.update()

# ====== CONSTRAINTS
print('Generating constraints...')

# --- number of eliminated columns is upper bounded by user provided constant
model.addConstr(quicksum(K[m] for m in range(numMutations)) <= args.maxMutationsToEliminate)


# --- Enforce three gametes rule
for i in range(numCells):
    for p in range(numMutations):
        for q in range(numMutations):
            model.addConstr(Y[i,p] + Y[i,q] - B[p,q,1,1] <= 1)
            model.addConstr(-Y[i,p] + Y[i,q] - B[p,q,0,1] <= 0)
            model.addConstr(Y[i,p] - Y[i,q] - B[p,q,1,0] <= 0)

# --- Null mutation present in each cell
for p in range(numMutations+1):
    model.addConstr(B[p,numMutations, 1, 0] == 0)


# --- Forbid conflict between columns (three gametes rule)
for p in range(numMutations):
    for q in range(numMutations):
        model.addConstr(B[p,q,0,1] + B[p,q,1,0] + B[p,q,1,1] <= 2 + K[p] + K[q])
    

# === Constraints for integrating VAF obtained from bulk data into the model
if usingBulk:
    validateSampleIDsConcordance(bulkMutations)
    sampleIDs = bulkMutations[0].getSampleIDs()
    bulkMutations.append(generateArtificialNullMutation(bulkMutations[0]))
    for p in range(numMutations):
        for q in range(p+1, numMutations):
            model.addConstr(A[p,q] + A[q,p] <= 1-K[p])
            model.addConstr(A[p,q] + A[q,p] <= 1-K[q])
            model.addConstr(A[p,q] + A[q,p] >= B[p,q,1,1] - K[p] - K[q])
    for p in range(numMutations+1):
        model.addConstr(A[p,p] == 0)
        for q in range(numMutations+1):
            model.addConstr(A[p, q] <= 1 - K[p])
            model.addConstr(A[p, q] <= 1 - K[q])
    
            if p<q:
                model.addConstr(A[p,q] + A[q,p] <= 1)

            model.addConstr(A[p,q] + B[p,q,0,1] <= 1 + K[p] + K[q])
            model.addConstr(B[p,q,1,0] + B[p,q,1,1] - A[p,q] <= 1 + K[p] + K[q])

            for sampleID in sampleIDs:
                VAF_p = bulkMutations[p].getVAF(sampleID)
                VAF_q = bulkMutations[q].getVAF(sampleID)
                if isTrueVAF:
                    VAF_p = bulkMutations[p].getTrueVAF(sampleID)
                    VAF_q = bulkMutations[q].getTrueVAF(sampleID)

                model.addConstr(A[p, q] * VAF_p * (1 + delta) >= A[p, q] * VAF_q)
                #'''
                for r in range(numMutations+1):
                    if r == q:
                        continue
                    VAF_r = bulkMutations[r].getVAF(sampleID)
                    if isTrueVAF:
                        VAF_r = bulkMutations[r].getTrueVAF(sampleID)
                    # Constraint 2
                    model.addConstr(
                                        VAF_p * (1 + delta) >= 
                                        VAF_q * (A[p, q] - A[r, q] - A[q, r]) + 
                                        VAF_r * (A[p, r] - A[r, q] - A[q, r])
                                    )
                #'''
            for r in range(numMutations+1):
                if r == q:
                    continue
                # Constraint 1.d
                model.addConstr(A[p, r] >= A[p, q] + A[q, r] - 1)
            

        candidateAncestors = [i for i in range(numMutations+1)]
        candidateAncestors.remove(p)

        if p < numMutations:
            model.addConstr(quicksum(A[s,p] for s in candidateAncestors) >= 1 - K[p])
        elif p == numMutations:
            model.addConstr(quicksum(A[s,p] for s in candidateAncestors) == 0)
        else:
            print("p index out of range. Exiting")
            sys.exit(2)


# --- Defining the objective function
objective = 0

for j in range(numMutations):
    numZeros = 0
    numOnes  = 0
    for i in range(numCells):
        if I[i][j] == 0:
            numZeros += 1
            objective += np.log(beta/(1-alpha)) * Y[i,j]
        elif I[i][j] == 1:
            numOnes += 1
            objective += np.log((1-beta)/alpha) * Y[i,j]
        
    objective += numZeros * np.log(1-alpha)
    objective += numOnes * np.log(alpha)
    objective -= K[j] * (numZeros * np.log(1-alpha) + numOnes * (np.log(alpha) + np.log((1-beta)/alpha)))

model.setObjective(objective, GRB.MAXIMIZE)
time_to_model = datetime.now() - start_model

# --- Optimize 
start_optimize = datetime.now()
model.optimize()



# ====== POST OPTIMIZATION
if model.status == GRB.Status.INFEASIBLE:
    print('The model is unfeasible.')
    exit(0)

time_to_opt = datetime.now() - start_optimize
time_to_run = datetime.now() - start_model



optimal_solution = model.ObjVal
print('Optimal solution: %f' % optimal_solution)

removedMutsIDs = []
sol_K = []
for j in range(numMutations):
    sol_K.append(nearestInt(float(K[j].X)))
    if sol_K[j] == 1:
        removedMutsIDs.append(mutIDs[j])

sol_Y = []
for i in range(numCells):
    sol_Y.append([nearestInt(float(Y[i,j].X)) for j in range(numMutations)])

conflictFreeMatrix = open("{}.CFMatrix".format(outfile), "w")
conflictFreeMatrix.write("cellID/mutID")
for j in range(numMutations):
    if sol_K[j] == 0:
        conflictFreeMatrix.write("\t" + mutIDs[j])

    
conflictFreeMatrix.write("\n")
for i in range(numCells):
    conflictFreeMatrix.write(cellIDs[i])
    for j in range(numMutations):
        if sol_K[j] == 0:
            conflictFreeMatrix.write("\t" + str(sol_Y[i][j]))
    conflictFreeMatrix.write("\n")
conflictFreeMatrix.close()


flips_0_1 = 0
flips_1_0 = 0
flips_3_0 = 0
flips_3_1 = 0
for i in range(numCells):
    for j in range(numMutations):
        if sol_K[j] == 0:
            if I[i][j] == 0 and sol_Y[i][j] == 1:
                flips_0_1 += 1
            elif I[i][j] == 1 and sol_Y[i][j] == 0:
                flips_1_0 += 1
            elif I[i][j] == 3 and sol_Y[i][j] == 0:
                flips_3_0 += 1
            elif I[i][j] == 3 and sol_Y[i][j] == 1:
                flips_3_1 += 1


# check if matrix is conflict free
conflictFree = True
for p in range(numMutations):
    if sol_K[p] == 1:
        continue
    for q in range(p + 1, numMutations):
        if sol_K[q] == 1:
            continue
        oneone = False
        zeroone = False
        onezero = False

        for r in range(numCells):
            if sol_Y[r][p] == 1 and sol_Y[r][q] == 1:
                oneone = True
            if sol_Y[r][p] == 0 and sol_Y[r][q] == 1:
                zeroone = True
            if sol_Y[r][p] == 1 and sol_Y[r][q] == 0:
                onezero = True

        if oneone and zeroone and onezero:
            conflictFree = False
            print('ERROR!!! Conflict in output matrix in columns (%d, %d)' % (p, q))


log = open('{}.log'.format(outfile), 'w+')    
# --- Input info
log.write('COMMAND: "{0}"\n'.format(' '.join(sys.argv)))
log.write('NUM_CELLS(ROWS): {0}\n'.format(str(numCells)))
log.write('NUM_MUTATIONS(COLUMNS): {0}\n'.format(str(numMutations)))
log.write('FN_WEIGHT: {0}\n'.format(str(beta)))
log.write('FP_WEIGHT: {0}\n'.format(str(alpha)))
log.write('COLUMN_ELIMINATION_WEIGHT: {0}\n'.format(str(args.colEliminationWeight)))
log.write('NUM_THREADS: {0}\n'.format(str(args.threads)))
log.write('MODEL_SOLVING_TIME_SECONDS: {0:.3f}\n'.format(time_to_opt.total_seconds()))
log.write('RUNNING_TIME_SECONDS: {0:.3f}\n'.format(time_to_run.total_seconds()))
if conflictFree:
    conflictFree = 'YES'
else:
    conflictFree = 'NO'
log.write('IS_CONFLICT_FREE: {0}\n'.format(conflictFree))
log.write('LIKELIHOOD: {0}\n'.format(str(optimal_solution)))
log.write('MIP_Gap_Value: %f\n' % model.MIPGap)
log.write('TOTAL_FLIPS_REPORTED: {0}\n'.format(str(flips_0_1 + flips_1_0 + flips_3_0 + flips_3_1)))
log.write('0_1_FLIPS_REPORTED: {0}\n'.format(str(flips_0_1)))
log.write('1_0_FLIPS_REPORTED: {0}\n'.format(str(flips_1_0)))
log.write('?_0_FLIPS_REPORTED: {0}\n'.format(str(flips_3_0)))
log.write('?_1_FLIPS_REPORTED: {0}\n'.format(str(flips_3_1)))
log.write('MUTATIONS_REMOVED_UPPER_BOUND: {0}\n'.format(str(args.maxMutationsToEliminate)))
log.write('MUTATIONS_REMOVED_NUM: {0}\n'. format(str(sum(sol_K))))
print('MUTATIONS_REMOVED_ID: {}\n'.format('.'.join(removedMutsIDs)))
log.write('MUTATIONS_REMOVED_ID: {}\n'.format(','.join(removedMutsIDs)))

log.write("-----------------------------------\n\n")
if usingBulk:
    for i in range(numMutations):
        for j in range(numMutations):
            log.write(bulkMutations[i].getID() + "\t" + bulkMutations[j].getID() + "\t" + str(nearestInt(float(A[i,j].X))))
            if sol_K[i] == 1:
                log.write("\t" + bulkMutations[i].getID() + " is eliminated.")
            if sol_K[j] == 1:
                log.write("\t" + bulkMutations[j].getID() + " is eliminated.")
            log.write("\n")
    for i in range(numMutations):
        log.write("NULL" + "\t" + bulkMutations[i].getID() + "\t" + str(nearestInt(float(A[numMutations,i].X))))
        if sol_K[i] == 1:
            log.write("\t" + bulkMutations[i].getID() + " is eliminated.")
        log.write("\n")
        log.write(bulkMutations[i].getID() + "\t" + "NULL" + "\t" + str(nearestInt(float(A[i,numMutations].X))))
        if sol_K[i] == 1:
            log.write("\t" + bulkMutations[i].getID() + " is eliminated.")
        log.write("\n")
log.close()

if args.drawTree:
    draw_tree("{}.CFMatrix".format(outfile), usingBulk, args.bulkFile)
if args.drawFarid:
    draw_farid("{}.CFMatrix".format(outfile), usingBulk, args.bulkFile)
