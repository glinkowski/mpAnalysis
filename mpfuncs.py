#import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
from os import listdir
from os.path import isfile
from os import rename




######## ######## ######## ########
# Function: Read in the Knowledge Graph
# Returns:
#   Edges - (Nx4) matrix of char strings
#       each row is: node, node, edge weight, edge type
#   Nodes - dictionary of nodes in the edge list
#       key = node name
#       value = list of indices of the rows in
#           the edge matrix where node (key) appears
# matrix of Edges (N,4), set of Vertices
def readEdgeFile(datafile, delimiter) :

    # get the number of lines in the file
    nLines = sum( 1 for line in open(datafile, "rb") )

    # assign space for edge list
#    dt = np.dtype('a23')
    dt = np.dtype.str
#    Edges = np.empty( [nLines,4], dtype=dt)
    Edges = np.empty( (nLines,4), dtype=dt)

    # dictionary to hold Node indices
    Nodes = dict()

    # Start reading from the file
    df = open(datafile, "rb")

    i = 0
    for line in df:
        # extract the data from the file
        line = line.rstrip()
        lv = line.split(delimiter)

        # insert into the edge list
        Edges[i] = lv

        # add node locations to dict
        if (lv[0] in Nodes.keys()) :
            Nodes[lv[0]].append(i)
        else :
            Nodes[lv[0]] = list()
            Nodes[lv[0]].append(i)
        #end if
        if (lv[1] in Nodes.keys()) :
            Nodes[lv[1]].append(i)
        else :
            Nodes[lv[1]] = list()
            Nodes[lv[1]].append(i)
        #end if

        i += 1
    #end loop

    # close the data file
    df.close()

#    print "  file contained {:,} lines".format(nLines)
    return Edges, Nodes
#end def ######## ######## ########

######## ######## ######## ########
# Function: Read in an edge file with corresponding
#       node dictionary file
# Returns:
#   Edges - (Nx4) matrix of char strings
#       each row is: node, node, edge weight, edge type
#   Nodes - dictionary of nodes in the edge list
#       key = node name
#       value = list of indices of the rows in
#           the edge matrix where node (key) appears
def readEdgeFilePlus(edgefile, nodefile, delimiter) :

    # get the number of lines in the file
    nLines = sum( 1 for line in open(datafile, "rb") )

    # assign space for edge list
#    dt = np.dtype('a23')
    dt = np.dtype.str
    Edges = np.empty( (nLines,4), dtype=dt)
    # Start reading from the file
    ef = open(edgefile, "rb")
    i = 0
    for line in ef:
        # extract the data from the file
        line = line.rstrip()
        lv = line.split(delimiter)
        # insert into the edge list
        Edges[i] = lv
        i += 1
    #end loop
    # close the data file
    df.close()

    # dictionary to hold Node indices
    Nodes = dict()
    # Start reading from file
    nf = open(nodefile, "rb")
    for line in nf:
        # extract the data
        line = line.rstrip()
        lv = line.split(delimiter)
        # Create the list to insert
        li = lv[1].split(',')
        # Add to the dict
        Nodes[lv[0]] = list()
        Nodes[lv[0]].extend(li)
    #end loop

    return Edges, Nodes
#end def ######## ######## ########
######## ######## ######## ########
# Function: Write the Knowledge Graph to a file
# Returns: nothing
def writeEdgeFile(network, datafile, delimiter) :

    f = open(datafile, 'wb')

    first = True
    for row in network :

        if first :
            f.write('{1}{0}{2}{0}{3}{0}{4}'.format(
                delimiter, row[0], row[1], row[2], row[3]))
            first = False
            continue
        #end if


        f.write('\n{1}{0}{2}{0}{3}{0}{4}'.format(
            delimiter, row[0], row[1], row[2], row[3]))
    #end loop

    f.close()

#end def ######## ######## ########
######## ######## ######## ########
# Function: Write the Knowledge Graph to a file
# Returns: nothing
def writeEdgeFilePlus(network, nodeDict, edgefile, dictfile, delimiter) :

    # Save the network (edge list) to a file
    f = open(edgefile, 'wb')

    first = True
    for row in network :

        if first :
            f.write('{1}{0}{2}{0}{3}{0}{4}'.format(
                delimiter, row[0], row[1], row[2], row[3]))
            first = False
            continue
        #end if


        f.write('\n{1}{0}{2}{0}{3}{0}{4}'.format(
            delimiter, row[0], row[1], row[2], row[3]))
    #end loop

    f.close()

#TODO: write the node dictionary to file

    # Save the dictionary to a file
    allKeys = nodeDict.keys()
    allKeys.sort()

    g = open(dictfile, 'wb')
    first = True
    for key in allKeys :
        # Don't print empty values from dictionary
        if not (nodeDict[key]) :
            continue
        #end if

        # Insert newline in all but first row
        if first :
            first = False
        else :
            g.write("\n")
        #end if


        # Print each (key, value), where value is a list
        # first print the key
        g.write("{}{}".format(key, delim))

        # then print the list of indeces
        run = 0
        templist = np.unique(nodeDict[key])
        for item in templist :
            if run == 0 :
                run = 1
                g.write("{}".format(item))
            else :
                g.write("{}{}".format(",",item))
            #end if
        #end loop
        g.write("\n")

    #end loop
    g.close()

#end def ######## ######## ########






#-# 
#-# ######## ######## ######## ########
#-# # Function: Read in the Knowledge Graph
#-# # Returns: matrix of Edges (N,4), set of Vertices
#-# def readGraphData(datafile, delimiter) :
#-# 
#-# #TODO: np.fromfile should be useful, if I could
#-# #   get it to work ...
#-# #    df = open(datafile, "rb")
#-# #    dt = np.dtype('a21')
#-# #    Edges = np.fromfile(df, dtype=dt, sep=" ", count=3)
#-# ##    Edges = np.fromfile(df, dtype=str, sep=" ", count=-1)
#-# ##    Edges = Edges.reshape( (nLines,4) )
#-# #    print Edges
#-# #    df.close()
#-# 
#-#     # get the number of lines in the file
#-#     nLines = sum( 1 for line in open(datafile, "rb") )
#-# 
#-#     # assign space for edge list
#-#     dt = np.dtype('a23')
#-#     Edges = np.empty( [nLines,4], dtype=dt)
#-# 
#-#     # use a set for faster lookup
#-#     Verts = set()
#-# 
#-#     # Start reading from the file
#-#     df = open(datafile, "rb")
#-# 
#-#     i = 0
#-#     for line in df:
#-#         # extract the data from the file
#-#         line = line.rstrip()
#-#         lv = line.split(delimiter)
#-# 
#-#         # insert into the edge list
#-#         Edges[i] = lv
#-# 
#-# #        print line
#-# #        print lv
#-# 
#-#         # insert into the vertice list
#-#         if lv[0] not in Verts:
#-#             Verts.add(lv[0])
#-#         if lv[1] not in Verts:
#-#             Verts.add(lv[1])
#-# 
#-#         i += 1
#-#     #end loop
#-# 
#-#     # close the data file
#-#     df.close()
#-# 
#-#     print "  file contained {:,} lines".format(nLines)
#-#     return Edges, Verts
#-# #end def ######## ######## ########

######## ######## ######## ########
# Function: Read in the Knowledge Graph
# Returns: set of all human & mouse genes in edge file
def readOnlyGenes(datafile, delimiter) :

    # use a set for faster lookup
    Verts = set()

    # Start reading from the file
    df = open(datafile, "rb")

    for line in df:
        # extract the data from the file
        line = line.rstrip()
        lv = line.split(delimiter)

        # insert into the vertice list
        if (lv[0][0] == 'E') :
            if lv[0] not in Verts:
              Verts.add(lv[0])
        if (lv[1][0] == 'E') :
            if lv[1] not in Verts:
              Verts.add(lv[1])
    #end loop

    # close the data file
    df.close()

    return Verts
#end def ######## ######## ########
# Function: Read in the Knowledge Graph
# Returns: set of all genes in the graph
#   can read from edge files with different gene names
#   based on the first two letters of gene designation
def readOnlyGenesFlex(datafile, delimiter, letters, numlet) :
    # use a set for faster lookup
    Verts = set()
    # Start reading from the file
    df = open(datafile, "rb")
    for line in df:
        # extract the data from the file
        line = line.rstrip()
        lv = line.split(delimiter)
        # insert into the vertice list
        if (lv[0][0:(numlet)] == letters) :
            if lv[0] not in Verts:
              Verts.add(lv[0])
        if (lv[1][0:(numlet)] == letters) :
            if lv[1] not in Verts:
              Verts.add(lv[1])
    #end loop
    # close the data file
    df.close()
    return Verts
#end def ######## ######## ########



######## ######## ######## ########
# Function: Read in the provided dataset (two files)
def readTwoSamples(path, file1, file2) :
    sampleNodes = list()
    # There are two files: _UP.txt and _DN.txt
    sf1 = open(path + file1, "rb")
    for line in sf1 :
        # remove any \n or whitespace to right
        sampleNodes.append( line.rstrip() )
    #end loop
    sf1.close()
    # read file 2
    sf2 = open(path + file2, "rb")
    for line in sf2 :
        sampleNodes.append( line.rstrip() )
    #end loop
    sf2.close()
    # remove the up/dn suffix from sample name
    sampleName = file1.rstrip("_UP.txt")
    sampleName = sampleName.rstrip("_DN.txt")
    # number of genes specified in this dataset
    sampleSize = len(sampleNodes)
    return sampleName, sampleNodes, sampleSize
#end def ######## ######## ########

######## ######## ######## ########
# Function: Read in the provided dataset (single file)
def readOneSample(path, samplefile) :
    sampleNodes = list()
    # Using one of the two files: _UP.txt and _DN.txt
    sf = open(path + samplefile, "rb")
    for line in sf :
        # remove any \n or whitespace to right
        sampleNodes.append( line.rstrip() )
    #end loop
    sf.close()
    # remove the .txt suffix from sample name
    sampleName = samplefile.rstrip(".txt")
    # number of genes specified in this dataset
    sampleSize = len(sampleNodes)
    return sampleName, sampleNodes, sampleSize
#end def ######## ######## ########




#TODO: When randomly selecting nodes,
#    Should I ensure the random selections are unique?
#    requires storage and computation
######## ######## ######## ########
# Function: Randomly select N vertices
def selectRandomNodes(N, VertSet) :

    # make copy: Do not delete from original set
    testSet = set(VertSet)
    # a set allows faster membership checks

    # randomly remove items from testSet
    while (len(testSet) > N) :
        testSet.remove(random.choice( list(testSet) ))

    return testSet
#end def ######## ######## ########




#-# ######## ######## ######## ########
#-# # Function: Count the paths in the test set
#-# def countPathsInSet(testSet, Edges, maxPathLen) :
#-# 
#-#     # convert to list so we can iterate
#-#     testList = list(testSet)
#-#     # will use both set and list below
#-# 
#-#     # Build sub-graph from test set
#-#     G = nx.Graph()
#-#     # add the nodes
#-#     G.add_nodes_from(testList)
#-# 
#-#     # step through edge list
#-#     for e in range(0,Edges.shape[0]) :
#-#         # if both nodes are in the test set
#-#         if (Edges[e,0] in testSet) and (Edges[e,1] in testSet) :
#-#             # then add that edge to the graph
#-#             G.add_edge( Edges[e,0], Edges[e,1], weight=Edges[e,2], type=Edges[e,3] )
#-#         #end if
#-#     #end loop
#-# 
#-#     # Count number simple paths in G
#-#     P = 0
#-# 
#-#     # find paths b/t each pair of vertices
#-#     # Note: edges are undirected
#-#     for i in range( 0, (len(testList)-1) ) :
#-#         for j in range( (i+1), len(testList) ) :
#-#             paths = nx.all_simple_paths(G, source=testList[i], target=testList[j], cutoff=maxPathLen)
#-#             # P is the running total count
#-#             P += len(list(paths))
#-#         #end loop
#-#     #end loop
#-# 
#-#     G.clear()
#-# 
#-#     return P
#-# #end def ######## ######## ########
#-# 
#-# ######## ######## ######## ########
#-# # Function: Count the paths in a gene-only set
#-# def countPathsInOnlyGenes(testSet, Edges, maxPathLen) :
#-# 
#-#     # convert to list so we can iterate
#-#     testList = list(testSet)
#-#     # will use both set and list below
#-# 
#-#     # Build sub-graph from test set
#-#     G = nx.Graph()
#-#     # add the nodes
#-#     G.add_nodes_from(testList)
#-# 
#-#     # step through edge list
#-#     for e in range(0,Edges.shape[0]) :
#-#         # if both nodes are in the test set
#-#         if (Edges[e,0] in testSet) and (Edges[e,1] in testSet) :
#-#             # then add that edge to the graph
#-#             G.add_edge( Edges[e,0], Edges[e,1], weight=Edges[e,2], type=Edges[e,3] )
#-#         # if an edge other than gene-gene
#-#         elif (Edges[e,1] in testSet) and (Edges[e,0][0] != 'E') :
#-#             G.add_edge( Edges[e,0], Edges[e,1], weight=Edges[e,2], type=Edges[e,3] )
#-#         #end if
#-#     #end loop
#-# 
#-#     # Count number simple paths in G
#-#     P = 0
#-# 
#-#     # find paths b/t each pair of vertices
#-#     # Note: edges are undirected
#-#     for i in range( 0, (len(testList)-1) ) :
#-#         for j in range( (i+1), len(testList) ) :
#-#             paths = nx.all_simple_paths(G, source=testList[i], target=testList[j], cutoff=maxPathLen)
#-#             # P is the running total count
#-#             P += len(list(paths))
#-#         #end loop
#-#     #end loop
#-# 
#-#     G.clear()
#-# 
#-#     return P
#-# #end def ######## ######## ########
#-# 
#-# 
#-# 
#-# 
#-# ######## ######## ######## ########
#-# # Function: Draw image of the entire graph
#-# def drawWholeGraph(Edges, Verts, toShow, toFile) :
#-#     A = nx.Graph()
#-#     # add the nodes
#-#     A.add_nodes_from(list(Verts))
#-#     # step through edge list
#-#     for e in range(0,Edges.shape[0]) :
#-#         # add each edge to the graph
#-#         A.add_edge( Edges[e,0], Edges[e,1], weight=Edges[e,2], type=Edges[e,3] )
#-#     #end loop
#-#     # draw the graph
#-#     nx.draw(A)
#-#     if toFile :
#-#         plt.savefig(toFile)
#-#     if toShow :
#-#         plt.show()
#-#     # clear graph
#-#     A.clear()
#-# #end def ######## ######## ########




######## ######## ######## ########
# Function: Pull data from set and move the files
#    Pull from both up/dn if available
#    Then move the file(s) from origpath to newpath
def getSampleBothFiles(origpath, newpath) :

    # Get (sorted) list of everything in directory
    dList = listdir(origpath)
    dList.sort()

    # Keep only the file names
    fNames = list()
    for f in range(0, len(dList)):
        if isfile(origpath + dList[f]):
            fNames.append(dList[f])
        #end if
    #end loop

    # Error: empty directory
    if (len(fNames) < 1) :
        print "The directory is empty."
        sName = ""
        sNodes = list()
        sSize = -1
    # There is only one file left
    elif (len(fNames) == 1) :
        sName, sNodes, sSize = readOneSample(origpath, fNames[0])
        rename(origpath + fNames[0], newpath + fNames[0])
    else :
        
        # get the name of the experiment
        # (remove the up/dn suffix from file name)
        sName = fNames[0].rstrip("_UP.txt")
        sName = sName.rstrip("_DN.txt")

        # get the name of the next file
        sName2 = fNames[1].rstrip("_UP.txt")
        sName2 = sName2.rstrip("_DN.txt")

        # Both are from same sample
        if (sName == sName2) :
            sName, sNodes, sSize = readTwoSamples(origpath, fNames[0], fNames[1])
            rename(origpath + fNames[0], newpath + fNames[0])
            rename(origpath + fNames[1], newpath + fNames[1])
        # This sample has only one file
        else :
            sName, sNodes, sSize = readOneSample(origpath, fNames[0])
            rename(origpath + fNames[0], newpath + fNames[0])
        #end if
    #end if

    return sName, sNodes, sSize
#end def ######## ######## ########

######## ######## ######## ########
# Function: Pull data from a single file in set
#    Pull only from either up or dn at a time
#    Then move the file(s) from origpath to newpath
def getSampleSingleFile(origpath, newpath) :

    # Get (sorted) list of everything in directory
    dList = listdir(origpath)
    dList.sort()

    # Keep only the file names
    fNames = list()
    for f in range(0, len(dList)):
        if isfile(origpath + dList[f]):
            fNames.append(dList[f])
        #end if
    #end loop

    # Error: empty directory
    if (len(fNames) < 1) :
        print "The directory is empty."
        sName = ""
        sNodes = list()
        sSize = -1
#    # There is only one file left
#    elif (len(fNames) == 1) :
#        sName, sNodes, sSize = readOneSample(origpath, fNames[0])
#        rename(origpath + fNames[0], newpath + fNames[0])
    else :
        
        # get the name of the experiment
        # (remove the .txt suffix from file name)
        sName = fNames[0].rstrip(".txt")

#        # get the name of the next file
#        sName2 = fNames[1].rstrip("_UP.txt")
#        sName2 = sName2.rstrip("_DN.txt")

#        # Both are from same sample
#        if (sName == sName2) :
#            sName, sNodes, sSize = readTwoSamples(origpath, fNames[0], fNames[1])
#            rename(origpath + fNames[0], newpath + fNames[0])
#            rename(origpath + fNames[1], newpath + fNames[1])
        # This sample has only one file
#        else :
        sName, sNodes, sSize = readOneSample(origpath, fNames[0])
        rename(origpath + fNames[0], newpath + fNames[0])
#        #end if
    #end if

    return sName, sNodes, sSize
#end def ######## ######## ########


######## ######## ######## ########
# Function: Get sorted list of files in directory
# Return: List of file names
def getFileListSorted(where) :
    # Get (sorted) list of everything in directory
    dList = listdir(where)
    dList.sort()
    # Keep only the file names
    fNames = list()
    for f in range(0, len(dList)) :
        if isfile(where + dList[f]) :
            fNames.append(dList[f])
        #end if
    #end loop
    fNames.sort()

    return fNames
#end def ######## ######## ########
