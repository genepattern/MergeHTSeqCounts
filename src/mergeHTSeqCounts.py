import os
import re
import sys
from optparse import OptionParser

sampleNames = []

def main():
    usage =  "Usage: %s --input=<input file> --output=<output prefix>"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--input",
            action="store", type="string", dest="inputFileList",
            help="A file containing a list of paths to HTSeq read counts files")
    parser.add_option("-o", "--output",
            action="store", type="string", dest="outputPrefix",
            help="The base name of the output file")

    (options, args) = parser.parse_args()

    inputFiles = readFileList(options.inputFileList)
    countTable = loadCountData(inputFiles)
    mergeCountFiles(countTable, options.outputPrefix)

def readFileList(fileList):
    files =[]

    fileReader = open(fileList, 'r')
    for line in fileReader:
        line = line.rstrip('\n')
        line = line.strip()
        result = line.split('\t')
        if len(result) >= 1:
            filepath = result[0]
            if not os.path.exists(filepath):
                sys.stderr.write('File does not exist: ' + filepath + '\n')
                exit(1)

            #if it is a directory then get all the files in that directory
            if os.path.isdir(filepath):
                subFiles = [os.path.join(filepath,fName) for fName in next(os.walk(filepath))[2]]
                subFiles.sort()
                files += subFiles
            else:
                files.append(filepath)
    return files

#files is a list of file paths
def loadCountData(files):
    maxRows = 0
    fileTables = []
    for i in range(len(files)):
        global sampleNames
        sampleName = os.path.basename(files[i])
        filebasename, file_extension = os.path.splitext(sampleName)
        sampleName = filebasename

        #specific to the VIB workshop data - strip out "all_counts"
        #string from the sample name
        sampleName = re.sub(r'_all_counts', "", sampleName)

        sampleNames.append(sampleName)
        lines = []
        with open(files[i]) as textFile:
            for line in textFile:
                if not line.startswith('__'):
                    lines.append(line.split())
        if i == 0:
            maxRows = len(lines)

        if maxRows != len(lines):
            sys.stderr.write('Number of rows differs in file ' + files[i] + '\n')
            exit(1)
        fileTables.append(lines)

    return fileTables

def mergeCountFiles(countTable, outputPrefix):
    #a list of 2d arrays for each sample
    if len(countTable) == 0:
        raise ValueError('Count table is empty\n')

    if len(outputPrefix) == 0:
        raise ValueError('No prefix given for the output file name.\n')

    outputFileName = outputPrefix + ".gct"

    finished = False
    output = None

    try:
        output = open(outputFileName, "wb")
        maxRows = len(countTable[0])

        if maxRows == 0:
            raise ValueError('No rows found in data.\n')
        output.write('#1.2\n')
        output.write('%s' % (maxRows))
        output.write('\t')
        output.write('%s' % (len(sampleNames)))
        output.write('\n')
        output.write('Name')
        output.write('\t')
        output.write('Description')
        output.write('\t')
        output.write('\t'.join(sampleNames))

        for r in range (maxRows):
            firstSampleRowName = countTable[0][r][0]

            for i in range (len(countTable)):
                sampleTable = countTable[i]

                if(len(sampleTable[r]) != 2):
                    raise ValueError('Expecting two columns found %s. Please check that the input is valid.\n' % (len(sampleTable[r])))

                currentSampleRowName = sampleTable[r][0]

                #check that the feature names are the same in each file
                if firstSampleRowName != currentSampleRowName:
                    raise ValueError('Error: row names differ in row %s' % (r+1) + '. ' +  firstSampleRowName + ' ' + currentSampleRowName + '\n')

                if i == 0:
                    output.write('\n')
                    output.write(firstSampleRowName)
                    #put the rowName also in the description column
                    output.write('\t')
                    output.write(firstSampleRowName)

                output.write('\t')

                readCount = sampleTable[r][1]
                output.write(readCount)
                finished = True
    except ValueError, err:
        sys.stderr.write(str(err))
        exit(1)
    finally:
        if output != None:
            output.close()
        if finished != True:
            #delete any partially created gct file
            os.remove(outputFileName)

main()



