import os
import re
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
            assert os.path.exists(filepath), 'File does not exist: ' + filepath

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

        assert maxRows == len(lines), '\nNumber of rows differs in file ' + files[i]
        fileTables.append(lines)

    return fileTables

def mergeCountFiles(countTable, outputPrefix):

    #a list of 2d arrays for each sample
    assert len(countTable) > 0, 'Count table is empty'

    assert len(outputPrefix) > 0, 'No prefix given for the output file name'
    outputFileName = outputPrefix + ".gct"
    output = open(outputFileName, "wb")
    maxRows = len(countTable[0])

    assert maxRows > 0, 'No rows found in data'
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
            assert len(sampleTable[r]) == 2, 'Expecting two columns found ' + len(sampleTable[r])

            currentSampleRowName = sampleTable[r][0]

            #check that the feature names are the same in each file
            assert firstSampleRowName == currentSampleRowName, 'Error: row names differ in row %s' % (r+1) + '. ' +  firstSampleRowName + ' ' + currentSampleRowName

            if i == 0:
                output.write('\n')
                output.write(firstSampleRowName)
                #put the rowName also in the description column
                output.write('\t')
                output.write(firstSampleRowName)

            output.write('\t')

            readCount = sampleTable[r][1]
            output.write(readCount)
    output.close()

main()



