import os
import re
import csv
from optparse import OptionParser

sampleNames = []

def main():
    usage =  "Usage: %s --input=<input file> --output=<output prefix>  --sampleInfoFile=<sampleinfo.file> --filenamesColumn=<filenames.column> --class.division.column=<class.division.column>"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--input",
            action="store", type="string", dest="inputFileList",
            help="A file containing a list of paths to HTSeq read counts files")
    parser.add_option("-o", "--output",
            action="store", type="string", dest="outputPrefix",
            help="The base name of the output file")
    parser.add_option("-s", "--sampleInfoFile",
            action="store", type="string", dest="inputSampleInfoFile",
            help="A sample info file containing HTSeq filenames and a column to base a cls file on")
    parser.add_option("-f", "--filenamesColumn",
            action="store", type="string", dest="inputFilenameColumn",
            help="The column in the sample info file that specifies a phenotype to use to assign classes and create a class file for the input samples. This is only relevant if a sample info file is provided.")
    parser.add_option("-c", "--classDivisionColumn",
            action="store", type="string", dest="inputClassColumn",
            help="The column in the sample info file that contains the filenames.  This is used to match up the class division column up to the appropriate input file. This can be either a column index (starting at 0) or a string matching the header of a column in the sample info file.")

    parser.add_option("-n", "--sampleNameColumn",
                      action="store", type="string", dest="sampleNameColumn",
                      help="The column in the sample info file that contains the sample name or the desired column name in the generated gct file.")

    (options, args) = parser.parse_args()

    inputFiles = readFileList(options.inputFileList)
    countTable = loadCountData(inputFiles)
    
    sampleInfo = None
    sampleInfoFilename = nonEmptyString(options.inputSampleInfoFile)
    inputFilenameColumn = nonEmptyString(options.inputFilenameColumn)
    inputClassColumn = nonEmptyString(options.inputClassColumn)
    sampleNameColumn = nonEmptyString(options.sampleNameColumn)

    # do some input validation here if a sample info file is provided
    # if either a sampleNameColumn or a inputClassColun is provided, we als MUST have an inputFilenameColumn
    if (sampleInfoFilename != None):
        if inputClassColumn is not None:
            assert inputFilenameColumn is not None, "An input filename column is needed to create a cls file"
        if sampleNameColumn is not None:
            assert inputFilenameColumn is not None, "An input filename column is needed to remap sample identifiers"

        sampleInfo = readSampleInfo(sampleInfoFilename,inputFilenameColumn, inputClassColumn, sampleNameColumn)
	
	
    mergeCountFiles(countTable, options.outputPrefix, sampleInfo)


def getValidIndex(headers, columnNameOrIndex):
    try:
        # assume its a column name
        index = headers.index(columnNameOrIndex)
    except ValueError:
        # ok, now try to see if its an integer
        try:
            index = int(columnNameOrIndex)
        except ValueError:
            raise ValueError(
                "The column index, " + columnNameOrIndex + ", is invalid or missing from the sample info file.")
    assert index < len(headers), 'The column index is higher than the number of columns in the sample info file.'
    return index

def nonEmptyString(str):
    if (str is not None):
        if len(str.strip()) > 0:
            return str.strip()
    return None

def readSampleInfo(filename, inputFilenameColumn, inputClassColumn, sampleNameColumn):
    # if there is no input filename, we can't do anything with the sample info so bail out even though
    # a file was provided but write a message to stderr
    if inputFilenameColumn is None:
        print("A sampleInfo file was provided but no input filename column was provided so the sampleInfo file will be ignored")
        return None

    sampleInfo = {}
    fileClassMap = {}
    headers = []
    classes = set()
    sampleInfo['writeClassFile'] = False
    sampleNameMap = None

    with open(filename,'rU') as infile:
        reader = csv.reader(infile, delimiter="\t")
        headers = next(reader)
        if inputClassColumn is not None:
            clsIdx = sampleInfo['clsColumnIndex'] = getValidIndex(headers, inputClassColumn)
            sampleInfo['writeClassFile'] = True

        fileIdx = sampleInfo['filenameIndex'] = getValidIndex(headers, inputFilenameColumn)
        sampleNameIdx = -1

        if (sampleNameColumn != None):
            sampleNameIdx = sampleInfo['sampleNameIndex'] = getValidIndex(headers, sampleNameColumn)
            sampleNameMap = {}

        for row in reader:
            # skip blank lines
            if (len(row) == len(headers)):
                if inputClassColumn is not None:
                    className = row[clsIdx].strip().replace(' ', '_')
                    fileClassMap[row[fileIdx]] = className
                    # count the number of unique classes defined in the approriate column
                    classes.add(className)
                # keep track of the sample name if different from the filename
                if (sampleNameColumn != None):
                    sampleNameMap[row[fileIdx]] = row[sampleNameIdx]

    sampleInfo['fileClassMap'] = fileClassMap
    #we want an ordered list of the classes (not a set) so that we can write the cls file properly
    sampleInfo['classes'] = list(classes)
    sampleInfo['numClasses'] = len(classes)
    sampleInfo['sampleNameMap'] = sampleNameMap

    return sampleInfo

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
    fileTables = {}
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
        fileTables[sampleName] =lines

    return fileTables

def mergeCountFiles(countTable, outputPrefix, sampleInfo=None):
    global sampleNames
    if (sampleInfo != None):
        writeClsFile = sampleInfo['writeClassFile']
        # now sort the sample names so that the output writes samples of the same class adjacent to each other
        print(sampleNames)
        print(sampleInfo['fileClassMap'])
        sampleNames = sorted(sampleNames, key=lambda sample: sampleInfo['fileClassMap'][sample])
        print(sampleNames)


    else:
        writeClsFile = False

    #a list of 2d arrays for each sample
    assert len(countTable) > 0, 'Count table is empty'

    assert len(outputPrefix) > 0, 'No prefix given for the output file name'
    outputFileName = outputPrefix + ".gct"

    output = open(outputFileName, "wb")

    maxRows = len(countTable[sampleNames[0]])

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

    if not writeClsFile:
        output.write('\t'.join(sampleNames))
    elif (sampleInfo['sampleNameMap'] is None):
        output.write('\t'.join(sampleNames))
    else:
        sampleNamesMapped = list(map(lambda x: sampleInfo['sampleNameMap'][x], sampleNames))
        output.write('\t'.join(sampleNamesMapped))

    if (writeClsFile ):
        clsFileName = outputPrefix + ".cls"
        clsFileOutput = open(clsFileName, "wb")
        clsFileOutput.write('%s' % (len(sampleNames)));
        clsFileOutput.write('\t')
        clsFileOutput.write('%s' % sampleInfo['numClasses'])
        clsFileOutput.write('\t1\n#\t')
        clsFileOutput.write('\t'.join(sampleInfo['classes']))
        clsFileOutput.write('\n')
        cls = list(map(lambda x: str(sampleInfo['classes'].index(sampleInfo['fileClassMap'][x])), sampleNames))
        clsFileOutput.write('\t'.join(cls))

    for r in range (maxRows):
        firstSampleRowName = countTable[sampleNames[0]][r][0]

        for i in range (len(countTable)):
            sampleTable = countTable[sampleNames[i]]
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



