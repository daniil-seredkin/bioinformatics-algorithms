def ReadGenome(genome):
    with open(genome + '.txt', 'r') as file:
        return file.read().replace('\n', '')

genome = ReadGenome('vibrio_cholerae')

# Use some sort of Knuth-Morris-Pratt algorithm to find frequency of the pattern
def PatternCount(text, pattern):
    count = 0
    patternLength = len(pattern)
    last = len(text) - patternLength + 1
    start = 0
    while start < last:
        lastScanIndex = start + patternLength
        matchedCount = 0
        for textIndex, patternIndex in zip(range(start, lastScanIndex), range(patternLength)):
            if (text[textIndex] is pattern[patternIndex]):
                matchedCount = matchedCount + 1
            else:
                break
        if matchedCount is patternLength:
            count = count + 1
            # Take a bigger step to protect us from scanning the same nucleotides more than one time
            start = start + patternLength
        else:
            start = start + 1
    return count

pattern = "TT"
print("%s is met in the vibrio cholerae %s time(-s)" % (pattern, PatternCount(genome, "TT")))

# Finding k-mers, inefficient algorithm
def FrequentWords(text, k):
    frequentPatterns = []
    counts = [None] * (len(text) - k + 1)
    for i in range(len(text) - k):
        pattern = text[i: i + k]
        counts[i] = PatternCount(text, pattern)
    maxCount = max(counts)

    for i in range(len(text) - k):
        if counts[i] is maxCount:
            frequentPatterns[text[i:i + k]] = counts[i]
    return frequentPatterns

# This algorithm has O(N^2) complexity and the genome has more than 1 million nucleotides...
# If you have enough time, then uncomment this
# print("The most frequent k-mer patterns are %s", FrequentWords(genome, 3))

def FrequencyTable(text, k):
    frequencyMap = {}
    for i in range(len(text) - k):
        pattern = text[i: i + k]
        if pattern in frequencyMap:
            frequencyMap[pattern] += 1
        else:
            frequencyMap[pattern] = 1
    return frequencyMap

print ("Frequency map of the genome is %s" % str(FrequencyTable(genome, 3)))

def MaxMap(frequencyMap):
    max = list(frequencyMap.values())[0]
    for key, value in frequencyMap.items():
        if value > max:
            max = value
    return max

def BetterFrequentWords(text, k):
    frequentPatters = {}
    frequencyMap = FrequencyTable(text, k)
    max = MaxMap(frequencyMap)

    for key, value in frequencyMap.items():
        if value == max:
            frequentPatters[key] = value
    return frequentPatters

print ("Frequency map of patterns with maximum occurencies in the genome is %s" % str(BetterFrequentWords(genome, 3)))

def Reverse(pattern):
    complementNucleotideMap = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversedPattern = ''
    for nucleotide in pattern:
        reversedPattern += complementNucleotideMap[nucleotide]
    return reversedPattern

nucleotidePattern = "ATGGCTATA"
print ("Reversed nucleotide pattern for %s is %s" % (nucleotidePattern, Reverse(nucleotidePattern)))


def FindOccurencies(text, pattern):
    patternLength = len(pattern)
    last = len(text) - patternLength + 1
    start = 0
    positions = []
    while start < last:
        lastScanIndex = start + patternLength
        matchedCount = 0
        for textIndex, patternIndex in zip(range(start, lastScanIndex), range(patternLength)):
            if (text[textIndex] is pattern[patternIndex]):
                matchedCount = matchedCount + 1
            else:
                break
        if matchedCount is patternLength:
            positions.append(start)
        start += 1
    return positions

nucleotidePattern = "ATGATCAAG"
genome = ReadGenome('vibrio_cholerae')

print("Pattern %s is found on positions %s in the genome" % (nucleotidePattern, FindOccurencies(genome, nucleotidePattern)))

# This is an optimized version of the algorithms represented in the book
def FindClumps(text, k, L, t):
    patterns = set()
    # It is enough to prepare the frequency map once and then avoid scanning the same positions again and again
    frequencyMap = FrequencyTable(text[0: L + k], k)
    for key, value in frequencyMap.items():
        if value >= t:
            patterns.add(key)

    for start in range(k, len(text) - L):
        firstPatternInPreviousWindow = text[start - 1: start + k - 1]

        frequencyMap[firstPatternInPreviousWindow] -= 1

        if frequencyMap[firstPatternInPreviousWindow] is 0:
            del frequencyMap[firstPatternInPreviousWindow]

        lastPatternInCurrentWindow = text[start + L - k: start + L]

        if lastPatternInCurrentWindow in frequencyMap:
            frequencyMap[lastPatternInCurrentWindow] += 1
        else:
            frequencyMap[lastPatternInCurrentWindow] = 1

        if frequencyMap[lastPatternInCurrentWindow] >= t:
            patterns.add(lastPatternInCurrentWindow)
    return patterns

genome = ReadGenome("E_coli")
print ("E. coli contains the following clumps: %s" % FindClumps(genome, 9, 500, 3))

