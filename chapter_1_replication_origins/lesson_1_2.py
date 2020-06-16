# Use Knuth-Morris-Pratt algorithm to find frequency of the pattern
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
            start = start + patternLength
        else:
            start = start + 1
    return count

print(PatternCount("ACTGAAGGGGTTATATAT", "TT"))