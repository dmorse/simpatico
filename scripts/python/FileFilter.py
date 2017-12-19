import os
import re
from os.path import *
from string  import *
from Directory import *

class FileFilter:

   def __init__(self):
      self.nYesFilter = 0
      self.yesFilters = []
      self.yesResults = []
      self.nNoFilter = 0
      self.noFilters = []
      self.noResults = []

   def addYesFilter(self, filter):
      self.yesFilters.append(re.compile(filter))
      self.yesResults.append(0)
      self.nYesFilter += 1

   def addNoFilter(self, filter):
      self.noFilters.append(re.compile(filter))
      self.noResults.append(0)
      self.nNoFilter += 1

   def clear(self):
      self.nYesFilter = 0
      self.yesFilters = []
      self.yesResults = []

   def clearResults(self):
      for i in range(self.nYesFilter):
         self.yesResults[i] = 0
      for i in range(self.nNoFilter):
         self.noResults[i] = 0

   def matchFile(self, filename):
    
      # Open and read the file 
      file = open(filename, 'r')
      lines = file.readlines()
      n  = len(lines)
      file.close()

      # Check for lines that match all filters
      self.clearResults() 
      for line in lines:
         for i in range(self.nYesFilter):
            if (self.yesFilters[i].search(line)):
               self.yesResults[i] += 1
         for i in range(self.nNoFilter):
            if (self.noFilters[i].search(line)):
               self.noResults[i] += 1
   
      # Success iff all yes-filters match and no no-filters 
      success = True 
      for i in range(self.nYesFilter):
         if (self.yesResults[i] == 0):
            success = False
      for i in range(self.nNoFilter):
         if (self.noResults[i] > 0):
            success = False

      return success 

   def matchFiles(self, dirName, pattern):
      dir = Directory(dirName)
      filenames = dir.filenames(pattern)

      matches = []
      for filename in filenames:
         if isfile(filename): 
            success = self.matchFile(filename)
            if (success):
               matches.append(filename)
      return matches
