# Library with paths to the analysis datasets
import os

# A function to find the main path of the project
def findMainPath():
    thisFilePath = os.path.abspath(__file__) # Aboslute path of this file
    thisFilePath = os.path.dirname(thisFilePath) # Directory of this file
    thisFilePath = os.path.dirname(thisFilePath) # One directory up = main path
    return thisFilePath

# Z->tautau datasets
v26Paths = {
"f25503me": ['/Users/user/Documents/AnalysisFW/data']
}

# Set up the example path for the analysis depending on the user machine name.
# If the user is not in the list, add the default path.
username = os.environ['USER']
if username in v26Paths:
    paths = v26Paths[username]
    paths = [findMainPath()+'/data/'] + paths
    #print(paths)
    v26Paths[username] = paths
else:
    v26Paths[username] = [findMainPath()+'/data/']

if __name__ == "__main__":
    print("This file is not meant to be executed --- it is a library of paths for the analysis datasets.")