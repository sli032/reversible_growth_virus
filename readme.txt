Version installed:

python = 3.5.6
pandas = 0.19.0
hoomd = 2.3.1


Run Code:
1.Generate parameter files:
  python2 param.py parent_param.yaml
2.Change runscript.sh file to run specific paramfiles
3.run shell script: 
  ./runscript.sh (all/allnew/filename)
  #if permission denied, run: chmod +x runscript.sh
  all: run all files regardless it is ran before or not
  allnew: only run new created files
  filename: explicitely point out the filename
