import pandas as pd

import uproot

from matplotlib import pyplot as plt

import numpy as np

import os

import sys

import argparse

#from string import split

json_2018 = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
runfill = lambda x :  float(x.split(":")[0])
fillrun = lambda x :  float(x.split(":")[1])

parser = argparse.ArgumentParser()
parser.add_argument('--lumi',  type=str, default=None,
                    help="Lumi section table path.")
parser.add_argument('--data',  type=str, default="./data.root",
                    help="ROOT file to be loeaded name.")
parser.add_argument('--dir',   type=str, default=None,
                    help="Directory, in ROOT file, where the tree is stored. If none, leave it blank.")
parser.add_argument('--tree',  type=str, default="tree",
                    help="Tree to be loaded.")
parser.add_argument('--run',  type=str, default="run",
                    help="Run leaf/column name (default = \'run\'")
parser.add_argument('--keys',  nargs="+",type=str, default=None,
                    help="Run leaf/column name (default = \'run\'")
parser.add_argument('--brilcalc', action='store_true',
                    help="Running brilcalc. If not a lumi file should be provided."))
parser.add_argument('--json',type =str, default=None)
args = parser.parse_args()

RUN = args.run
lumi_path = args.lumi

DUPLICATES = args.keys
print(DUPLICATES)


#############################

if not args.brilcalc and args.lumi is None:
    print("Select a luminosity table produced with brilcalc \n
           (with --lumi) or run --brilcalc to run the brilcalc \n
           tool.")
    sys.exit(1)

#############################

if args.brilcalc:
    print("> Running brilcalc:")
    if args.json is None:
        print(">> No Lumisection json file selected. Downloading default one (2018 Golden).")
        os.system("curl -o lumi.json " + json_2018)
        JSON_FILE = "lumi.json"
        os.system("brilcalc lumi -i lumi.json --output-style csv -o runs.csv")
    else:
        os.system("brilcalc lumi -i " + args.json + " --output-style csv -o runs.csv")
    print(">> The input lumisection table (" + args.lumi +") is being \n
              overwritten by the newly calculated one (runs.json")

    lumi_path = "runs.csv"

#############################

print("> Loading BrilCalc Lumi Table")

runs = pd.read_csv(args.lumi,delimiter=",")
runs[RUN] = runs["run:fill"].apply(runfill)
runs["fill"] = runs["run:fill"].apply(fillrun)
runs = runs.drop(["run:fill"],axis=1)

print("> Loading ROOT file into Pandas table")
tree = uproot.open(args.data)
if args.dir is not None:
    tree = tree[args.dir][args.tree]
else:
    tree = tree[args.tree]

tree = tree.pandas.df([RUN])


plt.figure(figsize=(12,9))

min_run = tree[RUN].min()
max_run = tree[RUN].max()

n,b = np.histogram(tree[RUN].values.astype(float),bins=1000,range=(min_run-200,max_run+200))
bw = (b[1]-b[0])*0.5
plt.errorbar(b[:-1] + bw,n,c="red",yerr=np.sqrt(n));

ax = plt.gca()
ax.margins(x=0,tight=False)
plt.ylim(0,)
ax.tick_params(axis = 'both', which = 'major', labelsize = 16)
plt.grid(True,axis="y")
plt.xlabel("Run Number",fontsize=17,fontstyle="italic")
plt.ylabel("Entries$",fontsize=17,fontstyle="italic")
#legend = plt.legend(title="Eras",fontsize=18,loc=2)
#plt.setp(legend.get_title(),fontsize=18,fontweight="bold")
#plt.barh(xnumbers,new_mean[99])
plt.title("Run Counts",fontsize=18,fontweight="bold")

plt.savefig('run_count.png')


counts_df = tree[RUN].value_counts().astype(float)
counts_df = pd.DataFrame(np.array([counts_df.index.astype(int),counts_df.values]).transpose(1,0),columns=[RUN,"counts"])
counts_df = counts_df[counts_df[RUN]>1]


runs = pd.merge(counts_df, runs, on=RUN)
runs["ratio"] = runs["counts"]*1000.0/runs["recorded"] # from mb to µb

xnames = runs[RUN]

plt.figure(figsize=(12,9))
ax = plt.gca()

plt.plot(xnames,runs["ratio"].values,"o",color="blue");
ax.tick_params(axis = 'both', which = 'major', labelsize = 16)
plt.xlabel("Run Number",fontsize=17,fontstyle="italic")
plt.ylabel("Entries / $(µb)^{-1}$",fontsize=17,fontstyle="italic")
plt.ylim(0,)
#legend = plt.legend(title="",fontsize=18,loc=2)
#plt.setp(legend.get_title(),fontsize=18,fontweight="bold")
#plt.barh(xnumbers,new_mean[99])
plt.title("Candidates Count",fontsize=22,fontstyle="italic")

plt.savefig('cand_count.png')
