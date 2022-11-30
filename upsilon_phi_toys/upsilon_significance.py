#!/usr/bin/env python
# coding: utf-8

# In[101]:


from goofit import *
import matplotlib.pyplot as plt
from scipy.special import binom
from scipy import integrate
import numpy as np
import contextlib
import os
from datetime import datetime
from tqdm import tqdm
import argparse
import sys


# In[99]:


def is_notebook() -> bool:
    print("here")
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


# In[81]:


start = np.floor(datetime.now().timestamp())


# In[98]:


if(is_notebook()):
    print("This is a notebook! Values will be fixed to ")
    
    class dummy():
        loops = 0
        events = 0
        goodnll = 0

    args = dummy()
    args.loops = 1000
    args.events = 18243+509
    args.goodnll = 20
    
    print("----> n loops   : %d"%(args.loops))
    print("----> n events  : %d"%(args.events))
    print("----> good nll  : %d"%(args.goodnll))
    
else:
    
    print("This is not a notebook! The inputs are important then!")
    
    parser = argparse.ArgumentParser(description='Toys for statistical significance')
    parser.add_argument('-l', '--loops',default=1000) 
    parser.add_argument('-e', '--events',default=18243+509) 
    parser.add_argument('-g', '--goodnll',default=20) 
    args = parser.parse_args()


# ### Variables

# In[6]:


aS = []
aS.append(Variable("a0",0.00214,0.00214-3*0.00080,0.00214+3*0.00080)) #0.00080
aS.append(Variable("a1",0.0115,0.0115-3*0.0042,0.0115+3*0.0042)) #0.0042
aS.append(Variable("a2",0.0094,0.0094-3*0.0035,0.0094+3*0.0035)) #0.0035
aS.append(Variable("a4",0.0127,0.0127-3*0.0046,0.0127+3*0.0046)) #0.0046
mean = Variable("mean",10.51225,10.51225-0.00021*3,10.51311+0.00021*3) #0.0046
sigma = Variable("sigma",0.0017,0.0017-0.00023*3,0.0017+0.00023*3) #0.0046


# In[ ]:


xmin = 10.48
xmax = 10.55
x_points = np.linspace(xmin,xmax,1000)


# ### Generation functions and parameters

# In[8]:


def bern_func(x,i=0,n=2):
    
    r = binom(n,i)*(x**i)*((1-x)**(n-i))
    #val = val + r
    return r

def bern_func_gen(x,lowX,uppX,i=0,n=2):
    xNorm = (x-lowX)/(uppX - lowX)
    r = bern_func(xNorm,i=i,n=n)
    return r

def gen_bern(x,lowX,uppX):
    r = 0
    xNorm = (x-lowX)/(uppX - lowX)
    for j in range(0,deg+1):
        r = r + aS_gen[j] * bern_func_gen(x,xmin,xmax,i=j,n=deg)
    return r

aS_gen = []
aS_gen.append(0.00214) #0.00080
aS_gen.append(0.0115) #0.0042
aS_gen.append(0.0094) #0.0035
aS_gen.append(0.0127) #0.0046

deg=3

y_integral = integrate.romberg(gen_bern, xmin, xmax, args=(xmin,xmax), show=False)

y_points = gen_bern(x_points,xmin,xmax)
y_max = y_points.max()
y_points = y_points/y_integral

binw = (xmax - xmin)/1000.0


# ### P.d.f.s definitions

# In[10]:


xvar = Observable("x",xmin,xmax)
bern = BernsteinPdf("bern", xvar,aS,0)
signal = GaussianPdf("gaus",xvar,mean,sigma)


# ## The MC Toys Loops

# In[14]:


print("Starting the MC Toys Run")
print("> n loops   : %d"%(args.loops))
print("> n events  : %d"%(args.events))
print("> good nll  : %d"%(args.goodnll))


# In[106]:


n_events = int(args.events)
for i in tqdm(range(args.loops)):
    
    try:
   
 
        ##this avoids to printout the fit
        #with open(os.devnull, 'w') as null, contextlib.redirect_stdout(null), contextlib.redirect_stderr(null):
        with open(os.devnull, "w", encoding='utf-8') as target:
            sys.stdout = target
            sys.stderr = target
            sys.stdin = target
            dats = []

            for i in range(n_events):
                this_x = np.random.uniform(xmin,xmax)
                this_y = gen_bern(this_x,xmin,xmax)
                this_roll = np.random.uniform(0.0,y_max*1.01)
                if this_y > this_roll:
                    dats.append(this_x)


            now = np.floor(datetime.now().timestamp())

            np.random.seed()

            data = UnbinnedDataSet(xvar)
            data.from_numpy([dats],filter=True)

            bern.setData(data)

            fitter = FitManager(bern)
            fitter.setVerbosity(-1)
            nll_bkg = fitter.fit().Fval()

            data = UnbinnedDataSet(xvar)
            data.from_numpy([dats],filter=True)

            sigfrac = Variable("sigFrac", 0.01, 0.0, 1.00)
            total = AddPdf("total", [sigfrac], [signal, bern])

            total.setData(data)

            fitter = FitManager(total)
            nll_sig = fitter.fit().Fval()
            fitter.setVerbosity(-1)
            
            delta_nll = -(nll_sig-nll_bkg)

            with open("delta_nlls%d.txt"%(start),"a") as F:
                    np.savetxt(F,[delta_nll])

            if(delta_nll)>args.goodnll:

                with open("data/%.2f_%d.txt"%(delta_nll,now),"wb") as F:
                    np.savetxt(F,dats)

                numbins = 200
                binw = (xmax-xmin)/numbins

                figure = plt.figure(figsize=(12,15))
                plt.hist(dats, bins=numbins, label='Data', density=False,histtype="step",lw="2",color="black")

                grid = bern.makeGrid()
                bern.setData(grid)

                main = np.array(bern.getCompProbsAtDataPoints()[0])
                main = main*float(len(dats))*binw

                xvals = grid.to_matrix().flatten()

                plt.plot(xvals, main, label='H0 Fit',lw=3,color="blue")

                total.setData(grid)
                main = np.array(total.getCompProbsAtDataPoints()[0])
                main = main*float(len(dats))*binw

                plt.plot(xvals, main, "--",label='H1 Fit',lw=3,color="red")

                xvals = grid.to_matrix().flatten()

                plt.legend(fontsize=18,loc="lower right")

                ax = plt.gca()
                ax.margins(x=0,tight=False)
                plt.ylim(0,)
                ax.tick_params(axis = 'both', which = 'major', labelsize = 16)

                plt.title("Generated Data + Goofit Fit",fontsize=18);
                plt.savefig("plots/%.2f_%d.png"%(delta_nll,now))
    
    finally:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        sys.stdin = sys.__stdin__


# In[ ]:




