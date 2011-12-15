#
#  shortcuts.py
#  General utilities that I like to use.
#
#  Created by Gorm Bruun Andresen on 24/11/2011.
#  Copyright (c) 2011 Department of Engineering, University of Aarhus. All rights reserved.
#

#Standard modules
from pylab import *
from scipy import *
import os

def get_positive(x):
    """Replaces all negative values with zeros."""
    return x*(x>0.)  #Possibly it has to be x>1e-10.    
    
def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
    """Wraper for savefig(). Saves figure to path+filename and prints a meassage to the screen."""
    
    figure(fignumber)
    savefig(path + figname, dpi=dpi)
    print 'Saved figure:',path + figname
    sys.stdout.flush()

