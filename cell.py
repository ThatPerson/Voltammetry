import math
import logging
import sys
import copy
import re
import numpy

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class RateConstant:
    """ Dependent variable rate constant calculation
            k0 is the rate constant at a dependent concentration of 0
            dep is the dependent concentration - for PFS this is 
                    the urea concentration
        Methods;
            k(v) - returns the rate constant for concentration v
    """
    k0 = 0
    dep = 0
    k = 0
    def __init__(self, k0, dep):
        self.k0 = k0
        self.dep = dep
        #self.k = k0
        #self.calculate_k(0)
    def calculate_k(self, v):
        self.k = self.k0 * math.exp(self.dep * v)
        return self.k
class Cell:
    """ On initialising takes in variables
        contents - array of the names of the components
            eg {"D": 100, "I": 0, "N": 0}
        reactions - array of the reactions which can occur
            eg [{"reactants": [{"symbol": "D", "V": 1}], 
            "products": [{"symbol": "I", "V": 1}], 
            "k": RateConstant(23000, -1.68)},...]
            In these V is the stoichiometry.
            I chose to put the reactions in this manner because while
             A + B -> C + D 
            may be more readable, the computer would have to store it
            anyway and I think this is a nice balance.
        v - factor exponentially impacting rate constant - eg [urea]
        
        Methods:
            iterate(dt) - performs one iteration.
            output(headers, frac) - returns an output string
            compare(other) - works out the difference in total quantity
                    between each pair.
    """
    contents = {}
    reactions = {}
    v = 0
    def __init__(self, components, reactions, v=0):
        self.contents = (components)
        self.reactions = (reactions)
        self.v = v
        for i in self.reactions:
            i['k'].calculate_k(self.v)

    
    def iterate(self, dt):
        """ 
            Performs one iteration of the model, moving forward a time
            interval given by dt. Generally dt should be set to around 1e-6. 
            """
        current = dict(self.contents) # Current is a copy of the contents so 
                                    #that the new values don't affect the old.
        q = 0
        for i in self.reactions:
            q = i["k"].k
            # The form of these is V k [A]^vA [B]^vB.
            # As the [A]^vA... bit is the same for reactants and products we 
            #precalculate this as q.
            for p in i["reactants"]: 
                try:
                    q = q * math.pow(self.contents[p["symbol"]], p["V"]) 
                except OverflowError:
                    # This can happen if the values in above are far too large.
                    logger.debug(self.contents)
                    logger.debug("%f %f\n", self.contents[p["symbol"]],p["V"])
                    logger.info("Overflow Error - infinity detected. \
                        Try a lower time step.")
                    exit()
                # This represents the [A]^vA [B]^vB... portion
            for p in i["reactants"]:
                # d[n]/dt = -vn k q
                dndt = -p["V"] * q
                current[p["symbol"]] = current[p["symbol"]] + dt * dndt
            for p in i["products"]:
                dndt = p["V"] * q
                current[p["symbol"]] = current[p["symbol"]] + dt * dndt
                
        t_diff = 0
        sel = self.contents
        for sd in current:
            t_diff += abs(current[sd] - sel[sd])
            if (t_diff > 0.00001):
                self.contents = dict(current)
                return 0
        self.contents = dict(current)
        return 1
    
    def compare(self, other, tol):
        """ Takes the last time step's output and compares to see how 
        different they are. """
        t_diff = 0
        sel = self.contents
        for sd in sel:
            t_diff += abs(sel[sd] - other[sd]) 
            if (t_diff > tol):
                return 0
        return 1
    
    def output(self, headers=1, frac = 1):
        """
            headers = 1 prints out headers, =0 doesn't.
            frac = 1 prints out fractions (eg 0.2, 0.4, 0.4), 
                   0 prints out actual concentrations
            Returns string.
            """
            
        return_str = ["v"]
        head = ''
        
        if (headers == 1):
            return_str += [i for i in self.contents]
            head = ', \t'.join(return_str) + '\n'
        pri = [str(self.v)]
        divisor = 1
        if (frac == 1):
            divisor = 0
            for i in self.contents:
                divisor = divisor + self.contents[i]
            if (divisor == 0):
                divisor = math.inf
        pri += [str(self.contents[i]/divisor) for i in self.contents]
        
        return head + ',\t'.join(pri) + "\n"
    

    


def read_file(fl):
    mode = 0
    contents = {}
    reactions = []
    e_reactions = []
    peake = 0
    scanrate = 0
    output = 1
    diff_coeff = 1/3.
    cells = 1
    with open(fl, 'r') as fi:
        for i in fi:
            
            if (i.strip() == "CONTENTS"):
                mode = 1
                #print(mode)
            if (i.strip() == "ECHEM"):
                mode = 3
            if (i.strip() == "REACTIONS"):
                mode = 2
            if (i.strip() == "EREACTIONS"):
                mode = 4
            if (mode == 1):
                kl = i.strip().split('=')
                print(kl)
                
                if (len(kl) == 2):
                    contents[kl[0].strip()] = float(kl[1])
            if (mode == 2):
                i = i.strip()
                try: 
                    k = float(re.search('\(([0-9.]+?)\)', i).group(1))
                except AttributeError:
                    continue
                print(k)
                i = i.replace(' ', '')
                preact = i.split('>')
                if (len(preact) == 2):
                    reactants = []
                    products = []
                    #reactants = preact[0].split('+')
                    for r in preact[0].split('+'):
                        reactants.append({'symbol':r, 'V':1})
                    preact[1] = preact[1][:preact[1].find('(')]
                    #products = preact[1].split('+')
                    for p in preact[1].split('+'):
                        products.append({'symbol':p, 'V':1})
                    reactions.append({'reactants':reactants, 'products': products, 'k': RateConstant(k, 0)})
            if (mode == 3):
                if (i[:len("PEAKE")] == "peake"):
                    peake = float(i[len("PEAKE")+1:])
                elif (i[:len("SCANRATE")] == "scanrate"):
                    scanrate = float(i[len("SCANRATE")+1:])
                elif (i[:len("CELLS")] == "cells"):
                    cells = int(i[len("CELLS")+1:])
                elif (i[:len("OUTPUT")] == "output"):
                    output = int(i[len("OUTPUT")+1:])
                elif (i[:len("diffusion")] == "diffusion"):
                    diff_coeff = float(i[len("DiFfusion")+1:])
            if (mode == 4):
                i = i.strip()
                l = i.split(",")
                if (len(l) < 2):
                    continue
                ox = l[0].strip()
                red = l[1].strip()
                pot = float(l[2])
                e_reactions.append({"oxid": ox, "red": red, "pot": pot})

    print(contents)
    print(reactions)
    return [contents, reactions, [peake, scanrate, cells, output, diff_coeff], e_reactions]

def FileReactor(c, r, s, t, f, q, e):
    ''' c and r define the contents and reactions respectively.
    s gives the time that the model will run, t gives the timestep.
    f gives the output file name, and q gives the electrochemical
    information ([peake, scanrate, cells, output]) '''
    fractions = []
    E = q[0]
    peake = q[0]
    scanrate = -q[1]
    output = q[3]
    d_c = q[4]
    I = 0
    for i in range(0, q[2]):
        fractions.append(Cell(c, r, 0)) # The cell doesn't know anything about EC
    with open(f, "w+") as f:
        
        rea = []
        for i in c:
            rea.append(i)
        f.write("t, n, E, I, " + ','.join(rea)+"\n")
        
        
        for i in range(0, int(s / t)): 
            if (E <= -peake):
                scanrate = abs(scanrate)
            elif (E >= peake):
                scanrate = -abs(scanrate)
            E = E + scanrate
            if (i % (1/t) == 0):
                logger.info("Done %d iterations..." % (i))
            
            
            # Here we need to;
            # 1. Set concentration of Ox and Red in cell 0
            # 2. Run Kinetic Model on all cells
            # 3. Diffuse between cells
            
            # Set concentration of Ox and Red in Cell 0
            I = 0
            for e_re in e:
                # e_re is of the form {"oxid":, "red":, "pot":}
                p = numpy.exp(E - e_re["pot"])
                print(fractions[0].contents)
                q = fractions[0].contents[e_re["oxid"]] + fractions[0].contents[e_re["red"]]
                new_ox = q - (q / (p+1))
                new_red = q / (p+1)
                I = I + fractions[0].contents[e_re["red"]] - new_red
                fractions[0].contents[e_re["oxid"]] = new_ox
                fractions[0].contents[e_re["red"]] = new_red
            
            # Run Kinetic Model
            for cell in range(0, len(fractions)):
                fractions[cell].iterate(t)

                if (cell % output == 0 and cell == 0):
                    f.write("%d, %d, %f, %f, %s" % (i, cell, E, I, fractions[cell].output(0, 0)))
            
            # Diffusion
            frac = copy.deepcopy(fractions)
            for cell_n in range(0, len(fractions)):
                l_n = cell_n - 1
                p_n = cell_n + 1
                if (l_n < 0):
                    l_n = 0
                if (p_n >= len(fractions)):
                    p_n = len(fractions) - 1
                for p in fractions[0].contents:
                    # Can't finetune extent of diffusion yet TODO
                    fractions[cell_n].contents[p] = (d_c * (frac[l_n].contents[p] +
                                                        frac[p_n].contents[p]) + 
                                                    (1 - 2*d_c) * frac[cell_n].contents[p])
                 #   fractions[cell_n].contents[p] = (frac[cell_n].contents[p] + 
                  #                                   frac[l_n].contents[p] + 
                   #                                  frac[p_n].contents[p])/3
                 
                    
                
                
                
    return 1 

def main():
    """ The main code, we only run this if we are being run directly
        so the program can be imported """
        
        
    if (len(sys.argv) < 2):
        logger.error("Requires at least one argument")
        exit()
        
    if (sys.argv[1] == "help"):
        print("This program takes arguments in the form;")
        print("  python3 cell.py MODEL OPTS")

        exit()
        
    elif (sys.argv[1] == "file"):
        [c, r, q, e] = read_file(sys.argv[2])
        seconds = float(sys.argv[3])
        timestep = float(sys.argv[4])
        FileReactor(c, r, seconds, timestep, sys.argv[5], q, e)
        # Take length and seconds as other arguments
        #
if __name__ == "__main__":
    main()
